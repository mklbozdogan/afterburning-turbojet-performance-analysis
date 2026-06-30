function turbojet()

    % === Inputs ===
    altitude = input('Altitude[m] (ex: 3048): ');
    Tt4 = input('Tt4 [K] (ex: 1400): ');
    Tt7 = input('Tt7 [K] (ex: 2000): ');
    Mmin = input('Initial Mach (ex: 0.5): ');
    Mmax = input('Final Mach (ex: 2.0): ');
    Mstep = input('Mach step (ex: 0.5): ');
    pic_start = input('Initial π_c (ex: 8): ');
    pic_end = input('Final π_c (ex: 24): ');
    pic_step = input('π_c step (ex: 4): ');

    Mvec = Mmin:Mstep:Mmax;
    pic_range = pic_start:pic_step:pic_end;

    % === Constants ===
    gamma = 1.4;
    cp = 1004.83;
    hPR = 42800000;
    R = cp * (1 - 1/gamma);

    % === Load atmosphere data ===
    atm_data = readtable('atmosphere.xlsx');
    h_vec = atm_data.h;
    T_vec = atm_data.T;
    P_vec = atm_data.P;

    [h_vec, ia, ~] = unique(h_vec, 'stable');
    T_vec = T_vec(ia);
    P_vec = P_vec(ia);

    theta_vec = T_vec / 288.15;
    delta_vec = P_vec / 101325;

    % === Interpolate for altitude ===
    theta = interp1(h_vec, theta_vec, altitude, 'linear', 'extrap');
    delta = interp1(h_vec, delta_vec, altitude, 'linear', 'extrap');

    T0 = theta * 288.15;
    a0 = sqrt(gamma * R * T0);

    % === Preallocate ===
    [X, ~] = meshgrid(pic_range, Mvec);

    Fm0 = zeros(size(X));
    f1 = zeros(size(X));
    etaP = zeros(size(X));
    etaT = zeros(size(X));
    etaO = zeros(size(X));
    TSFC = zeros(size(X));

    % === Main loop ===
    for i = 1:numel(Mvec)

        M0 = Mvec(i);

        tau_r = 1 + (gamma - 1)/2 * M0^2;
        tau_lambda = Tt4 / T0;
        tau_lambdaAB = Tt7 / T0;

        for j = 1:numel(pic_range)

            piC = pic_range(j);

            tau_c = piC^((gamma - 1)/gamma);

            numerator = tau_lambda / (tau_r * tau_c);
            denominator = tau_lambda - tau_r * (tau_c - 1);

            V9_a0_sq = (2 / (gamma - 1)) * tau_lambdaAB * ...
                       (1 - numerator / denominator);

            V9_a0_sq = max(V9_a0_sq, 0);
            V9_a0 = sqrt(V9_a0_sq);

            % Specific thrust
            Fm0(i,j) = a0 * (V9_a0 - M0);

            % Main combustor fuel-air ratio
            f1(i,j) = (cp * T0 / hPR) * ...
                      (tau_lambda - tau_r * tau_c);

            % Total fuel-air ratio with afterburner
            f_total = (cp * T0 / hPR) * ...
                      (tau_lambdaAB - tau_r);

            % Propulsive efficiency
            etaP(i,j) = (2 * M0) / (V9_a0 + M0);

            % Thermal efficiency
            etaT(i,j) = ((gamma - 1) * cp * T0 * ...
                        (V9_a0_sq - M0^2)) / ...
                        (2 * f_total * hPR);

            % Overall efficiency
            etaO(i,j) = etaP(i,j) * etaT(i,j);

            % TSFC
            TSFC(i,j) = f_total / Fm0(i,j);

        end
    end

    % === PLOTS ===

    figure;
    plot(pic_range, Fm0', 'LineWidth', 1.5);
    title('F/ṁ₀ vs \pi_c');
    xlabel('\pi_c');
    ylabel('F/ṁ₀ [m/s]');
    legend(arrayfun(@(x) sprintf('M = %.1f', x), Mvec, ...
        'UniformOutput', false), 'Location', 'best');
    grid on;

    figure;
    plot(pic_range, TSFC' * 1e6, 'LineWidth', 1.5);
    title('Thrust Specific Fuel Consumption (TSFC) vs \pi_c');
    xlabel('\pi_c');
    ylabel('TSFC [mg/(N·s)]');
    legend(arrayfun(@(x) sprintf('M = %.1f', x), Mvec, ...
        'UniformOutput', false), 'Location', 'best');
    grid on;

    figure;
    plot(pic_range, f1', 'LineWidth', 1.5);
    title('Fuel/Air Ratio f vs \pi_c');
    xlabel('\pi_c');
    ylabel('f');
    legend(arrayfun(@(x) sprintf('M = %.1f', x), Mvec, ...
        'UniformOutput', false), 'Location', 'best');
    grid on;

 figure;
plot(Mvec, f1, 'LineWidth', 1.5);
title('Fuel/Air Ratio f vs Mach Number');
xlabel('Mach Number');
ylabel('f');
legend(arrayfun(@(x) sprintf('\\pi_c = %.2f', x), pic_range, ...
    'UniformOutput', false), 'Location', 'best');
grid on;

    figure;
    plot(pic_range, etaT' * 100, 'LineWidth', 1.5);
    title('Thermal Efficiency \eta_T vs \pi_c');
    xlabel('\pi_c');
    ylabel('\eta_T (%)');
    legend(arrayfun(@(x) sprintf('M = %.1f', x), Mvec, ...
        'UniformOutput', false), 'Location', 'best');
    grid on;

        figure;
    plot(pic_range, etaP' * 100, 'LineWidth', 1.5);
    title('Propulsive Efficiency \eta_P vs \pi_c');
    xlabel('\pi_c');
    ylabel('\eta_P (%)');
    legend(arrayfun(@(x) sprintf('M = %.1f', x), Mvec, ...
        'UniformOutput', false), 'Location', 'best');
    grid on;

    figure;
    plot(pic_range, etaO' * 100, 'LineWidth', 1.5);
    title('Overall Efficiency \eta_O vs \pi_c');
    xlabel('\pi_c');
    ylabel('\eta_O (%)');
    legend(arrayfun(@(x) sprintf('M = %.1f', x), Mvec, ...
        'UniformOutput', false), 'Location', 'best');
    grid on;

end