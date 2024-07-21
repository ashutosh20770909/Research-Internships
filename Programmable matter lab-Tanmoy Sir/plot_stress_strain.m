function plot_stress_strain()
    % Constants  I have taken them by myself
    Rl = 0.0002;
    Q0 = 0;
    E = 1200000000;
    I = 0.0002;
    A0 = 0.0002;
    
    % Vary N0 such that sigma_eff/E ranges from 0 to 0.08 as we are trying
    % to match graph 5b
    sigma_eff_target = linspace(0, 0.08, 100) * E;
    N0_values = sigma_eff_target * A0;
    
    % Different values of theta_l to plot
    theta_l_values = [90, 120, 150, 180];
    colors = lines(length(theta_l_values)); 
    % Plotting
    figure;
    hold on;
    
    for j = 1:length(theta_l_values)
        theta_l = deg2rad(theta_l_values(j));
        
       
        sigma_eff_values = zeros(size(N0_values));
        epsilon_eff_values = zeros(size(N0_values));
        
        %  calculating different sigma_eff and epsilon_eff
        for i = 1:length(N0_values)
            N0 = N0_values(i);

            % Solve for beta_I and beta_II
            H = @(beta_I, beta) H_function(beta_I, beta, N0, Q0, E, I, A0, Rl);
            d_beta_dS = @(beta_I, beta) H(beta_I, beta).^2;
            integral_eq = @(beta_I, beta_II) integral(@(beta) sign(d_beta_dS(beta_I, beta)) ./ H(beta_I, beta), beta_I, beta_II) + theta_l * Rl;

            options = optimset('Display', 'off'); % Suppress fsolve output
            solution = fsolve(@(beta) integral_eq(beta(1), beta(2)), [pi, -pi], options);
            beta_I_sol = solution(1);
            beta_II_sol = solution(2);

            % Calculate xC and yC
            integrand_xC = @(beta) sign(d_beta_dS(beta_I_sol, beta)) .* (1 ./ H(beta_I_sol, beta)) .* (1 + (N0 * cos(beta) + Q0 * sin(beta)) / (E * A0)) .* cos(beta);
            integrand_yC = @(beta) sign(d_beta_dS(beta_I_sol, beta)) .* (1 ./ H(beta_I_sol, beta)) .* (1 + (N0 * cos(beta) + Q0 * sin(beta)) / (E * A0)) .* sin(beta);

            xC = integral(integrand_xC, beta_I_sol, beta_II_sol);
            yC = integral(integrand_yC, beta_I_sol, beta_II_sol);

            % Calculate xB and yB
            xB = xC - Rl * cos(pi/2 + beta_II_sol);
            yB = yC - Rl * sin(pi/2 - beta_II_sol);

            % Calculate sigma_eff and epsilon_eff
            sigma_eff_values(i) = N0 / A0;
            epsilon_eff_values(i) = abs(xB / (0.001) - 1);
        end
        
        % Plot current theta_l values
        plot(epsilon_eff_values*100, sigma_eff_values / E, 'o-', 'Color', colors(j, :), 'DisplayName', ['\theta_l = ', num2str(theta_l_values(j)), '^\circ']);
    end
    
    xlabel('Strain \epsilon_{eff} (%)');
    ylabel('Normalized stress \sigma_{eff}/E_S');
    title('Normalized Stress vs. Strain');
    legend('show');
    hold off;
end

function H = H_function(beta_I, beta, N0, Q0, E, I, A0, Rl)
    term1 = (N0 * (cos(beta_I) - cos(beta)) + Q0 * (sin(beta_I) - sin(beta))) / (E * I);
    term2 = 2 + (N0 * (cos(beta_I) + cos(beta)) + Q0 * (sin(beta_I) + sin(beta))) / (E * A0);
    H_beta_squared = term1 .* term2 + (1 / Rl)^2;
    H = sqrt(H_beta_squared);
end