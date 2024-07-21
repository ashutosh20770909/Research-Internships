%please use constants by yourself 
%constants I have used are just guesses or some i found in text,I was not
%able to extract them well


function H = H_function(beta_I, beta, N0, Q0, E, I, A0, Rl)
    
    term1 = (N0 * (cos(beta_I) - cos(beta)) + Q0 * (sin(beta_I) - sin(beta))) / (E * I);
    
    term2 = 2 + (N0 * (cos(beta_I) + cos(beta)) + Q0 * (sin(beta_I) + sin(beta))) / (E * A0);
    
    H_beta_squared = term1 .* term2 + (1 / Rl)^2;
    
     H = sqrt(H_beta_squared);
end


% PLEASE VERIFY THESE AS THESE ARE MOSTLY GUESSED,
% TO CHECK PUT CORRECT VALUES BY YOURSELF 
% I WAS NOT ABLE TO EXTARCT THESE CONSTANTS FROM PAPER AS THEY WERE MOSTLY
% DIMENSIONLESS
theta_l = pi/2;
Rl = 0.0002;
N0 = 150000; %normal stress
Q0 = 50000; %shear stress
E = 1200000000;%young's modulus
I = 0.0002;  %area moment of inertia
A0 = 0.0002;


H = @(beta_I, beta) H_function(beta_I, beta, N0, Q0, E, I, A0, Rl);

% Derivative d_beta_dS based = H(beta)^2
d_beta_dS = @(beta_I, beta) H(beta_I, beta).^2;

% Equation (15)
integral_eq = @(beta_I, beta_II) integral(@(beta) sign(d_beta_dS(beta_I, beta)) ./ H(beta_I, beta), beta_I, beta_II) + theta_l * Rl;

% Initial guesses for beta_I and beta_II as fsolve function requires
% initial guesses to start and finds closest value so function becomes 0
%i have taken guess values from
beta_I_sol = pi/2; 
beta_II_sol = -pi/2; %from figure4(c) page no.7,type accordingly for correct results


options = optimset('Display', 'iter');
solution = fsolve(@(beta) integral_eq(beta(1), beta(2)), [beta_I_sol, beta_II_sol], options);
beta_I_sol = solution(1);
beta_II_sol = solution(2);


integrand_xC = @(beta) sign(d_beta_dS(beta_I_sol, beta)) .* (1 ./ H(beta_I_sol, beta)) .* (1 + (N0 * cos(beta) + Q0 * sin(beta)) / (E * A0)) .* cos(beta);
integrand_yC = @(beta) sign(d_beta_dS(beta_I_sol, beta)) .* (1 ./ H(beta_I_sol, beta)) .* (1 + (N0 * cos(beta) + Q0 * sin(beta)) / (E * A0)) .* sin(beta);


xC = integral(integrand_xC, beta_I_sol, beta_II_sol);
yC = integral(integrand_yC, beta_I_sol, beta_II_sol);


xB = xC - Rl * cos(pi/2 + beta_II_sol);
yB = yC - Rl * sin(pi/2 - beta_II_sol);


fprintf('beta_I: %f\n', beta_I_sol);
fprintf('beta_II: %f\n', beta_II_sol);
fprintf('xC: %f\n', xC);
fprintf('yC: %f\n', yC);
fprintf('xB: %f\n', xB);
fprintf('yB: %f\n', yB);