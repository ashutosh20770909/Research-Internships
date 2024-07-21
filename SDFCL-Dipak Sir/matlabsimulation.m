% Constants and initial conditions
mass = 50; % Mass of the satellite in kg
radius = 5; % Radius of the satellite in meters
I = (2/5) * mass * radius^2 * eye(3); % Moments of inertia for a sphere

omega0 = 0.0011667; % Orbital angular velocity magnitude (rad/s)
q = [0; 0; 0; 1]; % Initial quaternion (no rotation)
omega = [0; 0; 0]; % Initial angular velocity (rad/s)

% Inclination and true anomaly
i = 40 * pi / 180; % Inclination in radians

% Magnetic field base strength
B0 = 7.812e15;

% Charge density (example values)
sigma_values = [1, 0.1, 0.01, 0.001, 0.0001]; % Different surface charge densities

% Orbital period (approximate)
T = 5400; % Orbital period in seconds (90 minutes)

% Number of orbits
num_orbits = 5; % Number of orbits to simulate

% Time span
tspan = [0, num_orbits * T]; % Simulation time for the specified number of orbits

% ODE solver options
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);

% Function to calculate A(q)
function A_q = calc_A(q)
    q1 = q(1); q2 = q(2); q3 = q(3); q4 = q(4);
    A_q = [q4^2 + q1^2 - q2^2 - q3^2, 2*(q1*q2 + q4*q3), 2*(q1*q3 - q4*q2);
           2*(q1*q2 - q4*q3), q4^2 - q1^2 + q2^2 - q3^2, 2*(q2*q3 + q4*q1);
           2*(q1*q3 + q4*q2), 2*(q2*q3 - q4*q1), q4^2 - q1^2 - q2^2 + q3^2];
end

% Function to calculate orbital position vector r
function r = calc_orbital_position(t, omega0, i)
    % Semi-major axis (assuming circular orbit)
    a = 6500000; % Semi-major axis in meters (adjust as needed)
    
    % Orbital angular position (true anomaly)
    theta = omega0 * t; % True anomaly

    % Calculate position in orbital plane
    r_orbital = a * [cos(theta); sin(theta); 0];

    % Rotation matrix to account for inclination
    R_i = [1, 0, 0;
           0, cos(i), sin(i);
           0, -sin(i), cos(i)];

    % Position vector in inertial frame
    r = R_i * r_orbital;
end

% Define the ODE function
function dydt = satellite_dynamics(t, y, I, omega0, B0, sigma, i)
    % Extract quaternion and angular velocity from state vector
    q = y(1:4);
    omega = y(5:7);
    
    % Convert quaternion to rotation matrix
    R_q = calc_A(q);
    
    % Earth's rotation rate in rad/s
    omega_E = 7.2921e-5; 
    tau = omega0 * t; % Orbital angular position (radians)
    
    % Calculate relative velocity
    V_c = [6500000 * omega_E * sin(i) * cos(omega0 * t); 65000000 * (omega0 - omega_E * cos(i)); 0];
    
    % Lorentz torque calculation
    q_charge = sigma * 1; % Assuming area A = 1 for simplicity
    rho_0 = [2; 2; 2]; % Assumed radius vector (example values)
    
    % Calculate position vector r in the orbital frame
    r = calc_orbital_position(t, omega0, i); % Calculate changing r vector
    r_norm = norm(r);
    
    % Calculate the magnetic field B_o
    N_hat = [0; 0; -1]; % Example dipole direction, adjust as needed
    B_o = B0 / r_norm^5 * (3 * (dot(N_hat, r) * r) - r_norm^2 * N_hat); % Magnetic field vector
    
    % Compute V_c cross B_o
    Vc_cross_Bo = cross(V_c, B_o);
    
    % Transform V_c cross B_o using rotation matrix
    A_cross_result = R_q * Vc_cross_Bo;
    
    % Calculate Coulomb torque
    T_coul = q_charge * cross(rho_0, A_cross_result);
    
    % Get the current e_zb (third column of the rotation matrix)
    e_zb = R_q(:, 3);
    
    % Euler's equation for rotational dynamics
    S_omega = [0, -omega(3), omega(2);
               omega(3), 0, -omega(1);
               -omega(2), omega(1), 0];
    T_gg = 3 * omega0^2 * S_omega * e_zb; % Gravity gradient torque
    
    % Calculate angular acceleration
    I_inv = inv(I);
    domega_dt = I_inv * (S_omega * I * omega + T_gg + T_coul);
    
    % Calculate quaternion rate
    dq_dt = 0.5 * [-q(2), -q(3), -q(4);
                   q(1), -q(4),  q(3);
                   q(4),  q(1), -q(2);
                  -q(3),  q(2),  q(1)] * omega;
    
    % Concatenate the derivatives
    dydt = [dq_dt; domega_dt];
end

% Initial state vector (quaternion and angular velocity)
y0 = [q; omega];

% Run simulations for different sigma values
for k = 1:length(sigma_values)
    sigma = sigma_values(k);
    [t, y] = ode15s(@(t, y) satellite_dynamics(t, y, I, omega0, B0, sigma, i), tspan, y0, options);
    
    % Extract quaternion and angular velocity results
    q = y(:, 1:4);
    omega = y(:, 5:7);
    
    % Plot the results for quaternions
    figure;
    subplot(2,1,1);
    plot(t, q);
    title(['Quaternions (sigma = ' num2str(sigma_values(k)) ')']);
    xlabel('Time (s)');
    ylabel('Quaternion components');
    legend('q1', 'q2', 'q3', 'q4');
    
    % Plot the results for angular velocity
    subplot(2,1,2);
    plot(t, omega);
    title(['Angular Velocity (sigma = ' num2str(sigma_values(k)) ')']);
    xlabel('Time (s)');
    ylabel('Angular velocity (rad/s)');
    legend('omega_x', 'omega_y', 'omega_z');
end
