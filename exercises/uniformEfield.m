% Parameters
q_proton = 1.6e-19;  % Charge of a proton (C)
m_proton = 1.67e-27; % Mass of a proton (kg)
Bz = 31000e-9;       % Magnetic field in z-direction (T)
Ey = 0.31;           % Electric field in y-direction (V/m)

% Initial conditions
v_perp = 1e5;        % Perpendicular velocity (m/s)
v_par = 0;           % Parallel velocity (m/s)
v0 = [v_perp, 0, v_par]; % Initial velocity vector
x0 = [0, 0, 0];      % Initial position vector

% Time parameters
T_g = 2 * pi * m_proton / (q_proton * Bz); % Gyration period
dt = T_g / 100;      % Time step (s)
t_end = 10 * T_g;    % Gyrate for 10 periods.
t = 0:dt:t_end;      
n_steps = length(t); % Number of time steps

% Magnetic and electric fields
B = [0, 0, Bz];      % Magnetic field vector
E = [0, Ey, 0];      % Electric field vector

% Arrays to store position and velocity
x = zeros(n_steps, 3); % Position array
v = zeros(n_steps, 3); % Velocity array

% Set initial conditions
x(1, :) = x0;
v(1, :) = v0;

% Leapfrog-Boris integration
for i = 1:n_steps-1
    x_mid = x(i, :) + 0.5 * dt * v(i, :);
    
    % Boris rotation step
    t_b = (q_proton / m_proton) * 0.5 * dt * B;
    v_minus = v(i, :) + dt * 0.5 * q_proton * E / m_proton; 
    
    % Rotation due to magnetic field
    v_prime = v_minus + cross(v_minus, t_b); 
    v_plus = v_minus + 2 / (1 + norm(t_b)^2) * cross(v_prime, t_b);
    
    v(i+1, :) = v_plus + 0.5 * dt * q_proton * E / m_proton;
    x(i+1, :) = x_mid + 0.5 * dt * v(i+1, :);
end

% Energy Conservation
energy = 0.5 * m_proton * sum(v.^2, 2); % Kinetic energy (J)

gyro_radius = v_perp * m_proton / (q_proton * Bz);

v_drift = Ey / Bz;

% Plot trajectory
figure;
plot3(x(:, 1), x(:, 2), x(:, 3), 'b', 'LineWidth', 1.5);
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');
title('3D Trajectory of a Proton in E \times B Field');
grid on;
view(3);

% Plot energy conservation
figure;
plot(t, energy, 'r', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Kinetic Energy (J)');
title('Energy Conservation');
grid on;

% Display results
fprintf('Gyro-radius: %.3e m\n', gyro_radius);
fprintf('E × B drift velocity: %.3e m/s\n', v_drift);
fprintf('Initial energy: %.3e J\n', energy(1));
fprintf('Final energy: %.3e J\n', energy(end));