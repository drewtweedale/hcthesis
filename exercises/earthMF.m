% Parameters
q_proton = 1.6e-19;       % (C)
m_proton = 1.67e-27;      % (kg)
B0 = 31000e-9;            % Magnetic field strength at equator (T)
RE = 6371e3;              % Earth radius (m)

% Initial conditions
r0 = [8 * RE, 0, 0];      % Starting position at 8RE
v0 = [9.5e6, 0, 0];       % (9,500 km/s)
% The fancy case:
v0 = [9.5e6, 9.5e6, 0];

% Time parameters
T_g = 2 * pi * m_proton / (q_proton * B0); % Gyration period
dt = T_g / 100;           % Time step (s)
t_end = 100000  * T_g;     % End time (simulate for 100000 gyrations)
t = 0:dt:t_end;          
n_steps = length(t);      % Number of time steps for iteration

% Arrays to store position and velocity
x = zeros(n_steps, 3);   
v = zeros(n_steps, 3);  

% Initial conditions
x(1, :) = r0;
v(1, :) = v0;

% Function to calculate dipole magnetic field
function B = dipole_field(r, B0, RE)
    x = r(1); y = r(2); z = r(3);
    r_mag = sqrt(x^2 + y^2 + z^2); 
    scale = (r_mag / RE)^3 * r_mag^2; % Scaling factor
    
    Bx = (-3 * B0 * x * z) / scale;
    By = (-3 * B0 * y * z) / scale;
    Bz = (B0 * (x^2 + y^2 - 2 * z^2)) / scale;
    
    B = [Bx, By, Bz];
end

% Leapfrog-Boris integration
for i = 1:n_steps-1
    % Half-step position update
    x_mid = x(i, :) + 0.5 * dt * v(i, :);
    
    % Calculate magnetic field at midpoint
    B = dipole_field(x_mid, B0, RE);
    
    % Boris rotation step
    t_b = (q_proton / m_proton) * 0.5 * dt * B;
    v_minus = v(i, :) + dt * 0.5 * q_proton * [0, 0, 0] / m_proton; % No electric field
    
    % Rotation due to magnetic field
    v_prime = v_minus + cross(v_minus, t_b); 
    v_plus = v_minus + 2 / (1 + norm(t_b)^2) * cross(v_prime, t_b);
    
    v(i+1, :) = v_plus + 0.5 * dt * q_proton * [0, 0, 0] / m_proton; % No electric field
    x(i+1, :) = x_mid + 0.5 * dt * v(i+1, :);
end

% Normalize positions in terms of Earth radii
x_RE = x / RE;

% Plot trajectory
figure;
plot3(x_RE(:, 1), x_RE(:, 2), x_RE(:, 3), 'b', 'LineWidth', 1.5);
hold on;
[x_earth, y_earth, z_earth] = sphere;
surf(x_earth, y_earth, z_earth, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', 'green');
xlabel('x (R_E)');
ylabel('y (R_E)');
zlabel('z (R_E)');
title('Trajectory of a Proton in Dipole Magnetic Field');
grid on;
view(3);
axis equal;

% Set axis limits
axis([-10 10 -10 10 -10 10]);

% Display results
fprintf('Gyration period: %.3e s\n', T_g);
fprintf('Simulation time: %.3e s\n', t_end);