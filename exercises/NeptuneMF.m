% Parameters
q_proton = 1.6e-19;       % (C)
m_proton = 1.67e-27;      % (kg)
B0 = 1.42e-5;             % Magnetic field strength at Neptune's equator (T)
RN = 24622e3;             % Neptune radius (m)
eta = 0;                  % Tilt angle 

% Initial conditions
r0 = [8 * RN, 0, 0];      % Starting position at 8RN
v0 = [9.5e6, 9.5e6, 0];   % Initial velocity (m/s)

% Time parameters
T_g = 2 * pi * m_proton / (q_proton * B0); % Gyration period
dt = T_g / 100;           % Time step (s)
t_end = 100000 * T_g;     % End time (simulate for 100000 gyrations)
t = 0:dt:t_end;          
n_steps = length(t);      % Number of time steps for iteration

% Arrays to store position and velocity
x = zeros(n_steps, 3);   
v = zeros(n_steps, 3);  

% Initial conditions
x(1, :) = r0;
v(1, :) = v0;

% Function to calculate dipole magnetic field
function B = dipole_field(r, B0, RN)
    x = r(1); y = r(2); z = r(3);
    r_mag = sqrt(x^2 + y^2 + z^2); 
    scale = (r_mag / RN)^3 * r_mag^2; % Scaling factor
    
    Bx = (-3 * B0 * x * z) / scale;
    By = (-3 * B0 * y * z) / scale;
    Bz = (B0 * (x^2 + y^2 - 2 * z^2)) / scale;
    
    B = [Bx, By, Bz];
end

% Function to calculate quadrupole magnetic field with eta
function B = quadrupole_field(r, RN, eta)
    x = r(1); y = r(2); z = r(3);
    Me = 3.12e-5; % Magnetic moment (scaled)
    Re = RN; % Neptune radius
    angle = eta * pi/2; % Tilt angle
    
    % Quadrupole field formulas with eta
    Bx = -(2.5 * x * (z^2 - y^2) * Me * Re^4) / (sqrt(x^2 + y^2 + z^2))^7 * sin(angle) + ...
          ((x^2 * y - 1.5 * y^3 + 3.5 * y * z^2) * Me * Re^4) / (sqrt(x^2 + y^2 + z^2))^7 * cos(angle);
    By = -((x^2 * y - 1.5 * y^3 + 3.5 * y * z^2) * Me * Re^4) / (sqrt(x^2 + y^2 + z^2))^7 * sin(angle) + ...
          (2.5 * x * (z^2 - y^2) * Me * Re^4) / (sqrt(x^2 + y^2 + z^2))^7 * cos(angle);
    Bz = -((-x^2 * z + 1.5 * z^3 - 3.5 * z * y^2) * Me * Re^4) / (sqrt(x^2 + y^2 + z^2))^7;
    
    B = [Bx, By, Bz];
end

% Function to calculate combined magnetic field (dipole + quadrupole)
function B = combined_field(r, B0, RN, eta)
    B_dipole = dipole_field(r, B0, RN);
    B_quadrupole = quadrupole_field(r, RN, eta);
    B = B_dipole + B_quadrupole; % Combine dipole and quadrupole fields
end

% Initialize particle hit counter
particle_hit_count = 0;

% Boris integration
for i = 1:n_steps-1
    % Half-step position update
    x_mid = x(i, :) + 0.5 * dt * v(i, :);
    
    % Calculate magnetic field at midpoint (combined dipole and quadrupole)
    B = combined_field(x_mid, B0, RN, eta);
    
    % Boris rotation step
    t_b = (q_proton / m_proton) * 0.5 * dt * B;
    v_minus = v(i, :) + dt * 0.5 * q_proton * [0, 0, 0] / m_proton; % No electric field
    
    % Rotation due to magnetic field
    v_prime = v_minus + cross(v_minus, t_b); 
    v_plus = v_minus + 2 / (1 + norm(t_b)^2) * cross(v_prime, t_b);
    v(i+1, :) = v_plus + 0.5 * dt * q_proton * [0, 0, 0] / m_proton; % No electric field
    x(i+1, :) = x_mid + 0.5 * dt * v(i+1, :);
    
    % Check if the particle is within 1 RN of Neptune's center
    if norm(x(i+1, :)) <= RN
        particle_hit_count = particle_hit_count + 1;
    end
end

% Normalize positions in terms of Neptune radii
x_RN = x / RN;

% Create a grid for magnetic field visualization
[x_grid, y_grid, z_grid] = meshgrid(-20:1:20, -20:1:20, -20:1:20); % Grid in RN units
x_grid = x_grid * RN;
y_grid = y_grid * RN;
z_grid = z_grid * RN;

% Compute magnetic field at each grid point
Bx_grid = zeros(size(x_grid));
By_grid = zeros(size(y_grid));
Bz_grid = zeros(size(z_grid));
for i = 1:size(x_grid, 1)
    for j = 1:size(x_grid, 2)
        for k = 1:size(x_grid, 3)
            r = [x_grid(i, j, k), y_grid(i, j, k), z_grid(i, j, k)];
            B = combined_field(r, B0, RN, eta); % Use combined field
            Bx_grid(i, j, k) = B(1);
            By_grid(i, j, k) = B(2);
            Bz_grid(i, j, k) = B(3);
        end
    end
end

% Plot trajectory
figure;
plot3(x_RN(:, 1), x_RN(:, 2), x_RN(:, 3), 'b', 'LineWidth', 2); % Particle trajectory
hold on;

% Plot Neptune
[x_neptune, y_neptune, z_neptune] = sphere;
surf(x_neptune, y_neptune, z_neptune, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', 'blue');

% Define multiple starting points for streamlines in both hemispheres
theta = linspace(0, 2*pi, 16); % 16 starting points around the equator
startx = 10 * RN * cos(theta); % Starting points in the x-y plane
starty = 10 * RN * sin(theta);
startz = RN * ones(size(theta)); 

% Plot magnetic field lines using streamline
for i = 1:length(theta)
    % Positive charges
    positiveq = streamline(x_grid / RN, y_grid / RN, z_grid / RN, Bx_grid, By_grid, Bz_grid, startx(i) / RN, starty(i) / RN, startz(i) / RN);
    set(positiveq, 'Color', 'r', 'LineWidth', 1);
    
    % Negative charges to get full magnetic field lines
    negativeq = streamline(x_grid / RN, y_grid / RN, z_grid / RN, -1 * Bx_grid, -1 * By_grid, -1 * Bz_grid, startx(i) / RN, starty(i) / RN, startz(i) / RN);
    set(negativeq, 'Color', 'r', 'LineWidth', 1);
end

% Add labels and grid
xlabel('x (R_N)');
ylabel('y (R_N)');
zlabel('z (R_N)');
title('Proton Trajectory and Combined Magnetic Field Lines Around Neptune');
grid on;
view(3);
axis equal;

% Set axis limits
axis([-15 15 -15 15 -15 15]);

% Display results
fprintf('Gyration period: %.3e s\n', T_g);
fprintf('Simulation time: %.3e s\n', t_end);
fprintf('Total particle collisions with Neptune: %d\n', particle_hit_count);