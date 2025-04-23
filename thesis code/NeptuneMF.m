% Load Voyager 2 data
data = load('../magnetic_field_data/COMPREHE.ASC');
INVALID_VALUE = 9999.99;
invalid_rows = any(abs(data(:,10:12)) >= INVALID_VALUE, 2);
clean_data = data(~invalid_rows,:);

% Parameters
q_proton = 1.6e-19;       % (C)
m_proton = 1.67e-27;      % (kg)
B0 = 1.42e-5;             % Magnetic field strength at Neptune's equator (T)
RN = 24765e3;             % Neptune radius (m) from PDS file
g10 =  0.09732;
g11 =  0.03220;
h11 = -0.09889;
g20 =  0.07448;
g21 =  0.00664; 
h21 =  0.11230;
g22 =  0.04499;
h22 = -0.00070;

% Neptune pole orientation (from JPL physical constants)
ra = deg2rad(298.90);     % Right ascension
dec = deg2rad(42.84);     % Declination

% Rotation matrix (RyRz)
R = [cos(ra)*cos(dec), -sin(ra), -cos(ra)*sin(dec);
     sin(ra)*cos(dec),  cos(ra), -sin(ra)*sin(dec);
     sin(dec),          0,        cos(dec)];

% Helper function to convert Cartesian B to spherical
function [Br, Btheta, Bphi] = cart2sph_field(B, theta_deg, phi_deg)
    theta = deg2rad(theta_deg);
    phi = deg2rad(phi_deg);
    Br = B(1)*sin(theta)*cos(phi) + B(2)*sin(theta)*sin(phi) + B(3)*cos(theta);
    Btheta = B(1)*cos(theta)*cos(phi) + B(2)*cos(theta)*sin(phi) - B(3)*sin(theta);
    Bphi = -B(1)*sin(phi) + B(2)*cos(phi);
end
    
% Extract cleaned columns
year = 1900 + clean_data(:,1);  % Convert from YY format
doy = clean_data(:,2);
hour = clean_data(:,3);
minute = clean_data(:,4);
second = clean_data(:,5);
millisecond = clean_data(:,6);
range_RN = clean_data(:,7);
latitude = clean_data(:,8);
longitude = clean_data(:,9);   
Br_meas = clean_data(:,10);
Btheta_meas = clean_data(:,11); % Southward component
Bphi_meas = clean_data(:,12);   % Azimuthal component

% Create datetime array
time_datetime = datetime(year, 1, doy) + hours(hour) + minutes(minute) + seconds(second) + milliseconds(millisecond);

% Apply longitude system correction (lambda_NLS = 360 - phi)
phi = mod(360 - longitude, 360); % Convert to standard spherical phi

% Convert to Cartesian coordinates with proper orientation
theta = 90 - latitude; % Convert to co-latitude
[x,y,z] = sph2cart(deg2rad(phi), deg2rad(theta), range_RN*RN);
positions = [x,y,z] * R'; 

% Calculate measured field magnitude
B_meas = sqrt(Br_meas.^2 + Btheta_meas.^2 + Bphi_meas.^2);

% Calculate model predictions at Voyager 2 positions
B_model = zeros(length(range_RN),3);
B_model_mag = zeros(length(range_RN),1);
    
for i = 1:length(range_RN)
    r = positions(i,:); % Use rotated position vector
    B = combined_field(r, B0, RN, g10, g11, h11, g20, g21, h21, g22, h22);
    [Br_mod, Btheta_mod, Bphi_mod] = cart2sph_field(B, theta(i), phi(i));
    B_model(i,:) = [Br_mod, Btheta_mod, Bphi_mod];
    B_model_mag(i) = norm(B);
end

% =========================================================================
% PLOTTING FIGURES
% =========================================================================

% Figure 1: Voyager 2 Trajectory with Rotation Axis
figure(1);
plot3(positions(:,1)/RN, positions(:,2)/RN, positions(:,3)/RN, 'b-', 'LineWidth', 2);
hold on;
[x_neptune, y_neptune, z_neptune] = sphere(50);
surf(x_neptune, y_neptune, z_neptune, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', 'blue');

% Add rotation axis (dark green line)
plot3([0 0], [0 0], [-3 3], 'Color', [0 0.5 0], 'LineWidth', 1);

axis([-12 12 -12 12 -12 12]);
xlabel('x (R_N)'); ylabel('y (R_N)'); zlabel('z (R_N)');
title('Voyager 2 Trajectory in Neptune-Fixed Coordinates');
grid on; view(45,30); axis equal;

% Time formatting for plots
time_hours = hours(time_datetime - time_datetime(1)); % Hours since start
start_time = floor(min(time_hours));
end_time = ceil(max(time_hours));

% Figure 2: Field Magnitude Comparison
figure(2);
plot(time_hours, B_meas, 'b-', 'LineWidth', 2);
hold on;
plot(time_hours, B_model_mag*1e9, 'r--', 'LineWidth', 2);

% Mark closest approach
[~, idx_min] = min(range_RN);
ca_time = time_hours(idx_min);
xline(ca_time, '--g', sprintf('CA: %.1f hours', ca_time), 'LabelVerticalAlignment', 'top');

xlabel('Hours Since 1989-08-24T18:00:00');
ylabel('|B| (nT)');
title('Magnetic Field Strength Comparison');
legend('Voyager 2', 'Model', 'Location', 'best');
grid on;
xlim([start_time end_time]);

% Figure 3: Component-wise Comparison
figure(3);
subplot(3,1,1);
plot(time_hours, Br_meas, 'b-', 'LineWidth', 1.5);
hold on;
plot(time_hours, B_model(:,1)*1e9, 'r--', 'LineWidth', 1.5);
ylabel('B_r (nT)');
title('Radial Component');
grid on;
xlim([start_time end_time]);

subplot(3,1,2);
plot(time_hours, Btheta_meas, 'b-', 'LineWidth', 1.5);
hold on;
plot(time_hours, B_model(:,2)*1e9, 'r--', 'LineWidth', 1.5);
ylabel('B_{\theta} (nT)');
title('Southward Component');
grid on;
xlim([start_time end_time]);

subplot(3,1,3);
plot(time_hours, Bphi_meas, 'b-', 'LineWidth', 1.5);
hold on;
plot(time_hours, B_model(:,3)*1e9, 'r--', 'LineWidth', 1.5);
ylabel('B_{\phi} (nT)');
title('Azimuthal Component');
xlabel('Hours Since 1989-08-24T18:00:00');
grid on;
xlim([start_time end_time]);

% Mark closest approach on all subplots
for i = 1:3
    subplot(3,1,i);
    xline(ca_time, '--g', 'LineWidth', 1);
end

% ====================================================================================
% BORIS ALGORITHM FOR PARTICLE TRACKING
% ====================================================================================

% Initial conditions
r0 = [8 * RN, 0, 0];      % Starting position at 8RN
v0 = [9.5e6, 9.5e6, 9.5e6];   % Initial velocity (m/s)

% Time parameters
T_g = 2 * pi * m_proton / (q_proton * B0); % Gyration period
dt = T_g / 100;           % Time step (s)
t_end = 300000 * T_g;     % End time (simulate for 100000 gyrations)
t = 0:dt:t_end;          
n_steps = length(t);      % Number of time steps for iteration

% Arrays to store position and velocity
x = zeros(n_steps, 3);   
v = zeros(n_steps, 3);  

% Initial conditions
x(1, :) = r0;
v(1, :) = v0;

% Initialize particle hit counter
particle_hit_count = 0;

% Boris integration
for i = 1:n_steps-1
    % Half-step position update
    x_mid = x(i, :) + 0.5 * dt * v(i, :);
    
    % Calculate magnetic field at midpoint (combined dipole and quadrupole)
    B = combined_field(x_mid, B0, RN, g10, g11, h11, g20, g21, h21, g22, h22);
    
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

% ====================================================================================
% PLOTTING MAGNETIC FIELD LINES WITH ROTATION AXIS
% ====================================================================================

% Create a grid for magnetic field lines
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
            B = combined_field(r, B0, RN, g10, g11, h11, g20, g21, h21, g22, h22);
            Bx_grid(i, j, k) = B(1);
            By_grid(i, j, k) = B(2);
            Bz_grid(i, j, k) = B(3);
        end
    end
end

% Plot trajectory
figure(4);
hold on;

% Plot Neptune
[x_neptune, y_neptune, z_neptune] = sphere;
surf(x_neptune, y_neptune, z_neptune, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', 'blue');

% Add rotation axis (dark green line)
plot3([0 0], [0 0], [-3 3], 'Color', [0 0.5 0], 'LineWidth', 1);

% Define starting points for magnetic field lines
num_streamlines = 100;
theta = linspace(0, 2*pi, num_streamlines);
radius = 5 * RN;
startx = radius * cos(theta);
starty = radius * sin(theta);
startz = .5 * RN * ones(size(startx)); % Slightly above equator

% Plot magnetic field lines
streamline_color = 'r';

for i = 1:num_streamlines
    h = streamline(x_grid / RN, y_grid / RN, z_grid / RN, Bx_grid, By_grid, Bz_grid, startx(i) / RN, starty(i) / RN, startz(i) / RN);
    g = streamline(x_grid / RN, y_grid / RN, z_grid / RN, -1 * Bx_grid, -1 * By_grid, -1 * Bz_grid, startx(i) / RN, starty(i) / RN, startz(i) / RN);
    hnegz = streamline(x_grid / RN, y_grid / RN, z_grid / RN, Bx_grid, By_grid, Bz_grid, startx(i) / RN, starty(i) / RN, -1 * startz(i) / RN);
    gnegz = streamline(x_grid / RN, y_grid / RN, z_grid / RN, -1 * Bx_grid, -1 * By_grid, -1 * Bz_grid, startx(i) / RN, starty(i) / RN, -1 * startz(i) / RN);
    
    % Set properties for the streamlines
    set(h, 'Color', streamline_color, 'LineWidth', 1);
    set(g, 'Color', streamline_color, 'LineWidth', 1);
    set(hnegz, 'Color', streamline_color, 'LineWidth', 1);
    set(gnegz, 'Color', streamline_color, 'LineWidth', 1);
end

% Add labels and grid
xlabel('x (R_N)');
ylabel('y (R_N)');
zlabel('z (R_N)');
title('Neptune Dipole + Quadrupole Field with Rotation Axis');
grid on;
view(3);
axis equal;

% Set axis limits
axis([-15 15 -15 15 -15 15]);

% Display results
fprintf('Gyration period: %.3e s\n', T_g);
fprintf('Simulation time: %.3e s\n', t_end);
fprintf('Total particle collisions with Neptune: %d\n', particle_hit_count);