% Load and clean Voyager 2 data
data = load('../magnetic_field_data/COMPREHE.ASC');
INVALID_VALUE = 9999.99;
invalid_rows = any(abs(data(:,10:12)) >= INVALID_VALUE, 2);
clean_data = data(~invalid_rows,:);

% Parameters
q_proton = 1.6e-19;       % (C)
m_proton = 1.67e-27;      % (kg)
B0 = 1.33e-5;             % Magnetic field strength at Neptune's equator (T)
RN = 24765e3;             % Neptune radius (m)
g10 =  0.09732 * 10^(-4); 
g11 =  0.03220 * 10^(-4); 
h11 = -0.09889 * 10^(-4);
g20 =  0.07448 * 10^(-4); 
g21 =  0.00664 * 10^(-4); 
h21 =  0.11230 * 10^(-4);
g22 =  0.04499 * 10^(-4); 
h22 = -0.00070 * 10^(-4);

% Extract cleaned data
year = 1900 + clean_data(:,1);  
doy = clean_data(:,2);
hour = clean_data(:,3); minute = clean_data(:,4); 
second = clean_data(:,5); millisecond = clean_data(:,6);
range_RN = clean_data(:,7); latitude = clean_data(:,8); 
longitude = clean_data(:,9); Br_meas = clean_data(:,10);
Btheta_meas = clean_data(:,11); Bphi_meas = clean_data(:,12);   

% Create datetime array
time_datetime = datetime(year, 1, doy) + hours(hour) + minutes(minute) + ...
                seconds(second) + milliseconds(millisecond);

% Longitude (W), converted to phi component.
phi = mod(360 - longitude, 360); 
phi_rad = deg2rad(phi);

% Latitude converted to theta (measured from Z axis)
theta = 90-latitude; 
theta_rad = deg2rad(theta);

[x,y,z] = sph2cartcoord(phi_rad, theta_rad, range_RN*RN);
positions = [x,y,z]; 

% Field comparison
B_meas = sqrt(Br_meas.^2 + Btheta_meas.^2 + Bphi_meas.^2);
B_model = zeros(length(range_RN),3);
for i = 1:length(range_RN)
    % Given position
    r = positions(i,:)';

    % Models: uncomment the desired model.

    % Dipole + quadrupole code derived from Victor's B-field derivations.
    % [Br_mod, Btheta_mod, Bphi_mod] = new_comb(r, theta_rad(i), phi_rad(i), RN);
    
    % Drew's dipole + quadrupole implementation.
    % [Br_mod, Btheta_mod, Bphi_mod] = combined_field(r, theta, phi, RN, g10, g11, h11, g20, g21, h21, g22, h22);
    
    % OTD model.
    [Br_mod, Btheta_mod, Bphi_mod] = dipole_field(r, theta_rad(i), phi_rad(i), RN, B0);
    B_model(i,:) = [Br_mod, Btheta_mod, Bphi_mod];
end
B_model_mag = sqrt(sum(B_model.^2,2));

% Plotting
time_hours = hours(time_datetime - time_datetime(1)); 
start_time = floor(min(time_hours));
end_time = ceil(max(time_hours));

% Figure 1: Voyager 2 Trajectory
figure(1);
plot3(positions(:,1)/RN, positions(:,2)/RN, positions(:,3)/RN, 'b-', 'LineWidth', 2);
hold on;
[x_neptune, y_neptune, z_neptune] = sphere(50);
surf(x_neptune, y_neptune, z_neptune, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', 'blue');
axis([-12 12 -12 12 -12 12]);
xlabel('x (R_N)'); ylabel('y (R_N)'); zlabel('z (R_N)');
title('Voyager 2 Trajectory in Neptune-Fixed Coordinates');
grid on; view(45,30); axis equal;

% Figure 2: Field Magnitude Comparison
% Figure 2: Field Magnitude Comparison (log scale)
figure(2);
semilogy(time_hours, B_meas, 'b-', 'LineWidth', 2);
hold on;
semilogy(time_hours, B_model_mag*1e9, 'r--', 'LineWidth', 2);

xlabel('Hours Since 1989-08-24T18:00:00');
ylabel('|B| (nT)');
title('Magnetic Field Strength Comparison');
legend('Voyager 2', 'Model', 'Location', 'best');
grid on;
xlim([start_time 20]);

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
xlabel('Hours Since 1989-08-24T18:00:00');
grid on;
xlim([start_time end_time]);

% Mark closest approach on all subplots
for i = 1:3
    subplot(3,1,i);
    xline(ca_time, '--g', 'LineWidth', 1);
end

% Field lines and Boris Algorithm
figure(4);
mfield_lines(B0, RN, g10, g11, h11, g20, g21, h21, g22, h22);
% boris_algo(B0, RN, q_proton, m_proton, g10, g11, h11, g20, g21, h21, g22, h22);