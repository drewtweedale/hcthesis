% This file is merely to graphically display the solutions for
% assignment 1.

% Parameters
q = 1.6e-19;      % Charge of the proton (Coulombs)
m = 1.67e-27;     % Mass of the proton (kg)
Bz = 1;           % Magnetic field strength (Tesla)
v_x0 = 1e6;       % Initial velocity in x-direction (m/s) (Arbitrary)
v_z0 = 0.5e6;     % Initial velocity in z-direction (m/s) (Arbitrary)
omega = q * Bz / m; % Cyclotron frequency (rad/s)

% Time vector
t = linspace(0, 1e-6, 1000); % Time from 0 to 1 microsecond

% Equations of motion
v_x = v_x0 * cos(omega * t);
v_y = -v_x0 * sin(omega * t);
v_z = v_z0 * ones(size(t));

x = (v_x0 / omega) * sin(omega * t);
y = -(v_x0 / omega) * cos(omega * t) + (v_x0 / omega);
z = v_z0 * t;

% Plot the trajectory in 3D
figure;
plot3(x, y, z, 'b', 'LineWidth', 1.5);
hold on;
grid on;
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');
title('Gyro Motion of a Charged Particle in a Magnetic Field');
legend('Particle Trajectory');
view(3);

% Plot projections onto coordinate planes
figure;

% x-y plane
subplot(2, 2, 1);
plot(x, y, 'r', 'LineWidth', 1.5);
xlabel('x (m)');
ylabel('y (m)');
title('Projection onto x-y Plane');
grid on;

% x-z plane
subplot(2, 2, 2);
plot(x, z, 'g', 'LineWidth', 1.5);
xlabel('x (m)');
ylabel('z (m)');
title('Projection onto x-z Plane');
grid on;

% y-z plane
subplot(2, 2, 3);
plot(y, z, 'b', 'LineWidth', 1.5);
xlabel('y (m)');
ylabel('z (m)');
title('Projection onto y-z Plane');
grid on;

% 3D trajectory for reference
subplot(2, 2, 4);
plot3(x, y, z, 'k', 'LineWidth', 1.5);
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');
title('3D Trajectory');
grid on;
view(3);
