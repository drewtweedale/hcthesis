% Parameters
q_proton = 1.6e-19;  
m_proton = 1.67e-27; 
q_ion = q_proton;    
m_ion = 4 * m_proton; 
Bz = 1;              

% Choose particle type
q = q_proton;        
m = m_proton;        

% Initial conditions
v_perp = 1e6;        
v_par = 0.5e6;       
v0 = [v_perp, 0, v_par]; 
x0 = [0, 0, 0];      

% Time parameters
dt = 1e-10;          
t_end = 2e-6;        
t = 0:dt:t_end;      
n_steps = length(t); 

% Magnetic field
B = [0, 0, Bz];      

% Arrays to store position and velocity
x = zeros(n_steps, 3); 
v = zeros(n_steps, 3); 

% Set initial conditions
x(1, :) = x0;
v(1, :) = v0;

% Leapfrog-Boris integration
for i = 1:n_steps-1
    % Half-step position update
    x_mid = x(i, :) + 0.5 * dt * v(i, :);
    
    % Electric field (zero in this case)
    E = [0, 0, 0];
    
    % Magnetic field force
    t_b = (q / m) * 0.5 * dt * B;
    v_minus = v(i, :) + dt * 0.5 * q * E / m; 
    
    % Rotation due to magnetic field
    v_prime = v_minus + cross(v_minus, t_b); 
    v_plus = v_prime / (1 + norm(t_b)^2);
    
    % Full velocity update
    v(i+1, :) = v_plus + 0.5 * dt * q * E / m;
    
    % Full position update
    x(i+1, :) = x_mid + 0.5 * dt * v(i+1, :);
end

% Calculate energy conservation
energy = 0.5 * m * sum(v.^2, 2); % Kinetic energy (J)

% Calculate gyroradius
gyro_radius = v_perp * m / (q * Bz);

% Plot trajectory
figure;
plot3(x(:, 1), x(:, 2), x(:, 3), 'b', 'LineWidth', 1.5);
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');
title('3D Helical Trajectory of a Charged Particle');
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
fprintf('Initial energy: %.3e J\n', energy(1));
fprintf('Final energy: %.3e J\n', energy(end));
