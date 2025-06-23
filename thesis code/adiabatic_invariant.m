% Parameters (shared with OD_sim.m)
q_proton = 1.6e-19;       % (C)
m_proton = 1.67e-27;      % (kg)
B0 = 1.42e-5;            % (T)
RN = 24765e3;            % (m)
offset_z = 0 * RN;      % (m)
KE_total = .1 * 1e6 * q_proton; % (eV)

% Test case: Choose a radial distance and pitch angle to simulate
test_R = 6 * RN;          % Example: 4 R_N
test_pitch_angle = 35;    % Example: 30 degrees

% Run the simulation and store intermediate states
[particle_lost, final_angle, time_steps, mu_values, energy_ratios] = runsim_with_checks(...
    test_pitch_angle, test_R, q_proton, m_proton, B0, RN, offset_z, KE_total);

% Plot the first adiabatic invariant over time
figure(1); clf;
subplot(2,1,1);
plot(time_steps, mu_values, 'b', 'LineWidth', 2);
xlabel('Time Step');
ylabel('Magnetic Moment \mu (J/T)');
title('First Adiabatic Invariant (\mu) Over Time');
grid on;

% Plot energy conservation (ratio of current energy to initial energy)
subplot(2,1,2);
plot(time_steps, energy_ratios, 'r', 'LineWidth', 2);
xlabel('Time Step');
ylabel('Energy Ratio (E/E_0)');
title('Energy Conservation Check');
yline(1, '--k', 'LineWidth', 1.5, 'DisplayName', 'Ideal Conservation');
grid on;
legend('Location', 'best');

% Calculate relative changes
mu_initial = mu_values(1);
mu_final = mu_values(end);
mu_relative_change = abs(mu_final - mu_initial) / mu_initial;
fprintf('Relative change in magnetic moment: %.6f\n', mu_relative_change);

energy_relative_change = abs(energy_ratios(end) - 1);
fprintf('Relative energy deviation from initial: %.6f\n', energy_relative_change);

function [particle_lost, final_pitch_angle, time_steps, mu_values, energy_ratios] = ...
         runsim_with_checks(pitch_angle_deg, r0_mag, q, m, B0, RN, offset_z, KE_total)
    
    % Convert to radians and calculate velocity
    pitch_angle = deg2rad(pitch_angle_deg);
    v_mag = sqrt(2 * KE_total / m);
    v_parallel = v_mag * cos(pitch_angle) * [0, 0, 1];
    v_perp = v_mag * sin(pitch_angle) * [1, 0, 0];
    v0 = v_parallel + v_perp;
    
    % Initial position
    r0 = [r0_mag, 0, 0];
    
    % Time parameters
    B_init = basic_dipole(r0, B0, RN);
    T_g = 2 * pi * m / (q * norm(B_init));
    dt = T_g / 100;  % Reduced steps for testing (increase if needed)
    L = norm(r0) / RN;
    drift_prefactor = (2 * pi * q * B0 * RN^2) / (3 * L * m * v_mag^2);
    drift_sin_factor = 0.35 + 0.15 * sin(pitch_angle);
    t_end = drift_prefactor * (1 / drift_sin_factor);
    n_steps = ceil(t_end / dt);
    
    % Initialize tracking
    x = r0;
    v = v0;
    particle_lost = false;
    final_pitch_angle = NaN;
    escape_radius = 2 * r0_mag;
    E = [0, 0, 0];
    
    % Preallocate arrays
    time_steps = zeros(n_steps, 1);
    mu_values = zeros(n_steps, 1);
    energy_ratios = zeros(n_steps, 1);
    initial_energy = 0.5 * m * norm(v0)^2;
    
    for i = 1:n_steps
        t = i * dt;
        time_steps(i) = i;
        
        % Boris algorithm (FIXED)
        B = basic_dipole(x, B0, RN);
        B_mag = norm(B);
        
        % Half-step electric field acceleration
        v_minus = v + (q ./ m) .* E .* 0.5 .* dt;
        
        % Rotation
        t_vec = (q ./ m) .* B .* (0.5 .* dt);
        t_mag = norm(t_vec);
        s_vec = 2 .* t_vec ./ (1 + t_mag^2);
        
        v_prime = v_minus + cross(v_minus, t_vec);
        v_plus = v_minus + cross(v_prime, s_vec);
        
        % Final half-step acceleration
        v = v_plus + (q ./ m) .* E .* 0.5 .* dt;
        x = x + v .* dt;
        
        % Calculate μ and energy
        v_parallel_current = dot(v, B) / B_mag;
        v_perp_current = sqrt(norm(v)^2 - v_parallel_current^2);
        mu_values(i) = (0.5 * m * v_perp_current^2) / B_mag;
        energy_ratios(i) = (0.5 * m * norm(v)^2) / initial_energy;
        
        if norm(x - [0, 0, offset_z]) <= RN
            B_collision = basic_dipole(x, B0, RN);
            pitch_angle_collision = atan2d(norm(cross(v, B_collision)), dot(v, B_collision));
            particle_lost = true;
            final_pitch_angle = pitch_angle_collision;
            break;
        end
    end
    
    % Trim unused preallocated space
    time_steps = time_steps(1:i);
    mu_values = mu_values(1:i);
    energy_ratios = energy_ratios(1:i);
end