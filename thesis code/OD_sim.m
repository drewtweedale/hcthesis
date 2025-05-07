% Parameters
q_proton = 1.6e-19;       % (C)
m_proton = 1.67e-27;      % (kg)
B0 = 1.42e-5;            % (T)
RN = 24765e3;            % (m)
offset_z = 0.5 * RN;      % (m)
KE_total = 1e6 * q_proton; % (eV)

% Radial distances to test
radial_distances = [2, 4, 6, 8, 10] * RN;
colors = lines(length(radial_distances)); % Different color for each distance

% Initialize results storage
results = struct('radius', num2cell(radial_distances), ...
                'eq_pitch_angles', [], ...
                'final_loss_angles', []);

% Main simulation loop over radial distances
for r_idx = 1:length(radial_distances)
    current_R = radial_distances(r_idx);
    fprintf('\n=== Testing R = %.1f R_N ===\n', current_R/RN);
    
    % Initialize storage for this radial distance
    eq_angles = [];
    final_angles = [];
    first_collision_printed = false;
    
    % First pass: coarse 15° steps from 0° to 90°
    stable_angle_found = false;
    for pitch_angle_deg = 0:15:90
        fprintf('Testing pitch angle %d° at R=%.1f R_N... ', pitch_angle_deg, current_R/RN);
        
        % Run simulation with original Boris algorithm
        [lost, final_angle] = runsim(pitch_angle_deg, current_R, ...
                                           q_proton, m_proton, B0, RN, offset_z, KE_total);
        
        if lost
            fprintf('Lost at %.1f°\n', final_angle);
            eq_angles = [eq_angles, pitch_angle_deg];
            final_angles = [final_angles, final_angle];
        else
            fprintf('Stable\n');
            stable_angle_found = true;
            stable_coarse_angle = pitch_angle_deg;
            break;
        end
    end
    
    % If we found a stable angle, do fine scan downward
    if stable_angle_found
        fprintf('Starting fine scan from %d° down to 0°...\n', stable_coarse_angle);
        
        for fine_angle = (stable_coarse_angle-1):-1:0
            fprintf('  Testing pitch angle %d°... ', fine_angle);
            [lost, fine_final] = runsim(fine_angle, current_R, ...
                                             q_proton, m_proton, B0, RN, offset_z, KE_total);
            
            if lost
                fprintf('Lost at %.1f°\n', fine_final);
                
                % Print the first collision found during fine scan
                if ~first_collision_printed
                    fprintf('FIRST COLLISION AT: %.1f° (lost at %.1f°)\n', ...
                           fine_angle, fine_final);
                    first_collision_printed = true;
                end
                
                eq_angles = [eq_angles, fine_angle];
                final_angles = [final_angles, fine_final];
            else
                fprintf('Stable\n');
            end
        end
    end
    
    % Store results for this radial distance
    results(r_idx).eq_pitch_angles = eq_angles;
    results(r_idx).final_loss_angles = final_angles;
end

% Plot results
figure(1); clf; hold on;
for r_idx = 1:length(radial_distances)
    if ~isempty(results(r_idx).eq_pitch_angles)
        % Extract unique data points
        [unique_eq_angles, idx] = unique(results(r_idx).eq_pitch_angles);
        unique_final_angles = results(r_idx).final_loss_angles(idx);
        
        % Plot with lines connecting unique points
        plot(unique_eq_angles, unique_final_angles, '-o', ...
            'Color', colors(r_idx,:), 'MarkerSize', 5, 'MarkerFaceColor', colors(r_idx,:), ...
            'DisplayName', sprintf('%.1f R_N', radial_distances(r_idx)/RN));
    else
        % Create dummy plot for legend (with the same style as the actual plots)
        plot(NaN, NaN, '-o', 'Color', colors(r_idx,:), 'MarkerSize', 5, 'MarkerFaceColor', colors(r_idx,:), ...
             'DisplayName', sprintf('%.1f R_N', radial_distances(r_idx)/RN));
    end
end

xlabel('Equatorial Pitch Angle (degrees)');
ylabel('Final Loss Pitch Angle (degrees)');
title('Particle Loss vs Equatorial Pitch Angle for .5 RN Dipole Offset');
legend('Location', 'best');
grid on;
axis([0 90 0 90]);
xticks(0:15:90);
yticks(0:15:90);

% Simulation function using  Boris algorithm
function [particle_lost, final_pitch_angle] = ...
         runsim(pitch_angle_deg, r0_mag, q, m, B0, RN, offset_z, KE_total)
    
    % Convert to radians and calculate velocity
    pitch_angle = deg2rad(pitch_angle_deg);
    v_mag = sqrt(2 * KE_total / m);
    v_parallel = v_mag * cos(pitch_angle) * [0, 0, 1];
    v_perp = v_mag * sin(pitch_angle) * [1, 0, 0];
    v0 = v_parallel + v_perp;
    
    % Initial position
    r0 = [r0_mag, 0, 0];
    
    % Time parameters
    T_g = 2 * pi * m / (q * B0);
    dt = T_g / 50;
    t_end = 450000 * T_g;
    n_steps = ceil(t_end/dt);
    
    % Initialize tracking
    x_current = r0;
    v_current = v0;
    particle_lost = false;
    final_pitch_angle = NaN;
    escape_radius = 2 * r0_mag;
    
    for i = 1:n_steps
        t = i*dt;
        
        % Update time step every 50 steps based on local B-field
        if mod(i,50) == 0
            B_current = basic_dipole(x_current, B0, RN);
            Tgyro = 2*pi*m/(q*norm(B_current));
            dt = Tgyro/50;
        end
        
        % Boris algorithm
        x_mid = x_current + 0.5 * dt * v_current;
        B = basic_dipole(x_mid, B0, RN);
        
        t_b = (q / m) * 0.5 * dt * B;
        v_minus = v_current + 0.5 * dt * q * [0, 0, 0] / m; % No E-field
        v_prime = v_minus + cross(v_minus, t_b);
        v_plus = v_minus + 2 / (1 + norm(t_b)^2) * cross(v_prime, t_b);
        v_next = v_plus + 0.5 * dt * q * [0, 0, 0] / m;
        x_next = x_mid + 0.5 * dt * v_next;
        
        % Escape check
        if norm(x_next) > escape_radius
            particle_lost = true;
            final_pitch_angle = pitch_angle_deg;
            return;
        end
        
        % Collision check
        if norm(x_next - [0, 0, offset_z]) <= RN
            % Get magnetic field at collision point
            B_collision = basic_dipole(x_next, B0, RN);
            B_col_dot = dot(v_next, B_collision);
            v_next_sq = v_next(1)^2 + v_next(2)^2 + v_next(3)^2;
            B_par = (B_col_dot / v_next_sq) * v_next;
            B_perp = B_collision - B_par;
        
            % Calculate pitch angle at collision (in degrees)
            pitch_angle_collision = atan2d(norm(B_par), norm(B_perp));
            particle_lost = true;
            final_pitch_angle = pitch_angle_collision;
            return;
        end
        
        x_current = x_next;
        v_current = v_next;
    end
    
    % If we get here, particle was stably trapped
    particle_lost = false;
end