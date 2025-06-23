% Parameters
q_proton = 1.6e-19;       % (C)
m_proton = 1.67e-27;      % (kg)
B0 = 1.42e-5;            % (T)
RN = 24765e3;            % (m)
offset_z = .5 * RN;       % (m)

% Energy levels to test (in MeV)
energy_levels = [1,2,3,4,5]; % MeV
colors = lines(length(energy_levels)); % Different color for each energy level

% Radial distances to test
radial_distances = [1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8] * RN;

% Initialize results storage
results = struct('energy', num2cell(energy_levels), ...
                'radial_distances', [], ...
                'eq_pitch_angles', []);

% Main simulation loop over energy levels
for e_idx = 1:length(energy_levels)
    current_energy = energy_levels(e_idx) * 1e6 * q_proton;
    fprintf('\n=== Testing Energy = %d MeV ===\n', energy_levels(e_idx));
    
    current_radial_distances = [];
    current_eq_angles = [];
    
    for r_idx = 1:length(radial_distances)
        current_R = radial_distances(r_idx);
        fprintf('  Testing R = %.1f R_N...\n', current_R/RN);
        
        % Step 1: Coarse search (1° steps) to find first stable angle
        coarse_angle = 0;
        stable_found = false;
        while coarse_angle <= 90 && ~stable_found
            fprintf('    Coarse testing pitch angle %.1f°... ', coarse_angle);
            [lost, ~] = runsim(coarse_angle, current_R, q_proton, m_proton, B0, RN, offset_z, current_energy);
            
            if lost
                fprintf('Lost\n');
                coarse_angle = coarse_angle + 1;
            else
                fprintf('Stable\n');
                stable_found = true;
            end
        end
        
        if stable_found
            % Step 2: Fine search (0.1° steps downward) to find loss boundary
            fine_angle = coarse_angle - 0.1;
            while fine_angle >= 0
                fprintf('      Fine testing pitch angle %.2f°... ', fine_angle);
                [lost, ~] = runsim(fine_angle, current_R, q_proton, m_proton, B0, RN, offset_z, current_energy);
                
                if lost
                    fprintf('Lost\n');
                    % Step 3: Ultra-fine search (0.01° steps upward) to pinpoint boundary
                    ultra_fine_angle = fine_angle + 0.01;
                    while ultra_fine_angle <= coarse_angle
                        fprintf('        Ultra-fine testing pitch angle %.2f°... ', ultra_fine_angle);
                        [lost, ~] = runsim(ultra_fine_angle, current_R, q_proton, m_proton, B0, RN, offset_z, current_energy);
                        
                        if ~lost
                            fprintf('Stable\n');
                            % Found the boundary - return previous angle
                            loss_angle = ultra_fine_angle - 0.01;
                            fprintf('      Boundary found at %.2f°\n', loss_angle);
                            current_radial_distances = [current_radial_distances, current_R/RN];
                            current_eq_angles = [current_eq_angles, loss_angle];
                            break;
                        else
                            fprintf('Lost\n');
                            ultra_fine_angle = ultra_fine_angle + 0.01;
                        end
                    end
                    break;
                else
                    fprintf('Stable\n');
                    coarse_angle = fine_angle; % Update our stable reference
                    fine_angle = fine_angle - 0.1;
                end
            end
        end
    end
    
    results(e_idx).radial_distances = current_radial_distances;
    results(e_idx).eq_pitch_angles = current_eq_angles;
end

% Calculate loss cone angles
L_values = linspace(1, 10, 100); % L-shell values from 1 to 10
valid_L = [];
theta_loss = [];
for L = L_values
    sin2_theta_loss = (4*L^6 - 3*L^5)^(-1/2);
    if sin2_theta_loss <= 1 && sin2_theta_loss >= 0
        valid_L = [valid_L, L];
        theta_loss = [theta_loss, asind(sqrt(sin2_theta_loss))];
    end
end

% Plot results
figure(1); clf; hold on;

% Plot theoretical loss cone
plot(valid_L, theta_loss, '-', 'LineWidth', 2, 'Color', 'k', ...
     'DisplayName', 'Theoretical Loss Cone');

% Plot simulation results for each energy level
for e_idx = 1:length(energy_levels)
    if ~isempty(results(e_idx).eq_pitch_angles) && length(results(e_idx).eq_pitch_angles) >= 4
        % Sort the data by radial distance for proper fitting
        [sorted_R, sort_idx] = sort(results(e_idx).radial_distances);
        sorted_angles = results(e_idx).eq_pitch_angles(sort_idx);
        
        % Create fit (using cubic polynomial)
        f = fit(sorted_R', sorted_angles', 'poly4');
        
        % Plot the data points (without legend entry)
        plot(sorted_R, sorted_angles, 'o', ...
             'Color', colors(e_idx,:), 'MarkerSize', 3, 'MarkerFaceColor', colors(e_idx,:), ...
             'HandleVisibility', 'off');
        
        % Plot the fit curve using predicted values
        x_fit = linspace(min(sorted_R), max(sorted_R), 100);
        y_fit = f(x_fit);
        plot(x_fit, y_fit, '--', 'Color', colors(e_idx,:), 'LineWidth', 1, ...
             'DisplayName', sprintf('%d MeV', energy_levels(e_idx)));
    else
        % Create dummy plot for legend
        plot(NaN, NaN, '--', 'Color', colors(e_idx,:), 'LineWidth', 1, ...
             'DisplayName', sprintf('%d MeV', energy_levels(e_idx)));
    end
end

xlabel('Radial Distance (R_N)');
ylabel('Equatorial Pitch Angle (degrees)');
title('Loss Cone Width vs Radial Distance for Different Energy Levels');
legend('Location', 'best');
grid on;
axis([1 8 0 40]);
xticks(0:1:8);
yticks(0:5:40);

% Print results to standard output
fprintf('\n=== Simulation Results ===\n');
for e_idx = 1:length(energy_levels)
    fprintf('\nEnergy: %d MeV\n', energy_levels(e_idx));
    if ~isempty(results(e_idx).eq_pitch_angles)
        % Sort results by radial distance for clean output
        [sorted_R, sort_idx] = sort(results(e_idx).radial_distances);
        sorted_angles = results(e_idx).eq_pitch_angles(sort_idx);
        
        fprintf('R (R_N)\tAngle (deg)\n');
        fprintf('-----------------\n');
        for i = 1:length(sorted_R)
            fprintf('%.1f\t%.2f\n', sorted_R(i), sorted_angles(i));
        end
    else
        fprintf('No stable trapping found for any pitch angle\n');
    end
end

% Simulation function using Boris algorithm (unchanged from original)
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
    T_g = 2 * pi * m / (q * norm(basic_dipole(r0,B0, RN))); % Gyroperiod
    dt = T_g / 1000; % Initial time step (small fraction of gyroperiod)
    
    % Calculate drift period (t_drift)
    L = r0_mag / RN; % L-shell value
    drift_prefactor = (2 * pi * q * B0 * RN^2) / (3 * L * m * v_mag^2);
    drift_sin_factor = 0.35 + 0.15 * sin(pitch_angle);
    t_drift = drift_prefactor * (1 / drift_sin_factor);
    
    % Ensure simulation runs for at least one full drift period
    t_end = t_drift;
    n_steps = ceil(t_end / dt);
    maxn_steps = 4000000;
    
    % Initialize tracking
    x = r0;
    v = v0;

    E = [0,0,0];
    final_pitch_angle = NaN;
    
    for i = 1:n_steps-1  
        if i == maxn_steps-1
            break
        end
        %==========================%
        B = basic_dipole(x, B0, RN);
        t_vec = (q ./ m) .* B .* 0.5 .* dt;
        s = 2 .* t_vec ./ (1 + norm(t_vec)^2);
        v_minus = v + (q ./ m) .* E .* 0.5 .* dt;
        v_prime = v_minus + cross(v_minus,t_vec);
        v_plus = v_minus + cross(v_prime,s);
        v = v_plus + (q ./ m) .* E .* 0.5 .* dt;
        x = x + v .* dt;
        %==========================%
        
        % Collision check
        if norm(x - [0, 0, offset_z]) <= RN
            B_collision = basic_dipole(x, B0, RN);
            pitch_angle_collision = atan2d(norm(cross(v, B_collision)), dot(v, B_collision));
            particle_lost = true;
            final_pitch_angle = pitch_angle_collision;
            return;
        end
    end
    
    % If we get here, particle was stably trapped for one full drift
    % period, or it flew out of Neptune's orbit.
    particle_lost = false;
end