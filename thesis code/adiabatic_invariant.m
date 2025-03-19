% Load data from NeptuneMF.m
run('NeptuneMF.m'); % This runs the NeptuneMF script to get x, v, and other variables

% Initialize arrays to store the first adiabatic invariant
mu = zeros(n_steps, 1);

% Calculate the first adiabatic invariant at each time step
for i = 1:n_steps
    % Position and velocity at current time step
    r = x(i, :);
    v_current = v(i, :);
    
    % Calculate magnetic field at current position
    B = combined_field(r, B0, RN, eta);
    B_mag = norm(B); % Magnetic field strength
    
    % Calculate perpendicular velocity component
    v_parallel = dot(v_current, B) / B_mag; 
    v_perp = sqrt(norm(v_current)^2 - v_parallel^2); 
    
    % Calculate magnetic moment (first adiabatic invariant)
    mu(i) = (0.5 * m_proton * v_perp^2) / B_mag;
end

% Plot the first adiabatic invariant over time
figure;
plot(t / T_g, mu, 'b', 'LineWidth', 2); % Time in units of gyration period
xlabel('Time (Gyration Periods)');
ylabel('Magnetic Moment \mu (J/T)');
title('First Adiabatic Invariant (\mu) Over Time');
grid on;

% Zoom out the y-axis scale
mu_mean = mean(mu); 
mu_range = max(mu) - min(mu); 
y_min = mu_mean - 10 * mu_range;
y_max = mu_mean + 10 * mu_range;
ylim([y_min, y_max]); % Set y-axis limits

% Display the relative change in mu
mu_initial = mu(1);
mu_final = mu(end);
relative_change = abs(mu_final - mu_initial) / mu_initial;
fprintf('Relative change in magnetic moment: %.6f\n', relative_change);