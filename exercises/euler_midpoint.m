% Exact solution of dy/dt = -y
exact_solution = @(t) exp(-t);

% Initial condition and parameters
t0 = 0;         % Initial time
y0 = 1;         % Initial condition
time_steps = [0.1, 0.25, 0.5]; % Different step sizes
t_end = 5;      % Final time

% Initialize results storage
errors_euler = {};  
errors_midpoint = {}; 

% Loop over different time step sizes
for idx = 1:length(time_steps)
    h = time_steps(idx); % Current time step
    t = t0:h:t_end;      % Time vector
    n_steps = length(t); % Number of steps
    
    % Forward Euler method
    y_euler = zeros(1, n_steps);
    y_euler(1) = y0; % Initial condition
    for i = 1:n_steps-1
        y_euler(i+1) = y_euler(i) + h * (-y_euler(i));
    end
    
    % Midpoint method
    y_midpoint = zeros(1, n_steps);
    y_midpoint(1) = y0; % Initial condition
    for i = 1:n_steps-1
        k1 = -y_midpoint(i);
        k2 = -(y_midpoint(i) + h/2 * k1);
        y_midpoint(i+1) = y_midpoint(i) + h * k2;
    end
    
    % Exact solution
    y_exact = exact_solution(t);
    
    % Calculate errors
    euler_error = abs(y_euler - y_exact);
    midpoint_error = abs(y_midpoint - y_exact);
    errors_euler{idx} = euler_error; % Store in cell array
    errors_midpoint{idx} = midpoint_error; % Store in cell array
    
    % Plot the solutions
    figure;
    plot(t, y_exact, 'k', 'LineWidth', 1.5); hold on;
    plot(t, y_euler, 'r--', 'LineWidth', 1.2);
    plot(t, y_midpoint, 'b-.', 'LineWidth', 1.2);
    title(['Solutions with h = ', num2str(h)]);
    xlabel('Time (t)');
    ylabel('y(t)');
    legend('Exact Solution', 'Euler Method', 'Midpoint Method');
    grid on;
end

% Plot error analysis
figure;
for idx = 1:length(time_steps)
    h = time_steps(idx);
    t = t0:time_steps(idx):t_end;
    subplot(length(time_steps), 1, idx);
    plot(t, errors_euler{idx}, 'r--', 'LineWidth', 1.2); hold on;
    plot(t, errors_midpoint{idx}, 'b-.', 'LineWidth', 1.2);
    title(['Error Analysis for h = ', num2str(h)]);
    xlabel('Time (t)');
    ylabel('Error');
    legend('Euler Error', 'Midpoint Error');
    grid on;
end
