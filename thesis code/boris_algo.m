function boris_algo(B0, RN, q, m, g10, g11, h11, g20, g21, h21, g22, h22)
    % Initial conditions
    r0 = [8*RN, 0, 0]; 
    v0 = [9.5e6, 9.5e6, 9.5e6];
    T_g = 2*pi*m/(q*B0); 
    dt = T_g/1000; 
    t_end = 300000*T_g;
    t = 0:dt:t_end; 
    n_steps = length(t);
    
    % Arrays
    x = zeros(n_steps, 3); 
    v = zeros(n_steps, 3);
    E = [0,0,0];
    x(1,:) = r0; 
    v(1,:) = v0;
    hits = 0;
    
    % Integration
    for i = 1:n_steps-1
        x_mid = x(i, :) + 0.5 * dt * v(i, :);
        
        % Field calculation
        [Br, Btheta, Bphi] = combined_field(x_mid, B0, RN, g10, g11, h11, g20, g21, h21, g22, h22);
        B = sph2cart_field(Br, Btheta, Bphi, x_mid);

        % Boris push
        t_b = (q./m) .* 0.5 .* dt .* B;
        s = 2 .* t_b ./ (1 + norm(t_b)^2);
        v_minus = v(i,:) + (q ./ m) .* E .* 0.5 .* dt;
        v_prime = v_minus + cross(v_minus,t_b);
        v_plus = v_minus + cross(v_prime,s);
        v(i+1,:) = v_plus + (q ./ m) .* E .* 0.5 .* dt;
        x(i+1,:) = x(i,:) + v(i+1,:) .* dt;
        
        % Hit detection
        if norm(x(i+1,:)) <= RN 
            hits = hits + 1; 
        end
    end
    
    % Plot
    figure(4); 
    hold on;
    [x_neptune, y_neptune, z_neptune] = sphere;
    surf(x_neptune, y_neptune, z_neptune, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', 'blue');
    plot3(x(:,1)/RN, x(:,2)/RN, x(:,3)/RN, 'r-', 'LineWidth', 1.5);
    plot3([0 0], [0 0], [-3 3], 'Color', [0 0.5 0], 'LineWidth', 1);
    xlabel('x (R_N)'); 
    ylabel('y (R_N)'); 
    zlabel('z (R_N)');
    title('Particle Trajectory'); 
    axis equal; 
    view(3); 
    grid on;
    fprintf('Particle collisions: %d\n', hits);
end