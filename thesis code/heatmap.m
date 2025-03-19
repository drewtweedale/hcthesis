function heatmap(B0, RN, eta)
    grid_step = 0.01;
    [y_grid, z_grid] = meshgrid(-5:grid_step:5, -5:grid_step:5);
    y_grid = y_grid * RN;
    z_grid = z_grid * RN;
    x_grid = zeros(size(y_grid));

    % Compute magnetic field magnitude at each grid point
    B_magnitude = zeros(size(y_grid));
    for i = 1:size(y_grid, 1)
        for j = 1:size(y_grid, 2)
            r = [x_grid(i, j), y_grid(i, j), z_grid(i, j)];
            B = combined_field(r, B0, RN, eta);
            B_magnitude(i, j) = norm(B); 
        end
    end

    y = 0.0001; 
    B_scaled = B_magnitude.^y;

    % Plot the heatmap
    figure;
    pcolor(y_grid / RN, z_grid / RN, B_scaled);
    shading interp;
    colorbar;
    xlabel('y (R_N)');
    ylabel('z (R_N)');
    title('Quadrupole Field Strength in the y-z Plane');
    axis equal;
    axis([-5 5 -5 5]);
end