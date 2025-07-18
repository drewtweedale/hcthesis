function [Br, Btheta, Bphi] = dipole_field(r, theta, phi, RN, B0)
    % Convert position to spherical coordinates
    offset = [0.17, 0.46, -.24] * RN;
    x = r(1) - offset(1); y = r(2) - offset(2); z = r(3) - offset(3);
    B0 = -B0;
        
    % Rotation angles
    tilt_y = -deg2rad(46.8);     % Counterclockwise rotation about y-axis.
    tilt_z = -deg2rad(259.5);     % Clockwise rotation about z-axis

    % Rotation matrices (global frame)
    Ry = [cos(tilt_y), 0, sin(tilt_y);
         0, 1, 0;
         -sin(tilt_y), 0, cos(tilt_y)];

    Rz = [cos(tilt_z), -sin(tilt_z), 0;
          sin(tilt_z),  cos(tilt_z), 0;
          0, 0, 1];

    % Combine into a single rotation matrix (applied to positions)
    R = Rz * Ry;  

    % Apply inverse rotation to position (rotate r into dipole frame)
    r_rot = R' * [x; y; z];  

    % Now compute B-field in dipole-aligned frame
    x_r = r_rot(1);
    y_r = r_rot(2);
    z_r = r_rot(3);
    r_mag_r = norm(r_rot);
    scale = (r_mag_r / RN)^3 * r_mag_r^2;

    Bx = (-3 * B0 * x_r * z_r) / scale;
    By = (-3 * B0 * y_r * z_r) / scale;
    Bz = (B0 * (x_r^2 + y_r^2 - 2 * z_r^2)) / scale;

    % Rotate B back into original (rotated dipole) frame
    B_rot = R * [Bx; By; Bz];

    % Convert to spherical field components
    B = cart2sph_field(B_rot(1), B_rot(2), B_rot(3), theta, phi);
    Br = B(1);
    Btheta = B(2);
    Bphi = B(3);
end