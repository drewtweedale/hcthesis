function [Br, Btheta, Bphi] = OTD(r, B0, RN)
    % Convert position to column vector
    r = r(:);
    
    offset = [0.05; 0.49; 0.0] * RN;
    
    % Apply offset to position vector
    r_shifted = r - offset;
    
    % Compute magnitude of shifted position
    r_mag = norm(r_shifted);
    
    % Define dipole tilt (46.9° about y-axis, toward -x)
    tilt_angle = deg2rad(46.9);
    Ry = [cos(tilt_angle), 0, -sin(tilt_angle);
          0,               1,  0;
          sin(tilt_angle), 0,  cos(tilt_angle)];
    
    % Compute dipole moment vector (tilted)
    m_vec = Ry * [0; 0; 1];
    
    % Compute Cartesian coordinates for the original position
    r_orig_mag = norm(r);
    theta = acos(r(3) / r_orig_mag);
    phi = atan2(r(2), r(1));
    
    % Define spherical unit vectors at the observation point
    e_r = [sin(theta)*cos(phi); sin(theta)*sin(phi); cos(theta)];
    e_theta = [cos(theta)*cos(phi); cos(theta)*sin(phi); -sin(theta)];
    e_phi = [-sin(phi); cos(phi); 0];
    
    % Compute the field in Cartesian coordinates
    r_hat = r_shifted / r_mag;
    m_dot_r = dot(m_vec, r_hat);
    
    scaling_factor = (B0) * (RN/r_mag)^3;
    
    % Calculate field vector in Cartesian coordinates
    B_cart = scaling_factor * (3 * m_dot_r * r_hat - m_vec);
    
    % Project onto spherical coordinate system
    Br = dot(B_cart, e_r);
    Btheta = dot(B_cart, e_theta);
    Bphi = dot(B_cart, e_phi);
end