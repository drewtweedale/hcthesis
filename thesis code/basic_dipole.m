function [Br, Btheta, Bphi] = OTD(r, B0, RN)
    x = r(1); 
    y = r(2);
    z = r(3);
    
    r_mag = norm([x, y, z]);  % Magnitude in meters
    theta = acos(z/r_mag);    % Polar angle (from +z axis)
    phi = atan2(y,x);         % Azimuthal angle
    
    % Dipole field components (standard equations)
    Br = 2*B0*(RN^3/r_mag^3)*cos(theta);
    Btheta = B0*(RN^3/r_mag^3)*sin(theta);
    Bphi = 0;
end