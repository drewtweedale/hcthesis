function [Br, Btheta, Bphi] = quadrupole_field(r, theta, phi, RN, g20, g21, h21, g22, h22)
    r_mag = norm(r);
    P20 = 0.5*(3*cos(theta).^2 - 1);
    P21 = -3*sin(theta).*cos(theta);
    P22 = 3*sin(theta).^2;
    
    dP20_dtheta = -3*sin(theta).*cos(theta);
    dP21_dtheta = -3*(cos(theta).^2 - sin(theta).^2);
    dP22_dtheta = 6*sin(theta).*cos(theta);
    
    % Corrected quadrupole terms with proper signs
    Br = 3*(RN^4/r_mag^4)*(g20*P20 + (g21*cos(phi) + h21*sin(phi)).*P21 + (g22*cos(2*phi) + h22*sin(2*phi)).*P22);
    
    Btheta = -(RN^4/r_mag^4)*(g20*dP20_dtheta + (g21*cos(phi) + h21*sin(phi)).*dP21_dtheta + (g22*cos(2*phi) + h22*sin(2*phi)).*dP22_dtheta);
    
    Bphi = -(RN^4/r_mag^4)*((-g21*sin(phi) + h21*cos(phi)).*P21/sin(theta) + 2*(-g22*sin(2*phi) + h22*cos(2*phi)).*P22/sin(theta));
end