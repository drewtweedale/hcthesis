function [Br, Btheta, Bphi] = quadrupole_field(r, RN, B0, varargin)
    x = r(1); y = r(2); z = r(3);
    r_mag = norm(r);
    theta = acos(z/r_mag);
    phi = atan2(y,x);

    if nargin == 3
        % Simple axial quadrupole case
        g20 = 1.0;  % Adjusted scaling
        g21 = 0;
        h21 = 0;
        g22 = 0;
        h22 = 0;
    else
        % Full Schmidt coefficients case
        g20 = varargin{1};
        g21 = varargin{2};
        h21 = varargin{3};
        g22 = varargin{4};
        h22 = varargin{5};
    end
    
    P20 = 0.5*(3*cos(theta)^2 - 1);
    P21 = -3*sin(theta)*cos(theta);
    P22 = 3*sin(theta)^2;
    
    dP20_dtheta = -3*sin(theta)*cos(theta);
    dP21_dtheta = -3*(cos(theta)^2 - sin(theta)^2);
    dP22_dtheta = 6*sin(theta)*cos(theta);
    
    % Corrected quadrupole terms with proper signs
    Br = 3*B0*(RN^4/r_mag^4)*(g20*P20 + (g21*cos(phi) + h21*sin(phi))*P21 + (g22*cos(2*phi) + h22*sin(2*phi))*P22);
    
    Btheta = -B0*(RN^4/r_mag^4)*(g20*dP20_dtheta + (g21*cos(phi) + h21*sin(phi))*dP21_dtheta + (g22*cos(2*phi) + h22*sin(2*phi))*dP22_dtheta);
    
    Bphi = -B0*(RN^4/r_mag^4)*((-g21*sin(phi) + h21*cos(phi))*P21/sin(theta) + 2*(-g22*sin(2*phi) + h22*cos(2*phi))*P22/sin(theta));
end