function B = dipole_field(r, RN, varargin)
    x = r(1); y = r(2); z = r(3);
    r_mag = sqrt(x^2 + y^2 + z^2);
    theta = acos(z/r_mag);
    phi = atan2(y,x);
    
    % Parse inputs and compute moments
    if nargin == 3
        B0 = varargin{1};
        Br = 2*B0*(RN^3/r_mag^3)*cos(theta);
        Btheta = B0*(RN^3/r_mag^3)*sin(theta);
        Bphi = 0;
    else
        % Schmidt coefficients case
        B0 = varargin{1};
        g10 = varargin{2};
        g11 = varargin{3};
        h11 = varargin{4};
        
        % Spherical components
        Br = 2*B0*(RN/r_mag)^3*(g10*cos(theta) - (g11*cos(phi) + h11*sin(phi))*sin(theta));
        Btheta = B0*(RN/r_mag)^3*(g10*sin(theta) + (g11*cos(phi) + h11*sin(phi))*cos(theta));
        Bphi = B0*(RN/r_mag)^3*(g11*sin(phi) - h11*cos(phi));
    end
    
    % Convert back to cartesian
    Bx = (Br*sin(theta) + Btheta*cos(theta))*cos(phi) + Bphi*sin(phi);  % Changed - to +
    By = (Br*sin(theta) + Btheta*cos(theta))*sin(phi) - Bphi*cos(phi);  % Changed + to -
    Bz = Br*cos(theta) - Btheta*sin(theta);
    
    B = [Bx, By, Bz];
end