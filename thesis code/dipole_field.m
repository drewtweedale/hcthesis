function [Br, Btheta, Bphi] = dipole_field(r, RN, varargin)
    % Convert position to spherical coordinates
    x = r(1); y = r(2); z = r(3);
    r_mag = norm(r);
    theta = acos(z/r_mag);
    phi = atan2(y,x);
    
    if nargin == 3
        % Simple axial dipole case
        B0 = varargin{1};
        
        % Standard dipole field equations (Jackson Classical Electrodynamics)
        Br = 2*B0*(RN^3/r_mag^3)*cos(theta);
        Btheta = B0*(RN^3/r_mag^3)*sin(theta);
        Bphi = 0;
    else
        % Tilted dipole using Schmidt-normalized coefficients
        B0 = varargin{1};
        g10 = varargin{2};
        g11 = varargin{3};
        h11 = varargin{4};
        
        % Proper spherical harmonic expansion for dipole
        Br = B0*(RN^3/r_mag^3)*(2*g10*cos(theta) + 2*(g11*cos(phi) + h11*sin(phi))*sin(theta));
        Btheta = B0*(RN^3/r_mag^3)*(g10*sin(theta) - (g11*cos(phi) + h11*sin(phi))*cos(theta));
        Bphi = B0*(RN^3/r_mag^3)*(g11*sin(phi) - h11*cos(phi));
    end
end