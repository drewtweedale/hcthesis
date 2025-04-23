function B = quadrupole_field(r, RN, B0, varargin)
    % Convert Cartesian to spherical coordinates
    x = r(1); y = r(2); z = r(3);
    r_mag = sqrt(x^2 + y^2 + z^2);
    theta = acos(z/r_mag);
    phi = atan2(y,x);

    % Parse inputs
    if nargin == 3
        % Simple axial quadrupole case
        g20 = 1 / (2 * RN^4); % Pulled from Vogt paper, but unsure about the scaling.
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
    
    % Associated Legendre Polynomials (l=2)
    P20 = 0.5 * (3*cos(theta)^2 - 1);   % m=0
    P21 = -3 * sin(theta) * cos(theta); % m=1 
    P22 = 3 * sin(theta)^2;             % m=2
    
    % Derivatives dP/dθ 
    dP20_dtheta = -3 * sin(theta) * cos(theta);
    dP21_dtheta = -3 * (cos(theta)^2 - sin(theta)^2);
    dP22_dtheta = 6 * sin(theta) * cos(theta);
    
    % Radial component 
    Br = 3 * B0 * (RN/r_mag)^4 * ( ...
        g20 * P20 + ...
        (g21 * cos(phi) + h21 * sin(phi)) * P21 + ...
        (g22 * cos(2*phi) + h22 * sin(2*phi)) * P22 );
    
    % Polar component 
    Btheta = -1 * B0 * (RN/r_mag)^4 * ( ...
        g20 * dP20_dtheta + ...
        (g21 * cos(phi) + h21 * sin(phi)) * dP21_dtheta + ...
        (g22 * cos(2*phi) + h22 * sin(2*phi)) * dP22_dtheta );
    
    % Azimuthal component
    Bphi = -1 * B0 *(RN/r_mag)^4 * ( ...
        ( -g21 * sin(phi) + h21 * cos(phi)) * P21 / sin(theta) + ...
        2 * (-g22 * sin(2*phi) + h22 * cos(2*phi)) * P22 / sin(theta) );
    
    % Convert to Cartesian 
    Bx_sph = Br*sin(theta)*cos(phi) + Btheta*cos(theta)*cos(phi) + Bphi*sin(phi);
    By_sph = Br*sin(theta)*sin(phi) + Btheta*cos(theta)*sin(phi) - Bphi*cos(phi);
    Bz_sph = Br*cos(theta) - Btheta*sin(theta);
    
    B = [Bx_sph, By_sph, Bz_sph];
end