function [Br, Btheta, Bphi] = dip_plus_quad(r, theta, phi, RN)

    % Victor's implementation 
    r = norm(r);
    % Schmidt coefficients (Tesla)
    g10 =  0.09732 * 10^(-4); 
    g11 =  0.03220 * 10^(-4); 
    h11 = -0.09889 * 10^(-4);
    g20 =  0.07448 * 10^(-4); 
    g21 =  0.00664 * 10^(-4); 
    h21 =  0.11230 * 10^(-4);
    g22 =  0.04499 * 10^(-4); 
    h22 = -0.00070 * 10^(-4);

    % Get Schmidt-normalized polynomials
    [P, dP] = schmidt_assc_legendre(theta);
    
    % Initialize components
    Br = 0.0;
    Btheta = 0.0;
    Bphi = 0.0;

    %% Dipole terms (n=1)
    n = 1;
    a_r_n1 = (RN/r)^(n+2);
    
    % Br components
    Br = Br + (n+1)*a_r_n1 * (g10*P(2,1) + ...
               (g11*cos(phi) + h11*sin(phi))*P(2,2));
    
    % Btheta components
    if sin(theta) > 1e-10
        Btheta = Btheta + a_r_n1*sin(theta) * (g10*dP(2,1) + ...
                         (g11*cos(phi) + h11*sin(phi))*dP(2,2));
    end
    
    % Bphi components
    if sin(theta) > 1e-10
        Bphi = Bphi + a_r_n1/sin(theta) * (g11*sin(phi) - h11*cos(phi))*P(2,2);
    end

    %% Quadrupole terms (n=2)
    n = 2;
    a_r_n2 = (RN/r)^(n+2);
    
    % Br components
    Br = Br + (n+1)*a_r_n2 * (g20*P(3,1) + ...
               (g21*cos(phi) + h21*sin(phi))*P(3,2) + ...
               (g22*cos(2*phi) + h22*sin(2*phi))*P(3,3));
    
    % Btheta components
    if sin(theta) > 1e-10
        Btheta = Btheta + a_r_n2*sin(theta) * (g20*dP(3,1) + ...
                         (g21*cos(phi) + h21*sin(phi))*dP(3,2) + ...
                         (g22*cos(2*phi) + h22*sin(2*phi))*dP(3,3));
    end
    
    % Bphi components
    if sin(theta) > 1e-10
        Bphi = Bphi + a_r_n2/sin(theta) * (1*(g21*sin(phi) - h21*cos(phi))*P(3,2) + ...
                        2*(g22*sin(2*phi) - h22*cos(2*phi))*P(3,3));
    end
end