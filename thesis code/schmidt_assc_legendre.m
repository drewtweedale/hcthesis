function [P, dP] = schmidt_assc_legendre(theta)
    % Function that outputs a matrix for Associated Legendre polynomials
    % and their derivatives.
    P = zeros(3,3);
    dP = zeros(3,3);

    % Associated Legendre polynomials
    P(1,1) = 1.0;
    P(2,1) = cos(theta);                    
    P(2,2) = -sin(theta);                 
    P(3,1) = 0.5*(3*cos(theta)^2 - 1);    
    P(3,2) = -3*cos(theta)*sin(theta); 
    P(3,3) = 3*sin(theta)^2;                
    
    % Derivatives (handling poles)
    if sin(theta) > 1e-10
        dP(2,1) = 1.0;
        dP(2,2) = cos(theta)/sin(theta);
        dP(3,1) = 3*cos(theta);
        dP(3,2) = 3*(2*cos(theta)^2 - 1)/sin(theta);
        dP(3,3) = -6*cos(theta);
    else
        dP(2,1) = 0;
        dP(2,2) = 0;
        dP(3,1) = 0;
        dP(3,2) = 0;
        dP(3,3) = 0;
    end
    
    % Schmidt normalization for m>0
    P(2,2) = P(2,2) * 1.0;
    dP(2,2) = dP(2,2) * 1.0;
    
    P(3,2) = P(3,2) / sqrt(3.0);
    dP(3,2) = dP(3,2) / sqrt(3.0);
    
    P(3,3) = P(3,3) / (2*sqrt(3.0));
    dP(3,3) = dP(3,3) / (2*sqrt(3.0));
end