% Function to calculate combined magnetic field (dipole + quadrupole)
function [Br, Btheta, Bphi] = new_comb(r, theta, phi, RN)
        % Trying a new implementation of the dipole + quadrupole model
        % (courtesy of Victor).
        % Define spherical variables
        r_mag = norm(r);
        
        % Schmidt Coefficients
        max_mul = 2;
        sf = 10^(-4);
        g = zeros(3, 3);
        h = zeros(3, 3);
        
        g(2,1) =  0.09732 * sf; 
        g(2,2) =  0.03220 * sf; 
        h(2,2) = -0.09889 * sf;
        g(3,1) =  0.07448 * sf; 
        g(3,2) =  0.00664 * sf; 
        h(3,2) =  0.11230 * sf;
        g(3,3) =  0.04499 * sf; 
        h(3,3) = -0.00070 * sf;

        % Associated Legendre Polynomials and their derivatives with respect to cos(theta):
        P = zeros(3, 3);
        dP = zeros(3, 3);

        P(1,1) = 1;
        P(2,1) = cos(theta);
        P(2,2) = -sin(theta);
        P(3,1) = .5*(3*cos(theta)^2-1);
        P(3,2) = -3*cos(theta)*sin(theta);
        P(3,3) = 3*sin(theta)^2;

        dP(1,1) = 0;
        dP(2,1) = 1;
        dP(2,2) = cot(theta);
        dP(3,1) = 3*cos(theta);
        dP(3,2) = (-3*(2*cos(theta)^2 -1))/sin(theta);
        dP(3,3) = -6*cos(theta);

        % Schmidt Normalization Coefficients
        S = zeros(3, 3);

        S(1,1) = 1;
        S(2,1) = 1;
        S(2,2) = 1;
        S(3,1) = 1;
        S(3,2) = 1/sqrt(3);
        S(3,3) = 1/(2*sqrt(3));

        Br = 0;
        Btheta = 0;
        Bphi = 0;
        
        for n = 1:max_mul
            Br_temp = 0;
            Btheta_temp = 0;
            Bphi_temp = 0;
            for m = 0:n
                Br_temp = Br_temp + S(n+1,m+1) * P(n+1,m+1) * (g(n+1,m+1)*cos(m*phi) + h(n+1,m+1)*sin(m*phi));
                Btheta_temp = Btheta_temp + sin(theta) * S(n+1,m+1) * dP(n+1,m+1) * (g(n+1,m+1)*cos(m*phi) + h(n+1,m+1)*sin(m*phi));
                Bphi_temp = Bphi_temp + (1/sin(theta)) * S(n+1,m+1) * P(n+1,m+1) * m * (g(n+1,m+1)*sin(m*phi) - h(n+1,m+1)*cos(m*phi));
            end 
            Br = Br + (n+1)*(RN/r_mag)^(n+2)*Br_temp;
            Btheta = Btheta + (RN/r_mag)^(n+2)*Btheta_temp;
            Bphi = Bphi + (RN/r_mag)^(n+2)*Bphi_temp;
        end


        


end