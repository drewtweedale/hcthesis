% Function to calculate combined magnetic field (dipole + quadrupole)
function [Br, Btheta, Bphi] = combined_field(r, theta, phi, RN, g10, g11, h11, g20, g21, h21, g22, h22)
    % Function that combines the individual component fields
    [Br_dipole, Btheta_dipole, Bphi_dipole] = dipole_field_sch(r, theta, phi, RN, g10, g11, h11);
    [Br_quadrupole, Btheta_quadrupole, Bphi_quadrupole] = quadrupole_field(r, theta, phi, RN, g20, g21, h21, g22, h22);
    
    Br =  Br_dipole + Br_quadrupole;
    Btheta = Btheta_dipole + Btheta_quadrupole;
    Bphi = Bphi_dipole + Bphi_quadrupole;

end