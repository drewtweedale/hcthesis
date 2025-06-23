% Function to calculate combined magnetic field (dipole + quadrupole)
function [Br, Btheta, Bphi] = combined_field(r, B0, RN, g10, g11, h11, g20, g21, h21, g22, h22)
    % USING SCHMIDT COEFFICIENTS
    [Br_dipole, Btheta_dipole, Bphi_dipole] = dipole_field(r, RN, B0, g10, g11, h11);
    [Br_quadrupole, Btheta_quadrupole, Bphi_quadrupole] = quadrupole_field(r, RN, B0, g20, g21, h21, g22, h22);

    % SYMMETRIC CASE
    % [Br_dipole, Btheta_dipole, Bphi_dipole] = dipole_field(r, RN, B0);
    % [Br_quadrupole, Btheta_quadrupole, Bphi_quadrupole] = quadrupole_field(r, RN, B0);

    Br =  Br_dipole + Br_quadrupole;
    Btheta = Btheta_dipole + Btheta_quadrupole;
    Bphi = Bphi_dipole + Bphi_quadrupole;

end