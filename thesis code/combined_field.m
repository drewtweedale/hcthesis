% Function to calculate combined magnetic field (dipole + quadrupole)
function B = combined_field(r, B0, RN, eta, g10, g11, h11, g20, g21, h21, g22, h22)
    B_dipole = dipole_field(r, RN, g10, g11, h11);
    % B_dipole = dipole_field(r, RN, B0);
    B_quadrupole = quadrupole_field(r, RN, g20, g21, h21, g22, h22);
    % B_quadrupole = quadrupole_field(r, RN);
    B = B_dipole + B_quadrupole;
end