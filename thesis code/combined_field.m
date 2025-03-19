% Function to calculate combined magnetic field (dipole + quadrupole)
function B = combined_field(r, B0, RN, eta)
    B_dipole = dipole_field(r, B0, RN);
    B_quadrupole = quadrupole_field(r, RN, eta);
    B = B_dipole + B_quadrupole;
end