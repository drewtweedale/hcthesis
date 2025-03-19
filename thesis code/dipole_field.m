function B = dipole_field(r, B0, RN)
    x = r(1); y = r(2); z = r(3);
    r_mag = sqrt(x^2 + y^2 + z^2); 
    scale = (r_mag / RN)^3 * r_mag^2; % Scaling factor
   
    Bx = (-3 * B0 * x * z) / scale;
    By = (-3 * B0 * y * z) / scale;
    Bz = (B0 * (x^2 + y^2 - 2 * z^2)) / scale;
   
    B = [Bx, By, Bz];
end