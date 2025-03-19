function B = quadrupole_field(r, RN, eta)
    x = r(1); y = r(2); z = r(3);
    Me = 3.12e-5; % Magnetic moment (scaled)
    Re = RN; % Neptune radius
    angle = eta * pi/2; % Tilt angle
   
    % Quadrupole field formulas with eta
    Bx = -(2.5 * x * (z^2 - y^2) * Me * Re^4) / (sqrt(x^2 + y^2 + z^2))^7 * sin(angle) + ...
          ((x^2 * y - 1.5 * y^3 + 3.5 * y * z^2) * Me * Re^4) / (sqrt(x^2 + y^2 + z^2))^7 * cos(angle);
   
    By = -((x^2 * y - 1.5 * y^3 + 3.5 * y * z^2) * Me * Re^4) / (sqrt(x^2 + y^2 + z^2))^7 * sin(angle) + ...
          (2.5 * x * (z^2 - y^2) * Me * Re^4) / (sqrt(x^2 + y^2 + z^2))^7 * cos(angle);
   
    Bz = -((-x^2 * z + 1.5 * z^3 - 3.5 * z * y^2) * Me * Re^4) / (sqrt(x^2 + y^2 + z^2))^7;
   
    B = [Bx, By, Bz];
end