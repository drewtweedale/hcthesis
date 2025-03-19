%%Quadrupole Geometry eta=0 Helper Function for Particle Tracker Model.
function [Bx,By,Bz]=quadrupoleeta0(x,y,z,t)
Me=(3.12e-5)*1;
Re=6371000;
w=7.272e-5;
Bx=-(-3.*x.*(x.^2+y.^2-4.*z.^2).*Me.*Re.^4)./(4.*(sqrt(x.^2+y.^2+z.^2)).^7);
By=-(-3.*y.*(x.^2+y.^2-4.*z.^2).*Me.*Re.^4)./(4.*(sqrt(x.^2+y.^2+z.^2)).^7);
Bz=-(-3.*z.*(3.*x.^2+3.*y.^2-2.*z.^2).*Me.*Re.^4)./(4.*(sqrt(x.^2+y.^2+z.^2)).^7);
end