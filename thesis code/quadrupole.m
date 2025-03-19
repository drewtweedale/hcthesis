%%Quadrupole Geometry Helper Function for Quadrupole Model.
function [Bx,By,Bz]=quadrupole(x,y,z,t)
Me=(3.12e-5)*1;
Re=6371000;
w=7.272e-5;
angle=pi/2;
Bx=-(2.5.*x.*(z.^2-y.^2).*Me.*Re.^4)./((sqrt(x.^2+y.^2+z.^2)).^7)*sin(angle)+((x.^2.*y-1.5.*y.^3+3.5.*y.*z.^2).*Me.*Re.^4)./((sqrt(x.^2+y.^2+z.^2)).^7)*cos(angle);
By=-((x.^2.*y-1.5.*y.^3+3.5.*y.*z.^2).*Me.*Re.^4)./((sqrt(x.^2+y.^2+z.^2)).^7)*sin(angle)+(2.5.*x.*(z.^2-y.^2).*Me.*Re.^4)./((sqrt(x.^2+y.^2+z.^2)).^7)*cos(angle);
Bz=-((-x.^2.*z+1.5.*z.^3-3.5.*z.*y.^2).*Me.*Re.^4)./((sqrt(x.^2+y.^2+z.^2)).^7);
% Bx=(x.*(3.*x.^2-2.*y.^2-7.*z.^2))./(2.*sqrt(x.^2+y.^2+z.^2).^7);
% By=(5.*y.*(x.^2-z.*2))./(2.*sqrt(x.^2+y.^2+z.^2).^7);
% Bz=(z.*(7.*x.^2+2.*y.^2-3.*z.^2))./(2.*sqrt(x.^2+y.^2+z.^2).^7);
end