function bdy = cupBdy(a, b, N, Rt)

N = N/2; 
s = ((1:N)-0.5)/N * pi;  % note half-offset, needed for easy reflection abt z
r = Rt* (1 - a*erf((s-pi/2)/a));  % radius: starts at 1+a, ends at 1-a
c = a; %*(1-b/pi);  % is theta rounding scale
sabs = @(x) exp(-(x/c).^2)*c/sqrt(pi)+x.*erf(x/c); % c-smoothed absval
th = b-a + 2*(1-(b-a)/pi)*sabs(s-pi/2);
%th = b + sqrt(log(exp((2*(1-b/pi)*(s-pi/2)).^2) + exp(c^2))); % theta

% coords in (rho,z) plane:
rho = r.*sin(th); z = r.*cos(th);  % theta down from z axis as in 3D cyl coords

Z = [rho -rho(end:-1:1)] + 1i*[z z(end:-1:1)]; % complex coords of full curve
% (appropriate for half-integer offset
Zp = perispecdiff(Z);    % take derivative to get Z' at each node

Zn = (Zp/1i)./abs(Zp);   % unit normals

bdy.x = real(Z); 
bdy.y = imag(Z); 
bdy.nx = real(Zn); 
bdy.ny = imag(Zn); 
bdy.px = real(Zp)./abs(Zp); 
bdy.py = imag(Zp)./abs(Zp);

bdy.x = bdy.x'; 
bdy.y = bdy.y'; 
bdy.px = bdy.px'; 
bdy.py = bdy.py'; 
bdy.nx = bdy.nx'; 
bdy.ny = bdy.ny';

end
