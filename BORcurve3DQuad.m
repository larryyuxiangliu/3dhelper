function t = BORcurve3DQuad(f, f_t, N, RR)
% BORcurve3DQuad(f, N, RR) returns a structure that contains the locations
% of the boundary pts in the r-z plane; 
%
% Larry Liu, 06/11/2014 

tt=((1:N)-1/2)/N*pi-pi/2;
%tt = (1:N)/N * pi - pi/2;  
tt = reshape(tt, N, 1); 
t.r=RR * f(tt).*cos(tt); 
t.z=RR * f(tt).*sin(tt);
t.rb = -f(tt).*sin(tt) + f_t (tt).* cos(tt); 
t.zb = f(tt).* cos(tt) + f_t (tt).* sin(tt); 
t.rn = t.zb ./ sqrt (t.rb.^2 + t.zb.^2); 
t.zn = -t.rb./sqrt (t.rb.^2 + t.zb.^2); 

