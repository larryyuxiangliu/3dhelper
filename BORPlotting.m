function BORPlotting()
clear all;
close all;

% ============= Plot w = 4 ================================================
N = 50;
M = 1.2 * N;    % the number of target points on the objects along 1 direction; so totally M1^2 pts;
Rt = 1;    % the scaling factor of the object inside of the box
hs = 0.1;
k = 10;
a = 0.3; w = 4; p = 0;
f = @(s) 1+a*cos(w*(s-p));
f_t = @(s) -a*w*sin(w*(s-p));

targetBOR = BORcurve3DQuad(f, f_t, M, Rt);           % boundary pts in the r-z plane;
sourceBOR = BORsource3DQuad(f, f_t, N, Rt, hs);  % MFS pts in the r-z plane;

plot(targetBOR.r, targetBOR.z, '.', sourceBOR.r, sourceBOR.z, 'r*');
hold on;
set(gca, 'FontSize', 16);
h1 = xlabel('$\rho$');
h2 = ylabel('$z$');
set(h1,'Interpreter','latex');
set(h2,'Interpreter','latex');
%h3 = title('the bdy pts and src pts in the $\rho-z$ plane');
%set(h3,'Interpreter','latex');
legend('boundary points', 'source points');
axis equal;

saveas(gcf, '/Users/Yuxiang/Documents/AAResearch/Paper1_acoustic periodic scattering/plotw4.eps', 'epsc2');


% ============= Plot w = 8 ================================================
w = 8;
hs = 0.045;
f = @(s) 1+a*cos(w*(s-p));
f_t = @(s) -a*w*sin(w*(s-p));
targetBOR = BORcurve3DQuad(f, f_t, M, Rt);           % boundary pts in the r-z plane;
sourceBOR = BORsource3DQuad(f, f_t, N, Rt, hs);  % MFS pts in the r-z plane;
figure;
plot(targetBOR.r, targetBOR.z, '.', sourceBOR.r, sourceBOR.z, 'r*');
set(gca, 'FontSize', 16);
h1 = xlabel('$\rho$');
h2 = ylabel('$z$');
set(h1,'Interpreter','latex');
set(h2,'Interpreter','latex');
%h3 = title('the bdy pts and src pts in the $\rho-z$ plane');
%set(h3,'Interpreter','latex');
legend('boundary points', 'source points');
axis equal;

saveas(gcf, '/Users/Yuxiang/Documents/AAResearch/Paper1_acoustic periodic scattering/plotw8.eps', 'epsc2');


% ============= Plot the cup shape=========================================
aa = 0.2;  % thickness: cannot be >1/3 otherwise not smooth to emach
bb = pi/6;  % controls approx opening angle in radians (keep small for resonant)
h = 0.08;
bdyBOR = cupBdyBOR(aa, bb, M, Rt);
srcBOR = srcBdyBOR(aa, bb, N, h, Rt);
figure;
plot(bdyBOR.r, bdyBOR.z, '.', srcBOR.r, srcBOR.z, 'r*');
set(gca, 'FontSize', 16);
h1 = xlabel('$\rho$');
h2 = ylabel('$z$');
set(h1,'Interpreter','latex');
set(h2,'Interpreter','latex');
%h3 = title('the bdy pts and src pts in the $\rho-z$ plane');
%set(h3,'Interpreter','latex');
legend('boundary points', 'source points');
axis equal;

saveas(gcf, '/Users/Yuxiang/Documents/AAResearch/Paper1_acoustic periodic scattering/plotcup.eps', 'epsc2');
end


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

end



function sourceBOR = BORsource3DQuad(f, f_t, N, RR, hs)
% BORsource3DQuad(f, f_t, N,RR, hs) returns a structure that contains the locations
% of the MFS source pts in the r-z plane;
%
% Larry Liu, 06/11/2014
if 0,    % doing porportional to local speed
    t=((1:N)-1/2)/N*pi-pi/2;
    %t = (1:N)/N *pi - pi/2;
    %t = t*1.03;
    
    t = reshape(t, N,1);
    
    ss.r = RR * f(t).*cos(t);
    ss.z = RR * f(t).*sin(t);
    
    ss.tr = (-f(t).*sin(t) + f_t(t).*cos(t));
    ss.tz = (f(t).*cos(t) + f_t(t).*sin(t));
    
    
    ll = sqrt(ss.tr.^2 +ss.tz.^2);
    
    ss.nr = ss.tz./ll;
    ss.nz = -ss.tr./ll;
    
    h = hs * ll;
    
    
    sourceBOR.r = ss.r - h.* ss.nr;
    sourceBOR.z = ss.z - h.* ss.nz;
else   % using complexification
    if 0,
        
        t=((1:2*N)-1/2)/(2*N) *2*pi - pi/2;
        %t = (1:N)/N *pi - pi/2;
        %t = t*1.03;
        
        t = reshape(t, 2*N,1);
        
        ss.r = RR * f(t).*cos(t);
        ss.z = RR * f(t).*sin(t);
        
        z = ss.r + 1i * ss.z;
        zz = imagshiftcurve(z,hs);
        sourceBOR.r = real(zz(1:N));
        sourceBOR.z = imag(zz(1:N));
    else
        
        t = 1i * hs + ((1:N)-1/2)/N*pi-pi/2;
        t = reshape(t, N, 1);
        t = RR * exp(1i*t).*f(t);
        sourceBOR.r = real(t);
        sourceBOR.z = imag(t);
        
    end
end
end




function srcBOR = srcBdyBOR(a, b, N, h, Rt)

bb = cupBdy(a, b, 2*N, Rt);
z = bb.x + 1i * bb.y;
zz = imagshiftcurve(z,h);
srcBOR.r = real(zz(1:N));
srcBOR.z = imag(zz(1:N));

end


function bdyBOR = cupBdyBOR(a, b, N, Rt)




bdy = cupBdy(a, b, 2*N, Rt);


bdyBOR.r = bdy.x(1:N);
bdyBOR.z = bdy.y(1:N);
bdyBOR.rn = bdy.nx(1:N);
bdyBOR.zn = bdy.ny(1:N);
bdyBOR.rb = bdy.px(1:N);
bdyBOR.zb = bdy.py(1:N);

end


