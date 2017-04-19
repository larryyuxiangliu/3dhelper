% This code deals with 3D acoustic scattering problem (Helmholtz equation)
% from periodic strucutres where the periodicity only appears in 2D.
% We need to impose Rayleigh-Bloch Expansion on the top and bottom and also
% quasi-periodic condition on the left and right as well as front and back
% of the box;
%
% Notes: This code uses the GMRES-FMM technique, details about this
% technique can be seen in Hel_3D_Periodic_14X/notes.pdf;
%
% Notes: This code uses BOR technique, which is use FFT to solve the linear
% system in each Fourier Mode;
%
% By Larry Liu 07/14/2014
close all;
clear all;
%addpath('/home/larry/MATLAB/fmm');
addpath('/home/larry/MATLAB/fmmlib3d-1.2/matlab');  % maxlevel = 3;
%addpath('/home/larry/MATLAB/mex_lowrankC');

tic;
%%%% =========== setup of the parameters ==================================
N1 = 60;    % the number of boundary pts on the object in the r-z plane;
M1 = 1.2 * N1;    % the number of MFS pts in the r-z plane;
P = 60;     % number of Fourier Modes;
q = 220; % for each number in the Fourier Modes;
PP = 70;
N0 = 19;     % N0 is related to the number of the RB expansion coefficients
M2 = ceil(sqrt(pi) * N0);
L = 1;       % this is half the periodicity in both directions
R = 3.5;     % the radius of the proxy points;
Rt = 0.5;    % the scaling factor of the object inside of the box
hs = 0.1;
nei = 1;     % 0 or 1
ex = 2*L;  ey = 2*L;  ez = 0;    % periodic unit vector
k = 10;  theta = -pi/4;   phi = pi/3;
N_iter = 50;
type = 'd';
acc = 1e-11;
RP = 0.75;
NP = 25;
RQ = 0.75;
NQ = 25;
% theta and phi are the direction of the incident plane wave traveling in
% polar coordinate system
kx = k * cos(theta) * cos(phi);  ky = k * cos(theta) * sin(phi);  kz = k * sin(theta);
% the x, y, z component of the wavenumber vector k;
NN = getRBdegree(N0, kx, ky, ex, ey);     % NN is the number of coefficients for RB coefficients at top / bottom, so totally 2 * NN.

alpha = exp(1i*kx*ex);       % the quasi-periodic phase in x direction
beta = exp(1i*ky*ey);        % the quasi-periodic phase in y direction
%%%% ======================= END ==========================================

%%%% ================= Get the target and MFS points ======================
%a = 0.3;  w = 4; p = 0;
a = 0.3;   w = 4;   p = 0;
%a = 0.2; w = 2; p = pi/2;
f = @(s) 1+a*cos(w*(s-p));
f_t = @(s) -a*w*sin(w*(s-p));
%aa = 0.2;  % thickness: cannot be >1/3 otherwise not smooth to emach
%bb = pi/6;  % controls approx opening angle in radians (keep small for resonant)
[TT,BO,LL,RR,FF,BA] = squarequad(M2, L);   % box target points;

bdy3D = curve3DQuad(f, f_t, M1, P, Rt);           % curve target points and its normal derivatives;
src3D = source3DQuad(f, f_t, N1, P, Rt, hs);      % curve MFS points, propotional to local speed;

% bdyBOR = cupBdyBOR(aa, bb, M1, Rt);
% srcBOR = srcBdyBOR(aa, bb, N1, h, Rt);
% src3D = srcBdy3D(aa, bb, N1, h, P, Rt);
% bdy3D = cupBdy3D(aa, bb, M1, P, Rt);

xs = src3D.x;  ys = src3D.y;  zs = src3D.z;
xt = bdy3D.x;  yt = bdy3D.y;  zt = bdy3D.z;


% ========================== END ==========================================

bb = rand(N1*P,1);

%tic; uu = applyAelse(xs,ys,zs,bb, k, bdy3D,nei, alpha, beta,ex,ey, type); toc;

tic;

[U, K, V, IS, IT, Rk] = getSkeleton(bdy3D, src3D, k, acc, RP, NP, RQ, NQ,...
    nei,alpha, beta, ex, ey, type);
toc;
tic;
u = applyAelse_matvec(U, K, V, IS, IT, bb, nei, bdy3D, type);
toc;
Rk
%norm(u - uu)/sqrt (length(u))
figure;   % real souce to the proxy sphere;
axes1 = axes('Parent',figure,'PlotBoxAspectRatio',[434 342.3 273.84],...
    'DataAspectRatio',[1 1 1]);
view(axes1,[0.5 0]);
hold(axes1,'all');
axis equal;
[xp, yp, zp] = squaresourcequad(NP, RP);
source = src3D;
i = -1;   j = 0;
xs = source.x + i*ex;  ys = source.y + j*ey ;  zs = source.z;
hold on;
plot3(xs, ys, zs, 'r.');
plot3(xt, yt, zt, 'g.');
plot3(xp, yp, zp, 'b*')
legend('srcs', 'targets', 'proxy target')
xlabel('x')
ylabel('y')
zlabel('z')
saveas(gcf, '/Users/Yuxiang/Documents/AAResearch/ID_15S/getv.eps', 'epsc2');


figure;   % real souce to the proxy sphere;
axes1 = axes('Parent',figure,'PlotBoxAspectRatio',[434 342.3 273.84],...
    'DataAspectRatio',[1 1 1]);
view(axes1,[0.5 0]);
hold(axes1,'all');
axis equal;
[xp, yp, zp] = squaresourcequad(NQ, RQ);
i = -1; j = 0;
xp = xp + i * ex;   yp = yp + j * ey;    zp = zp;
xs = source.x + i*ex;  ys = source.y + j*ey ;  zs = source.z;
target = bdy3D;
xt = target.x;  yt = target.y;  zt = target.z;
hold on;
plot3(xs, ys, zs, 'g.');
plot3(xt, yt, zt, 'b.');
plot3(xp, yp, zp, 'r*')
legend('srcs', 'targets', 'proxy src')
xlabel('x')
ylabel('y')
zlabel('z')
saveas(gcf, '/Users/Yuxiang/Documents/AAResearch/ID_15S/getu.eps', 'epsc2');



ii = 1;
figure;
for i = -nei:nei
    for j = - nei:nei
        if (i~=0) || (j~=0)
            hold on;
            rs = Rk(ii,1);
            xs = source.x + i*ex;  ys = source.y + j*ey ;  zs = source.z;
            xsf = xs(IS{ii}(1:rs));    ysf = ys(IS{ii}(1:rs));  zsf = zs(IS{ii}(1:rs));
            rt = Rk(ii,2);
            xtf = xt(IT{ii}(1:rt));   ytf = yt(IT{ii}(1:rt));     ztf = zt(IT{ii}(1:rt));
            plot3(xsf, ysf, zsf, 'r*', xtf, ytf, ztf, 'b.');
            ii = ii + 1;
            
        end
    end
end
axis equal;
xlabel('x')
ylabel('y')
zlabel('z')
saveas(gcf, '/Users/Yuxiang/Documents/AAResearch/ID_15S/skeleton_all.eps', 'epsc2');



ii = 1;
figure;
for i = -nei:nei
    for j = - nei:nei
        if (i==-1) && (j==0)
            hold on;
            rs = Rk(ii,1);
            xs = source.x + i*ex;  ys = source.y + j*ey ;  zs = source.z;
            xsf = xs(IS{ii}(1:rs));    ysf = ys(IS{ii}(1:rs));  zsf = zs(IS{ii}(1:rs));
            rt = Rk(ii,2);
            xtf = xt(IT{ii}(1:rt));   ytf = yt(IT{ii}(1:rt));     ztf = zt(IT{ii}(1:rt));
            plot3(xsf, ysf, zsf, 'r*', xtf, ytf, ztf, 'b.');
            ii = ii + 1;
            
        end
    end
end
axis equal;
xlabel('x')
ylabel('y')
zlabel('z')

saveas(gcf, '/Users/Yuxiang/Documents/AAResearch/ID_15S/skeleto_one.eps', 'epsc2');