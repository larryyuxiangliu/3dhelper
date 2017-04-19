% This is a test code for 3D acoustic scattering from axisymmetric objects;
%
% Using this code to pick a good value of M, N, P, hs for different shape
% and different wavenumber k;
%
% Notes: This is to test: use BOR technique (FFT) to test a non-periodic
% case
%
% Notes: This is to deal with transmission scattering problem;
%
% By Larry Liu 06/19/2014

close all;
clear all;
addpath('/home/larry/MATLAB/fmm');
filename = 'trans_nonper_results.txt';
fid = fopen(filename, 'a');
tic;
%%%% =========== setup of the parameters ==================================

N = 400;    % the number of target points on the objects along 1 direction; so totally M1^2 pts;
M = 1.2*N;    % number of MFS points inside of the objects along 1 direction; so totally N1^2 pts;
P = 200;
Q = 1200;
Rt = 0.5;    % the scaling factor of the object inside of the box
hp = 0.046;
hm = -0.046;
kp = 30;  theta = -pi/4;   phi = pi/3;
km = 40;
% theta and phi are the direction of the incident plane wave traveling in
% polar coordinate system
kx = kp * cos(theta) * cos(phi);  ky = kp * cos(theta) * sin(phi);  kz = kp * sin(theta);

%%%% ================= Get the target and MFS points ======================
a = 0.3; w = 8; p = 0;
f = @(s) 1+a*cos(w*(s-p));
f_t = @(s) -a*w*sin(w*(s-p));
fprintf(fid, 'The parameters are M = %d, N = %d, P = %d, Q = %d, Rt = %2.2f, hp = %2.2f, hm = %2.2f, kp = %d, km = %d, w = %d \n', ...
    M, N, P, Q, Rt, hp, hm, kp,km,w);

target = curve3DQuad(f, f_t, M, P, Rt);           % curve target points and its normal derivatives;
xt = target.x;  yt = target.y;  zt = target.z;
source_p = source3DQuad(f, f_t, N, Q, Rt, hp);      % curve MFS points, propotional to local speed;
xsp = source_p.x;  ysp = source_p.y;  zsp = source_p.z;
source_m = source3DQuad(f, f_t, N, Q, Rt, hm);      % curve MFS points, propotional to local speed;
xsm = source_m.x;  ysm = source_m.y;  zsm = source_m.z;

targetBOR = BORcurve3DQuad(f, f_t, M, Rt);           % boundary pts in the r-z plane;
sourceBOR_p = BORsource3DQuad(f, f_t, N, Rt, hp);  % MFS pts in the r-z plane;
sourceBOR_m = BORsource3DQuad(f, f_t, N, Rt, hm);
An = EvalMFSHelm3DKernalBOR_trans(targetBOR, sourceBOR_p, sourceBOR_m, kp,km, P, Q);
[f,g] = QuasiPerVecFilling3DBOR_trans(targetBOR, kx, ky, kz, P);
fn = fft(f, [], 2)/P;
gn = fft(g, [], 2)/P;
bn = zeros(2*M, P);
c = zeros(2*N, P);
for n=1:P
    bn (1:M, n) = fn(:,n);
    bn (M + 1: 2*M, n) = gn(:,n);
    c(:,n) = An(:,:,n)\bn(:,n);
end

norm(c)
toc;
fprintf(fid, 'Doing matrix solving takes %4.4f seconds. \n', toc);
error = zeros(P, 1);
for n=1:P
    error(n) = norm(An(:,:,n)*c(:,n)-bn(:,n))/sqrt(length(bn(:,n)));
end
norm(error)/sqrt(P)

%cc = coeffBORto3D(c, Q);

c1 = c(1:N, :);
c2 = c(N+1:2*N, :);
cc1 = coeffBORto3D(c1, Q);
cc2 = coeffBORto3D(c2, Q);

%%% ============================Get the matrix and vector =================


if 1,
    tic;
    M_test = 128; P_test = 128;
    f = @(s) 1+a*cos(w*(s-p));
    f_t = @(s) -a*w*sin(w*(s-p));
    target_test = curve3DQuad(f, f_t, M_test, P_test, Rt);
    [u_exact, un_exact] = QuasiPerVecFilling3D_trans(target_test, kx, ky, kz);
    
    x_test = target_test.x;  y_test = target_test.y;  z_test = target_test.z;
    [up, Ep] = fmm3DHelpoteval(xsp, ysp, zsp, cc1, kp, x_test, y_test, z_test);
    [um, Em] = fmm3DHelpoteval(xsm, ysm, zsm, cc2, km, x_test, y_test, z_test);
    
    u = up - um;
    u = u.';
    
    Ep = -Ep.'; % change from a row vector to a column vector;  negative sign does not know why;
    Em = -Em.';
    
    unp = Ep(:,1).*target_test.nx + Ep(:,2).*target_test.ny + Ep(:,3).*target_test.nz;
    unm = Em(:,1).*target_test.nx + Em(:,2).*target_test.ny + Em(:,3).*target_test.nz;
    
    un = unp - unm;
    err1 = norm(u - u_exact)/sqrt (length(u))
    err2 = norm(un - un_exact)/sqrt (length(un))
    toc;
    fprintf(fid, 'Doing err evalution takes %4.4f seconds. \n', toc);
    fprintf(fid, 'boundary matching error is %e, normal derivative matching error is %e \n \n', err1, err2);
end
