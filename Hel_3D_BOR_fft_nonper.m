% This is a test code for 3D acoustic scattering from axisymmetric objects;
%
% Using this code to pick a good value of M, N, P, hs for different shape
% and different wavenumber k;
%
% Notes: This is to test: use BOR technique (FFT) to test a non-periodic
% case and find out the best parameters;
%
% By Larry Liu 05/16/2014

close all;
clear all;
tic;

%%%% =========== setup of the parameters ==================================
N = 100; 
M = 1.2 * N; 
P = 100;    % number of Fourier Modes; 
Q = 100;    % number of trap rule quadrature pts;  
Rt = 0.65;     % the scaling factor of the object inside of the box
hs = 0.2;  % the sources are inside of the boundary, so it is a minus; 
k = 20;  theta = -pi/4;   phi = pi/3;
% theta and phi are the direction of the incident plane wave traveling in
% polar coordinate system
kx = k * cos(theta) * cos(phi);  ky = k * cos(theta) * sin(phi);  kz = k * sin(theta);

%%%% ================= Get the target and MFS points ======================
%a = 0.2; w = 2; p = pi/2;
%a = 0.3; w = 4; p = 0;
a = 0;  w = 0; p = 0; 
f = @(s) 1+a*cos(w*(s-p));
f_t = @(s) -a*w*sin(w*(s-p));
% std a=0.3, w = 4, p = 0;

source = source3DQuad(f, f_t, N, Q, Rt, hs);      % the MFS source points in 3D space; 
xs = source.x;  ys = source.y;  zs = source.z;

targetBOR = BORcurve3DQuad(f, f_t, M, Rt);        % boundary pts in the rho-z plane;
sourceBOR = BORsource3DQuad(f, f_t, N,Rt, hs);    % MFS pts in the rho-z plane;
BORrs = sourceBOR.r;     BORzs = sourceBOR.z; 
An = EvalMFSHelm3DKernalBOR(targetBOR, BORrs,BORzs, k, P, Q);
f = QuasiPerVecFilling3DBOR(targetBOR, kx, ky, kz, P);
fn = fft(f, [], 2)/P;
c = zeros(N,P);
for n=1:P
    c(:,n) = An(:,:,n)\fn(:,n);
end
norm(c)
error = zeros(P,1);
for n=1:P
    error(n) = norm(An(:,:,n)*c(:,n)-fn(:,n));
end
norm(error)

cc = coeffBORto3D(c, Q);

%%% ============================Get the matrix and vector =================
if 1,
    M_test = 120;
    f = @(s) 1+a*cos(w*(s-p));
    f_t = @(s) -a*w*sin(w*(s-p));
    target_test = curve3DQuad(f, f_t, M_test, P, Rt);
    b_exact = QuasiPerVecFilling3D(target_test, kx, ky, kz); 
    %b_mfs = Dir3DHelPotEval(target_test, xs, ys, zs, cc, k);
    x_test = target_test.x;  y_test = target_test.y;  z_test = target_test.z;
    b_mfs = fmm3DHelpoteval(xs, ys, zs, cc, k, x_test, y_test, z_test);
    b_mfs = b_mfs.'; % change from a row vector to a column vector; 
    norm(b_exact-b_mfs)/sqrt(length(b_exact))
end

toc;
