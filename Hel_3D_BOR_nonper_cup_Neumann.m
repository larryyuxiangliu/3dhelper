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

N = 250;    % number of MFS points inside of the objects along 1 direction; so totally N1^2 pts;
M = 1.2 * N;    % the number of target points on the objects along 1 direction; so totally M1^2 pts;
P = 120;
Q = 250;
Rt = 0.5;    % the scaling factor of the object inside of the box
h = 0.1;
k = 30;  theta = -pi/4;   phi = pi/3;
% theta and phi are the direction of the incident plane wave traveling in
% polar coordinate system
kx = k * cos(theta) * cos(phi);  ky = k * cos(theta) * sin(phi);  kz = k * sin(theta);

%%%% ================= Get the target and MFS points ======================
aa = 0.2;  % thickness: cannot be >1/3 otherwise not smooth to emach
bb = pi/6;  % controls approx opening angle in radians (keep small for resonant)
bdyBOR = cupBdyBOR(aa, bb, M, Rt); 
srcBOR = srcBdyBOR(aa, bb, N, h, Rt); 
src3D = srcBdy3D(aa, bb, N, h, Q, Rt); 
bdy3D = cupBdy3D(aa, bb, M, P, Rt); 

xs = src3D.x;  ys = src3D.y;  zs = src3D.z;

An = EvalMFSHelm3DKernalBOR_Neumann(bdyBOR, srcBOR, k, P, Q);
[~,g] = QuasiPerVecFilling3DBOR_trans(bdyBOR, kx, ky, kz, P);
gn = fft(g, [], 2)/P;
c = zeros(N,P);
for n = 1:P
    [AU(:,:,n), AS, AV(:,:,n)] = svd(An(:,:,n), 0);
    ss = diag(AS);
    nA = length(ss);
    rA = length(find (ss >1e-10));   % choose a cut, 1e-10.
    iA(:,n) =  [1./ss(1:rA); zeros(nA-rA,1)];
end    

for n=1:P
    %c(:,n) = An(:,:,n)\gn(:,n);
    c(:,n) = AV(:,:,n) * (iA(:,n) .* (AU(:,:,n)' * gn(:,n) ));
end
norm(c)
error = zeros(P,1);
for n=1:P
    error(n) = norm(An(:,:,n)*c(:,n)-gn(:,n));
end
norm(error)

cc = coeffBORto3D(c, Q);

%%% ============================Get the matrix and vector =================
if 1,
    M_test = 128;  P_test = 128;
    
    target_test = cupBdy3D(aa, bb, M_test, P_test, Rt);
    [~,bn] = QuasiPerVecFilling3D_trans(target_test, kx, ky, kz);
    x_test = target_test.x;  y_test = target_test.y;  z_test = target_test.z;
    [u, E] = fmm3DHelpoteval(xs, ys, zs, cc, k, x_test, y_test, z_test);
    E = -E.';  % change from a row vector to a column vector;
    un = E(:,1).*target_test.nx + E(:,2).*target_test.ny + E(:,3).*target_test.nz;
    err = norm(un-bn)/sqrt(length(bn))
    
end

toc;
