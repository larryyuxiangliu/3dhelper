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
addpath('/home/larry/MATLAB/fmmlib3d-1.2/matlab');  %maxlevel = 2; 
tic;
filename = 'dirichilet_per_results_ID.txt';
fid = fopen(filename, 'a');
%%%% =========== setup of the parameters ==================================
N1 = 100;    % the number of boundary pts on the object in the r-z plane;
M1 = 1.2 * N1;    % the number of MFS pts in the r-z plane;
P = 100;     % number of Fourier Modes;
q = 100;
PP = 20;
N0 = 9;     % N0 is related to the number of the RB expansion coefficients
M2 = ceil(sqrt(pi) * N0); 
L = 1;       % this is half the periodicity in both directions
R = 3.5;     % the radius of the proxy points;
Rt = 0.65;    % the scaling factor of the object inside of the box
hs = 0.2;
nei = 1;     % 0 or 1
ex = 2*L;  ey = 2*L;  ez = 0;    % periodic unit vector
k = 3;  theta = -pi/4;   phi = pi/3;
N_iter = 50;
type = 'd';            % menas Dirichilet problem;
% theta and phi are the direction of the incident plane wave traveling in
% polar coordinate system
acc = 1e-11; 
RP = 0.75; 
NP = 25; 
RQ = 0.75; 
NQ = 25;

kx = k * cos(theta) * cos(phi);  ky = k * cos(theta) * sin(phi);  kz = k * sin(theta);
% the x, y, z component of the wavenumber vector k;
NN = getRBdegree(N0, kx, ky, ex, ey);     % NN is the number of coefficients for RB coefficients at top / bottom, so totally 2 * NN.

alpha = exp(1i*kx*ex);       % the quasi-periodic phase in x direction
beta = exp(1i*ky*ey);        % the quasi-periodic phase in y direction
%%%% ======================= END ==========================================

%%%% ================= Get the target and MFS points ======================
%a = 0.3;  w = 4; p = 0;
a = 0;   w = 0;   p = 0; 
%a = 0.2; w = 2; p = pi/2;
f = @(s) 1+a*cos(w*(s-p));
f_t = @(s) -a*w*sin(w*(s-p));
[TT,BO,LL,RR,FF,BA] = squarequad(M2, L);   % box target points;
fprintf(fid, '---------------- Now maxlevel = 2 ----------------------------------- \n');
fprintf(fid, 'This is a special test for q = P when considering far field evalution. \n');
fprintf(fid, 'The parameters are M1 = %d, N1 = %d, P = %d, q = %d, M2 = %d, PP = %d, N0 = %d, Rt = %2.2f, hs = %2.2f, k = %d, N_iter = %d, w = %d \n', ...
    M1, N1, P, q, M2, PP, N0, Rt, hs,k, N_iter, w);

target = curve3DQuad(f, f_t, M1, P, Rt);           % curve target points and its normal derivatives;
targetBOR = BORcurve3DQuad(f, f_t, M1, Rt);           % boundary pts in the r-z plane;
sourceBOR = BORsource3DQuad(f, f_t, N1,Rt, hs);  % MFS pts in the r-z plane;
xt = target.x;  yt = target.y;  zt = target.z;
BORrs = sourceBOR.r;
BORzs = sourceBOR.z;

source = source3DQuad(f, f_t, N1, P, Rt, hs);      % curve MFS points, propotional to local speed;
xs= source.x;  ys = source.y;  zs = source.z;
% ========================== END ==========================================

% ============================= to fill the matrices ======================
A_0 = EvalMFSHelm3DKernalBOR(targetBOR, BORrs, BORzs, k, P, q);

B = QuasiPerB3D(target, PP,NN, k);
C = QuasiPerC3D(TT,BO,LL,RR,FF,BA, xs,ys,zs, k, ex, ey, nei, alpha, beta);
Q = QuasiPerQ3D(TT,BO,LL,RR,FF,BA, M2, PP, N0,NN, k, kx, ky, ex,ey, alpha, beta);

[U, K, V, IS, IT] = getSkeleton(target, source, k, acc, RP, NP, RQ, NQ,...
    nei,alpha, beta, ex, ey, type);


bb = QuasiPerVecFilling3D(target,kx, ky, kz);
toc;
fprintf(fid, 'filling the matrix takes %4.4f seconds. \n', toc);

% ============================== END ======================================


% ====================== doing SVD or QR on each block of A_0 and Q =======
tic;

for n=1:P
    
    [AU(:,:,n), AS, AV(:,:,n)] = svd(A_0(:,:,n), 0);
    ss = diag(AS);
    nA = length(ss);
    rA = length(find (ss >1e-10));   % choose a cut, 1e-10.
    iA(:,n) =  [1./ss(1:rA); zeros(nA-rA,1)];
    
end

[QU, QS, QV] = svd(Q, 0);
ss = diag(QS);
nQ = length(ss);
rQ = length(find (ss >1e-10));      % choose a cut, 1e-10.
iQ =  [1./ss(1:rQ); zeros(nQ-rQ,1)];

toc; 
% ==================== get the right hand side A_0^{\dag}*f ===============

fprintf(fid, 'Doing SVD on the matrix takes %4.4f seconds. \n', toc);
% ================================= END ===================================


% ========================== do the GMRES iteration =======================
tic;

%ff = @(x)  apply_GMRES_Matrix_BOR_SVD_test(target, AU, iA, AV, QU, iQ,QV, B, C, ...
%    M1, P, P, xs, ys, zs, [],[],[],x, k, [], nei, alpha, beta, ex, ey, type);

  ff = @(x)  apply_GMRES_Matrix_BOR_SVD_ID(target,AU,iA,AV,QU,iQ,QV, ...
        B, C, U, K,V,M1, P, P, x, nei,type); 

[xx, ~, ~, it, res] = gmres(ff, bb, [], 1e-10, N_iter);


XX = reshape(xx, P, M1);
XX = XX.';
xn = fft(XX, [],2)/P;
for n=1:P
    Axn(:,n) = AV(:,:,n) * (iA(:,n) .* (AU(:,:,n)' * xn(:,n) ));
end
eta = coeffBORto3D(Axn, P);   % this is used to do the far-field evalution;


y = C * eta;

xi = -QV * (iQ .* (QU' * y));

toc;


fprintf(fid, 'GMRES iteration takes %4.4f seconds and has %d iterations. \n', toc, it(2));
% ============================= END =======================================

tic;

if 1,
    M_test = 120;   P_test = 120;
    %a = 0.3;  w = 4; p = 0;
    a = 0;  w = 0;  p = 0; 
    f = @(s) 1+a*cos(w*(s-p));
    f_t = @(s) -a*w*sin(w*(s-p));
    error = getBdyError(M_test, P_test, f, f_t, Rt,hs,[],xs,ys,zs, k,[], kx, ky, kz, N1,P,q, PP, eta, xi, nei ,alpha, beta, ex, ey, type)
end



if 1,
    Mtest = 40;
    [err1, err2] = getPerError (Mtest, L,  xs, ys, zs, eta, xi, PP, k, nei ,alpha, beta, ex, ey)
    flux = getFlux(xi, PP, N0, NN, k, kx, ky, kz, ex, ey);
    err_f = abs(flux - abs(kz))
end

toc;
fprintf(fid, 'Doing error evaluation take %4.4f seconds \n', toc);
fprintf(fid, 'boundary matching error is %e, flux matching error is %e \n', error, err_f);
fprintf(fid, 'x-direction periodicity matching error is %e, y-direction periodicity matching error is %e.\n', err1, err2);

fprintf(fid, '------------------------------------------------------------ \n\n');
fclose(fid); 


%Periodic3DPlotting(targetBOR, ex, ey, kx, ky, kz, N1,P,PP, N0, NN, xs, ys, zs, k,nei, alpha,beta, eta, xi);
