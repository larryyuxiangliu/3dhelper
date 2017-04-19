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
%
% Modified into a simplifed version on 10/14/2014, which uses SVD,
% spherical harmonics and building matrcies blocks without leaving the
% options;

close all;
clear all;
%addpath('/home/larry/MATLAB/fmm');  % the old default one; 
addpath('/home/larry/MATLAB/fmmlib3d-1.2/matlab');  % maxlevel = 3;  
tic;
filename = 'trans_per_results_maxlevel3.txt';
fid = fopen(filename, 'a');
type = 't';
%%%% =========== setup of the parameters ==================================
N1 = 400;    % the number of boundary pts on the object in the r-z plane;
M1 = 1.2 * N1;    % the number of MFS pts in the r-z plane;
P = 200;     % number of Fourier Modes;
q = 1200;
PP = 86;
N0 = 25;     % N0 is related to the number of the RB expansion coefficients
M2 = ceil(sqrt(pi) * N0); 
L = 1;       % this is half the periodicity in both directions
R = 3.5;     % the radius of the proxy points;
Rt = 0.5;    % the scaling factor of the object inside of the box
hp = 0.045;
hm = -0.045;
nei = 1;     % 0 or 1
ex = 2*L;  ey = 2*L;  ez = 0;    % periodic unit vector
kp = 40;  theta = -pi/4;   phi = pi/3;
km = 60;
N_iter = 140;

% theta and phi are the direction of the incident plane wave traveling in
% polar coordinate system
kx = kp * cos(theta) * cos(phi);  ky = kp * cos(theta) * sin(phi);  kz = kp * sin(theta);
% the x, y, z component of the wavenumber vector k;
NN = getRBdegree(N0, kx, ky, ex, ey);     % NN is the number of coefficients for RB coefficients at top / bottom, so totally 2 * NN.

alpha = exp(1i*kx*ex);       % the quasi-periodic phase in x direction
beta = exp(1i*ky*ey);        % the quasi-periodic phase in y direction
%%%% ======================= END ==========================================

%%%% ================= Get the target and MFS points ======================
a = 0.3;  w = 8; p = 0;
f = @(s) 1+a*cos(w*(s-p));
f_t = @(s) -a*w*sin(w*(s-p));
[TT,BO,LL,RR,FF,BA] = squarequad(M2, L);            % box target points;
fprintf(fid, '================= A_else does not contain source points xsm, ysm, zsm =============================\n');
fprintf(fid, 'This is to test the case when q = P for far field evalution\n');

fprintf(fid, 'The parameters are M1 = %d, N1 = %d, P = %d, q = %d, M2 = %d, PP = %d, N0 = %d, Rt = %3.3f, hp = %3.3f, hm = %3.3f, kp = %d, km = %d, N_iter = %d, w = %d \n', ...
    M1, N1, P, q, M2, PP, N0, Rt, hp, hm, kp, km, N_iter, w);


bdy3D = curve3DQuad(f, f_t, M1, P, Rt);             % curve target points and its normal derivatives;
src3D_p = source3DQuad(f, f_t, N1, P, Rt, hp);      % curve MFS points, propotional to local speed;
src3D_m = source3DQuad(f, f_t, N1, P, Rt, hm);

xsp = src3D_p.x;  ysp = src3D_p.y;  zsp = src3D_p.z;
xsm = src3D_m.x;  ysm = src3D_m.y;  zsm = src3D_m.z;

src3D_pf = source3DQuad(f, f_t, N1, P, Rt, hp);
src3D_mf = source3DQuad(f, f_t, N1, P, Rt, hm);

xspf = src3D_pf.x;  yspf = src3D_pf.y;  zspf = src3D_pf.z;
xsmf = src3D_mf.x;  ysmf = src3D_mf.y;  zsmf = src3D_mf.z;
% these points are for far-field evalution, we don't need that big q;

xt = bdy3D.x;  yt = bdy3D.y;  zt = bdy3D.z;

bdyBOR = BORcurve3DQuad(f, f_t, M1, Rt);           % boundary pts in the r-z plane;
srcBOR_p = BORsource3DQuad(f, f_t, N1,Rt, hp);     % MFS pts in the r-z plane;
srcBOR_m = BORsource3DQuad(f, f_t, N1,Rt, hm);     % MFS pts in the r-z plane;
% ========================== END ==========================================

% ============================= to fill the matrices ======================
A_0 = EvalMFSHelm3DKernalBOR_trans(bdyBOR, srcBOR_p, srcBOR_m, kp,km, P, q);

[B, Bn] = QuasiPerB3D(bdy3D, PP, NN, kp);
B = [B; Bn];

C1 = QuasiPerC3D(TT,BO,LL,RR,FF,BA, xsp, ysp, zsp, kp, ex, ey, nei, alpha, beta);
C = [C1, zeros(size(C1,1), size(C1,2))];

Q = QuasiPerQ3D(TT,BO,LL,RR,FF,BA, M2, PP, N0, NN, kp, kx, ky, ex,ey, alpha, beta);

[b, bn] = QuasiPerVecFilling3D_trans(bdy3D, kx, ky, kz) ;
bb = [b; bn];
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
fprintf(fid, 'Doing SVD on the matrix takes %4.4f seconds. \n', toc);
% ================================= END ===================================


% ========================== do the GMRES iteration =======================

if 1,
    tic;
    ff = @(x)  apply_GMRES_Matrix_BOR_SVD_test(bdy3D, AU, iA, AV,QU,iQ,QV,...
        B, C, M1, P, P,  xsp, ysp, zsp, x, kp, ...
        nei, alpha, beta, ex, ey, type);
    
    [xx, ~, ~, it, res] = gmres(ff, bb, [], 1e-11, N_iter);
    
    XX = reshape(xx, P, 2*M1);
    XX = XX.';
    xn = fft(XX, [],2)/P;
    for n=1:P
        Axn(:,n) = AV(:,:,n) * (iA(:,n) .* (AU(:,:,n)' * xn(:,n) ));
    end
    eta = coeffBORto3D(Axn, P);   % this is used to do the far-field evalution;
    y = C * eta;
    
    xi = -QV * (iQ .* (QU' * y));
    
    toc;
end
fprintf(fid, 'GMRES iteration takes %4.4f seconds and has %d iterations. \n', toc, it(2));
% ============================= END =======================================
if 1,
    tic;
    M_test = 120;   P_test = 120;
    error = getBdyError(M_test, P_test, f, f_t,Rt,hp, hm, xsp,ysp,zsp, kp,km, kx, ky, kz,N1, P,q, PP,eta, xi, nei ,alpha, beta, ex, ey, type);
    
    
    Mtest = 40;
    [err1, err2] = getPerError (M_test, L,  xsp, ysp, zsp, eta(1:end/2), xi, PP, kp, nei ,alpha, beta, ex, ey);
    flux = getFlux(xi, PP, N0,NN, kp, kx, ky, kz, ex, ey);
    err_f = abs(flux - abs(kz));
    
    toc;
end
fprintf(fid, 'Doing error evaluation take %4.4f seconds \n', toc);
fprintf(fid, 'boundary matching error is %e. \n', error);
fprintf(fid, 'x-direction periodicity matching error is %e, y-direction periodicity matching error is %e.\n', err1, err2);
fprintf(fid, 'Periodicity matching condition is %e.\n', (err1 + err2)/2 );
fprintf(fid, 'flux matching error is %e. \n', err_f);



a1 = whos('A_0');  
a2 = whos('AU'); 
a3 = whos('AS'); 
a4 = whos('AV'); 
a5 = whos('B'); 
a6 = whos('C'); 
a7 = whos('Q'); 
a8 = whos('QU'); 
a9 = whos('QS'); 
a10 = whos('QV'); 
byte = a1.bytes + a2.bytes + a3.bytes + a4.bytes + a5.bytes + a6.bytes + ...
    a7.bytes + a8.bytes + a9.bytes + a10.bytes; 

%save('byte_trans_low', 'byte'); 
%save('eta_trans_low', 'eta'); 
%save('xi_trans_low', 'xi');


fprintf(fid, 'The whole program takes %4.2f GB of RAM. \n', byte/(1e9));

fprintf(fid, '=============================================================================\n\n');
fclose(fid);