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

curpath = pwd;
addpath(sprintf('%s/fmmlib3d-1.2/matlab',curpath))
% clear
clear
%addpath('/home/larry/MATLAB/fmm');  % the old default one;
addpath('/home/larry/MATLAB/fmmlib3d-1.2/matlab');  % maxlevel = 3;
tic;
filename = 'neumann_per_results_maxlevel3.txt';
fid = fopen(filename, 'a');
type = 'n';
%%%% =========== setup of the parameters ==================================
N1 = 360;    % the number of boundary pts on the object in the r-z plane;
M1 = 1.2 * N1;    % the number of MFS pts in the r-z plane;
P = 180;     % number of Fourier Modes;
q = 600;
PP = 86;
N0 = 25;     % N0 is related to the number of the RB expansion coefficients
M2 = ceil(sqrt(pi) * N0);
L = 1;       % this is half the periodicity in both directions
R = 3.5;     % the radius of the proxy points;
Rt = 0.5;    % the scaling factor of the object inside of the box
h = 0.08;
nei = 1;     % 0 or 1
ex = 2*L;  ey = 2*L;  ez = 0;    % periodic unit vector
k = 40;  theta = -pi/4;   phi = pi/3;
N_iter = 100;
fprintf(fid, '=============================================================================\n');
fprintf(fid, 'This is to test the case when q = P for far field evalution\n');
fprintf(fid, 'The parameters are M1 = %d, N1 = %d, P = %d, q = %d, M2 = %d, PP = %d, N0 = %d, Rt = %2.2f, h = %2.2f, k = %d, N_iter = %d \n', ...
    M1, N1, P, q, M2, PP, N0, Rt, h, k, N_iter);

% theta and phi are the direction of the incident plane wave traveling in
% polar coordinate system
kx = k * cos(theta) * cos(phi);  ky = k * cos(theta) * sin(phi);  kz = k * sin(theta);
% the x, y, z component of the wavenumber vector k;
NN = getRBdegree(N0, kx, ky, ex, ey);     % NN is the number of coefficients for RB coefficients at top / bottom, so totally 2 * NN.

alpha = exp(1i*kx*ex);       % the quasi-periodic phase in x direction
beta = exp(1i*ky*ey);        % the quasi-periodic phase in y direction
%%%% ======================= END ==========================================

%%%% ================= Get the target and MFS points ======================
aa = 0.2;  % thickness: cannot be >1/3 otherwise not smooth to emach
bb = pi/6;  % controls approx opening angle in radians (keep small for resonant)

[TT,BO,LL,RR,FF,BA] = squarequad(M2, L);            % box target points;

bdyBOR = cupBdyBOR(aa, bb, M1, Rt);
srcBOR = srcBdyBOR(aa, bb, N1, h, Rt);
src3D = srcBdy3D(aa, bb, N1, h, P, Rt);
bdy3D = cupBdy3D(aa, bb, M1, P, Rt);

xs = src3D.x;  ys = src3D.y;  zs = src3D.z;
xt = bdy3D.x;  yt = bdy3D.y;  zt = bdy3D.z;


% ========================== END ==========================================

% ============================= to fill the matrices ======================
A_0 = EvalMFSHelm3DKernalBOR_Neumann(bdyBOR, srcBOR, k, P, q);

[~, Bn] = QuasiPerB3D(bdy3D, PP, NN, k);
B = Bn;

C = QuasiPerC3D(TT,BO,LL,RR,FF,BA, xs, ys, zs, k, ex, ey, nei, alpha, beta);

Q = QuasiPerQ3D(TT,BO,LL,RR,FF,BA, M2, PP, N0, NN, k, kx, ky, ex,ey, alpha, beta);

[~, bb] = QuasiPerVecFilling3D_trans(bdy3D, kx, ky, kz) ;

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
    ff = @(x)  apply_GMRES_Matrix_BOR_SVD_test(bdy3D, AU,iA,AV,QU,iQ,QV,...
        B, C, M1, P, P,  xs, ys, zs, x, k, ...
        nei, alpha, beta, ex, ey, type);
    
    [xx, ~, ~, it, res] = gmres(ff, bb, [], 1e-11, N_iter);
    
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
end
fprintf(fid, 'GMRES iteration takes %4.4f seconds and has %d iterations. \n', toc, it(2));
% ============================= END =======================================
if 1,
    tic;
    M_test = 120;   P_test = 120;
    aa = 0.2;  % thickness: cannot be >1/3 otherwise not smooth to emach
    bb = pi/6;  % controls approx opening angle in radians (keep small for resonant)
    
    error = getBdyError(M_test, P_test, aa, bb,Rt,h, [],xs,ys,zs, k,0, kx, ky, kz, N1, P, q,PP,eta, xi, nei ,alpha, beta, ex, ey, type);
    
    
    Mtest = 40;
    [err1, err2] = getPerError (M_test, L,  xs, ys, zs, eta, xi, PP, k, nei ,alpha, beta, ex, ey);
    flux = getFlux(xi, PP, N0,NN, k, kx, ky, kz, ex, ey);
    err_f = abs(flux - abs(kz));
    
    toc;
end
fprintf(fid, 'Doing error evaluation take %4.4f seconds \n', toc);
fprintf(fid, 'boundary matching error is %e, flux matching error is %e \n', error, err_f);
fprintf(fid, 'x-direction periodicity matching error is %e, y-direction periodicity matching error is %e.\n', err1, err2);
fprintf(fid, 'Periodicity matching condition is %e.\n', (err1 + err2)/2 );
fprintf(fid, '=============================================================================\n\n');
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

%save('byte_neu', 'byte');
%save('eta_neu', 'eta');
%save('xi_neu', 'xi');
fprintf(fid, 'The whole program takes %4.2f GB of RAM. \n', byte/(1e9));

fprintf(fid, '=============================================================================\n\n');
fclose(fid);



