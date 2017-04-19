% test for fast direct solver for

close all;
clear all;
addpath('/home/larry/MATLAB/fmmlib3d-1.2/matlab');  % maxlevel = 3;
%%%% =========== setup of the parameters ==================================
N1 = 50;    % the number of boundary pts on the object in the r-z plane;
M1 = 1.2 * N1;    % the number of MFS pts in the r-z plane;
P = 60;     % number of Fourier Modes;
q = P;
L = 1;       % this is half the periodicity in both directions
R = 3.5;     % the radius of the proxy points;
Rt = 0.5;    % the scaling factor of the object inside of the box
hs = 0.3;
nei = 1;     % 0 or 1
ex = 2*L;  ey = 2*L;  ez = 0;    % periodic unit vector
k = 3;  theta = -pi/4;   phi = pi/3;
% theta and phi are the direction of the incident plane wave traveling in
% polar coordinate system
kx = k * cos(theta) * cos(phi);  ky = k * cos(theta) * sin(phi);  kz = k * sin(theta);
% the x, y, z component of the wavenumber vector k;

alpha = exp(1i*kx*ex);       % the quasi-periodic phase in x direction
beta = exp(1i*ky*ey);        % the quasi-periodic phase in y direction

NP = 50;
RP = 0.75;
NQ = 50;
RQ = 0.75;
acc = 1e-10;
%%%% ======================= END ==========================================

%%%% ================= Get the target and MFS points ======================
a = 0;  w = 0; p = 0;
%a = 0.2; w = 2; p = pi/2;
f = @(s) 1+a*cos(w*(s-p));
f_t = @(s) -a*w*sin(w*(s-p));

target = curve3DQuad(f, f_t, M1, P, Rt);           % curve target points and its normal derivatives;
xt = target.x;  yt = target.y;  zt = target.z;

source = source3DQuad(f, f_t, N1, P, Rt, hs);      % curve MFS points, propotional to local speed;
xs= source.x;  ys = source.y;  zs = source.z;


A = EvalMFSHelm3DKernal (target, xs, ys, zs, k);
b = QuasiPerVecFilling3D(target,kx, ky, kz);
tic;
c = A\b;
toc;



if 0,
    tic;
    %%%%%% real source to target proxy sphere =================================
    [xp, yp, zp] = squaresourcequad(NP, RP);      % box proxy points, points on the big sphere;
    t.x = xp;      t.y = yp;    t.z = zp;
    as = EvalMFSHelm3DKernal (t, xs, ys, zs, k);
    [TS, IS] = id_decomp(as,acc);
    
    rs = size(TS, 1)
    
    V = [eye(rs), TS];
    V(:,IS) =  V(:, :);    % change the indices;
    xsf = xs(IS(1:rs));    ysf = ys(IS(1:rs));  zsf = zs(IS(1:rs));
    
    %===== proxy source to real target=========================================
    [xp, yp, zp] = squaresourcequad(NQ, RQ);
    
    at = EvalMFSHelm3DKernal (target, xp, yp, zp, k);
    at = at';
    [TT, IT] = id_decomp(at,acc);
    
    rt = size(TT,1)
    
    U = [eye(rt), TT];
    U = U';
    
    U(IT,:) = U(:, :);  % change the indices;
    xtf = xt(IT(1:rt));   ytf = yt(IT(1:rt));     ztf = zt(IT(1:rt));
    
    
    
    % ============== Build the K matrix ===========================
    ske.x = xtf;  ske.y = ytf;   ske.z = ztf;
    K = EvalMFSHelm3DKernal (ske, xsf, ysf, zsf, k);
    toc;
end

tic; 
[U, K, V] = getSelfSkeleton(target, source, k, acc, RP, NP, RQ, NQ);
toc; 

norm(U*(K*(V*c)) - b)
norm(A*c-b)