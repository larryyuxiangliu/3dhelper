function Q = QuasiPerQ3D_impedance(t, bo, l, r, f, ba,  M,  P, N0, NN, k,kx, ky, ex,ey, alpha, beta)
%
% function Q =  QuasiPerQ3D(t, bo, l, r, f, ba,  M,  P,NN, k,kx, ky, ex,ey, alpha, beta) returns the matrix Q;
%
% By Larry Liu, 05/18/2014


m = M^2; % number of points on each surface of the square;
N = (P+1)^2;  % number of target points;
%%%% ========== This is the Z in the formula ==========================
rscale = 1;
center = [0;0;0];
nterms = P;
ztarg = [l.x.'; l.y.'; l.z.'];
znor = [l.nx.'; l.ny.'; l.nz.'];
[Al, Anl] = h3dtamatrix(k,rscale,center,nterms,ztarg,znor);

ztarg = [r.x.'; r.y.'; r.z.'];
znor = [r.nx.'; r.ny.'; r.nz.'];
[Ar, Anr] = h3dtamatrix(k,rscale,center,nterms,ztarg,znor);

ztarg = [f.x.'; f.y.'; f.z.'];
znor = [f.nx.'; f.ny.'; f.nz.'];
[Af, Anf] = h3dtamatrix(k,rscale,center,nterms,ztarg,znor);

ztarg = [ba.x.'; ba.y.'; ba.z.'];
znor = [ba.nx.'; ba.ny.'; ba.nz.'];
[Aba, Anba] = h3dtamatrix(k,rscale,center,nterms,ztarg,znor);


ztarg = [t.x.'; t.y.'; t.z.'];
znor = [t.nx.'; t.ny.'; t.nz.'];
[At, Ant] = h3dtamatrix(k,rscale,center,nterms,ztarg,znor);


ztarg = [bo.x.'; bo.y.'; bo.z.'];
znor = [bo.nx.'; bo.ny.'; bo.nz.'];
[Abo, Anbo] = h3dtamatrix(k,rscale,center,nterms,ztarg,znor);
Q = zeros(7*m, N + NN);


Q(1:m,1:N) = Ar - alpha * Al;
Q(m+1:2*m, 1:N ) = Anr - alpha * Anl;
Q(2*m+1:3*m,1:N) = Aba - beta * Af;
Q(3*m+1:4*m, 1:N ) = Anba - beta * Anf;
Q(1:4*m, N+1:N + NN) = 0;
%%%========================================================================

%%%%======== This is the second block in the formula, the 0 matrix ========

% ============= This is to build up the V and W matrix in the formula =====


Q(4*m+1:5*m, 1:N) = At;
Q(5*m+1:6*m, 1:N) = Ant;

lambda = 1; 
Q(6*m+1:7*m, 1:N) = 1i * lambda * Abo + Anbo;
W = RBW3D_impedance(t, bo, M, N0, NN, k, kx, ky, ex, ey);
Q(4*m+1:6*m, N+1 :N + NN) = -W;

%==========================================================================

