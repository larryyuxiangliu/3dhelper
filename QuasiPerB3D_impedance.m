function [B, Bn] = QuasiPerB3D_impedance(bdy3D, PP, NN, k)
%
% function [B, Bn] = QuasiPerB3D(bdy3D, PP, NN, k)
%
% Larry Liu, 10/14/2014 


M = length(bdy3D.x);

N = (PP+1)^2;  % number of degree for spherical harmonics;
B = zeros(M, N +  NN) ;
Bn = zeros(M, N + NN);

%%%% ========== This is the Z in the formula ==========================
rscale = 1;
center = [0;0;0];
nterms = PP;
ztarg = [bdy3D.x.'; bdy3D.y.'; bdy3D.z.'];
znor = [bdy3D.nx.'; bdy3D.ny.'; bdy3D.nz.'];
[B(:,1:N), Bn(:,1:N)] = h3dtamatrix(k,rscale,center,nterms,ztarg,znor);



