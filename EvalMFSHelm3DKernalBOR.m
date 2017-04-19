function u = EvalMFSHelm3DKernalBOR(targetBOR, rs,zs, k, P, Q)
%
% EvalMFSHelm3DKernalBOR (t, xs, ys, zs, k, P, NQ) returns the FFTed matrix kernel of the j-th
% MFS points at i-th location; 
% 
% Inputs: 
%       
%         targetBOR: the strucutre that contains the locations of the
%                    boundary points in the r-z plane 
%         rs, zs:    the locations of the MFS source pts in the r-z plane
%         k, P, Q:   as usual, the wavenumber, the number of Fourier Modes
%                    and the number of Trap rule quadature pts. 
% 
% Outputs: u, M x N x P matrix, contains the kernel matrix for different
% Fourier modes. 
% 
% Also see EvalMFSHel3DKernal.m
% 
% By Larry Liu 06/14/2014 

r = targetBOR.r; 
z = targetBOR.z; 

M = length(r);
N = length(rs);
u = zeros(M, N, P);

for i=1:M
    for j=1:N
        u(i,j,:) = ringkernelffthem(r(i),z(i),rs(j),zs(j),k,P,Q); 
    end
end


        



