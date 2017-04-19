function a = EvalMFSHelm3DKernalBOR_trans(targetBOR, sourceBOR_p, sourceBOR_m, kp,km, P, Q)
%
% EvalMFSHelm3DKernalBOR_trans(targetBOR, sourceBOR_p, sourceBOR_m, kp,km, P, Q)
% returns the FFTed matrix kernel of the j-th MFS points at i-th location; 
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
% Notes: Transmission problems; 
% 
% By Larry Liu 06/19/2014 

r = targetBOR.r; 
z = targetBOR.z; 
rn = targetBOR.rn; 
zn = targetBOR.zn; 
rsp = sourceBOR_p.r; 
zsp = sourceBOR_p.z; 
rsm = sourceBOR_m.r; 
zsm = sourceBOR_m.z; 

M = length(r);
N = length(rsp);
up = zeros(M, N, P); um = up; 
unp = up;  unm = up; 

for i=1:M
    for j=1:N
        [up(i,j,:),unp(i,j,:)] = ringkernelffthemtrans(r(i),z(i),rn(i), zn(i),rsp(j),zsp(j),kp,P,Q); 
        [um(i,j,:),unm(i,j,:)] = ringkernelffthemtrans(r(i),z(i),rn(i), zn(i),rsm(j),zsm(j),km,P,Q); 
    end
end

a = zeros(2*M, 2*N, P); 
a(1:M, 1:N, :) = up; 
a(1:M, N+1:2*N, :) = -um; 
a(M+1:2*M, 1:N, :) = unp; 
a(M+1:2*M, N+1:2*N, :) = -unm; 

        



