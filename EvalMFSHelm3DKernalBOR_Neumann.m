function a = EvalMFSHelm3DKernalBOR_Neumann(targetBOR, sourceBOR, k, P, Q)
%
% EvalMFSHelm3DKernalBOR_Neumann(targetBOR, sourceBOR, k, P, Q)
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
% Notes: Neumann problem; 
% 
% By Larry Liu 06/19/2014 

r = targetBOR.r; 
z = targetBOR.z; 
rn = targetBOR.rn; 
zn = targetBOR.zn; 
rs = sourceBOR.r; 
zs = sourceBOR.z; 

M = length(r);
N = length(rs);
un = zeros(M, N, P); 

for i=1:M
    for j=1:N
        [~,un(i,j,:)] = ringkernelffthemtrans(r(i),z(i), rn(i), zn(i),rs(j),zs(j),k,P,Q); 
    end
end

a = un; 
        



