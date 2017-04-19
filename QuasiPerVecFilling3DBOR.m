function f = QuasiPerVecFilling3DBOR(targetBOR, kx, ky, kz, P) 
% 
% f = QuasiPerVecFilling3DBOR(targetBOR, kx, ky, kz, P) returns the RHS
% of the linear system; 
% 
% f needs to be "FFTed" in the main scrpit 
% 
% Inputs: 
% 
%       targetBOR : the structure which contains the r and z coordinates of
%       the boundary points in the r-z plane. 
%       kx, ky, kz: the wavenumber vector 
%       P: the number of Fourier modes 
% 
% Outputs:
%
%       f (M X P): the value for each boundary point in different phi
%       values
% 
% Notes: This only applies to the BOR case, i.e., when we use FFT 
% 
% Larry Liu, 06/11/2014

r = targetBOR.r; 
z = targetBOR.z; 

M = length(r); 
f = zeros(M,P);   xx = f;  yy = f;  zz = f; 
for n=1:P
    for i=1:M
        xx(i,n) = r(i) * cos((n-1)/P*2*pi);
        yy(i,n) = r(i) * sin((n-1)/P*2*pi);
        zz(i,n) = z(i);         
    end
end

f= -exp (1i*(xx * kx + yy * ky + zz * kz)); 
