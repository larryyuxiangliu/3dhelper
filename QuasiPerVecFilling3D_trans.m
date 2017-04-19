function [b, bn] = QuasiPerVecFilling3D_trans(t, kx, ky, kz) 
% 
% [b, bn] = QuasiPerVecFilling3D_trans(t, kx, ky, kz) returns the RHS vector
% for 3D acoustic scattering over periodic structures; 
% 
% Notes: non-BOR, non-FFT, transmission problem 
%
% Larry Liu, 06/19/2014


[u, un] = pweval(t, kx, ky, kz); 
b = - u;  
bn = -un; 
