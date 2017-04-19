function un=ringkernelffthem(r,z,rs,zs,k,P,Q)
% RINGKERNELFFTHEM - output ring kernel for all Fourier modes
%
% u=ringkernelfft(r,z,rs,zs,k,P,Q) returns vector the ring kernel for all
% Fourier modes. 
%
% Inputs:    
%            P: the number of Fourier Modes.
%            Q: the number of trp rule points.
%            k: the wavenumer.
%            (r,z): location of boundary pt in the r-z plane. 
%            (rs, zs): location of MFS pt in the r-z plane. 
%
% Outputs:   
%            un(1 X P): the ring kernel for all Fourier modes for the pair
%            of pts (r, z) and (rs, zs). 
% 
% This is to creat a 3D ringkerl for BOR technique FFT...
% Also see: RINGKERNEL RINGKERNELFFT RINGKERNELFFTTEST
% By Alex, Spring 2013


j = 0:Q-1;
rk = sqrt((zs-z)^2 + (rs*cos(2*pi*j/Q)-r).^2+(rs*sin(2*pi*j/Q)).^2); 
u = (2*pi/Q)*exp(1i*k*rk)./rk;

un_q=fft(u);
%un_q is the frequency domain of the data . 


if Q>=P, 
	un = [un_q(1:P/2), un_q(end-P/2+1:end)];
%un = reshape(un, 1, P); 
%un is the Fourier Coefficients of the data u with P Fourier Modes.
	
else 
	un = zeros(1, P); 
	
	un(:,1:Q/2) = un_q(:, 1:Q/2);
	
	un(:,end-Q/2+1:end) = un_q(:,end-Q/2+1:end);   
	
		
end
