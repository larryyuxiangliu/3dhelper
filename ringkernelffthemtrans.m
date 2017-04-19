function [u,un]=ringkernelffthemtrans(r,z,rn,zn,rs,zs,k,P,Q)
% RINGKERNELFFTHEMTRANS - output ring kernel for all Fourier modes
%
% u=ringkernelfft(r,z,rn,zn,rs,zs,k,P,Q) returns vector u and un
% of values
% I:
% O:
%Q: the number of trp rule points.
%P: is the number of Fourier Modes.
%k: the frequency of the source.
%(r,z) are the boundary points.
%(rs,zs) are the source points.
%(rn,zn) are the normal unit vector of boundary points.

% This is to get the Fourier coefficients of data on the boundary point
% produced by a ring charge. 

% Also see: RINGKERNELFFTHEM



j = 0:Q-1;
rd = sqrt((zs-z)^2 + (rs*cos(2*pi*j/Q)-r).^2+(rs*sin(2*pi*j/Q)).^2);
%rd = sqrt((zs-z)^2 + (r*cos(2*pi*j/Q)-rs).^2+(r*sin(2*pi*j/Q)).^2);
% rd is the distance between the boundary points to the charge rings. 

f = (2*pi/Q).*exp(1i*k*rd)./rd;
un_f=fft(f);


nx=rn;
ny=0;
nz=zn;

Rx= rs*cos(2*pi*j/Q)-r;
Ry= rs*sin(2*pi*j/Q);
Rz= zs-z;
cosalpha= (nx.*Rx+ny.*Ry+nz*Rz)./rd;

g = (2*pi/Q).*cosalpha.*(1-1i*k*rd).*exp(1i*k*rd)./(rd.^2);
un_g = fft(g);

if Q>=P, 
	u = [un_f(1:P/2), un_f(end-P/2+1:end)];
	un = [un_g(1:P/2), un_g(end-P/2+1:end)];
	
else 
	u = zeros(1, P); 
	
	u(:,1:Q/2) = un_f(:, 1:Q/2);
	
	u(:,end-Q/2+1:end) = un_f(:,end-Q/2+1:end);   
	
	un = zeros(1, P);
	
	un(:,1:Q/2) = un_g(:, 1:Q/2);
	
	un(:,end-Q/2+1:end) = un_g(:,end-Q/2+1:end);   
		
end
