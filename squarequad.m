function [t,bo,l,r,f,ba] = squarequad(M, L)
% [t,bo,l,r,f,ba] = squarequad(M, L) returns 6 structures, each structure
% means a surface, top, bottom, left, right, front and back. each strueture
% contains x, y and z coordiante as well as the normal derivatives of the 
% surface points in a 3D cartesian coordinate system; 
% 
% M is the number of points along each direction on each surface, so there
% are totally 6 * M^2 points; 
% L is the half of the length of the square;  
%
% Notes: To make sure that the points on the top and bottom surface are
% uniformly distributed, while the points on the other 4 surfaces are
% Guassian quadrature points; 
% 
% By Larry Liu, May 16, 2014 

x0 = L * gauss(M);
%x0 = linspace (-L, L, M); 
[X,Y] = meshgrid(x0,x0); 
x1 = 2*((1:M)/M) - 1; 
x1 = L * x1; 
x1 = reshape(x1, M, 1); 
[X1, Y1] = meshgrid(x1, x1); 
x1 = reshape(X1, M^2, 1); 
y1 = reshape(Y1, M^2, 1); 
x0 = reshape(X, M^2,1); 
y0 = reshape(Y, M^2,1); 
LL = L*ones(M^2,1); 
m = M^2; 

%t.x = x0;  t.y = y0;  t.z = LL ;   
%t.nx = zeros(m,1); t.ny = zeros(m,1); t.nz = ones(m,1); 
%bo.x = x0; bo.y = y0; bo.z = -LL; 
%bo.nx =zeros(m,1); bo.ny = zeros(m,1); bo.nz = ones(m,1); 


t.x = x1;  t.y = y1;  t.z = LL ;   
t.nx = zeros(m,1); t.ny = zeros(m,1); t.nz = ones(m,1); 
bo.x = x1; bo.y = y1; bo.z = -LL; 
bo.nx =zeros(m,1); bo.ny = zeros(m,1); bo.nz = ones(m,1); 

l.x = -LL;  l.y = x0;  l.z = y0;   
l.nx = ones(m,1); l.ny = zeros(m,1); l.nz = zeros(m,1); 
r.x = LL; r.y = x0; r.z = y0;      
r.nx = ones(m,1); r.ny = zeros(m,1); r.nz = zeros(m,1); 

f.x = x0; f.y = -LL; f.z = y0;     
f.nx = zeros(m,1); f.ny = ones(m,1); f.nz = zeros(m,1); 
ba.x = x0;  ba.y = LL; ba.z = y0;  
ba.nx = zeros(m,1); ba.ny = ones(m,1); ba.nz = zeros(m,1); 



