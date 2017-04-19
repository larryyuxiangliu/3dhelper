function Periodic3DPlotting(targetBOR, ex, ey, kx, ky, kz, N1,P,PP, N0,NN, xsp, ysp, zsp, k, nei, alpha,beta, eta, xi)
%
% Periodic3DPlotting is to plot the field on a plane z = -3 and x = 3 given
% all the coefficients;
%
% Notes:
% 1. get the solution on the center box and apply phase to get the other
% copies;
%
% 2. for z > 1 or z < -1, first use RB coefficients to get the field at
% -1 < x < 1 and -1 < y < 1 and then apply phase to get the rest;
%
% 3. for z > 1 or z < -1, the field equals the sum of incident wave
% and the RB expansion;
%
% 4. for the center box, the field quals the sum of the proxy points, MFS
% points and the incident wave;
%
% Larry Liu, 07/14/2014


for ii=-1:1
    for jj=-1:1
        showsurffunc(targetBOR.r, targetBOR.z, ii,jj,ex, ey, kx, ky, kz, P);
        hold on;
    end
end
% get the values on the surface of the object;


x0 = -1:0.01:1;
dim = length(x0);
[X, Y] = meshgrid(x0, x0);
x_test = reshape(X, 1, dim^2);
y_test = reshape(Y, 1, dim^2);
z_test = -3 * ones(1, dim^2);
Z = -3 * ones(dim, dim);
u_test = exp(1i * (kx * x_test + ky * y_test + kz * z_test));

an = xi((PP+1)^2 + 1: (PP+1)^2 + NN);
bn = xi((PP+1)^2 + NN + 1 : end);
% get the RB coefficients an and bn;
ns = -2*N0:2*N0;   % to make sure that the circle with radius kappa contains all the grid points;
kappa = pi*N0;
ii = 1;
for m=ns,     % index of freq n in the list from 1:nd
    for n=ns,
        kappanx = kx + 2*pi*m/ex;   % NOT TRUE FOR GENERAL e DIRECTION - fix me
        kappany = ky + 2*pi*n/ey;
        if sqrt(kappanx^2 + kappany^2) <= kappa
            kn = sqrt(k^2 - kappanx^2 - kappany^2);   % vertical cmpnt of wavevector of this RB mode
            u_test = u_test + bn(ii)* exp(1i * (kappanx * x_test  +  kappany * y_test -  kn * (z_test+1)));
            
            ii = ii + 1;
        end
        
    end
end
U = reshape(u_test, dim, dim);
surf(X, Y, Z, real(U));
axis equal; shading interp; xlabel('x'); ylabel('y'); zlabel('z');
hold on;
U_b = U * beta;
surf(X, Y+ey, Z, real(U_b));
axis equal; shading interp; xlabel('x'); ylabel('y'); zlabel('z');
hold on;
U_f = U * beta^(-1);
surf(X, Y-ey, Z, real(U_f));
axis equal; shading interp; xlabel('x'); ylabel('y'); zlabel('z');
hold on;
U_l = U * alpha^(-1);
surf(X-ex, Y, Z, real(U_l));
axis equal; shading interp; xlabel('x'); ylabel('y'); zlabel('z');
hold on;
U_r = U * alpha;
surf(X+ex, Y, Z, real(U_r));
axis equal; shading interp; xlabel('x'); ylabel('y'); zlabel('z');
hold on;
U_lb = U * beta * alpha^(-1);
surf(X-ex, Y+ey, Z, real(U_lb));
axis equal; shading interp; xlabel('x'); ylabel('y'); zlabel('z');
hold on;
U_rb = U * beta * alpha;
surf(X+ex, Y+ey, Z, real(U_rb));
axis equal; shading interp; xlabel('x'); ylabel('y'); zlabel('z');
hold on;
U_lf = U * beta^(-1)* alpha^(-1);
surf(X-ex, Y-ey, Z, real(U_lf));
axis equal; shading interp; xlabel('x'); ylabel('y'); zlabel('z');
hold on;
U_rf = U * beta^(-1)* alpha;
surf(X+ex, Y-ey, Z, real(U_rf));
axis equal; shading interp; xlabel('x'); ylabel('y'); zlabel('z');
hold on;


y0 = -1:0.01:1;
dim = length(y0);
[Y, Z] = meshgrid(y0, y0);
y_test = reshape(Y, 1, dim^2);
z_test = reshape(Z, 1, dim^2);
x_test = 1 * ones(1, dim^2);
X = 3 * ones(dim, dim);
u_test2 = zeros(1, dim^2);
%u_test2 = fmm3DHelpoteval(xp, yp, zp, xi(1:N2^2), kp, x_test, y_test, z_test);
rscale = 1;
center = [0; 0; 0];
nterms = PP;
clm = xi(1:(PP+1)^2);
ztarg = [x_test; y_test; z_test];
znor = [x_test; zeros(1, dim^2); zeros(1,dim^2)];
u_test2 = h3dtaeval(k,rscale,center,nterms,clm,ztarg,znor).';
u_test2 = u_test2 + fmm3DPerHelpoteval(xsp, ysp, zsp, eta(1:N1*P), k, x_test, y_test, z_test, nei, alpha, beta, ex, ey);
u_test2 = u_test2 + exp(1i * (kx * x_test + ky * y_test + kz * z_test));
%u_test3 = fmm3DHelpoteval(xp, yp, zp, xi(1:N2^2), kp, x_test, y_test, z_test);
U1 = reshape(u_test2, dim, dim);
U1 = U1 * alpha;         % the value at the plane x = 3;
surf(X, Y, Z, real(U1));
axis equal; shading interp; xlabel('x'); ylabel('y'); zlabel('z');
hold on;
U = U1 * beta;
surf(X, Y+ey, Z, real(U));
axis equal; shading interp; xlabel('x'); ylabel('y'); zlabel('z');
hold on;
U = U1 * beta^(-1);
surf(X, Y-ey, Z, real(U));
axis equal; shading interp; xlabel('x'); ylabel('y'); zlabel('z');
hold on;


u_test2 = zeros(1,dim^2);
z_test = z_test + 2;
u_test2 = exp(1i * (kx * x_test + ky * y_test + kz * z_test)) ;
ns = -2*N0:2*N0;   % to make sure that the circle with radius kappa contains all the grid points;
kappa = pi*N0;
ii = 1;
for m=ns,    % index of freq n in the list from 1:nd
    for n=ns,
        kappanx = kx + 2*pi*m/ex;   % NOT TRUE FOR GENERAL e DIRECTION - fix me
        kappany = ky + 2*pi*n/ey;
        if sqrt(kappanx^2 + kappany^2) <= kappa
            kn = sqrt(k^2 - kappanx^2 - kappany^2);   % vertical cmpnt of wavevector of this RB mode
            
            
            u_test2 = u_test2 + an(ii)* exp(1i * ((kappanx * x_test ) +  kappany * y_test + kn * (z_test-1)));
            
            ii = ii + 1;
        end
    end
end


U1 = reshape(u_test2, dim, dim);

U1 = U1 * alpha;
surf(X, Y, Z+2, real(U1));
axis equal; shading interp; xlabel('x'); ylabel('y'); zlabel('z');
hold on;
U = U1 * beta;
surf(X, Y+ey, Z+2, real(U));
axis equal; shading interp; xlabel('x'); ylabel('y'); zlabel('z');
hold on;
U = U1 * beta^(-1);
surf(X, Y-ey, Z+2, real(U));
axis equal; shading interp; xlabel('x'); ylabel('y'); zlabel('z');
hold on;

u_test2 = zeros(1,dim^2);
z_test = z_test - 4;
u_test2 = exp(1i * (kx * x_test + ky * y_test + kz * z_test));
ns = -2*N0:2*N0;   % to make sure that the circle with radius kappa contains all the grid points;
kappa = pi*N0;
ii = 1;
for m=ns,    % index of freq n in the list from 1:nd
    for n=ns,
        kappanx = kx + 2*pi*m/ex;   % NOT TRUE FOR GENERAL e DIRECTION - fix me
        kappany = ky + 2*pi*n/ey;
        if sqrt(kappanx^2 + kappany^2) <= kappa
            kn = sqrt(k^2 - kappanx^2 - kappany^2);   % vertical cmpnt of wavevector of this RB mode
            
            
            u_test2 = u_test2 + bn(ii)* exp(1i * ((kappanx * x_test ) +  kappany * y_test - kn * (z_test+1)));
            
            ii = ii + 1;
        end
    end
end


U1 = reshape(u_test2, dim, dim);

U1 = U1 * alpha;
surf(X, Y, Z-2, real(U1));
axis equal; shading interp; xlabel('x'); ylabel('y'); zlabel('z');
hold on;
U = U1 * beta;
surf(X, Y+ey, Z-2, real(U));
axis equal; shading interp; xlabel('x'); ylabel('y'); zlabel('z');
hold on;
U = U1 * beta^(-1);
surf(X, Y-ey, Z-2, real(U));
axis equal; shading interp; xlabel('x'); ylabel('y'); zlabel('z');
hold on;


