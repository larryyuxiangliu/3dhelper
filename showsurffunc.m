function showsurffunc(r,h,ii,jj,ex, ey, kx, ky, kz, P)
% plot a 3d surf with func on it in color.
% explain inputs.

M = numel(r);

phi=(0:(P-1))/P * 2* pi;

x = zeros(M, P);  y = x;  z = x;  f = x;  
for i=1:M
    for j=1:P
        x(i,j)=r(i)*cos(phi(j)) + ii * ex;
        y(i,j)=r(i)*sin(phi(j)) + jj * ey;
        z(i,j)=h(i);  
        f(i,j)= -exp (1i*(x(i,j) * kx + y(i,j) * ky + z(i,j) * kz));
    end
end

surf(x,y,z, real(f));
axis equal; shading interp; xlabel('x'); ylabel('y'); zlabel('z');
%lightangle(-45, 0); %caxis(2*[-1 1]); 
colorbar

end

