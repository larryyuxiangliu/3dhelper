function [u, un] = EvalMFSHelm3DKernal (t, xs, ys, zs, k)
%
% EvalMFSHelm3DKernal (t, xs, ys, zs, k) returns the matrix kernel of the j-th
% MFS points at i-th location;
% Inputs: t, 1x1 structure:  t.x, t.y, t.z, t.nx, t.ny, t.nz;
%
% Outputs: u, M x N matrix:
% By Larry Liu 05/16/2014

M = length(t.x);
N = length(xs);

if 1,
    [XS, XT] = meshgrid (xs, t.x);
    Rx = XS - XT;
    clear XS;
    clear XT;
    
    [YS, YT] = meshgrid (ys, t.y);
    Ry = YS - YT;
    clear YS;
    clear YT:
    
    [ZS, ZT] = meshgrid (zs, t.z);
    Rz = ZS - ZT;
    clear ZS;
    clear ZT;
    
    R = sqrt(Rx.^2+ Ry.^2+ Rz.^2);
    u = exp(1i*k*R)./R;
    
    if nargout >1
        un = zeros(M, N);
        [XS, NX] = meshgrid (xs, t.nx);
        clear XS;
        [YS, NY] = meshgrid (ys, t.ny);
        clear YS;
        [ZS, NZ] = meshgrid (zs, t.nz);
        clear ZS;
        cosphi=(NX.*Rx + NY.*Ry + NZ.* Rz)./R;
        un = cosphi.* (1 - 1i* k * R) .* exp(1i * k * R)./(R.^2);
    end
    
    
else
    u = zeros(M,N);  un = u;
    for j = 1:N
        Rx = xs(j) - t.x;
        Ry = ys(j) - t.y;
        Rz = zs(j) - t.z;
        R = sqrt (Rx.^2 + Ry.^2 + Rz.^2);
        u(:,j) = exp(1i*k*R)./R;
        if nargout >1
            cosphi = (t.nx.*Rx + t.ny.*Ry + t.nz .*Rz)./R;
            un(:,j) = cosphi.* (1 - 1i* k * R) .* exp(1i * k * R)./(R.^2);
            
        end
    end
    
end


