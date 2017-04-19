function W = RBW3D(t, bo, M2, N0, NN, k, kx, ky, ex, ey)
% function W = RBW3D(t,bo, M2, N0, k, kx, ky, ex, ey) returns the matrix W
% of matrix Q for 3D acoustic scattering from periodic structure;
%
% By Larry Liu, 05/18/2014

M = M2^2;
ms = -2*N0:2*N0;   % to make sure that the circle with radius kappa contains all the grid points; 
ii = 1;
kappa = pi*N0; 
W = zeros(4*M, 2*NN); 
for m=ms,    % index of freq n in the list from 1:nd
    for n=ms,       
        kappanx = kx + 2*pi*m/ex;   % NOT TRUE FOR GENERAL e DIRECTION - fix me
        kappany = ky + 2*pi*n/ey; 
        if sqrt(kappanx^2 + kappany^2) <= kappa 
            kn = sqrt(k^2 - kappanx^2 - kappany^2);   % vertical cmpnt of wavevector of this RB mode
            W(1:M, ii) = exp( 1i * (kappanx * t.x + kappany * t.y)); % RB rep has exp(i*kn*(y-y0))
            W(M+1:2*M, ii) = 1i * kn * W(1:M, ii);
            W(2*M+1:3*M, NN + ii) = exp( 1i * (kappanx * bo.x + kappany * bo.y));
            W(3*M+1:4*M,  NN + ii) = -1i * kn * W(2*M+1:3*M, NN + ii);
            ii = ii + 1;
        end
    end
end

