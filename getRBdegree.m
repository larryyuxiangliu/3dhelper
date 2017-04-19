function NN = getRBdegree(N0, kx, ky, ex, ey)
% 
% getRBdegree - Get the number of RB coefficients at top/bottom given the
% wavenumber kx, ky, periodic unit length ex, ey and the parameter N0. N0
% is controlled by the incident wavenumber k. 
%
% getRBdegree(N0, kx, ky, ex, ey) 
% 
% Larry Liu, 08/03/2014

ii = 0;   % counter to count the number of unkowns for RM expansion. 
ms = -2*N0:2*N0;  % to make sure that the circle with radius kappa contains all the grid points; 
kappa = pi * N0; 
for m=ms,    % index of freq n in the list from 1:nd
    for n=ms,       
        kappanx = kx + 2*pi*m/ex;   % NOT TRUE FOR GENERAL e DIRECTION - fix me
        kappany = ky + 2*pi*n/ey; 
        if sqrt(kappanx^2 + kappany^2) <= kappa 
            ii = ii+1; 
        end
    end
end

NN = ii; 