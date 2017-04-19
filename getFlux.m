function flux = getFlux(xi, PP, N0, NN, k, kx, ky, kz, ex, ey)
%
% flux = getFlux(N2, N0, k, kx, ky, kz, ex, ey,xi) returns the flux for the
% scattering problem; for details, see notes Hel_3D_Periodic_14X/notes.pdf
%
% Inputs: N2, number of proxy points in each direction, so totally N2^2 pts
%         N0, degree of RB expansion in each direction in each side, so
%         totally 2 * (2 * N0+1)^2
%
%         k,kx,ky,kz incident wavenumber and its components;
%         ex, ey the periodic vector;
%
%         xi, the coefficients xi, made up from proxy points coefficients and
%             RB coefficients;
%
% Outputs: flux, the flux from the periodic objects;
%
% Larry Liu, 07/14/2014




an = xi((PP+1)^2 + 1: (PP+1)^2 + NN);
bn = xi((PP+1)^2 + NN + 1 : end);
% get the RB coefficients an and bn;
ns = -2*N0:2*N0;   % to make sure that the circle with radius kappa contains all the grid points; 
kappa = pi*N0;
ii = 1;
flux_t = 0;   % the flux from the top
for m=ns,     % index of freq n in the list from 1:nd
    for n=ns,
        kappanx = kx + 2*pi*m/ex;   % NOT TRUE FOR GENERAL e DIRECTION - fix me
        kappany = ky + 2*pi*n/ey;
        if sqrt(kappanx^2 + kappany^2) <= kappa
            kn = sqrt(k^2 - kappanx^2 - kappany^2);   % vertical cmpnt of wavevector of this RB mode
            
            flux_t = flux_t + kn * abs(an(ii))^2;  % see notes for details;
            
            ii = ii + 1;
        end
    end
end


ii = 1;
flux_b = 0;  % flux from the bottom;
for m=ns,    % index of freq n in the list from 1:nd
    for n=ns,
        kappanx = kx + 2*pi*m/ex;   % NOT TRUE FOR GENERAL e DIRECTION - fix me
        kappany = ky + 2*pi*n/ey;
        if sqrt(kappanx^2 + kappany^2) <= kappa
            kn = sqrt(k^2 - kappanx^2 - kappany^2);   % vertical cmpnt of wavevector of this RB mode
            if (m==0) && (n==0)
                flux_b = flux_b + kn * abs( bn(ii) + exp(-1i*kz*1))^2;
            else
                flux_b = flux_b + kn *abs( bn(ii))^2;
            end
            
            ii = ii + 1;
        end
    end
end

flux = flux_t + flux_b;
flux = real(flux);

