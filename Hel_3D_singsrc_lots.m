% This code deals with 3D acoustic scattering problem (Helmholtz equation)
% from periodic strucutres where the periodicity only appears in 2D.
% We need to impose Rayleigh-Bloch Expansion on the top and bottom and also
% quasi-periodic condition on the left and right as well as front and back
% of the box;
%
% Notes: In this case, the wave is coming from a singular source inside of
% the box;
%
% Notes: This is a test code to see how to pick the right periodic
% paramters; there is only one single source in each box;
%
% By Larry Liu 05/16/2014

close all;
clear all;
addpath('/home/larry/MATLAB/fmm');
filename = 'singsrc_per.txt';
fid = fopen(filename, 'a');
%%%% =========== setup of the parameters ==================================
fprintf(fid, '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
fprintf(fid, 'Now it is the new test for new way of choosing M2 based on N0');

PPs = 70;
N0s = 21:2:25;

ii = 1;


for PP = PPs   % degree of spherical harmonics, so totally (PP+1)^2 coefficients;
    for N0 = N0s
        
        tic;
        M = ceil(sqrt(pi) * N0);
        N = 0;
        
        flag = 1;    % computational flag, flag = 1 means to use spherical harmonics which is faster; flag = 0 means to use proxy points;
        L = 1;       % this is half the periodicity in both directions
        R = 3.5;     % the radius of the proxy points;
        ex = 2*L;  ey = 2*L;  ez = 0;    % periodic unit vector
        nei = 1;
        k = 30;
        theta = -pi/4;  phi = pi/3;
        kx = k * cos(theta) * cos(phi);  ky = k * sin(theta) * sin(phi);
        kz = k * sin(theta);
        %kx = k;  ky = k;
        NN = getRBdegree(N0, kx, ky, ex, ey);
        alpha = exp(1i*kx*ex);       % the quasi-periodic phase in x direction
        beta = exp(1i*ky*ey);        % the quasi-periodic phase in y direction
        
        x_sig = 0.91;  y_sig = 0.91;  z_sig = 0.91;
        fprintf(fid, 'The paramaters are M = %d, PP = %d, N0 = %d, k = %d, point is at (%2.2f, %2.2f, %2.2f) \n', M, PP, N0, k, x_sig, y_sig, z_sig);
        %%%% ======================= END ==========================================
        
        %%%% ================= Get the target and MFS points ======================
        [TT,BO,LL,RR,FF,BA] = squarequad(M, L);   % box target points;
        
        
        %%% ============================Get the matrix and vector =================
        if 1,
            
            Q = QuasiPerQ3D(TT,BO,LL,RR,FF,BA,M, PP, N0, NN, k, kx, ky, ex, ey, alpha, beta);
            
            b = -QuasiPerC3D(TT,BO,LL,RR,FF,BA, x_sig,y_sig,z_sig,  k, ex, ey, nei, alpha, beta);
            %[QU, QS, QV] = svd(Q, 0);
            %ss = diag(QS);
            %nQ = length(ss);
            %rQ = length(find (ss >1e-10));      % choose a cut, 1e-10.
            %iQ =  [1./ss(1:rQ); zeros(nQ-rQ,1)];
            %c = QV * (iQ .* (QU' * b));
            c = Q\b;
            norm(Q*c-b);
        end
        
        
        if 1,
            
            M_test = 40;
            
            
            [err1, err2] = getPerError (M_test, L,  x_sig, y_sig, z_sig, 1, c, PP, k, nei, alpha, beta, ex, ey);
            
            error = (err1 + err2)/2;
            fprintf(fid, 'x-direction periodicity matching is %e, y-direction periodicity matching is %e \n', err1, err2 );
            % results(ii,:) = strcat(num2str(M, '%02d'), '&', num2str(PP, '%02d'), '&', num2str(N0, '%02d'), '&', num2str(8*M^2, '%05d'),'-by-', num2str( (PP+1)^2+2*(2*N0+1)^2, '%05d') , '&', num2str(error, '%10.5e\n'), '\\ \hline');
            toc;
            fprintf(fid, 'the time for this test is %2.2f \n', toc);
            
        end
        ii = ii+1;
    end
end

fclose(fid);
