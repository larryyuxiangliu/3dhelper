
addpath('/home/larry/MATLAB/fmm');
tic;
filename = 'trans_nonper_results.txt';
fid = fopen(filename, 'a');


%%%% =========== setup of the parameters ==================================
Ns = 400:40:520; 
Qs = 800:200:1400;
hps = 0.044:0.002:0.044;
for N = Ns; 
for Q = Qs
for  hp = hps; 
            Rt = 0.5;    % the scaling factor of the object inside of the box
            P = 240; 
            hm = -hp; 
            kp = 40;  theta = -pi/4;   phi = pi/3;
            km = 60;
            % theta and phi are the direction of the incident plane wave traveling in
            % polar coordinate system
            kx = kp * cos(theta) * cos(phi);  ky = kp * cos(theta) * sin(phi);  kz = kp * sin(theta);
            
            %%%% ================= Get the target and MFS points ======================
            %a = 0.2; w = 2; p = pi/2;
            a = 0.3; w = 8; p = 0;
            f = @(s) 1+a*cos(w*(s-p));
            f_t = @(s) -a*w*sin(w*(s-p));
            % std a=0.3, w = 4, p = 0;
			fprintf(fid, '=====================================It ueses SVD in this program =================\n');
			fprintf(fid, 'non-periodic transmission problem with kp = %d, km = %d, w = %d, Rt = 0.5\n', kp, km, w); 
			M = 1.2 * N; 

            
            target = curve3DQuad(f, f_t, M, P, Rt);           % curve target points and its normal derivatives;
            xt = target.x;  yt = target.y;  zt = target.z;
            source_p = source3DQuad(f, f_t, N, Q, Rt, hp);      % curve MFS points, propotional to local speed;
            xsp = source_p.x;  ysp = source_p.y;  zsp = source_p.z;
            source_m = source3DQuad(f, f_t, N, Q, Rt, hm);      % curve MFS points, propotional to local speed;
            xsm = source_m.x;  ysm = source_m.y;  zsm = source_m.z;
            
            targetBOR = BORcurve3DQuad(f, f_t, M, Rt);           % boundary pts in the r-z plane;
            sourceBOR_p = BORsource3DQuad(f, f_t, N, Rt, hp);  % MFS pts in the r-z plane;
            sourceBOR_m = BORsource3DQuad(f, f_t, N, Rt, hm);
            An = EvalMFSHelm3DKernalBOR_trans(targetBOR, sourceBOR_p, sourceBOR_m, kp,km, P, Q);
  			%AQ = 0;   AR = 0;           
            for n=1:P
    
    [AU(:,:,n), AS, AV(:,:,n)] = svd(An(:,:,n), 0);
    ss = diag(AS);
    nA = length(ss);
    rA = length(find (ss >1e-11));   % choose a cut, 1e-10.
    iA(:,n) =  [1./ss(1:rA); zeros(nA-rA,1)];
 %  iA(:,:,n) = pinv(An(:,:,n));  
 % [AQ(:,:,n), AR(:,:,n)] = qr(An(:,:,n));
end
            
            [f,g] = QuasiPerVecFilling3DBOR_trans(targetBOR, kx, ky, kz, P);
            fn = fft(f, [], 2)/P;
            gn = fft(g, [], 2)/P;
            bn = zeros(2*M, P);
            c = zeros(2*N, P);
            for n=1:P
                bn (1:M, n) = fn(:,n);
                bn (M + 1: 2*M, n) = gn(:,n);
       			c(:,n) = An(:,:,n)\bn(:,n);
       			c(:,n) = AV(:,:,n) * (iA(:,n) .* (AU(:,:,n)' * bn(:,n) ));
            	%c(:,n) = An(:,:,n) \ bn(:,n); 
            	%opts.UT = true;
    
        		%c(:,n) =  linsolve( AR(:,:,n) , AQ(:,:,n)' * bn(:,n), opts);

            end
            normc = norm(c); 
            error = zeros(P, 1);
            for n=1:P
                error(n) = norm(An(:,:,n)*c(:,n)-bn(:,n))/sqrt(length(bn(:,n)));
            end
            res = norm(error)/sqrt(P);
            
            %cc = coeffBORto3D(c, Q);
            
            c1 = c(1:N, :);
            c2 = c(N+1:2*N, :);
            cc1 = coeffBORto3D(c1, Q);
            cc2 = coeffBORto3D(c2, Q);
            
            
            
            
            %%% ============================Get the matrix and vector =================
            
          
            
            if 1,
               
                M_test = 128; P_test = 128;
                f = @(s) 1+a*cos(w*(s-p));
                f_t = @(s) -a*w*sin(w*(s-p));
                target_test = curve3DQuad(f, f_t, M_test, P_test, Rt);
                [u_exact, un_exact] = QuasiPerVecFilling3D_trans(target_test, kx, ky, kz);
                
                x_test = target_test.x;  y_test = target_test.y;  z_test = target_test.z;
                [up, Ep] = fmm3DHelpoteval(xsp, ysp, zsp, cc1, kp, x_test, y_test, z_test);
                [um, Em] = fmm3DHelpoteval(xsm, ysm, zsm, cc2, km, x_test, y_test, z_test);
                
                u = up - um;
                u = u.';
                
                Ep = -Ep.'; % change from a row vector to a column vector;  negative sign does not know why;
                Em = -Em.';
                
                unp = Ep(:,1).*target_test.nx + Ep(:,2).*target_test.ny + Ep(:,3).*target_test.nz;
                unm = Em(:,1).*target_test.nx + Em(:,2).*target_test.ny + Em(:,3).*target_test.nz;
                
                un = unp - unm;
                err1 = norm(u - u_exact)/sqrt (length(u)); 
                err2 = norm(un - un_exact)/sqrt (length(un)); 
                
            end
 	        toc; 
            
                    fprintf(fid, 'M = %d, N = %d, P = %d, Q = %d, h is %4.3f, norm of c is %e, and residual is %e; \nerr1 is %e, err2 is %e, time for the code is %4.2f seconds \n\n ',...
                M, N, P, Q, hp, normc, res, err1, err2, toc);
            


end
end
end

fclose(fid);


