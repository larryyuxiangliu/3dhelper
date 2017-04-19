function error = getBdyError(M_test, P_test, f, f_t, Rt,hs, hm, ...
    xsp,ysp,zsp, kp,km, kx, ky, kz, N,P, q,...
    PP,eta, xi, nei ,alpha, beta, ex, ey, type)
%
% getBdyError returns the error on the boundary;
%
% Larry Liu, 07/31/2014

if type == 'd',
    xs = xsp;
    ys = ysp;
    zs = zsp;
    k = kp;
    target_test = curve3DQuad(f, f_t, M_test, P_test, Rt);
    b_test = QuasiPerVecFilling3D(target_test, kx, ky, kz);
    
    x_test = target_test.x;   y_test = target_test.y;  z_test = target_test.z;
    
    c = reshape(eta, P, N);
    c = c.';
    cc = fft(c, [],2)/(2*pi);
    ccc = coeffBORto3D(cc, q);
    
    source = source3DQuad(f, f_t, N, q, Rt, hs);
    xss = source.x;  yss = source.y;  zss = source.z;
    u2 = fmm3DHelpoteval(xss, yss, zss, ccc, k, x_test, y_test, z_test);
    u2 = u2.';
   
    u1 = applyAelse (xs,ys,zs,eta, kp, target_test,nei, alpha, beta,ex,ey, type); 
    b_mfs =  u1 + u2;
    rscale = 1;
    center = [0; 0; 0];
    nterms = PP;
    clm = xi(1:(PP+1)^2);
    
    ztarg = [target_test.x.'; target_test.y.'; target_test.z.'];
    znor = [target_test.nx.'; target_test.ny.'; target_test.nz.'];
    u1 = h3dtaeval(k,rscale,center,nterms,clm,ztarg,znor);
    
    b_mfs = b_mfs + u1;
    error = norm(b_mfs - b_test)/sqrt(length(b_mfs));
    
elseif type =='n',  % when deals with Neumann scattering;
    xs = xsp;
    ys = ysp;
    zs = zsp;
    
    k = kp;
    if isnumeric(f),
        target_test = cupBdy3D(f, f_t, M_test, P_test, Rt);
        [~ , bn_test] = QuasiPerVecFilling3D_trans(target_test, kx, ky, kz);
        
        x_test = target_test.x;   y_test = target_test.y;  z_test = target_test.z;
        
        c = reshape(eta, P, N);
        c = c.';
        cc = fft(c, [],2)/(2*pi);
        ccc = coeffBORto3D(cc, q);
        source = srcBdy3D(f, f_t, N, hs, q, Rt);
    else
        target_test = curve3DQuad(f, f_t, M_test, P_test, Rt);
        [~ , bn_test] = QuasiPerVecFilling3D_trans(target_test, kx, ky, kz);
        
        x_test = target_test.x;   y_test = target_test.y;  z_test = target_test.z;
        
        c = reshape(eta, P, N);
        c = c.';
        cc = fft(c, [],2)/(2*pi);
        ccc = coeffBORto3D(cc, q);
        source = source3DQuad(f, f_t, N, q, Rt, hs);
    end
    xss = source.x;  yss = source.y;  zss = source.z;
    [~, E] = fmm3DHelpoteval(xss, yss, zss, ccc, k, x_test, y_test, z_test);
    E = - E.';
    u2 = E(:,1).*target_test.nx + E(:,2).*target_test.ny + E(:,3).*target_test.nz;
    
    u1 = applyAelse (xs,ys,zs,eta, kp,target_test,nei, alpha, beta,ex,ey, type); 
    bn_mfs =  u1 + u2;
    
    
    
    rscale = 1;
    center = [0; 0; 0];
    nterms = PP;
    clm = xi(1:(PP+1)^2);
    
    ztarg = [target_test.x.'; target_test.y.'; target_test.z.'];
    znor = [target_test.nx.'; target_test.ny.'; target_test.nz.'];
    [~, un_sh] = h3dtaeval(kp,rscale,center,nterms,clm,ztarg,znor);
    
    bn_mfs = bn_mfs + un_sh;
    error = norm(bn_mfs - bn_test)/sqrt(length(bn_mfs));
    
    
else      % when deals with transmission scattering;
    
    target_test = curve3DQuad(f, f_t, M_test, P_test, Rt);
    b_test = QuasiPerVecFilling3D_trans(target_test, kx, ky, kz);
    
    x_test = target_test.x;   y_test = target_test.y;  z_test = target_test.z;
    
    c = reshape(eta, P, 2*N);
    c = c.';
    cc = fft(c, [],2)/(2*pi);
    ccc = coeffBORto3D(cc, q);
    
    source = source3DQuad(f, f_t, N, q, Rt, hs);
    xss = source.x;  yss = source.y;  zss = source.z;
    u2p = fmm3DHelpoteval(xss, yss, zss, ccc(1:end/2), kp, x_test, y_test, z_test);
    
    
    source = source3DQuad(f, f_t, N, q, Rt, hm);
    xss = source.x;  yss = source.y;  zss = source.z;
    u2m = fmm3DHelpoteval(xss, yss, zss, ccc(end/2+1:end), km, x_test, y_test, z_test);
    
    
    u2 = u2p - u2m;
    u2 = u2.';
    

    u1 = applyAelse (xsp,ysp,zsp, eta(1:end/2), kp,target_test,nei, alpha, beta,ex,ey, type); 
    u1 = u1(1:end/2);  % get the only values, not derivatives; 
    b_mfs =  u1 + u2;
    
    
    
    rscale = 1;
    center = [0; 0; 0];
    nterms = PP;
    clm = xi(1:(PP+1)^2);
    
    ztarg = [target_test.x.'; target_test.y.'; target_test.z.'];
    znor = [target_test.nx.'; target_test.ny.'; target_test.nz.'];
    u1 = h3dtaeval(kp,rscale,center,nterms,clm,ztarg,znor);
    
    b_mfs = b_mfs + u1;
    error = norm(b_mfs - b_test)/sqrt(length(b_mfs));
    
    
    
end
