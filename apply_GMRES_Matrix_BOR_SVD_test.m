function uu = apply_GMRES_Matrix_BOR_SVD_test(bdy3D,AU,iA,AV,QU,iQ,QV, B, C, M1, P, q, xsp, ysp, zsp, x, kp, nei, alpha, beta, ex, ey, type)
% 
% Notes: This code with "test" means that in GMRES iteration, we use
% another conditioning, as Arvind suggested. 
% 
% apply_GMRES_Matrix_BOR is to build up the function ff(x), which takes an
% input vector x0 to its next GMRES step. Details can be seen in the notes
% Hel_3D_Periodic_14X/notes.pdf;
%
% Inputs: AU, AV: the U and V matrix of the svd output of each diagonal
%                 block of A_0; so they have dimensionality 3;
%             iA: the inverse of the singular values of each diagonal block
%                 of A_0; so it has dimensionality 2;
%
%         QU, QV: the U and V matrix of the svd output of Q, it has dim 2;
%             iQ: the inverse of the singular values of Q, it has dim 1 (vector);
%
%           B, C: the matrix B and C in the notes;
%
%     (xs,ys,zs): the location of source points in the 3D Cartesian
%                 coordinate system;
%     (xt,yt,zt): the boundary points in the 3D Cartesian coordinate
%                 system;
%     nei, alpha, beta, ex, ey:  are the periodic parameters;
%
%             M1: the number of boundary points in each Fourier mode;
%             N1: the number of MFS points in each Fourier mode;
%              P: the number of Fourier modes;
%
% Output:    uu: the output vector, which means uu = ff(x);
%
% Notes: This is the BOR case, which means we are using FFT and ultilize the
% axisymmetric property of the objects; This code is used in the
% Hel_3D_GMRES_FMM_BOR.m main code;
%
% Notes: This subroutine uses SVD to the matrix factorization.
%
% Also see: apply_GMRES_Matrix.m
%
% Larry Liu
% 07/14/2014


if type == 'd'    % when deals with Dirichilet scattering;
    %================= calculate the x1 = A_0^{\dag} * x ==================
    X = reshape(x, P, M1);
    X = X.';
    xn = fft(X, [],2)/P;
    for n=1:P
        Axn(:,n) = AV(:,:,n) * (iA(:,n) .* (AU(:,:,n)' * xn(:,n) ));
    end   
    x1 = coeffBORto3D(Axn, q);
    %======================================================================  
    
    u1 = C*x1;  % apply the C matrix  
    u2 = QV * (iQ .* (QU' * u1));  % apply the Q^{\dag} matrix;    
    u3 = B*u2;  % apply the B matrix;
    

    u5 = applyAelse (xsp,ysp,zsp,x1, kp,bdy3D,nei, alpha, beta,ex,ey, type); 
    % Using FMM to do the field evaluation. i.e., apply A_else matrix;
    
    uu = x + u5 - u3;
    % the output is:
    % uu = x +  A_else * A_0^{\dag} * x - B * Q^{\dag} * C * A_0^{\dag} * x;
    
    
    
elseif type =='n',       % when deals with Neumann scattering;
    
    %================= calculate the x1 = A_0^{\dag} * x ==================
    X = reshape(x, P, M1);
    X = X.';
    xn = fft(X, [],2)/P;
    for n=1:P
        Axn(:,n) = AV(:,:,n) * (iA(:,n) .* (AU(:,:,n)' * xn(:,n) ));
    end   
    x1 = coeffBORto3D(Axn, q);
    %======================================================================
    
    u1 = C*x1;  % apply the C matrix
    u2 = QV * (iQ .* (QU' * u1));  % apply the Q^{\dag} matrix;
    u3 = B*u2;  % apply the B matrix;

    u5 = applyAelse (xsp,ysp,zsp,x1, kp,bdy3D,nei, alpha, beta,ex,ey, type);   
    
    uu = x + u5 - u3;
    % the output is:
    % uu = x +  A_else * A_0^{\dag} * x - B * Q^{\dag} * C * A_0^{\dag} * x;
    
    
else        % when deals with transmission
    %================= calculate the x1 = A_0^{\dag} * x ==================
    X = reshape(x, P, 2*M1);
    X = X.';
    xn = fft(X, [],2)/P;
    for n=1:P
        Axn(:,n) = AV(:,:,n) * (iA(:,n) .* (AU(:,:,n)' * xn(:,n) ));
    end   
    x1 = coeffBORto3D(Axn, q);
    %======================================================================
    
    u1 = C*x1;  % apply the C matrix
    u2 = QV * (iQ .* (QU' * u1));  % apply the Q^{\dag} matrix;
    u3 = B*u2;  % apply the B matrix;
    
    u5 = applyAelse (xsp,ysp,zsp,x1(1:end/2), kp,bdy3D,nei, alpha, beta,ex,ey, type); 
    % Using FMM to do the field evaluation. i.e., apply A_else matrix;
    
    uu = x + u5 - u3;
    % the output is:
    % uu = x +  A_else * A_0^{\dag} * x - B * Q^{\dag} * C * A_0^{\dag} * x;    
    
    
end


