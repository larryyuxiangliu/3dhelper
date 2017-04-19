function uu = apply_GMRES_Matrix_BOR_SVD(bdy3D,AU,iA,AV,QU,iQ,QV, B, C, M1, P, xsp, ysp, zsp, x, kp , nei, alpha, beta, ex, ey, type)
% 
% Notes: This code without "test" means that in GMRES iteration, we use
% Barnett-Greengard scheme, not Arvind one. 
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

    
    u1 = C*x;  % apply the C matrix
    
    u2 = QV * (iQ .* (QU' * u1));  % apply the Q^{\dag} matrix;
    
    u3 = B*u2;  % apply the B matrix;
    
    
    bn = reshape(u3, P, M1);
    bn = bn.';
    fbn = fft(bn, [],2)/P;
    % this line and previous line are used to change the vector into its
    % Fourier mode and then apply matrix A_0^{\dag} in each Fourier mode;
    
    for n=1:P
        Afbn(:,n) = AV(:,:,n) * (iA(:,n) .* (AU(:,:,n)' * fbn(:,n)));
    end
    % apply matrix A_0^{\dag} in each Fourier mode;
    
    
    b = coeffBORto3D(Afbn, P);
    %Afbn (:,2:end) = Afbn(:, end:-1:2); % flip FOurier coeffs m to -m
    %b = zeros(N1*P,1);
    %for i=1:N1
    %    b((i-1)*P+1:i*P) = fft(Afbn(i,:))*2*pi/P;
    %end
    % the previous 5 lines of code is to transfer the Fourier mode vector into
    % its original 3D space vector;
    
    u4 = b;
    
    u5 = applyAelse (xsp,ysp,zsp,x, kp,bdy3D,nei, alpha, beta,ex,ey, type); 
    % Using FMM to do the field evaluation. i.e., apply A_else matrix;
    
    bn = reshape(u5, P, M1);
    bn = bn.';
    fbn = fft(bn, [],2)/P;
    % same, the previous 2 lines are used to change a vector into its Fourier
    % mode;
    
    for n=1:P
        Afbn(:,n) = AV(:,:,n) * (iA(:,n) .* (AU(:,:,n)' * fbn(:,n) ));
    end
    % apply A_0^{\dag} in each Fourier mode;
    
    b = coeffBORto3D(Afbn, P);
    
    %Afbn (:,2:end) = Afbn(:, end:-1:2); % flip FOurier coeffs m to -m
    %b = zeros(N1*P,1);
    %for i=1:N1
    %    b((i-1)*P+1:i*P) = fft(Afbn(i,:))*2*pi/P;
    %end
    % go back to its orginal 3D space;
    
    u6 = b;
    uu = x + u6 - u4;
    % the output is:
    % uu = x + A_0^{\dag} * A_else*x - A_0^{\dag} * B * Q^{\dag} * C * x;
    
    
    
elseif type =='n',       % when deals with Neumann scattering;
    
    
    u1 = C*x;  % apply the C matrix
    
    u2 = QV * (iQ .* (QU' * u1));  % apply the Q^{\dag} matrix;
    
    u3 = B*u2;  % apply the B matrix;
    
    
    
    bn = reshape(u3, P, M1);
    bn = bn.';
    fbn = fft(bn, [],2)/P;
    % this line and previous line are used to change the vector into its
    % Fourier mode and then apply matrix A_0^{\dag} in each Fourier mode;
    
    for n=1:P
        Afbn(:,n) = AV(:,:,n) * (iA(:,n) .* (AU(:,:,n)' * fbn(:,n)));
    end
    % apply matrix A_0^{\dag} in each Fourier mode;
    
    
    b = coeffBORto3D(Afbn, P);
    %Afbn (:,2:end) = Afbn(:, end:-1:2); % flip FOurier coeffs m to -m
    %b = zeros(N1*P,1);
    %for i=1:N1
    %    b((i-1)*P+1:i*P) = fft(Afbn(i,:))*2*pi/P;
    %end
    % the previous 5 lines of code is to transfer the Fourier mode vector into
    % its original 3D space vector;
    
    u4 = b;

    u5 = applyAelse (xsp,ysp,zsp,x, kp,bdy3D,nei, alpha, beta,ex,ey, type); 
    % Using FMM to do the field evaluation. i.e., apply A_else matrix;
    
    bn = reshape(u5, P, M1);
    bn = bn.';
    fbn = fft(bn, [],2)/P;
    % same, the previous 2 lines are used to change a vector into its Fourier
    % mode;
    
    for n=1:P
        Afbn(:,n) = AV(:,:,n) * (iA(:,n) .* (AU(:,:,n)' * fbn(:,n) ));
    end
    % apply A_0^{\dag} in each Fourier mode;
    
    b = coeffBORto3D(Afbn, P);
    
    u6 = b;
    uu = x + u6 - u4;
    % the output is:
    % uu = x + A_0^{\dag} * A_else*x - A_0^{\dag} * B * Q^{\dag} * C * x;
    
    
else        % when deals with transmission
    %AQ = AU; 
    %AR = iA; 
    u1 = C*x;  % apply the C matrix
    
    u2 = QV * (iQ .* (QU' * u1));  % apply the Q^{\dag} matrix;
    
    u3 = B*u2;  % apply the B matrix;
    
    
    
    bn = reshape(u3, P, 2*M1);
    bn = bn.';
    fbn = fft(bn, [],2)/P;
    % this line and previous line are used to change the vector into its
    % Fourier mode and then apply matrix A_0^{\dag} in each Fourier mode;
    
    for n=1:P
        Afbn(:,n) = AV(:,:,n) * (iA(:,n) .* (AU(:,:,n)' * fbn(:,n)));
        %Afbn(:,n) = iA(:,:,n) \ fbn(:,n); 
        %opts.UT = true;
    
    	%Afbn(:,n) =  linsolve( AR(:,:,n) , AQ(:,:,n)' * fbn(:,n), opts);
    end
    % apply matrix A_0^{\dag} in each Fourier mode;
    
    
    b = coeffBORto3D(Afbn, P);
    %Afbn (:,2:end) = Afbn(:, end:-1:2); % flip FOurier coeffs m to -m
    %b = zeros(N1*P,1);
    %for i=1:N1
    %    b((i-1)*P+1:i*P) = fft(Afbn(i,:))*2*pi/P;
    %end
    % the previous 5 lines of code is to transfer the Fourier mode vector into
    % its original 3D space vector;
    
    u4 = b;
    
     u5 = applyAelse (xsp,ysp,zsp,x, kp,bdy3D,nei, alpha, beta,ex,ey, type); 
    % Using FMM to do the field evaluation. i.e., apply A_else matrix;
    
    bn = reshape(u5, P, 2*M1);
    bn = bn.';
    fbn = fft(bn, [],2)/P;
    % same, the previous 2 lines are used to change a vector into its Fourier
    % mode;
    
    for n=1:P
        Afbn(:,n) = AV(:,:,n) * (iA(:,n) .* (AU(:,:,n)' * fbn(:,n) ));
    	%Afbn(:,n) = iA(:,:,n) \ fbn(:,n); 
    	%opts.UT = true;
    
    	%Afbn(:,n) =  linsolve( AR(:,:,n) , AQ(:,:,n)' * fbn(:,n), opts);
    end
    % apply A_0^{\dag} in each Fourier mode;
    
    b = coeffBORto3D(Afbn, P);
    
    u6 = b;
    uu = x + u6 - u4;
    % the output is:
    % uu = x + A_0^{\dag} * A_else*x - A_0^{\dag} * B * Q^{\dag} * C * x;
    
    
    
end


