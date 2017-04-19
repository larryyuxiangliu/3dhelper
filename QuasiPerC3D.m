function C = QuasiPerC3D(t, bo, l, r, f,ba, xs,ys,zs, k, ex,ey, nei, alpha, beta)
% 
% C = QuasiPerC(t, b, l, r, xs,ys, M1, N2,k, ex,ey, nei, alpha) returns the 
% matrix C in the notes; 
% 
% Notes: the N1 is the number of the source points in the r-z plane; NP is
% the number of source points in the phi direction 
% 
% By Larry Liu 05/18/2014

m = length(t.x); 
N = length(xs);  
C = zeros(8*m, N); 

% [u_l, un_l] = EvalMFSHelm3DKernal (l, xs, ys, zs, k); 
% [u_r, un_r] = EvalMFSHelm3DKernal (r, xs, ys, zs, k); 
% 
% [u_f, un_f] = EvalMFSHelm3DKernal (f, xs, ys, zs, k); 
% [u_ba, un_ba] = EvalMFSHelm3DKernal (ba, xs, ys, zs, k); 

% here I used cancellation; 

[u_bo, un_bo] = EvalMFSHelm3DKernalQuasi(bo, xs, ys, zs, k, ex, ey, nei, alpha, beta); 
[u_t, un_t] = EvalMFSHelm3DKernalQuasi(t, xs, ys, zs, k, ex, ey, nei, alpha, beta); 


[u_l, un_l] = EvalMFSHelm3DKernalQuasi(l, xs, ys, zs, k, ex, ey, nei, alpha, beta); 
[u_r, un_r] = EvalMFSHelm3DKernalQuasi(r, xs, ys, zs, k, ex, ey, nei, alpha, beta); 

[u_f, un_f] = EvalMFSHelm3DKernalQuasi(f, xs, ys, zs, k, ex, ey, nei, alpha, beta); 
[u_ba, un_ba] = EvalMFSHelm3DKernalQuasi(ba, xs, ys, zs, k, ex, ey, nei, alpha, beta); 

C(1:m, :) = u_r - alpha * u_l; 
C(m+1:2*m, :) = un_r - alpha * un_l; 
C(2*m+1:3*m, :) = u_ba - beta * u_f; 
C(3*m+1:4*m, :) = un_ba - beta * un_f; 

% C(1:m, :) = alpha^(-1) * u_r * (beta^(-1) + beta + 1) - alpha^2 * u_l *  (beta ^(-1) + 1 + beta); 
% C(m+1:2*m, :) = alpha^(-1) * un_r * (beta^(-1) + beta + 1) - alpha^2 * un_l *  (beta ^(-1) + 1 + beta); 
% 
% C(2*m+1:3*m, :) = beta^(-1) * u_ba * (alpha^(-1) + alpha + 1) - beta^2 * u_f * (alpha ^(-1) + 1 + alpha); 
% C(3*m+1:4*m, :) = beta^(-1) * un_ba * (alpha^(-1) + alpha + 1) - beta^2 * un_f * (alpha ^(-1) + 1 + alpha); 

C(4*m+1:5*m, :) = u_t; 
C(5*m+1:6*m, :) = un_t; 

C(6*m+1:7*m, :) = u_bo; 
C(7*m+1:8*m, :) = un_bo; 

