function [u, un] = EvalMFSHelm3DKernalQuasi(t, xs, ys, zs, k, ex, ey, nei, alpha, beta)
%
% EvalMFSHelm3DKernalQuasi(t, xs, ys, zs, k, ex, ey, nei, alpha) returns the
% matrix kenerl but also taking the number of neigbor copies of MFS pts into account,
% together with the phase difference alpha; so it is a direct sum over all
% the neighbor source points;
%
% Inputs: t, 1x1 structure:  t.x, t.y, t.z, t.nx, t.ny, t.nz;
%
% Outputs: u, M x N matrix
% 
% Also see: EvalMFSHel3DKernal.m
% 
% By Larry Liu 05/16/2014

if nargout == 1,
    u  = EvalMFSHelm3DKernal (t, xs, ys, zs, k);
    
    if nei>0
        M = length(t.x);     N = length(xs);
        u = zeros(M, N);
        for i = - nei : nei
            for j = -nei : nei
                u1 = EvalMFSHelm3DKernal(t, xs + i*ex, ys + j*ey, zs, k); 
                u = u + alpha^(i) * beta^(j) * u1; 
            end
        end
        
    end
    
    
else
    [u, un] = EvalMFSHelm3DKernal (t, xs, ys, zs, k);
    
    if nei>0
        M = length(t.x);    N = length(xs);
        u = zeros(M, N);    un = zeros(M,N);       
        for i = - nei : nei
            for j = -nei : nei
                [u1, un1] = EvalMFSHelm3DKernal(t, xs+i*ex, ys+j*ey, zs, k);
                u = u + alpha^(i) * beta^(j) * u1;
                un = un + alpha^(i) * beta^(j) * un1;
            end
        end
    end
    
end
end
