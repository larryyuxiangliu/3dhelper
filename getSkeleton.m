function [U, K, V, IS, IT, RR] = getSkeleton(target, source,  k, acc, RP, NP, RQ, NQ, nei,alpha, beta, ex, ey, type)

ii = 1;

for i = -nei:nei
    for j = -nei:nei
        if ((i ~= 0) || (j ~= 0))
            %%%%%% real source to target proxy sphere =================================
            [xp, yp, zp] = squaresourcequad(NP, RP);      % box proxy points, points on the big sphere;
            xs = source.x + i*ex;  ys = source.y + j*ey ;  zs = source.z;
            xt = xp;  yt = yp;  zt = zp;
            
            [XS, XT] = meshgrid (xs, xt);
            [YS, YT] = meshgrid (ys, yt);
            [ZS, ZT] = meshgrid (zs, zt);
            
            Rx = XS - XT;
            Ry = YS - YT;
            Rz = ZS - ZT;
            R = sqrt(Rx.^2+ Ry.^2+ Rz.^2);
            as = exp(1i*k*R)./R;
            
            [TS, IS{ii}] = id_decomp(as,acc);
            %[TS, IS{ii}] = comp_id_decomp(as,acc);
            rs = size(TS, 1);
            RR(ii, 1) = rs;
            V{ii} = [eye(rs), TS];
            V{ii}(:,IS{ii}) = V{ii} (:, :);    % change the indices; 
            xsf = xs(IS{ii}(1:rs));    ysf = ys(IS{ii}(1:rs));  zsf = zs(IS{ii}(1:rs));
            
            %==========================================================================
            
            %===== proxy source to real target=========================================
            [xp, yp, zp] = squaresourcequad(NQ, RQ);
            xs = xp + i * ex;   ys = yp + j * ey;    zs = zp;
            xt = target.x;  yt = target.y;  zt = target.z;
            if type == 'd'
                at = EvalMFSHelm3DKernal (target, xs, ys, zs, k);
                at = at';
                [TT, IT{ii}] = id_decomp(at,acc);
                %[TT, IT{ii}] = comp_id_decomp(at,acc);
                rt = size(TT,1);
                RR(ii,2) = rt;
                U{ii} = [eye(rt), TT];
                U{ii} = U{ii}';
                
                U{ii}(IT{ii},:) = U{ii} (:, :);  % change the indices; 
                
                xtf = xt(IT{ii}(1:rt));   ytf = yt(IT{ii}(1:rt));     ztf = zt(IT{ii}(1:rt));
            elseif type == 'n'
                [~, at] = EvalMFSHelm3DKernal (target, xs, ys, zs, k);
                at = at';
                [TT, IT{ii}] = id_decomp(at,acc);
                %[TT, IT{ii}] = comp_id_decomp(at,acc);
                rt = size(TT,1);
                RR(ii,2) = rt;
                U{ii} = [eye(rt), TT];
                U{ii} = U{ii}';
                U{ii}(IT{ii},:) = U{ii} (:, :);  % change the indices; 
                xtf = xt(IT{ii}(1:rt));   ytf = yt(IT{ii}(1:rt));     ztf = zt(IT{ii}(1:rt));
                xtn = target.nx(IT{ii}(1:rt));  ytn = target.ny(IT{ii}(1:rt)); ztn = target.nz(IT{ii}(1:rt));
            else
                [a1, a2] = EvalMFSHelm3DKernal (target, xs, ys, zs, k);
                % Here, I do two IDs, one for the potential values, another
                % for normal derivatives; 
                
                % === get the potential values =========
                at = a1';
                [TT, IT{ii,1}] = id_decomp(at,acc);
                %[TT, IT{ii,1}] = comp_id_decomp(at,acc);
                rt = size(TT,1);
                RR(ii,2) = rt;
                U{ii,1} = [eye(rt), TT];
                U{ii,1} = U{ii,1}';
                U{ii, 1}(IT{ii, 1},:) = U{ii, 1} (:, :);  % change the indices; 
                xtf = xt(IT{ii,1}(1:rt));   ytf = yt(IT{ii,1}(1:rt));     ztf = zt(IT{ii,1}(1:rt));
                xtn = target.nx(IT{ii,1}(1:rt));  ytn = target.ny(IT{ii,1}(1:rt)); ztn = target.nz(IT{ii,1}(1:rt));
                
                t1.x = xtf;  t1.y = ytf;  t1.z = ztf;
                t1.nx = xtn;  t1.ny = ytn;  t1.nz = ztn;
                
                % =========== get the normal derivatives ================
                at = a2';
                [TT, IT{ii,2}] = id_decomp(at,acc);
                %[TT, IT{ii,2}] = comp_id_decomp(at,acc);
                rt = size(TT,1);
                RR(ii,3) = rt;
                U{ii,2} = [eye(rt), TT];
                U{ii,2} = U{ii,2}';
                U{ii, 2}(IT{ii, 2},:) = U{ii, 2} (:, :);  % change the indices; 
                xtf = xt(IT{ii,2}(1:rt));   ytf = yt(IT{ii,2}(1:rt));     ztf = zt(IT{ii,2}(1:rt));
                xtn = target.nx(IT{ii,2}(1:rt));  ytn = target.ny(IT{ii,2}(1:rt)); ztn = target.nz(IT{ii,2}(1:rt));
                
                t2.x = xtf;  t2.y = ytf;  t2.z = ztf;
                t2.nx = xtn;  t2.ny = ytn;  t2.nz = ztn;
                
            end
            
            
            
            % ============== Build the K matrix ===========================
            
            
            if type == 'd'
                t.x = xtf;  t.y = ytf;  t.z = ztf;
                UU = EvalMFSHelm3DKernal (t, xsf, ysf, zsf, k);
                
                K{ii} = alpha^(i) * beta^(j) * UU;
            elseif type == 'n'
                
                t.x = xtf;  t.y = ytf;  t.z = ztf;
                t.nx = xtn;  t.ny = ytn;  t.nz = ztn;
                [~, UU] = EvalMFSHelm3DKernal (t, xsf, ysf, zsf, k);
                
                K{ii} = alpha^(i) * beta^(j) * UU;
            else
                
           
                [UU] = EvalMFSHelm3DKernal (t1, xsf, ysf, zsf, k);
                
                K{ii,1} = alpha^(i) * beta^(j) * UU;
                
                [~, UU] = EvalMFSHelm3DKernal (t2, xsf, ysf, zsf, k);
                
                K{ii,2} = alpha^(i) * beta^(j) * UU;
            end
            
            
            
            
            % =============================================================
            
            
            ii = ii + 1;
        end
        
    end
end


