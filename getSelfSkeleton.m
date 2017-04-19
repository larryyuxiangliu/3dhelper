function [U, K, V] = getSelfSkeleton(target, source,  k, acc, RP, NP, RQ, NQ)





%%%%%% real source to target proxy sphere =================================
[xp, yp, zp] = squaresourcequad(NP, RP);      % box proxy points, points on the big sphere;
xs = source.x;  ys = source.y;  zs = source.z;
xt = xp;  yt = yp;  zt = zp;

[XS, XT] = meshgrid (xs, xt);
[YS, YT] = meshgrid (ys, yt);
[ZS, ZT] = meshgrid (zs, zt);

Rx = XS - XT;
Ry = YS - YT;
Rz = ZS - ZT;
R = sqrt(Rx.^2+ Ry.^2+ Rz.^2);
as = exp(1i*k*R)./R;

[TS, IS] = id_decomp(as,acc);
%[TS, IS{ii}] = comp_id_decomp(as,acc);
rs = size(TS, 1);
V = [eye(rs), TS];
V(:,IS) = V(:, :);    % change the indices;
xsf = xs(IS(1:rs));    ysf = ys(IS(1:rs));  zsf = zs(IS(1:rs));

%==========================================================================

%===== proxy source to real target=========================================
[xp, yp, zp] = squaresourcequad(NQ, RQ);
xs = xp;   ys = yp;    zs = zp;
xt = target.x;  yt = target.y;  zt = target.z;

at = EvalMFSHelm3DKernal (target, xs, ys, zs, k);
at = at';
[TT, IT] = id_decomp(at,acc);
%[TT, IT{ii}] = comp_id_decomp(at,acc);
rt = size(TT,1);
U = [eye(rt), TT];
U = U';

U(IT,:) = U(:, :);  % change the indices;

xtf = xt(IT(1:rt));   ytf = yt(IT(1:rt));     ztf = zt(IT(1:rt));



% ============== Build the K matrix ===========================

t.x = xtf;  t.y = ytf;  t.z = ztf;
UU = EvalMFSHelm3DKernal (t, xsf, ysf, zsf, k);

K = UU;
% =============================================================





