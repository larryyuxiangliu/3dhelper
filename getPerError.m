function [err1, err2] = getPerError (M, L,  xs, ys, zs, eta, c, P, k, nei ,alpha, beta, ex, ey)
%
% getPerError evaluates the periodicity error on the box;
%
% Larry Liu, 07/31/2014

[~,~,LL,RR,FF,BA] = squarequad(M, L);


u_l = fmm3DPerHelpoteval(xs, ys, zs, eta, k, LL.x, LL.y, LL.z, nei, alpha, beta, ex, ey);
rscale = 1;
center = [0; 0; 0];
nterms = P;
clm = c(1:(P+1)^2);
l = LL;
r = RR;
f = FF;
ba = BA;
ztarg = [l.x.'; l.y.'; l.z.'];
znor = [l.nx.'; l.ny.'; l.nz.'];
u_l = u_l + h3dtaeval(k,rscale,center,nterms,clm,ztarg,znor).';

u_r = fmm3DPerHelpoteval(xs, ys, zs, eta, k, RR.x, RR.y, RR.z, nei, alpha, beta, ex, ey);
ztarg = [r.x.'; r.y.'; r.z.'];
znor = [r.nx.'; r.ny.'; r.nz.'];
u_r = u_r + h3dtaeval(k,rscale,center,nterms,clm,ztarg,znor).';

u_f = fmm3DPerHelpoteval(xs, ys, zs, eta, k, FF.x, FF.y, FF.z, nei, alpha, beta, ex, ey);
ztarg = [f.x.'; f.y.'; f.z.'];
znor = [f.nx.'; f.ny.'; f.nz.'];
u_f = u_f +  h3dtaeval(k,rscale,center,nterms,clm,ztarg,znor).';

u_ba = fmm3DPerHelpoteval(xs, ys, zs, eta, k, BA.x, BA.y, BA.z, nei, alpha, beta, ex, ey);
ztarg = [ba.x.'; ba.y.'; ba.z.'];
znor = [ba.nx.'; ba.ny.'; ba.nz.'];
u_ba = u_ba +  h3dtaeval(k,rscale,center,nterms,clm,ztarg,znor).';

err1 = norm(alpha * u_l - u_r)/sqrt(length(u_l));
err2 = norm(beta * u_f - u_ba)/sqrt(length(u_f));


