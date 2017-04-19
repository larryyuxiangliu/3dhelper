function u = applyAelse (xsp,ysp,zsp,c, k,bdy3D,nei, alpha, beta,ex,ey, type)
% applyAelse is a subroutine to use FMM to apply matrix A_else to a vector;
%
% Larry Liu, 07/14/2014

u_fmm = 0;   un_fmm = 0;
xt = bdy3D.x;  yt = bdy3D.y;  zt = bdy3D.z;
for i= -nei : nei
    for j= -nei : nei
        if ((i~=0) || (j~=0)),  % make sure it deos not include the center one, which is A_0
            [u, E] = fmm3DHelpoteval(xsp + i*ex, ysp + j*ey, zsp, c, k, xt, yt, zt);
            u_fmm = u_fmm + alpha^(i) * beta^(j) * (u.');
            E = -E.';
            un_fmm = un_fmm + (E(:,1).*bdy3D.nx + E(:,2).*bdy3D.ny + E(:,3).*bdy3D.nz)* alpha^(i) * beta^(j);
            
        end
    end
end


if type == 'd'
    u = u_fmm;
elseif type == 'n'
    u = un_fmm;
else
    u = [u_fmm; un_fmm];
end




% elseif type == 'n',    % when deals with Neumann scattering;
%
%     for i= -nei : nei
%         for j= -nei : nei
%             if ((i~=0) || (j~=0)),  % make sure it deos not include the center one, which is A_0
%                 [~, E1] = fmm3DHelpoteval(xsp + i*ex, ysp + j*ey, zsp, c, kp, xt, yt, zt);
%                 E1 = -E1.';
%                 un1_fmm = un1_fmm + (E1(:,1).*bdy3D.nx + E1(:,2).*bdy3D.ny + E1(:,3).*bdy3D.nz)* alpha^(i) * beta^(j);
%
%
%
%             end
%         end
%     end
%
%
%     uun = un1_fmm;
%     uu = 0;
%
%
% elseif type == 't'
%
%     for i= -nei : nei
%         for j= -nei : nei
%             if ((i~=0) || (j~=0)),  % make sure it deos not include the center one, which is A_0
%                 [u1, E1] = fmm3DHelpoteval(xsp + i*ex, ysp + j*ey, zsp, c(1:end/2), kp, xt, yt, zt);
%
%                 u1_fmm = u1_fmm + alpha^(i) * beta^(j) * (u1.');
%                 E1 = -E1.';
%                 un1_fmm = un1_fmm + (E1(:,1).*bdy3D.nx + E1(:,2).*bdy3D.ny + E1(:,3).*bdy3D.nz)* alpha^(i) * beta^(j);
%
%
%
%             end
%         end
%     end
%
%     uu = u1_fmm;
%     uun = un1_fmm;
%
% else         % when deals with transmission
%
%     for i= -nei : nei
%         for j= -nei : nei
%             if ((i~=0) || (j~=0)),  % make sure it deos not include the center one, which is A_0
%                 [u1, E1] = fmm3DHelpoteval(xsp + i*ex, ysp + j*ey, zsp, c(1:end/2), kp, xt, yt, zt);
%                 [u2, E2] = fmm3DHelpoteval(xsm + i*ex, ysm + j*ey, zsm, c(end/2 + 1:end), km, xt, yt, zt);
%                 u1_fmm = u1_fmm + alpha^(i) * beta^(j) * (u1.');
%                 E1 = -E1.';
%                 un1_fmm = un1_fmm + (E1(:,1).*bdy3D.nx + E1(:,2).*bdy3D.ny + E1(:,3).*bdy3D.nz)* alpha^(i) * beta^(j);
%
%                 u2_fmm = u2_fmm + alpha^(i) * beta^(j) * (u2.');
%                 E2 = -E2.';
%                 un2_fmm = un2_fmm + (E2(:,1).*bdy3D.nx + E2(:,2).*bdy3D.ny + E2(:,3).*bdy3D.nz)* alpha^(i) * beta^(j);
%
%             end
%         end
%     end
%
%     uu = u1_fmm - u2_fmm;
%     uun = un1_fmm - un2_fmm;
% end
