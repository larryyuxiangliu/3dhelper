function bdy3D = cupBdy3D(a, b, N, P, Rt)
 


p = (0:P-1)/P*2*pi;
p = p'; 
bdyBOR = cupBdyBOR(a, b, N, Rt); 
bdy3D.x = kron(bdyBOR.r, cos(p));
bdy3D.y = kron(bdyBOR.r, sin(p));
bdy3D.z = kron(bdyBOR.z, ones(P,1));

bdy3D.tx = kron(bdyBOR.rb, cos(p));
bdy3D.ty = kron(bdyBOR.rb, sin(p));
bdy3D.tz = kron(bdyBOR.zb, ones(P,1));

bdy3D.px = kron(bdyBOR.r, -sin(p));
bdy3D.py = kron(bdyBOR.r, cos(p));
bdy3D.pz = zeros(N*P,1);



aa= -cross ([bdy3D.tx, bdy3D.ty, bdy3D.tz], [bdy3D.px, bdy3D.py, bdy3D.pz], 2);  %#ok<*NODEF>
bdy3D.nx = aa(:,1);  
bdy3D.ny = aa(:,2); 
bdy3D.nz = aa(:,3); 
ll = sqrt(bdy3D.nx.^2 + bdy3D.ny.^2 + bdy3D.nz.^2); 
bdy3D.nx = bdy3D.nx./ll; 
bdy3D.ny = bdy3D.ny./ll; 
bdy3D.nz = bdy3D.nz./ll; 


end