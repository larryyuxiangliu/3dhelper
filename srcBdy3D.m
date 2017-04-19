function src3D = srcBdy3D(a, b, N, h, P, Rt)


p = (0:P-1)/P*2*pi;
p = p'; 

srcBOR = srcBdyBOR(a, b, N, h, Rt);

src3D.x = kron(srcBOR.r, cos(p));
src3D.y = kron(srcBOR.r, sin(p));
src3D.z = kron(srcBOR.z, ones(P,1));






end