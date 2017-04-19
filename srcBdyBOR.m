function srcBOR = srcBdyBOR(a, b, N, h, Rt)

bb = cupBdy(a, b, 2*N, Rt);
z = bb.x + 1i * bb.y;
zz = imagshiftcurve(z,h);
srcBOR.r = real(zz(1:N));
srcBOR.z = imag(zz(1:N));


end