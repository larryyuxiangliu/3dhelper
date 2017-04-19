function bdyBOR = cupBdyBOR(a, b, N, Rt)
 



bdy = cupBdy(a, b, 2*N, Rt); 


bdyBOR.r = bdy.x(1:N); 
bdyBOR.z = bdy.y(1:N); 
bdyBOR.rn = bdy.nx(1:N); 
bdyBOR.zn = bdy.ny(1:N); 
bdyBOR.rb = bdy.px(1:N); 
bdyBOR.zb = bdy.py(1:N); 

end