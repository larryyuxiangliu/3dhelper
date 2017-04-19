function [u] = applyAelse_matvec(Us, Ks, Vs, c, nei, target, type)

M = length(target.x);

if (type == 'd') || (type == 'n')
	u = zeros(M, 1);
else 
	u = zeros(2*M, 1); 
end

ii = 1;
c0 = c;
for i = -nei:nei
    for j = -nei:nei
        if ((i ~= 0) || (j ~= 0))
            if (type == 'd') || (type == 'n')
                U = Us{ii};
                K = Ks{ii};
                V = Vs{ii};

                
                
               
                u0 = U * (K * (V * c));
                

                
                u = u + u0;
                ii = ii + 1;
            else
                
                U = Us{ii, 1};
                K = Ks{ii, 1};
                V = Vs{ii};

                
                

                u0 = U * (K * (V * c));
                

                
                
                U = Us{ii, 2};
                K = Ks{ii, 2};

                u0n = U * (K * (V * c));

                u0 = [u0; u0n];
                u = u + u0;
                ii = ii + 1;
            end
            
            
            
        end
        
    end
end
