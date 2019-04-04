function par = mlePar(y,err)
par = sum(y./err.^2)./sum(1./err.^2); 
end