
function fun = zMcDonald(p,z)
    
    fun = abs(z).^(p).*besselk(p,abs(z));
    fun(z==0) = pi^(1/2)*gamma(2*p)/gamma(p+1/2)*2^(-p);

end