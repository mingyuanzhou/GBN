function prob = pmf_bessel(alpha,nu,u)
    PMF = log(besseli(nu,alpha,1))-nu.*log(alpha/2);
    prob = exp(u.*log(alpha.^2/4)-gammaln(u+nu+1)-gammaln(u+1) -alpha-PMF);
end