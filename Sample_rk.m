function [r_k,gamma0,c0]=Sample_rk(XTplusOne_sum,r_k,p_jTplusOne,gamma0,c0,IsNoSample,e0,f0,a0,b0)
if nargin<6
    IsNoSample=false;
end
if nargin<7
    e0=1;
    f0=1;
    a0=0.01;
    b0=0.01;
end

if size(XTplusOne_sum,2)>1
    XTplusOne_sum = full(sum(XTplusOne_sum,2));
end

KT=length(r_k);


if ~IsNoSample
    c0 = randg(e0+gamma0)/(f0+sum(r_k));
    sumlogpi = sum(log(max(1-p_jTplusOne,realmin)));
    p_prime = -sumlogpi./(c0-sumlogpi);
    % L_k = CRT_sum_mex_matrix_v1(sparse(Xt'),r_k')';
    %L_k = full(sum(XTplusOne,2));
    %XTplusOne_sum
    gamma0 = randg(a0 + CRT_sum_mex_v1(XTplusOne_sum,gamma0/KT))/(b0 - log(max(1-p_prime,realmin)));
    r_k = randg(gamma0/KT+XTplusOne_sum)./(c0-sumlogpi);
else
    c0 = (e0+gamma0)/(f0+sum(r_k));
    sumlogpi = sum(log(max(1-p_jTplusOne,realmin)));
    p_prime = -sumlogpi./(c0-sumlogpi);
    % L_k = CRT_sum_mex_matrix_v1(sparse(Xt'),r_k')';
    %L_k = full(sum(XTplusOne,2));
    %XTplusOne_sum;
    temp = gamma0/KT*sum(psi(XTplusOne_sum+gamma0/KT)-psi(gamma0/KT));
    gamma0 = (a0 + temp)/(b0 - log(max(1-p_prime,realmin)));
    r_k = (gamma0/KT+XTplusOne_sum)./(c0-sumlogpi);
end