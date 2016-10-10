function p_j = Calculate_pj(c_j,T)
p_j=cell(T+1,1);
N=length(c_j{2});
p_j{1}=(1-exp(-1))*ones(1,N);
p_j{2} = 1./(1+c_j{2});
for t=3:T+1
    temp = -log(max(1-p_j{t-1},realmin));
    p_j{t} = temp./(temp+c_j{t});
    %p_j{t} = 1./(1+c_j{t}./(-log(max(1-p_j{t-1},realmin))));
    if nnz(isnan(p_j{t}))
        warning('pj Nan');
        p_j{t}(isnan(p_j{t}))= eps ;
    end
end
