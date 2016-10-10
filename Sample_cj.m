function c_j=Sample_cj(Theta,r_k,T,IsNoSample,e0,f0)
if nargin<4
    IsNoSample=false;
end
if nargin<5
    e0=1;
    f0=1;
end
if ~IsNoSample
    N=size(Theta{1},2);
    c_j=cell(T+1,1);
    c_j{1}=ones(1,N);
    for t=2:T
        c_j{t} = randg(e0+sum(Theta{t},1))./(f0+sum(Theta{t-1},1));
    end
    c_j{T+1} = randg(e0+sum(r_k)*ones(1,N))./(f0+sum(Theta{T},1));
else
    N=size(Theta{1},2);
    c_j=cell(T+1,1);
    c_j{1}=ones(1,N);
    for t=2:T
        c_j{t} = (e0+sum(Theta{t},1))./(f0+sum(Theta{t-1},1));
    end
    c_j{T+1} = (e0+sum(r_k)*ones(1,N))./(f0+sum(Theta{T},1));
end