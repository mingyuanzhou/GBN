function x = truncated_Poisson_rnd(lambda)
% draw random samples from a truncated Poisson distribution
% P(x=k) = poisspdf(k,lambda)/(1-exp(-lambda))
% Coded by Mingyuan Zhou, 042014
% Using a rejection sampler described in "Mingyuan Zhou, Infinite edge partition models for overlapping community detection and link prediction," Artificial Intelligence and Statistics (AISTATS2015), JMLR W&CP, vol. 38, San Diego, CA, May 2015.
%%Demo:
%lambda=0.3;aa=hist(truncated_Poisson_rnd(lambda*ones(1,10000)),1:20);
%plot(1:20,aa/sum(aa),'r',1:20,poisspdf(1:20,lambda)./(1-exp(-lambda)));

lambda1=lambda(lambda>1);
lambda2=lambda(lambda<=1);
x=zeros(length(lambda),1);
x1 = zeros(length(lambda1),1);
x2 = zeros(length(lambda2),1);

while 1
    dex=find(x1==0);
    if isempty(dex)
        break
    else
        lambdadex=lambda1(dex);
        temp = poissrnd(lambdadex);
        idex = temp>0;
        x1(dex(idex))=temp(idex);
    end
end
x(lambda>1)=x1;

while 1
    dex=find(x2==0);
    if isempty(dex)
        break
    else
        lambdadex=lambda2(dex);
        temp = 1+poissrnd(lambdadex);
        idex = rand(size(temp))<1./(temp);
        x2(dex(idex))=temp(idex);
    end
end
x(lambda<=1)=x2;

