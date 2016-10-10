function x = bessel_rnd_v1(alpha,nu)
% Algorithms of Devroye, Luc. "Simulating Bessel random variables." 
% Statistics & probability letters 57.3 (2002): 249-257.
% Algorithm 1:  Page 7 of http://132.206.3.210/bessel.pdf
% Algorithm 2: Page 8 of http://132.206.3.210/bessel.pdf
% Algorithm 3: Yuan, Lin, and John D. Kalbfleisch. "On the Bessel distribution and related problems." Annals of the Institute of Statistical Mathematics 52.3 (2000): 438-447.

% Coded by Mingyuan Zhou, March, 2015
%
m = fix((sqrt(alpha.^2+nu.^2)-nu)/2);
A = alpha.^2+nu.^2;
B = A + 2*nu+1;
A = sqrt(A);
B = sqrt(B);
Q = 0.5*alpha.^2./(nu+A) + (1+alpha.^2.*(1+B-A)./(2*(nu+A).*(nu+1+B))).^2;
x=zeros(size(alpha));
dex=find(alpha>0);

q= min(1./(Q*sqrt(648)),1/3);


while 1
    
    Y=zeros(length(dex),1);
    if isempty(dex)
        break
    else
        U = rand(length(dex),1);
        W = rand(length(dex),1);
        S = (rand(length(dex),1)>0.5)*2 - 1;
        qq = q(dex);
        %idex = (U<=(1+2./q(dex))./(1+4./q(dex))));
        idex = U<= (0.5+0.5./(1+4./qq));
        
        %         size(V)
        %         size(Y)
        %         size(qq(idex));
        if nnz(idex)>0
            V = rand(nnz(idex),1);
            Y(idex)= V.*(1/2+1./qq(idex));
        end
        if nnz(~idex)>0
            %Y(~idex)=1/2+1./qq(~idex)+randg(ones(nnz(~idex),1))./qq(~idex);
            Y(~idex)=1/2+1./qq(~idex)+exp(1)./qq(~idex);
        end
        X = S.*round(Y);
        
        idex = m(dex)+X>=0;
        idex(idex) =  log(W(idex))+min(0,(1+qq(idex)/2-qq(idex).*Y(idex))) <=...
            (2*X(idex).*log(alpha(dex(idex))/2)+gammaln(m(dex(idex))+nu+1)+gammaln(m(dex(idex))+1)-gammaln(m(dex(idex))+X(idex)+nu+1)-gammaln(m(dex(idex))+X(idex)+1));
        %
        %        % idex = (W<0 | X+m(dex) (log(abs(W))+min(0,1+qq/2-qq.*Y)<= 2*X.*log(alpha(dex)/2)+gammaln(m(dex)+nu+1)+gammaln(m(dex)+1)-gammaln(m(dex)+X+nu+1)-gammaln(m(dex)+X+1)));
        %         idex = W.*min(1,exp(1+qq/2-qq.*Y)) <= ...
        %               (alpha(dex)/2).^(2*X).*gamma(m(dex)+nu+1)./gamma(m(dex)+X+nu+1).*gamma(m(dex)+1)./gamma(m(dex)+X+1);
        x(dex(idex))=X(idex)+m(dex(idex));
        dex(idex)=[];
    end
end
%Using Gaussian approximation
