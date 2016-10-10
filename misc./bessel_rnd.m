function x = bessel_rnd(alpha,nu,Algorithm,m0)
% Coded by Mingyuan Zhou, March, 2015
%
m = fix((sqrt(alpha.^2+nu.^2)-nu)/2);

if nargin<3
    Algorithm=1;
end
if nargin<4
    m0=20;
end

x=zeros(size(alpha));
% 
dex=find(alpha>0&m<=m0);
PMF = besseli(nu,alpha(dex),1).*(alpha(dex)/2).^(-nu);
u=0;
prob = exp(-gammaln(u+nu+1)-gammaln(u+1) -alpha(dex));
while u<=2*m0+1
    
    idex= prob./PMF >= rand(length(dex),1);
    x(dex(idex))=u;
    dex(idex)=[];
    PMF(idex)=[];
    if isempty(dex)
        break
    else
        PMF=PMF-prob(~idex);
    end
    u=u+1;
    prob = exp(u.*log(alpha(dex).^2/4)-gammaln(u+nu+1)-gammaln(u+1) -alpha(dex));
    %u
    %prob
    %PMF
end

dex = [find(m>m0);dex];

if Algorithm==1
    %% The Algorithm in page 7
%     pm=zeros(size(alpha));
%     w = zeros(size(alpha));
%     temp =zeros(size(alpha));
%     
%     pm(dex)=pmf_bessel(alpha(dex),nu,m(dex));
%     w(dex)=1+pm(dex)/2;
    
    ppm=pmf_bessel(alpha(dex),nu,m(dex));
    ww=1+ppm/2;
    temp = gammaln(m(dex)+nu+1)+gammaln(m(dex)+1);
    
    
  %  dex=find(x==-1);
    
    while ~isempty(dex)
        %  dex=find(x==-1);
        Y=zeros(length(dex),1);
        U = rand(length(dex),1);
        %     W = rand(length(dex),1);
        S = (rand(length(dex),1)>0.5)*2 - 1;
        %ww = w(dex);
        %ppm=pm(dex);
        mm=m(dex);
        if length(nu)>1
            nunu=nu(dex);
        end
        alpha_alpha=alpha(dex);
        %idex = (U<=(1+2./q(dex))./(1+4./q(dex))));
        idex = U<= ww./(1+ww);
        
        %         size(V)
        %         size(Y)
        %         size(qq(idex));
        if nnz(idex)>0
            V = rand(nnz(idex),1);
            Y(idex)= V.*ww(idex)./ppm(idex);
        end
        if nnz(~idex)>0
            %Y(~idex)=1/2+1./qq(~idex)+randg(ones(nnz(~idex),1))./qq(~idex);
            Y(~idex)=(ww(~idex)+exp(1))./ppm(~idex);
        end
        X = S.*round(Y);
        
        idex = m(dex)+X>=0;
        if nnz(idex)>0
            if length(nu)>1
                nu1=nunu(idex);
            else
                nu1=nu;
            end
            idex(idex) = log(rand(nnz(idex),1))+min(0,(ww(idex)-ppm(idex).*Y(idex))) <=...
                (2*X(idex).*log(alpha_alpha(idex)/2)+ temp(idex) -gammaln(mm(idex)+X(idex)+nu1+1)-gammaln(mm(idex)+X(idex)+1));
            x(dex(idex))=X(idex)+mm(idex);
            dex(idex)=[];
            ppm(idex)=[];
            ww(idex)=[];
            temp(idex)=[];
        end
        %
        %        % idex = (W<0 | X+m(dex) (log(abs(W))+min(0,1+qq/2-qq.*Y)<= 2*X.*log(alpha(dex)/2)+gammaln(m(dex)+nu+1)+gammaln(m(dex)+1)-gammaln(m(dex)+X+nu+1)-gammaln(m(dex)+X+1)));
        %         idex = W.*min(1,exp(1+qq/2-qq.*Y)) <= ...
        %               (alpha(dex)/2).^(2*X).*gamma(m(dex)+nu+1)./gamma(m(dex)+X+nu+1).*gamma(m(dex)+1)./gamma(m(dex)+X+1);
        
        
    end
%% The Algorithm in page 8
elseif Algorithm==2
    
    mu = alpha/2.*besseli(nu+1,alpha,1)./besseli(nu,alpha,1);
    var = alpha.^2/4.*besseli(nu+2,alpha,1)./besseli(nu,alpha,1) + mu - mu.^2;
    
    sigma = sqrt(var + (m-mu).^2);
    q= min(1./(sigma*sqrt(648)),1/3);
    
    
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
elseif Algorithm==3
    
    
    alpha1=alpha(dex);
    mu = alpha1/2.*besseli(nu+1,alpha1,1)./besseli(nu,alpha1,1);
    sigma = sqrt(alpha1.^2/4.*besseli(nu+2,alpha1,1)./besseli(nu,alpha1,1) + mu - mu.^2);
    temp =rand(size(mu));
    x(dex) = round(mu + sigma.*norminv(temp+(1-temp).*normcdf(-mu./sigma)));
end

