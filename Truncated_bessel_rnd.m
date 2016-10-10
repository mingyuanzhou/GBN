function x = Truncated_bessel_rnd(varargin)
% Coded by Mingyuan Zhou, March, 2015
% Modified in Oct, 2015
if nargin==1
    alpha = varargin{1};
    nu=-1;
    x=zeros(size(alpha));
    dex=alpha>0;
    
    alpha = alpha(dex);
    Mode = max(fix((sqrt(alpha.^2+nu.^2)-nu)/2),1);
    PMF = besseli(nu,alpha,1).*(alpha/2).^(-nu);
    
    x(dex) = Rand_Truncated_PMF_bessel(PMF,Mode,alpha);
else
    c=varargin{1};
    Xmask=varargin{2};
    Phi=varargin{3};
    Theta=varargin{4};
    
    iijj = Xmask>0;
    
    alpha = Phi*Theta;
    alpha = full(2*sqrt(c*Xmask(iijj).*alpha(iijj)));
    
    nu=-1;
    dex=alpha>0;
    alpha = alpha(dex);
    Mode = max(fix((sqrt(alpha.^2+nu.^2)-nu)/2),1);
    PMF = besseli(nu,alpha,1).*(alpha/2).^(-nu);
    x(dex) = Rand_Truncated_PMF_bessel(PMF,Mode,alpha);
end