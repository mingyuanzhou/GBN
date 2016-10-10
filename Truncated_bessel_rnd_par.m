function x = Truncated_bessel_rnd_par(alpha,NumBlock)

% Coded by Mingyuan Zhou, March, 2015
% Modified in Oct, 2015
if nargin<2
    NumBlock=12;
end

nu=-1;
x=zeros(1,length(alpha));
dex=alpha>0;
if nnz(dex)<numel(dex)
    alpha = alpha(dex);
end
alpha = alpha(:)';

TotalLen = length(alpha);
NumBlock = min(NumBlock,TotalLen);
BlockSize = ceil(TotalLen/NumBlock);
%Dim2 = BlockSize*ones(1,NumBlock);
%Dim2(end) = TotalLen - sum(Dim2(1:end-1));
alpha1 = cell(1,NumBlock);
for block=1:NumBlock
    blockdex = (block-1)*BlockSize+1 : min(block*BlockSize,TotalLen);
    alpha1{block} = alpha(blockdex);
end
%alpha = mat2cell(alpha,1,Dim2);
parfor block=1:NumBlock
    Mode = max(fix((sqrt(alpha1{block}.^2+nu.^2)-nu)/2),1);
    PMF = besseli(nu,alpha1{block},1).*(alpha1{block}/2).^(-nu);
    alpha1{block}= Rand_Truncated_PMF_bessel(PMF,Mode,alpha1{block});
end
x(dex)=[alpha1{:}];