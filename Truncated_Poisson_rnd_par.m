function x = Truncated_Poisson_rnd_par(alpha,NumBlock)

% Coded by Mingyuan Zhou, March, 2015
% Modified in Oct, 2015
if nargin<2
    NumBlock=12;
end
alpha = alpha(:)';
TotalLen = length(alpha);
NumBlock = min(NumBlock,TotalLen);
BlockSize = fix(TotalLen/NumBlock);
Dim2 = BlockSize*ones(1,NumBlock);
for i=1:mod(TotalLen,NumBlock)
    Dim2(i)=Dim2(i)+1;
end
% alpha1 = cell(1,NumBlock);
% for block=1:NumBlock
%     blockdex = (block-1)*BlockSize+1 : min(block*BlockSize,TotalLen);
%     alpha1{block} = alpha(blockdex);
% end
alpha = mat2cell(alpha,1,Dim2);
parfor block=1:NumBlock
    alpha{block}= (truncated_Poisson_rnd(alpha{block}))';
end
x=[alpha{:}];