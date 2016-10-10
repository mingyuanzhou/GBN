function [ZSDS,WSZS] = Multrnd_Matrix_mex_fast_v1_par(Xt,Phi,Theta,NumBlock)
if nargin<4
    NumBlock=12;
end
TotalLen = size(Theta,2);
NumBlock = min(NumBlock,TotalLen);
BlockSize = fix(TotalLen/NumBlock);
ZSDS=cell(1,NumBlock);

% Xt_temp = cell(1,NumBlock);
% for block=1:NumBlock
%     blockdex = (block-1)*BlockSize+1 : min(block*BlockSize,TotalLen);
%     C{block} = Theta(:,blockdex);
%     Xt_temp{block} = Xt(:,blockdex);
% end
% A = 0;

Dim2 = BlockSize*ones(1,NumBlock);
for i=1:mod(TotalLen,NumBlock)
    Dim2(i)=Dim2(i)+1;
end

Theta = mat2cell(Theta,size(Theta,1),Dim2);
Xt = mat2cell(Xt,size(Xt,1),Dim2);

WSZS=0;
parfor block=1:NumBlock
    %blockdex = (block-1)*BlockSize+1 : min(block*BlockSize,TotalLen);
    [ZSDS{block},temp] = Multrnd_Matrix_mex_fast_v1(Xt{block},Phi,Theta{block});
    WSZS = WSZS + temp;
end
ZSDS=[ZSDS{:}];