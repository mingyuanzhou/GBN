function R=randg_par(R,NumBlock)
if nargin<3
    NumBlock=12;
end
TotalLen = size(R,2);
NumBlock = min(NumBlock,TotalLen);
BlockSize = fix(TotalLen/NumBlock);
Dim2 = BlockSize*ones(1,NumBlock);
for i=1:mod(TotalLen,NumBlock)
    Dim2(i)=Dim2(i)+1;
end
R = mat2cell(R,size(R,1),Dim2);

parfor block=1:NumBlock
    %blockdex = (block-1)*BlockSize+1 : min(block*BlockSize,TotalLen);
    R{block} = randg(R{block});
end
R=[R{:}];