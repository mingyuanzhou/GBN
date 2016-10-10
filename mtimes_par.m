function B=mtimes_par(A,B,NumBlock)
if nargin<3
    NumBlock=12;
end
%C=zeros(size(A,1),size(B,2));
TotalLen = size(B,2);
NumBlock = min(NumBlock,TotalLen);

BlockSize = fix(TotalLen/NumBlock);

% C=cell(1,NumBlock);
% for block=1:NumBlock
%     blockdex = (block-1)*BlockSize+1 : min(block*BlockSize,TotalLen);
%     C{block} = B(:,blockdex);
% end

Dim2 = BlockSize*ones(1,NumBlock);
for i=1:mod(TotalLen,NumBlock)
    Dim2(i)=Dim2(i)+1;
end
%Dim2(end) = TotalLen - sum(Dim2(1:end-1));
B = mat2cell(B,size(B,1),Dim2);

parfor block=1:NumBlock
    %blockdex = (block-1)*BlockSize+1 : min(block*BlockSize,TotalLen);
    B{block} = mtimes(A,B{block});
end
B=[B{:}];