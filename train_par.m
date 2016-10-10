function CC=train_par(TrainY, sparse_TrainX_Transpose, CC,NumBlock)
if nargin<4
    NumBlock=4;
end
%
TotalLen = length(CC);
NumBlock = min(NumBlock,TotalLen);
BlockSize = fix(TotalLen/NumBlock);
Dim2 = BlockSize*ones(1,NumBlock);
for i=1:mod(TotalLen,NumBlock)
    Dim2(i)=Dim2(i)+1;
end
CC = mat2cell(CC(:)',1,Dim2);

parfor block=1:NumBlock
    %blockdex = (block-1)*BlockSize+1 : min(block*BlockSize,TotalLen);
    for i=1:length(CC{block})
        CC{block}(i)=train(double(TrainY), sparse_TrainX_Transpose, ['-s 0 -c ', num2str(CC{block}(i)), ' -v 5 -q ']);
    end
end
CC=[CC{:}];