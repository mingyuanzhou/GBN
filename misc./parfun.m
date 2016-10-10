function varargout=parfun(fun,NumBlock,slice,output_concatenate,output_add, varargin)

TotalLen = size(varargin{slice(1)},2);
NumBlock = min(NumBlock,TotalLen);
NumOutput = nargout;

varargout = cell(1,NumOutput);
BlockSize = ceil(TotalLen/NumBlock);

M = cell(1,NumBlock);
%varargin_block= varargin;
varargin_block= cell(1,NumBlock);
for block=1:NumBlock
    blockdex = (block-1)*BlockSize+1 : min(block*BlockSize,TotalLen);
    varargin_block{block} = varargin;
    for i=slice
        varargin_block{block}{i}=varargin{i}(:,blockdex);
    end 
end

parfor block=1:NumBlock
    out = cell(1,NumOutput);
    [out{1:NumOutput}] = feval(fun,varargin_block{block}{1:end});
    M{block}=out;
end

M = [M{1:end}];
for i=output_add
    varargout{i} = sum(reshape([M{i:NumOutput:end}],size(M{i},1),size(M{i},2),[]),3);
%     varargout{i} = M{i};
%     for j=(i+NumOutput):NumOutput:length(M)
%         varargout{i} = varargout{i} + M{j};
%     end
end

for i=1:output_concatenate
    varargout{i}= [M{i:NumOutput:end}];
end
