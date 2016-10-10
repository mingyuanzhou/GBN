function [Xtrain,Xtest,WS,DS,WordTrainS,DocTrainS] = PartitionX_v1(X,TrainPercentage)
[wi,di,cc]=find(X);
wdi=find(X>0);
ii=zeros(sum(cc),1);
jj=zeros(sum(cc),1);
kk=zeros(sum(cc),1);
count=0;
for i=1:length(wdi)
    ii(count+1:count+cc(i))=wi(i);
    jj(count+1:count+cc(i))=di(i);
    count=count+cc(i);
end
WS=ii;
DS=jj;
[P,N] = size(X);
TestLen = N;
dex = 1:N;
WordTrainS = true(size(DS));
DocTrainS = true(size(DS));
for ii = 1:TestLen
    blockdex = find(DS==dex(ii));
    if ~isempty(blockdex)
        DocTrainS(blockdex)=false;
        aa = randperm(length(blockdex));    
        TestStartLocation = max(round(length(blockdex)*TrainPercentage)+1,2);
        if TestStartLocation==1 
            TestStartLocation=2;        
        end
        WordTrainS(blockdex(aa(TestStartLocation:end))) = false;       
       % WordTrainS(blockdex(aa(ceil(length(blockdex)*0.2):end))) = false;       
    end   
end
Xtrain = sparse(WS(WordTrainS),DS(WordTrainS),1,P,N);
Xtest = sparse(WS(~WordTrainS),DS(~WordTrainS),1,P,N);
DocTrainS = WordTrainS;
