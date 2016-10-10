%% MNISTPreprocessing
% % % % % % % % MinibatchAll   =   floor(TrainSize/MinibatchSize)     ;
% % % % % % % % iset    =   mod( (1:MinibatchAll)-1 , 10 )  ;
% % % % % % % % ndx = []    ;
% % % % % % % % for i = iset
% % % % % % % %     tmp = find(train_label==i);
% % % % % % % %     ndx = [ndx; tmp(1:MinibatchSize)];
% % % % % % % % end
% % % % % % % % X       =   train_mnist(:,ndx);
% % % % % % % % Xlabel  =   train_label(ndx);
N   =   size(X_all,2)   ;
ndx     =   zeros(1,N)  ;
MinibatchAllMax     =   100000;
MinibatchEachNum    =   (MinibatchSize/10)  ;
for i = 0:9
    tmp = find(prepar.Y == i);
    MinibatchAllMax     =   min(MinibatchAllMax,floor(length(tmp) / MinibatchEachNum))    ;
    for ii  =   1:MinibatchAllMax
        ndx((ii-1)*MinibatchSize + i*MinibatchEachNum+ (1:MinibatchEachNum))    ...
            =  tmp((ii-1)*MinibatchEachNum +(1:MinibatchEachNum))   ;
    end
end
ndx     =   ndx(1:(MinibatchAllMax*MinibatchSize))  ;
X_all_for_Training   =   X_all(:,ndx)   ;



