parpool(2)
%%======== 20news-bydate   ============
SystemInRuningLinux=1;
DataType    =   'Count';
dataname    =   '20newsBydate';
IsBinaryClassificaiton = false;
if SystemInRuningLinux
    addpath('data/20news-bydate')
    addpath('liblinear-2.1/matlab/');
else
    addpath('E:\0Research\Data\20news-bydate')
end
load train.data     ;           load test.data  ;
test(:,1)=max(train(:,1))+test(:,1);
train_test = [train;test];
X_all =sparse(train_test(:,2),train_test(:,1),train_test(:,3));
load train.label    ;          load test.label  ;
prepar.trindx   =   1 : length(train)  ;
prepar.teindx   =    length(train)  + (1: length(test))  ;
prepar.Y        =   [train;test];
if IsBinaryClassificaiton
    if dataset==1
        %alt.atheism 1 vs talk.religion.misc 20
        train(train>1&train<20)=[];
        test(test>1&test<20)=[];
        prepar.trindx   =   1 : length(train)  ;
        prepar.teindx   =    length(train)  + (1: length(test))  ;
        dex = (prepar.Y>1) & (prepar.Y<20)  ;
        dataname    =   '20newsAtheismVsReligion';
    elseif dataset==2
        %talk talk.politics.guns 17 vs talk.politics.mideast 18
        train(train<17|train>18)=[];
        test(test<17|test>18)=[];
        prepar.trindx   =   1 : length(train)  ;
        prepar.teindx   =    length(train)  + (1: length(test))  ;
        dex = (prepar.Y<17) | (prepar.Y>18)  ;
        dataname    =   '20newsGunsVsMideast';
    elseif dataset==3
        % % %comp comp.sys.ibm.pc.hardware 4 vs comp.sys.mac.hardware 5
        train(train<4|train>5)=[];
        test(test<4|test>5)=[];
        prepar.trindx   =   1 : length(train)  ;
        prepar.teindx   =    length(train)  + (1: length(test))  ;
        dex = (prepar.Y<4) | (prepar.Y>5)  ;
        dataname    =   '20newsPcVsMac';
    elseif dataset==4
        %sci sci.electronics 13 vs sci.med 14
        train(train<13|train>14)=[];
        test(test<13|test>14)=[];
        prepar.trindx   =   1 : length(train)  ;
        prepar.teindx   =    length(train)  + (1: length(test))  ;
        dex = (prepar.Y<13) | (prepar.Y>14)  ;
        dataname    =   '20newsElecVsMed';
    end
    X_all(:,dex)=[];
    prepar.Y(dex)=[];
end
X_all_save=X_all;
ACC=zeros(1,8);
clear train test dex stopwords tmp train_test ;

for iter=1:8
    switch iter
        case {1,2}
            %Remove stopwords, Top2000
            fid=fopen('stop-word-list.txt');
            stopwords=textscan(fid, '%s');
            stopwords=stopwords{1};
            fclose(fid);
            fid = fopen('vocabulary.txt');
            WO = textscan(fid, '%s');
            fclose(fid);
            WO = WO{1};
            dex = true(length(WO),1);
            for i=1:length(WO)
                dex(i)=1-nnz(strcmp(WO(i),stopwords));
            end
            WO=WO(dex);
            X_all=X_all_save(dex,:);
            [~,dex]=sort(sum(X_all,2),'descend');
            X_all = X_all(dex(1:2000),:);
            if iter==2
                X_all = bsxfun(@rdivide, X_all,max(sum(X_all,1),realmin));
            end
        case {3,4}
            %Do not remove stopwords, Top2000
            X_all=X_all_save(dex,:);
            [~,dex]=sort(sum(X_all,2),'descend');
            X_all = X_all(dex(1:2000),:);
            if iter==4
                X_all = bsxfun(@rdivide, X_all,max(sum(X_all,1),realmin));
            end
        case {5,6}
            %No preprocessing, all terms
            X_all = X_all_save;
            if iter==6
                X_all = bsxfun(@rdivide, X_all,max(sum(X_all,1),realmin));
            end
        case {7,8}
            %Remove stopwords and terms appear less than five times
            fid=fopen('stop-word-list.txt');
            stopwords=textscan(fid, '%s');
            stopwords=stopwords{1};
            fclose(fid);
            fid = fopen('vocabulary.txt');
            WO = textscan(fid, '%s');
            fclose(fid);
            WO = WO{1};
            dex = true(length(WO),1);
            for i=1:length(WO)
                dex(i)=1-nnz(strcmp(WO(i),stopwords));
            end
            WO=WO(dex);
            X_all=X_all_save(dex,:);
            tmp     =   (sum(X_all,2)>=5)    ; %words appear at least 5 times
            %tmp     =   (sum(X_all>0,2)>=5)    ;%words appear in at lest 5 documents
            WO = WO(tmp);
            X_all=X_all(tmp,:);
            prepar.WO        =   WO  ;
            prepar.stopwords        =   stopwords  ;
            if iter==8
                X_all = bsxfun(@rdivide, X_all,max(sum(X_all,1),realmin));
            end
    end
    
    %Choose L2 regularized multiclass logisitic regression or svm for classication
    
    
    TrainX      =   X_all(:,prepar.trindx) + 0*eps;  TrainY  =   prepar.Y(prepar.trindx)     ;
    TestX       =   X_all(:,prepar.teindx)+ 0*eps;  TestY  =   prepar.Y(prepar.teindx)     ;
    %%==================================== NO PreprocessBfTest ===================================
    %     tic;
    %     BestModel = train(double(TrainY), sparse(TrainX'), ['-C -s 0']);
    %     CVtime  =   toc
    %     tic ;
    %     option = ['-s 0 -c ', num2str(BestModel(1)), ' -q'];
    %     model = train(double(TrainY), sparse(TrainX'), option);
    %     Traintime   =   toc
    %     tic
    %     [predicted_label_LR, Accuracy_LR, prob_estimates_LR] = predict(double(TestY), sparse(TestX'), model, ' -b 1');
    %     Testtime = toc
    %%==================================== NO PreprocessBfTest ===================================
    CC=2.^(-10:1:15);
    ModelOut=zeros(1,length(CC));
    parfor ij=1:length(CC)
        ModelOut(ij) = train(double(TrainY), sparse(TrainX'), ['-s 0 -c ', num2str(CC(ij)), ' -v 5 -q ']);
    end
    [~,maxdex]=max(ModelOut);
    option = ['-s 0 -c ', num2str(CC(maxdex)), ' -q'];
    model = train(double(TrainY), sparse(TrainX'), option);
    [predicted_label_LR, Accuracy_LR, prob_estimates_LR] = predict(double(TestY), sparse(TestX'), model, ' -b 1');
    nameaaa     =   ['Main_PretrainTrim_',dataname,'_feature',num2str(iter),'.mat']     ;
    
    save(nameaaa,'predicted_label_LR','Accuracy_LR','prob_estimates_LR')   ;
    
    
    ACC(iter)=Accuracy_LR(1)
    
end


