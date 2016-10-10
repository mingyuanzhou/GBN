switch ToBeAnalized
    case 1  
        %%============================ MNIST ======================================
        if SystemInRuningLinux
            load('data/mnist_gray.mat')     ;   % Load Data  
        else
            load data\mnist_gray.mat     ;   % Load Data   
        end
%         [~,IX]  =   sort(train_label,'ascend')  ;   train_mnist     =   train_mnist(:,IX)   ;    train_label  =   train_label(IX)  ;
%         [~,IX]  =   sort(test_label,'ascend')   ;   test_mnist      =   test_mnist(:,IX)   ;    test_label  =   test_label(IX)  ;
%         TrainSize      =   1e3   ;       TestSize        =   0.5e3   ;
%         MinibatchSize   =   250 ;
%         K   =   [128,64]     ;
        %======== Prepare Train Data  ===========
        ndx = []; mtrain = TrainSize / 10;
        if mtrain < 6000
            for ii = 0:1:9
                tmp = find(train_label==ii);
                ndx = [ndx; tmp(1:mtrain)];
            end
        else
            ndx = [1:60000];
        end
        X       =   train_mnist(:,ndx);
        Xlabel  =   train_label(ndx);
        rng(0,'twister');
        DataType    =   'Positive';
        dataname    =   'MNIST';
        %======== Prepare Test Data   ===========
        ndx = []; mtest = TestSize/10  ;
        if mtest < 1000
            for ii = 0:1:9
                tmp = find(test_label==ii);
                ndx = [ndx; tmp(1:mtest)];
            end
        else
            ndx=[1:10000];
        end
        Xtest       =   test_mnist(:,ndx);
        Xtestlabel  =   test_label(ndx);
        %======== Combine Train and Test ===========
        clear train_mnist train_label test_mnist test_label ;
        X_all   =   [X,Xtest]   ;
        prepar.trindx   =   1:length(Xlabel)    ;
        prepar.teindx   =   length(Xlabel) + (1:length(Xtestlabel))     ;
        prepar.Y        =   [Xlabel(:);Xtestlabel(:)]   ;
        prepar.Y(prepar.Y == 0)     =   10  ; 
%       MNISTPreprocessing  ;
        clear X Xlabel  ndx Xtest Xtestlabel;

    case 9 
        %%======== 20news-bydate   ============
        DataType    =   'Count';
        dataname    =   '20newsBydate';
        if SystemInRuningLinux
                addpath('data/20news-bydate')
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
        if 1
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
            X_all=X_all(dex,:);
        end
        tmp     =   (sum(X_all,2)>=5)    ; %words appear at least 5 times
        %tmp     =   (sum(X_all>0,2)>=5)    ;%words appear in at lest 5 documents
        WO = WO(tmp);
        X_all=X_all(tmp,:);
        
        if dataset==5
            [~,dex]=sort(sum(X_all,2),'descend');
            X_all = X_all(dex(1:2000),:);
            WO = WO(dex(1:2000));
            dataname    =   '20newsBydateTop2000';
        end
        
        prepar.WO        =   WO  ;
        prepar.stopwords        =   stopwords  ;
        
        clear train test dex stopwords tmp train_test WO ;
        
    
    otherwise
        error('Wrong "ToBeAnalized"')
end