%% ======== 20news-bydate = = ==========
Para.DataType = 'Count';
dataname = '20newsBydate';
addpath('liblinear-2.1/matlab/'); 
%first download and compile the liblinear software package, whose logistic
%regression model will be applied to the inferred feature vectors
addpath('data/');

dataset = 5; %[0,5,3,4]

%dataset = 0; Multiclass classfication, 11,269 documents for training, 7505
%for testing; Vocabulary size = 33420 after removing stop-words and
%infrequent words

%dataset = 5; Multiclass classfication, 11,269 documents for training, 7505
%for testing; Vocabulary size = 2000, including the most frequent 2000 words
%after removing stop-words

%dataset = 1, 2, 3, or 4; binary classificaiton of two news groups removing
%stop-words and infrequent words

%%

%Run the Demo Code with datset in [3,4], K0 in [25,50,100,200,400,600,800],
%and trial in [1,2,..,12] to reproduce the results used to plot Figure 9 of
%"Agumentable Gamma Belief Networks"

%Run the Demo Code with datset = 0, K0 in [50,100,200,400,600,800], and
%trial in [1,2,3,4,5] to reproduce the results used to plot Figure 10 of
%"Agumentable Gamma Belief Networks"

%Run the Demo Code with datset = 5, K0 in [32,64,128,256,512], and trial in
%[1,2,3,4,5] to reproduce the results used to plot Figure 11 of
%"Agumentable Gamma Belief Networks"

%Run the Demo code with dataset = 0, eta = 0.1, K0 = 800, T=5 to produce
%the results used to plot Figures 3-7 and 17-21 of "Agumentable Gamma
%Belief Networks"


%%

%% read the 20newsgroups data
if dataset == 0 || dataset == 5
    IsBinaryClassificaiton = false;
else
    IsBinaryClassificaiton = true;
end

addpath('data/20news-bydate')
load train.data;           
load test.data;
test(:,1)=max(train(:,1))+test(:,1);
train_test = [train;test];
X_all =sparse(train_test(:,2),train_test(:,1),train_test(:,3));
load train.label;          
load test.label;
Para.train_idx = 1 : length(train);
Para.test_idx = length(train)  + (1: length(test));
Para.Y =  [train;test];
if IsBinaryClassificaiton
    if dataset==1
        %alt.atheism 1 vs talk.religion.misc 20
        train(train>1&train<20)=[];
        test(test>1&test<20)=[];
        Para.train_idx = 1 : length(train);
        Para.test_idx = length(train)  + (1: length(test));
        dex = (Para.Y>1) & (Para.Y<20);
        dataname = '20newsAtheismVsReligion';
    elseif dataset==2
        %talk talk.politics.guns 17 vs talk.politics.mideast 18
        train(train<17|train>18)=[];
        test(test<17|test>18)=[];
        Para.train_idx = 1 : length(train);
        Para.test_idx = length(train)  + (1: length(test));
        dex = (Para.Y<17) | (Para.Y>18);
        dataname = '20newsGunsVsMideast';
    elseif dataset==3
        % % %comp comp.sys.ibm.pc.hardware 4 vs comp.sys.mac.hardware 5
        train(train<4|train>5)=[];
        test(test<4|test>5)=[];
        Para.train_idx = 1 : length(train);
        Para.test_idx = length(train)  + (1: length(test));
        dex = (Para.Y<4) | (Para.Y>5);
        dataname = '20newsPcVsMac';
    elseif dataset==4
        %sci sci.electronics 13 vs sci.med 14
        train(train<13|train>14)=[];
        test(test<13|test>14)=[];
        Para.train_idx = 1 : length(train);
        Para.test_idx = length(train)  + (1: length(test));
        dex = (Para.Y<13) | (Para.Y>14);
        dataname = '20newsElecVsMed';
    end
    X_all(:,dex)=[];
    Para.Y(dex)=[];
end

%% remove stop-words
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

%% remove infrequent words
tmp = (sum(X_all,2)>=5); %words appear at least 5 times
%tmp = (sum(X_all>0,2)>=5);%words appear in at lest 5 documents
WO = WO(tmp);
X_all=X_all(tmp,:);

%% choose top 2000 words if dataset = 5
if dataset==5
    [~,dex]=sort(sum(X_all,2),'descend');
    X_all = X_all(dex(1:2000),:);
    WO = WO(dex(1:2000));
    dataname = '20newsBydateTop2000';
end

Para.WO =  WO;
Para.stopwords =  stopwords;

clear train test dex stopwords tmp train_test WO;

%% set model parameters



if dataset==0 || dataset==5
    trial = 1; %select from [1,2,3,4,5]
    if dataset==0
        K0 = 400; %select from [50,100,200,400,600]
    else
        K0 = 128; %select from [32,64,128,256,512]
    end
    eta = 0.05; %[0.01,0.05,0.1];
    T = 5;
    if (dataset==0 && K0>400) || dataset == 5 %reduce the number of iterations to save computation
        Para.TrainBurnin = [1000,500*ones(1,T-1)];
    else
        Para.TrainBurnin = [1000,1000*ones(1,T-1)];
    end
    Para.TrainCollection = [500,500*ones(1,T-1)];
else
    trial = 1; %select from [1,2,3,4,5,6,7,8,9,10,11,12]
    K0 = 200; %select from [25,50,100,200,400,600,800]
    eta = 0.01; %[0.01,0.05,0.1];
    T = 8;
    Para.TrainBurnin = [1000,1000*ones(1,T-1)];       
    Para.TrainCollection = [1000,1000*ones(1,T-1)];
end

Para.TestBurnin = 500;       Para.TestCollection = 500;

K = ones(1,T)*K0;  

Para.percentage=101; %all word tokens will be used
Para.dataname = dataname;
Para.CollectionStep=5;

%% GBN Training

% Whether to perform testing inside the code for training
Para.IsTestInside = 1;

% =================== Settings = ==================
Para.IsDisplay = 0;          Para.FigureGap = 100;
Para.SaveDataGap = 100;
Para.ParallelProcessing = 0;  Para.NumBlockParallel = 8;

Para.AddBiasTerm = false; %wheter to add bias term to logistic regression
%Para.AddBiasTerm  = true;

%% =================== Training And Testing ===================
rng(trial,'twister');

Para.name_save = [dataname,'_K0_',num2str(K0),'_T_',num2str(T),'_eta',...
    num2str(round(eta*1000)),'_Trial',num2str(trial),'.mat'];

tic
[ParaGlobal,ParaLocal,Accuracy_all] = GBN_Layerwise(X_all,K,T,eta,Para);
TimeAll = toc;

if trial==1
    save(['results1/ZGlobal',Para.name_save],'ParaGlobal','Accuracy_all','Para','K','T','eta');
    %save(['results1/ZLocal',nameccc],'ParaLocal');
end
save(['results1/',Para.name_save],'Accuracy_all');

