%% ======== 20news-bydate = = ==========
Para.DataType = 'Binary';
dataname = '20newsBydate_Binary';
addpath('liblinear-2.1/matlab/'); 
%first download and compile the liblinear software package, whose logistic
%regression model will be applied to the inferred feature vectors
addpath('data/');

dataset = 0; 

%%

%Run the Demo Code with K0 in [50,100,200,400,600,800], and trial in
%[1,2,3,4,5] to reproduce the results used to plot Figure 14 of
%"Agumentable Gamma Belief Networks"


%%

%% read the 20newsgroups data
IsBinaryClassificaiton = false;
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

X_all = double(X_all>=1);

Para.WO =  WO;
Para.stopwords =  stopwords;

clear train test dex stopwords tmp train_test WO;

%% set model parameters
trial = 1; %select from [1,2,3,4,5]
K0 = 100; %select from [50,100,200,400,600]
eta = 0.05;
T = 5;
if (dataset==0 && K0>400) %reduce the number of iterations to save computation
    Para.TrainBurnin = [1000,500*ones(1,T-1)];
else
    Para.TrainBurnin = [1000,1000*ones(1,T-1)];
end
Para.TrainCollection = [500,500*ones(1,T-1)];
Para.TestBurnin = 500;       Para.TestCollection = 500;

K = ones(1,T)*K0;  

Para.percentage=101; %all word tokens will be used
Para.dataname = dataname;
Para.CollectionStep=5;

Para.c=1;

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

