%% ======== 20news-bydate = = ==========
Para.DataType = 'Positive';
dataname = 'MNIST_Positive';
addpath('liblinear-2.1/matlab/'); 
%first download and compile the liblinear software package, whose logistic
%regression model will be applied to the inferred feature vectors
addpath('data/');

%%

%Run the Demo Code with TrainSize = 60000, TestSize = 10000, K0 in
%[50,100,200,400], and trial in [1,2,3,4,5] to reproduce the results used
%to plot Figures 15 and 22-24 of "Agumentable Gamma Belief Networks"


%%

%% read the 20newsgroups data

load('data/mnist_gray.mat')     ;   % Load Data

TrainSize = 1000;
TestSize = 1000;

%TrainSize = 60000;
%TestSize = 10000;

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
Para.train_idx   =   1:length(Xlabel)    ;
Para.test_idx   =   length(Xlabel) + (1:length(Xtestlabel))     ;
Para.Y        =   [Xlabel(:);Xtestlabel(:)]   ;
Para.Y(Para.Y == 0)     =   10  ;
%         MNISTPreprocessing  ;
clear X Xlabel  ndx Xtest Xtestlabel;




%% set model parameters
trial = 1; %select from [1,2,3,4,5]
K0 = 50; %select from [50,100,200,400]
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

%Para.c = 1;
Para.c = ones(1,size(X_all,2));

%% GBN Training

% Whether to perform testing inside the code for training
Para.IsTestInside = 1;

% =================== Settings = ==================
Para.IsDisplay = 1;          Para.FigureGap = 100;
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
