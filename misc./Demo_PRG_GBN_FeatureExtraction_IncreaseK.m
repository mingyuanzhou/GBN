% clear all
i1=1;
i2=1;
i3=1;
i4=1;
%parpool(3)
% 


K0_all = [25,50,100,200,400,600,800,1000]
eta_all=[0.01,0.05,0.1];
K0 = K0_all(i1) %[50,100,200,400,600]
eta0=eta_all(i2)
trial = i3

%
%Settings.TrimStrategy  =   6 %increase by four times at most for T>=2

Settings.TrimStrategy  =   7 %adaptive tuncation for T>=2, bo upbound


Settings.PreTrainYes = 0

T   =   5

%eta0 = [eta0,0.1*ones(1,T-1)];
Settings.TrainBurnin  =   [750,750*ones(1,T-1)]/1     ;
Settings.TrainCollection  =  [250,250*ones(1,T-1)]/1;
Settings.TrainSampleSpace   =   1   ;

Settings.TestBurnin  =   500/1     ;
Settings.TestCollection  =   500/1     ;
Settings.TestSampleSpace   =   1   ;


Settings.ParallelProcessing      =   i4

POOL.NumWorkers = 24;
Settings.NumBlockParallel   =   POOL.NumWorkers*1;


%POOL = gcp; % If no pool, do not create new one.




Settings.IsDisplay   =  1   ;          Settings.FigureGap       =   10; %1e5   ;


%clear,clc,close all;
SystemInRuningLinux  =   1   ;   %   1:linux       ;;;     0 : windows



% if SystemInRuningLinux
addpath('liblinear-2.1/matlab/');
addpath('data/');
% else
%   addpath('C:\liblinear-2.1\matlab\');
% end

%% Data
 

Settings.maxIter_preTrain    =   100     ;
%Settings.TrainBurnin  =   [1500,1000*ones(1,T-1)]     ;       Settings.TrainCollection  =   [1000,500*ones(1,T-1)];            Settings.TrainSampleSpace   =   1   ;

%if SystemInRuningLinux
    load('data/mnist_gray.mat')     ;   % Load Data
%else
%    load E:\0Research\Data\MNIST\mnist_gray.mat     ;   % Load Data
%end
%         [~,IX]  =   sort(train_label,'ascend')  ;   train_mnist     =   train_mnist(:,IX)   ;    train_label  =   train_label(IX)  ;
%         [~,IX]  =   sort(test_label,'ascend')   ;   test_mnist      =   test_mnist(:,IX)   ;    test_label  =   test_label(IX)  ;
%         TrainSize      =   1e3   ;       TestSize        =   0.5e3   ;
%         MinibatchSize   =   250 ;
%         K   =   [128,64]     ;
TrainSize = 60000/1;
TestSize = 10000/1;
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
%         MNISTPreprocessing  ;
clear X Xlabel  ndx Xtest Xtestlabel;


Settings.c = 1;

%Settings.c = ones(size(X_all,1),1);

Settings.c = ones(1,size(X_all,2));

%% GBN Training

SuPara.Epsilon     =   0    ;
SuPara.ac  =   1   ;   SuPara.bc      =   1   ;
SuPara.a0pj    =   0.01   ;   SuPara.b0pj    =   0.01   ;
%SuPara.e0cj    =   5   ;   SuPara.f0cj    =   5   ;
SuPara.e0cj    =   1   ;   SuPara.f0cj    =   1   ;
SuPara.e0c0    =   1    ;   SuPara.f0c0    =   1    ;
SuPara.a0gamma  =   1   ;   SuPara.b0gamma      =   1   ;

super1or2   =   2   ;   SuPara.super1or2    =   super1or2   ;
switch super1or2
    case 1
        %%=======================   SuPara 1   ================================
        SuPara.eta  =   0.01     ;
    case 2
        %%=======================   SuPara 2   ================================
        SuPara.eta  =   eta0     ;
    case 3
        %%=======================   SuPara 3   ================================
        SuPara.eta  =   2.^(-7:1:-3)   ;
    case 4
        %%=======================   SuPara 3   ================================
        SuPara.eta  =   2.^(-8 + 0.2*(1:T))   ;
end
%%===================  Settings  ===================
Settings.SaveDataGap       =   200   ;   Settings.TestGap    =   1e5  ;
Settings.IsInferEta     =   zeros(1,T)   ;     Settings.IsTestInside   =    1 ;

Settings.AllForPhi  =   0       ;

   



%Settings.TrainBurnin_Layer1  =   1500     ;       Settings.TrainCollection_Layer1  =   1000     ;

Settings.ProcessMethod  =   [2]   ;   % 1 => Original Theta        % 2 => Theta / sum(Theta,1)        % 3 => log(Theta)      % 4 => (Theta).^(0.05)


%maxTrial    =   5   ;
%for K0  =   [50,100,200,400,600]

K = ones(1,T)*K0;   %     Accuracies = []     ;   Accuracies_CC     =   [];
%for trial   =    1:maxTrial

%%=================== Training And Testing ===================
rng(trial,'twister');
tic
%                 [ParaGlobal,ParaLocal,Accuracy_all]  =   GBN_LayerwisePretrain_trim_cong_V2_Combined(X_all,prepar,K,T,trial,DataType,dataname,SuPara,Settings)        ;
%                 [ParaGlobal,ParaLocal,Accuracy_all]  =   GBN_LayerwisePretrain_trim_cong_V3_Seperated(X_all,prepar,K,T,trial,DataType,dataname,SuPara,Settings)        ;
%                   [ParaGlobal,ParaLocal,Accuracy_all]  =   GBN_LayerwisePretrain_trim_cong_Journal(X_all,prepar,K,T,trial,DataType,dataname,SuPara,Settings)        ;

if length(Settings.c)==1
    nameccc     =   ['Mnist_K0_',num2str(K0),'_T_',num2str(T),'_eta0',num2str(round(eta0(1)*1000)),...
        '_PreTrain',num2str(Settings.PreTrainYes),'_AllForPhi',num2str(Settings.AllForPhi),'_TrimStrategy',num2str(Settings.TrimStrategy),'_Trial',num2str(trial),'_Data.mat']     ;
else
    nameccc     =   ['Mnist_CJ_K0_',num2str(K0),'_T_',num2str(T),'_eta0',num2str(round(eta0(1)*1000)),...
    '_PreTrain',num2str(Settings.PreTrainYes),'_AllForPhi',num2str(Settings.AllForPhi),'_TrimStrategy',num2str(Settings.TrimStrategy),'_Trial',num2str(trial),'_Data.mat']     ;
end

Settings.nameccc = nameccc;

[ParaGlobal,ParaLocal,Accuracy_all]  =   GBN_Layerwise_IncreaseK(X_all,prepar,K,T,trial,DataType,dataname,SuPara,Settings)        ;
TimeAll     =   toc     ;
%         Accuracies  =   [Accuracies,Accuracy_all]   ;

%         AccThetaFig  =   []  ;
%         for iii     =   1:T
%             AccThetaFig  =   [AccThetaFig,Accuracy_all{iii}.Accuracy_Theta]    ;
%         end, figure,plot(AccThetaFig)
%         AccThetaFreqFig  =   []  ;
%         for iii     =   1:T
%             AccThetaFreqFig  =   [AccThetaFreqFig,Accuracy_all{iii}.Accuracy_ThetaFreq]    ;
%         end, figure,plot(AccThetaFreqFig)
%%===================================  Save Data  ==========================================
%nameccc     =   ['Main_PretrainTrim2_',dataname,'_K0_',num2str(K0),'_T_',num2str(T),'_SuPara',num2str(super1or2),...
%    '_AllForPhi',num2str(Settings.AllForPhi),'_TrimStrategy',num2str(Settings.TrimStrategy),'_Trial',num2str(trial),'_Data.mat']     ;


%save(nameccc,'-v7.3')   ;
if trial==1
    save(['results2/ZGlobal',nameccc],'ParaGlobal','Accuracy_all','prepar','Settings','DataType','SuPara','super1or2')   ;
    save(['results2/ZLocal',nameccc],'ParaLocal')   ;
end
for t=1:T
    Accuracy_all{t}.DataSaving=[];
    Accuracy_all{t}.DataSavingThetaFreqAver=[];
end
save(['results2/',nameccc],'Accuracy_all');


%%===================================  Trial End  ==========================================
%end
%end
%% FiguresExhibition
% FiguresExhibition   =   1   ;
% if FiguresExhibition
%         for Tcurrent    =   T
%                 for tnow    =   1:Tcurrent
%                     if tnow == 1
%                         Phitmp  =   ParaGlobal{Tcurrent}.Phi{tnow} ;
%                     else
%                         Phitmp  =   Phitmp * ParaGlobal{Tcurrent}.Phi{tnow} ;
%                     end
% %                     figure(Tcurrent*10+tnow),
% %                     for ii  =   1:16
% %                         subplot(4,4,ii),plot(Phitmp(:,ii));
% %                     en
%                     if strcmp(DataType,'Binary')
%                         figure(Tcurrent*10+tnow),DispDictionaryImshowNegative(1-exp(-Phitmp)) ; title(['Layer',num2str(tnow)]);
%                     else
%                         figure(Tcurrent*10+tnow),DispDictionaryImshowNegative(Phitmp);title(['Layer',num2str(tnow)]);
%                     end
%                 end
%         end
% end

