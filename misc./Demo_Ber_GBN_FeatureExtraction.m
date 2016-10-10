clear all
%parpool(2)

i1=5;
i2=1;
i3=1;
i4=1;

T   =  5
K0_all = [25,50,100,200,400,500,600,800,1000]
eta_all=[0.01,0.05,0.1];
K0 = K0_all(i1) %[50,100,200,400,600]
eta0=eta_all(i2)
%eta0=[eta_all(i2),0.1*ones(1,T-1)];
trial = i3
Settings.TrimStrategy  =   5
Settings.PreTrainYes = 0


Settings.TrainBurnin  =   [400,400*ones(1,T-1)]/1     ;
Settings.TrainCollection  = [100,100*ones(1,T-1)]/1;
Settings.TrainSampleSpace   =   1   ;

Settings.TestBurnin  =  250/1    ;       
Settings.TestCollection  =250/1     ;
Settings.TestSampleSpace   =   1   ;


Settings.ParallelProcessing      =   i4    

POOL.NumWorkers = 12;
Settings.NumBlockParallel   =   POOL.NumWorkers*1;


%POOL = gcp; % If no pool, do not create new one.




Settings.IsDisplay   =   0   ;          Settings.FigureGap       =   10   ;


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


%%============================ CalTech_101_Silhouettes ======================================
DataType    =   'Binary';
DataType = 'Positive';

if 0
dataname    =   'CalTech_101_Silhouettes';
SmallPatch=false;
TargetOne = false;
DataOrdering = false;
if SmallPatch
    load('caltech101_silhouettes_16_split1.mat')  ;
else
    load('caltech101_silhouettes_28_split1.mat')  ;
end

TrainData       =   [(train_data).',val_data.']     ;       clear train_data val_data        ;
TrainLabel      =   [train_labels(:);val_labels(:)]  ;       clear train_labels val_labels      ;
TestData        =   (test_data).'     ;       clear test_data          ;
TestLabel       =   test_labels(:)           ;       clear test_labels        ;

if TargetOne
    TrainData       =   1-TrainData     ;    
    TestData        =   1-TestData       ;     
end

%%================ Ordering  ====================
if DataOrdering
    [TrainLabel,tmp]    =   sort(TrainLabel)    ;
    TrainData   =   TrainData(:,tmp)    ;
    [TestLabel,tmp]    =   sort(TestLabel)    ;
    TestData   =   TestData(:,tmp)    ;
end
%%================ Ordering end  ================

for ii  =   1:101
    SampleNumPerClassTrain(ii)   =   sum(TrainLabel==ii)    ;
    SampleNumPerClassTest(ii)   =   sum(TestLabel==ii)    ;
    %             SampleNumPerClassVal(ii)   =   sum(val_labels==ii)    ;
end

else
%     dataname    =   'OCR';
%     load OCR;
%     TrainData = X(:,1:42152);
%     TestData = X(:,42153:end);
%     TrainLabel = Y(1:42152);
%     TestLabel = Y(42153:end);
    
%     %%
%      dataname = 'MnistBinary';
%     %if SystemInRuningLinux
%     load('data/mnist_gray.mat')     ;   % Load Data
% %else
% %    load E:\0Research\Data\MNIST\mnist_gray.mat     ;   % Load Data
% %end
% %         [~,IX]  =   sort(train_label,'ascend')  ;   train_mnist     =   train_mnist(:,IX)   ;    train_label  =   train_label(IX)  ;
% %         [~,IX]  =   sort(test_label,'ascend')   ;   test_mnist      =   test_mnist(:,IX)   ;    test_label  =   test_label(IX)  ;
% %         TrainSize      =   1e3   ;       TestSize        =   0.5e3   ;
% %         MinibatchSize   =   250 ;
% %         K   =   [128,64]     ;
% TrainSize = 60000/10;
% TestSize = 10000/10;
% %======== Prepare Train Data  ===========
% ndx = []; mtrain = TrainSize / 10;
% if mtrain < 6000
%     for ii = 0:1:9
%         tmp = find(train_label==ii);
%         ndx = [ndx; tmp(1:mtrain)];
%     end
% else
%     ndx = [1:60000];
% end
% X       =   train_mnist(:,ndx);
% Xlabel  =   train_label(ndx);
% rng(0,'twister');
% DataType    =   'Positive';
% %DataType    =   'Binary';
% 
% %======== Prepare Test Data   ===========
% ndx = []; mtest = TestSize/10  ;
% if mtest < 1000
%     for ii = 0:1:9
%         tmp = find(test_label==ii);
%         ndx = [ndx; tmp(1:mtest)];
%     end
% else
%     ndx=[1:10000];
% end
% Xtest       =   test_mnist(:,ndx);
% Xtestlabel  =   test_label(ndx);
% %======== Combine Train and Test ===========
% clear train_mnist train_label test_mnist test_label ;
% X_all   =   [X,Xtest]  ;
% 
% X_all = double(X_all>rand(size(X_all)));
% 
% prepar.trindx   =   1:length(Xlabel)    ;
% prepar.teindx   =   length(Xlabel) + (1:length(Xtestlabel))     ;
% prepar.Y        =   [Xlabel(:);Xtestlabel(:)]   ;
% prepar.Y(prepar.Y == 0)     =   10  ;
% %         MNISTPreprocessing  ;
% clear X Xlabel  ndx Xtest Xtestlabel;

ToBeAnalized = 9;
dataset = 3
switch ToBeAnalized
    case 1          %%========= MNIST =================
        TrainSize      =   1e3   ;       TestSize        =   1e3   ;
    case 2        %%======== MIT_CBCL_FACE =============
        TrainSize      =   2429   ;
    case 3        %%========= PIE FACE =================
        TrainSize      =   680   ;       TestSize        =   680   ;
    case 4        %%========= GBN ToyData =============
        % K   =   [80,8]  ;
    case 5        %%========= HTTP ================
        TrainSize      =   600   ;       TestSize        =   2400   ;
        DataUsed    =  2   ;   %   1:HRRP_Data_1HRRP_MCA         2: HRRP_Data_250PowAverHRRP_MCA
    case 6        %%========== MSTAR 3 ==============
        TrainSize      =   698   ;       TestSize        =   1365   ;
        SmallPatch  =   1   ;
        NormalizationMethod  =   2   ;  %  0:  Nothing          1:divide_max_all        2:divide_max_respect      3: divide_energy
    case 7        %%========= CalTech_101_Silhouettes =======
        SmallPatch  =   0   ;               %    1: 16*16                         0: 28*28
        TargetOne   =   1   ;
        DataOrdering    =   1   ;
    case 8         %%======== MPEG7_CE-Shape-1_Part_B =======
        TrainSize      =   700   ;       TestSize        =   1400 - TrainSize    ;
    case 9
        %%======== 20news-bydate   ============
        switch dataset
            case {0,5}
                IsBinaryClassificaiton  =   false   ;
            otherwise
                IsBinaryClassificaiton  =   true   ;
                %dataset = 3; %%comp comp.sys.ibm.pc.hardware 4 vs comp.sys.mac.hardware 5
                %dataset = 4; %%sci sci.electronics 13 vs sci.med 14
        end
    case 10        %%========= CNAE-9  ===============
        TrainPersentage  =   0.8     ;
    otherwise
        error('Wrong "ToBeAnalized"')
end
LoadData_GBN    ;

dataname    =   '20newsbinary';

    DataType = 'Binary';
X_all = double(X_all>rand(size(X_all)));
end


% X_all   =   [TrainData , TestData]   ;
% prepar.trindx   =   1 : length(TrainLabel)  ;
% prepar.teindx   =    length(TrainLabel)  + (1: length(TestLabel))  ;
% prepar.Y        =   [TrainLabel(:) ; TestLabel(:)]   ;
% %prepar.classnames   =   classnames  ;
% TrainSize   =   length(prepar.trindx)   ;
% TestSize    =   length(prepar.teindx)   ;
% %         X_all_for_Training  =   X_all   ;
% clear TrainData  TrainLabel  TestData  TestLabel  classnames val_data val_labels   SampleNumPerClassTrain  SampleNumPerClassTest;



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

if length(Settings.c)==1
    nameccc     =   [dataname,DataType,'_K0_',num2str(K0),'_T_',num2str(T),'_eta0',num2str(round(eta0(1)*1000)),...
        '_PreTrain',num2str(Settings.PreTrainYes),'_AllForPhi',num2str(Settings.AllForPhi),'_TrimStrategy',num2str(Settings.TrimStrategy),'_Trial',num2str(trial),'_Data.mat']     ;
else
    nameccc     =   [dataname,DataType,'_CJ_K0_',num2str(K0),'_T_',num2str(T),'_eta0',num2str(round(eta0(1)*1000)),...
    '_PreTrain',num2str(Settings.PreTrainYes),'_AllForPhi',num2str(Settings.AllForPhi),'_TrimStrategy',num2str(Settings.TrimStrategy),'_Trial',num2str(trial),'_Data.mat']     ;
end

Settings.nameccc=nameccc;

tic
%                 [ParaGlobal,ParaLocal,Accuracy_all]  =   GBN_LayerwisePretrain_trim_cong_V2_Combined(X_all,prepar,K,T,trial,DataType,dataname,SuPara,Settings)        ;
%                 [ParaGlobal,ParaLocal,Accuracy_all]  =   GBN_LayerwisePretrain_trim_cong_V3_Seperated(X_all,prepar,K,T,trial,DataType,dataname,SuPara,Settings)        ;
%                   [ParaGlobal,ParaLocal,Accuracy_all]  =   GBN_LayerwisePretrain_trim_cong_Journal(X_all,prepar,K,T,trial,DataType,dataname,SuPara,Settings)        ;
[ParaGlobal,ParaLocal,Accuracy_all]  =   GBN_Layerwise(X_all,prepar,K,T,trial,DataType,dataname,SuPara,Settings)        ;
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
    save(['results1/ZGlobal',nameccc],'ParaGlobal','Accuracy_all','prepar','Settings','DataType','SuPara','super1or2')   ;
    save(['results1/ZLocal',nameccc],'ParaLocal')   ;
end
for t=1:T
    Accuracy_all{t}.DataSaving=[];
    Accuracy_all{t}.DataSavingThetaFreqAver=[];
end
save(['results1/',nameccc],'Accuracy_all');


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

