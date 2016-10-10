K0 = 200; %[25,50,100,200,400,600,800]
eta= 0.05; %[0.01,0.05,0.1];
trial = 1; %[1,2,3,4,5]
T = 5;
K = ones(1,T)*K0;

Para.TrainBurnin = [1000,500*ones(1,T-1)];
Para.TrainCollection = [500,500*ones(1,T-1)];

Para.TrainBurnin = sum(Para.TrainBurnin)+sum(Para.TrainCollection)-500;
Para.TrainCollection=500;
T = 1;
K = ones(1,T)*K0;

Para.percentage=0.3;
Para.dataname = [dataname,'_',num2str(Para.percentage*100)];
Para.CollectionStep=5;
Para.train_idx = 1:size(X_all,2);



ToBeAnalized    =   9   ;


K0_all = [25,50,100,200,400,600,800,1000]
eta_all = [0.01,0.05,0.1];
if dataset==5
    K0_all = [16,32,64,128,256,512,1024]
end

K0 = 50  %[50,100,200,400,600]
eta0=eta_all(i3)
trial = i4
Para.TrimStrategy  =   i5
Para.PreTrainYes =i6

if i1==0 || i1==5
    T   =   5
    Para.TrainBurnin  =   [1000,1000*ones(1,T-1)]     ;       Para.TrainCollection  =   [500,500*ones(1,T-1)];            Para.TrainSampleSpace   =   1   ;
else
    T   =   8
    Para.TrainBurnin  =   [1000,1000*ones(1,T-1)]     ;       Para.TrainCollection  =   [1000,1000*ones(1,T-1)];            Para.TrainSampleSpace   =   1   ;
end

if i1==0 && K0>500
    T   =   5
    Para.TrainBurnin  =   [1000,500*ones(1,T-1)]     ;       Para.TrainCollection  =   [500,500*ones(1,T-1)];            Para.TrainSampleSpace   =   1   ;
end

%clear,clc,close all;
SystemInRuningLinux  =   1   ;   %   1:linux       ;;;     0 : windows

% if SystemInRuningLinux
addpath('liblinear-2.1/matlab/');
addpath('data/');
% else
%   addpath('C:\liblinear-2.1\matlab\');
% end

%% Data
 

Para.maxIter_preTrain    =   100     ;
%Para.TrainBurnin  =   [1500,1000*ones(1,T-1)]     ;       Para.TrainCollection  =   [1000,500*ones(1,T-1)];            Para.TrainSampleSpace   =   1   ;

switch ToBeAnalized
    case 1          %%========= MNIST =================
        TrainSize      =   1e3   ;       TestSize        =   1e3   ;
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
    otherwise
        error('Wrong "ToBeAnalized"')
end
LoadData_GBN    ;

%% Direct SVM
%%Choose L2 regularized multiclass logisitic regression or svm for classication

% DirectSVMProcessing     =  1   ;
% if DirectSVMProcessing
%     TrainX      =   X_all(:,prepar.trindx) + eps;  TrainY  =   prepar.Y(prepar.trindx)     ;
%     TestX       =   X_all(:,prepar.teindx)+ eps;  TestY  =   prepar.Y(prepar.teindx)     ;
%     %%==================================== NO PreprocessBfTest ===================================
% %     tic;
% %     BestModel = train(double(TrainY), sparse(TrainX'), ['-C -s 0']);
% %     CVtime  =   toc
% %     tic ;
% %     option = ['-s 0 -c ', num2str(BestModel(1)), ' -q'];
% %     model = train(double(TrainY), sparse(TrainX'), option);
% %     Traintime   =   toc
% %     tic
% %     [predicted_label_DSVM, Accuracy_DSVM, prob_estimates_DSVM] = predict(double(TestY), sparse(TestX'), model, ' -b 1');
% %     Testtime = toc
%     %%==================================== NO PreprocessBfTest ===================================
%     CC=2.^(-14:2:10);
%     ModelOut=zeros(1,length(CC));
%     parfor ij=1:length(CC)
%         ModelOut(ij) = train(double(TrainY), sparse(TrainX'), ['-s 0 -c ', num2str(CC(ij)), ' -v 3 -q ']);
%     end
%     [~,maxdex]=max(ModelOut);
%     option = ['-s 0 -c ', num2str(CC(maxdex)), ' -q'];
%     model = train(double(TrainY), sparse(TrainX'), option);
%     [predicted_label_DSVM, Accuracy_DSVM, prob_estimates_DSVM] = predict(double(TestY), sparse(TestX'), model, ' -b 1');
%     nameaaa     =   ['Main_PretrainTrim_',dataname,'_DirectSVM_TrainSize_',num2str(TrainSize),'_TestSize_',num2str(TestSize),'_Data.mat']     ;
%     save(nameaaa,'predicted_label_DSVM','Accuracy_DSVM','prob_estimates_DSVM')   ;
% else
%     nameccc     =   ['Main_PretrainTrim_',dataname,'_DirectSVM_TrainSize_',num2str(TrainSize),'_TestSize_',num2str(TestSize),'_Data.mat']     ;
%     load(nameccc)   ;
% end

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
Para.IsDisplay   =   0   ;          Para.FigureGap       =   1000   ;
Para.SaveDataGap       =   200   ;   Para.TestGap    =   1e5  ;
Para.IsInferEta     =   zeros(1,T)   ;     Para.IsTestInside   =    1 ;

Para.AllForPhi  =   0       ;

Para.ParallelProcessing      =   0   ;     Para.NumBlockParallel   =   8   ;

%Para.TrainBurnin_Layer1  =   1500     ;       Para.TrainCollection_Layer1  =   1000     ;

Para.TestBurnin  =   500     ;       Para.TestCollection  =   500     ;       Para.TestSampleSpace   =   1   ;
Para.ProcessMethod  =   [2]   ;   % 1 => Original Theta        % 2 => Theta / sum(Theta,1)        % 3 => log(Theta)      % 4 => (Theta).^(0.05)


%maxTrial    =   5   ;
%for K0  =   [50,100,200,400,600]

K = ones(1,T)*K0;   %     Accuracies = []     ;   Accuracies_CC     =   [];
%for trial   =    1:maxTrial

%%=================== Training And Testing ===================
rng(trial,'twister');

nameccc     =   ['Trim_',dataname,'_K0_',num2str(K0),'_T_',num2str(T),'_eta0',num2str(round(eta0*1000)),...
    '_PreTrain',num2str(Para.PreTrainYes),'_AllForPhi',num2str(Para.AllForPhi),'_TrimStrategy',num2str(Para.TrimStrategy),'_Trial',num2str(trial),'_Data.mat']     ;
Para.nameccc=nameccc;

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
%    '_AllForPhi',num2str(Para.AllForPhi),'_TrimStrategy',num2str(Para.TrimStrategy),'_Trial',num2str(trial),'_Data.mat']     ;

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

