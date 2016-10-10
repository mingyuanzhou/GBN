% ToBeAnalized    =   9   ;
% T   =   5 ;
% 
% dataset = 3;
% K0 = 25; %[50,100,200,400,600]
% eta0 = 0.05;
% %Settings.TrainBurnin  =   [2500,1500*ones(1,T-1)]     ;       Settings.TrainCollection  =   [0,0*ones(1,T-1)];            Settings.TrainSampleSpace   =   1   ;
% Settings.TrimStrategy  =   4;
% trial=6
% %  1=>TrimAllAfBurnin      % 2=>TrimTcurrentAfBurnin     %   3=>TrimAllAfEnd      %  4=>TrimTcurrentAfEnd
% Settings.PreTrainYes     =   1   ;   
%Settings.TrainBurnin  =   [1000,1000*ones(1,T-1)]     ;       Settings.TrainCollection  =   [1000,1000*ones(1,T-1)];            Settings.TrainSampleSpace   =   1   ;


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

DataType    =   'Binary';
X_all = double(X_all>=1);

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
Settings.IsDisplay   =   0   ;          Settings.FigureGap       =   1000   ;
Settings.SaveDataGap       =   200   ;   Settings.TestGap    =   1e5  ;
Settings.IsInferEta     =   zeros(1,T)   ;     Settings.IsTestInside   =    1 ;

Settings.AllForPhi  =   0       ;

Settings.c=1;

%Settings.TrainBurnin_Layer1  =   1500     ;       Settings.TrainCollection_Layer1  =   1000     ;

Settings.TestBurnin  =   500     ;       Settings.TestCollection  =   500     ;       Settings.TestSampleSpace   =   1   ;
Settings.ProcessMethod  =   [2]   ;   % 1 => Original Theta        % 2 => Theta / sum(Theta,1)        % 3 => log(Theta)      % 4 => (Theta).^(0.05)


%maxTrial    =   5   ;
%for K0  =   [50,100,200,400,600]

K = ones(1,T)*K0;   %     Accuracies = []     ;   Accuracies_CC     =   [];
%for trial   =    1:maxTrial

%%=================== Training And Testing ===================
rng(trial,'twister');

nameccc     =   ['Trim_BerPo_',dataname,'_K0_',num2str(K0),'_T_',num2str(T),'_eta0',num2str(round(eta0*1000)),...
    '_PreTrain',num2str(Settings.PreTrainYes),'_AllForPhi',num2str(Settings.AllForPhi),'_TrimStrategy',num2str(Settings.TrimStrategy),'_Trial',num2str(trial),'_Data.mat']     ;
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

