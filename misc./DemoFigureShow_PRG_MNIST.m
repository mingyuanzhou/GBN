%close all;
 path = '/Users/zhoum/Box Sync/GBN_results/results_PRG_GBN/'
  path = '/Users/zhoum/Box Sync/GBN_results/results_PRG_GBN_results2/'
 T   =   5  ;



eta_all=[0.01,0.05,0.1];
%
% dataset = i1 %0,1,3,4
% K0 = K0_all(i2) %[50,100,200,400,600]
% eta0=eta_all(i3)
% %trial = i4
% %Settings.TrimStrategy  =   i5
% Settings.PreTrainYes =i6

%dataname = '20newsPcVsMac'
% dataname    =   '20newsAtheismVsReligion';
 %dataname    =   '20newsElecVsMed';
% dataname    =   'CNAE9';
super1or2    = 2 ;
Settings.AllForPhi  =   0       ;

eta0 = 0.05;
Settings.TrimStrategy  = 5;           %  1=>TrimAllAfBurnin      % 2=>TrimTcurrentAfBurnin     %   3=>TrimAllAfEnd      %  4=>TrimTcurrentAfEnd
Settings.PreTrainYes = 0;

%K0Set111  =   K0_all; %[50,100,200,400,600] ; %[50 100 200 400]  ;
%K0Set111  =  [50,100,200,400,600,800] %,800]%,600,800] %,1000]%,400]; %,400,600]; % [50,100,200,400,600]; %[50 100 200 400]  ;

 %,1024];%,1024];
 %,512,1024]
%K0Set111 =800
K0Set111  =  [25,50,100,200,400,600]
TrialSet111    =   1:5    ;

ClassificationAccuracyShow 	=   zeros(length(TrialSet111) , T , length(K0Set111))   ;
ClassificationAccuracyShowFreq 	=   zeros(length(TrialSet111) , T , length(K0Set111))   ;
ClassificationAccuracyShowFreqAver 	=   zeros(length(TrialSet111) , T , length(K0Set111))   ;




for trialii111     =   1:length(TrialSet111)
    trial       =   TrialSet111(trialii111)   ;
    for kii111     =   1:length(K0Set111)
        K0111      =   K0Set111(kii111)  ;
        %         nameccc     =   ['Main_PretrainTrim_',dataname,'_K0_',num2str(K0111),'_T_',num2str(T),'_SuPara',num2str(super1or2),...
        %                 '_AllForPhi',num2str(Settings.AllForPhi),'_TrimStrategy',num2str(Settings.TrimStrategy),'_Trial',num2str(trial),'_Data.mat']     ;
        %nameccc     =   ['Main_PretrainTrim_',dataname,'_K0_',num2str(K0111),'_T_',num2str(T),'_SuPara',num2str(super1or2),...
        %        '_AllForPhi',num2str(Settings.AllForPhi),'_TrimStrategy',num2str(Settings.TrimStrategy),'_Trial',num2str(trial),'_Data.mat']     ;
        if 1
        nameccc     =   [path,'Mnist_CJ_K0_',num2str(K0111),'_T_',num2str(T),'_eta0',num2str(round(eta0*1000)),...
            '_PreTrain',num2str(Settings.PreTrainYes),'_AllForPhi',num2str(Settings.AllForPhi),'_TrimStrategy',num2str(Settings.TrimStrategy),'_Trial',num2str(trial),'_Data.mat']     ;
        else
        nameccc     =   [path,'Mnist_K0_',num2str(K0111),'_T_',num2str(T),'_eta0',num2str(round(eta0*1000)),...
            '_PreTrain',num2str(Settings.PreTrainYes),'_AllForPhi',num2str(Settings.AllForPhi),'_TrimStrategy',num2str(Settings.TrimStrategy),'_Trial',num2str(trial),'_Data.mat']     ;
        end
        
        if exist(nameccc,'file')
             load(nameccc,'Accuracy_all')   ;
         else
             
             %fprintf('dataset=0;K0=%d;T=%d;eta0=%f;Settings.PreTrainYes=%d;Settings.TrimStrategy=%d;trial=%d;\n',K0111,T,eta0, Settings.PreTrainYes,Settings.TrimStrategy,trial)
             i1=5;
             i2 = find(K0_all==K0111);
             i3 = find(eta_all==eta0);
             i4 =trial;
             i5=Settings.TrimStrategy;
             i6=Settings.PreTrainYes;
             fprintf('matlab  -nodisplay -nosplash -nodesktop -r "Call_%s(%d,%d,%d,%d,%d,%d)" -logfile %s_%d_%d_%d_%d_%d_%d.txt\n', core,i1,i2,i3,i4,i5,i6, core,i1,i2,i3,i4,i5,i6);
             %warning(nameccc);
         end
        %load(nameccc,'Accuracy_all')   ;
        for t   =   1:T
            %             ClassificationAccuracyShow(trialii111,t,kii111)    =  Accuracy_all{t}.Accuracy_Theta   ;
%            ClassificationAccuracyShowFreq(trialii111,t,kii111)    =  Accuracy_all{t}.Accuracy_ThetaFreq   ;
            if  t<=5
                ClassificationAccuracyShowFreqAver(trialii111,t,kii111)    =  Accuracy_all{t}.Accuracy_ThetaFreqAver   ;
            else
                ClassificationAccuracyShowFreqAver(trialii111,t,kii111) =90;
            end
        end
    end
end
%fclose(fid);
% save('NEWSGROUP20_ClassificationAccuracySET.mat','ClassificationAccuracyShow')   ;  %   ClassificationAccuracyShow{Tcurrent}(trialii,kii)    =  accuracy(1)   ;
%%================================================  Save Data  End  ============================================================

%%================================================  Reshow Figure Begin  ===========================================================
% load('NEWSGROUP20_ClassificationAccuracy_Data.mat')   ;
LineSettings    =  {'--*','-.o',':s','-d','-->','-.h',':x','-p'}
%{':*','--o','-.s','-d',':>','--h','-.x','-p'}
%{':b*','--ko','-.ms','-rd',':>','--h','-.x','-p'}

%{'b:*','k--o','m-.s','r-d','c-<',':x','--+','-.h','-p',':','--','-.','-'}   ;
xlabelsize  =   16  ;   ylabelsize  =   16  ;   gcafontsize     =   16  ;
linewidthsize   =   1.5   ;   Markersize = 10 ;

%
% H = figure;set(H,'color','w'),hold on     ;
% for kii111     =   1:length(K0Set111)
%     plot(1:T,mean(ClassificationAccuracyShow(:,:,kii111),1),LineSettings{kii111},...
%         'LineWidth',linewidthsize,'MarkerSize',Markersize)   ;  %    ylim([77,85]);
% end
% xlabel('Number of layers T','FontSize',xlabelsize),ylabel('Classification accuracy','FontSize',ylabelsize),
% legend('K_{1max}=50','K_{1max}=100','K_{1max}=200','K_{1max}=400','K_{1max}=600'),set(gca,'fontsize',gcafontsize);%xlim([50 600])
%
%
% H = figure;set(H,'color','w'),hold on     ;
% for Tcurrent    =   1:T
%     plot(K0Set111,mean(squeeze(ClassificationAccuracyShow(:,Tcurrent,:)),1),LineSettings{Tcurrent},...
%         'LineWidth',linewidthsize,'MarkerSize',Markersize)   ;  %    ylim([77,85]);
% end
% xlabel('K_{1max}','FontSize',xlabelsize),ylabel('Classification accuracy','FontSize',ylabelsize),
% legend('T=1','T=2','T=3','T=4','T=5'),set(gca,'fontsize',gcafontsize);xlim([50 600])

%=============================
if 0
    ClassificationAccuracyShow  =   ClassificationAccuracyShowFreq  ;
    
    H = figure;set(H,'color','w'),hold on     ;
    for kii111     =   1:length(K0Set111)
        plot(1:T,mean(ClassificationAccuracyShow(:,:,kii111),1),LineSettings{kii111},...
            'LineWidth',linewidthsize,'MarkerSize',Markersize)   ;  %    ylim([77,85]);
    end
    xlabel('Number of layers T','FontSize',xlabelsize),ylabel('Classification accuracy','FontSize',ylabelsize),
    %legend('K_{1max}=50','K_{1max}=100','K_{1max}=200','K_{1max}=400','K_{1max}=600'),
    legend_text={};
    for KKK=K0Set111
        legend_text{end+1}=['K_{1max}',num2str(KKK)];
    end
    legend(legend_text)
    set(gca,'fontsize',gcafontsize);%xlim([50 600])
    
%     H = figure;set(H,'color','w'),hold on     ;
%     for kii111     =   min(5,max(K0Set111)) %1:length(K0Set111)
%         boxplot(squeeze(ClassificationAccuracyShow(:,:,kii111)));
%         xlabel('Number of layers T','FontSize',xlabelsize),
%         ylabel('Classification accuracy','FontSize',ylabelsize),
%     end
%     
    H = figure;set(H,'color','w'),hold on     ;
    for Tcurrent    =   1:T
        errorbar(K0Set111,mean(squeeze(ClassificationAccuracyShow(:,Tcurrent,:)),1),...
            std(squeeze(ClassificationAccuracyShow(:,Tcurrent,:)),1),LineSettings{Tcurrent},...
            'LineWidth',linewidthsize,'MarkerSize',Markersize)   ;  %    ylim([77,85]);
        %plot(K0Set111,mean(squeeze(ClassificationAccuracyShow(:,Tcurrent,:)),1),LineSettings{Tcurrent},...
          %  'LineWidth',linewidthsize,'MarkerSize',Markersize)   ;  %    ylim([77,85]);
    end
    xlabel('K_{1max}','FontSize',xlabelsize),ylabel('Classification accuracy','FontSize',ylabelsize),
    legend('T=1','T=2','T=3','T=4','T=5'),set(gca,'fontsize',gcafontsize);
    %xlim([50 600])
end
%=============================
ClassificationAccuracyShow  =   ClassificationAccuracyShowFreqAver  ;
%ClassificationAccuracyShow  =   ClassificationAccuracyShowFreq  ;
H = figure;set(H,'color','w'),hold on     ;
for kii111     =   1:length(K0Set111)
   % plot(1:T,mean(ClassificationAccuracyShow(:,:,kii111),1),LineSettings{kii111},...
   %     'LineWidth',linewidthsize,'MarkerSize',Markersize)   ;  %    ylim([77,85]);
    
    errorbar(1:T,mean(ClassificationAccuracyShow(:,:,kii111),1),...
        std(ClassificationAccuracyShow(:,:,kii111),1), LineSettings{kii111},...
        'LineWidth',linewidthsize,'MarkerSize',Markersize)   ;
    %set(gca,'xscale','log')
end
xlabel('Number of layers T')
%,'FontSize',xlabelsize),
ylabel('Classification accuracy')
%,'FontSize',ylabelsize),
%legend('K_{1max}=50','K_{1max}=100','K_{1max}=200','K_{1max}=400','K_{1max}=600'),set(gca,'fontsize',gcafontsize);%xlim([50 600])
legend_text={};
for KKK=K0Set111
    legend_text{end+1}=['K_{1max} = ',num2str(KKK)];
end
legend(legend_text)
%set(gca,'fontsize',gcafontsize);%xlim([50 600])
box on
title('(a)')
%xlim([0.8,7])
% 
% H = figure; %(100);%set(H,'color','w'),
% hold on     ;
% for kii111     =   min(5,length(K0Set111))  % 1:length(K0Set111)
%     boxplot(squeeze(ClassificationAccuracyShow(:,:,kii111)));hold on
%     xlabel('Number of layers T','FontSize',xlabelsize),
%     ylabel('Classification accuracy','FontSize',ylabelsize),
% end
% hold on;
% plot(mean(squeeze(ClassificationAccuracyShow(:,:,kii111)),1),'d-','Linewidth',2);hold on
% 
% H = figure; %(100);%set(H,'color','w'),
% hold on     ;
% for kii111     =   min(6,length(K0Set111))  % 1:length(K0Set111)
%     boxplot(squeeze(ClassificationAccuracyShow(:,:,kii111)));hold on
%     xlabel('Number of layers T','FontSize',xlabelsize),
%     ylabel('Classification accuracy','FontSize',ylabelsize),
% end
% hold on;
% plot(mean(squeeze(ClassificationAccuracyShow(:,:,kii111)),1),'d-','Linewidth',2);hold on
% 
% 
H = figure;set(H,'color','w'),hold on     ;
for Tcurrent    =   1:T
    %plot(K0Set111,mean(squeeze(ClassificationAccuracyShow(:,Tcurrent,:)),1),LineSettings{Tcurrent},...
    %    'LineWidth',linewidthsize,'MarkerSize',Markersize)   ;  %    ylim([77,85]);
    errorbar(K0Set111,mean(squeeze(ClassificationAccuracyShow(:,Tcurrent,:)),1),...
            std(squeeze(ClassificationAccuracyShow(:,Tcurrent,:)),1),LineSettings{Tcurrent},...
            'LineWidth',linewidthsize,'MarkerSize',Markersize)   ;  %    ylim([77,85]);
%set(gca,'xscale','log')
end
xlabel('K_{1max}')
%,'FontSize',xlabelsize),
ylabel('Classification accuracy')
%,'FontSize',ylabelsize),
legend_text={};
for t=1:T
    legend_text{end+1}=['T=',num2str(t)];
end
legend(legend_text)
%set(gca,'fontsize',gcafontsize);
%xlim([35,815])
% %%================================================  Reshow Figure End  ===========================================================
box on
title('(b)')
%if i1==1
    ax=gca;
set(ax,'XTick',K0Set111)
set(ax,'XTickLabel',K0Set111)
%end
% 
% H = figure;set(H,'color','w'),hold on     ;
% for Tcurrent    =   1:T
%     plot(K0Set111,(squeeze(ClassificationAccuracyShow(:,Tcurrent,:))).',LineSettings{Tcurrent},...
%         'LineWidth',linewidthsize,'MarkerSize',Markersize)   ;  %    ylim([77,85]);
% end
% xlabel('K_{1max}','FontSize',xlabelsize),ylabel('Classification accuracy','FontSize',ylabelsize),
% legend('T=1','T=2','T=3','T=4','T=5'),set(gca,'fontsize',gcafontsize);
% %xlim([50 600])




%%=====================================================================================================================================
%%=====================================================================================================================================
%%=====================================================================================================================================
%%=====================================================================================================================================




% % % % T   =   5   ;
% % % % K0Set111  =   [50 100 200 400 600]  ;
% % % % TrialSet111    =   1:5     ;
% % % %
% % % % ClassificationAccuracyShow1 	=   cell(T,1)   ;
% % % % for Tcurrent    =   1:T
% % % %     for kii111     =   1:length(K0Set111)
% % % %         K0111      =   K0Set111(kii111)  ;
% % % %         for trialii111     =   1:length(TrialSet111)
% % % %             trial       =   TrialSet111(trialii111)   ;
% % % %             nameccc     =   ['MultiClass_K_',num2str(K0111),'_T_',num2str(T),'_Tcurrent_',num2str(Tcurrent),'_Trial_',num2str(trial),'_Data_After_Classification.mat']     ;
% % % %             load(nameccc)   ;
% % % %             ClassificationAccuracyShow1{Tcurrent}(trialii111,kii111)    =  accuracy(1)   ;
% % % %         end
% % % %     end
% % % % end
% % % % % save('NEWSGROUP20_ClassificationAccuracy_Data_128_512.mat','ClassificationAccuracyShow1')   ;  %   ClassificationAccuracyShow{Tcurrent}(trialii,kii)    =  accuracy(1)   ;
% % % % % save('NEWSGROUP20_ClassificationAccuracy_Data.mat','ClassificationAccuracyShow1')   ;  %   ClassificationAccuracyShow{Tcurrent}(trialii,kii)    =  accuracy(1)   ;
% % % % %%================================================  Save Data  End  ============================================================
% % % %
% % % % %%================================================  Reshow Figure Begin  ===========================================================
% % % % % load('NEWSGROUP20_ClassificationAccuracy_Data.mat')   ;
% % % % LineSettings    =   {'b:*','k--o','m-.s','r-d','c-<'}   ;
% % % % xlabelsize  =   16  ;   ylabelsize  =   16  ;   gcafontsize     =   16  ;
% % % % linewidthsize   =   2   ;   Markersize = 15 ;
% % % %
% % % %
% % % % H = figure;set(H,'color','w'),hold on     ;
% % % % for Tcurrent    =   1:T
% % % %     plot(K0Set111,mean(ClassificationAccuracyShow1{Tcurrent},1),LineSettings{Tcurrent},...
% % % %         'LineWidth',linewidthsize,'MarkerSize',Markersize)   ;  %    ylim([77,85]);
% % % % end
% % % % xlabel('K_{1max}','FontSize',xlabelsize),ylabel('Classification accuracy','FontSize',ylabelsize),
% % % % legend('T=1','T=2','T=3','T=4','T=5'),set(gca,'fontsize',gcafontsize);xlim([50 600])
% % % % %%================================================  Reshow Figure End  ===========================================================





%%=====================================================================================================================================
%%=====================================================================================================================================
%%=====================================================================================================================================
%%=====================================================================================================================================






% % clear,clc,close all     ;
% % %%================================================  Save Data  Begin  ===========================================================
% % T   =   8   ;
% % K0Set  =   [50 100 200 400 600]  ;
% % TrialSet    =   1:5     ;
% %
% % ClassificationAccuracyShow 	=   cell(T,1)   ;
% % for i   =   1:T
% %     ClassificationAccuracyShow{i} 	=   zeros(length(TrialSet),length(K0Set))   ;
% % end
% % for Tcurrent    =   1:T
% %     for kii     =   1:length(K0Set)
% %         K0      =   K0Set(kii)  ;
% %         for trialii     =   1:length(TrialSet)
% %             trial       =   TrialSet(trialii)   ;
% %             nameccc     =   ['MultiClass_K_',num2str(K0),'_T_',num2str(T),'_Tcurrent_',num2str(Tcurrent),'_Trial_',num2str(trial),'_Data_Before_Classification.mat']     ;
% %             if exist(nameccc,'file')
% %                 load(nameccc)   ;
% %                 %%==========================  CLASSIFICATION  BEGIN  ===================================================
% %                 addpath('/home/xd1708/Documents/CYL/liblinear-1.96/matlab/')  ;
% %                 CC=2.^(-10:2:18);
% %                 ModelOut=zeros(1,length(CC));
% %                 parfor ij=1:length(CC)
% %                     ModelOut(ij) = train(double(TrainY), sparse(TrainX'), ['-s 0 -c ', num2str(CC(ij)), ' -v 3  -q ']);
% %                 end
% %                 [~,maxdex]=max(ModelOut);
% %                 num2str(CC(maxdex));
% %                 option = ['-s 0 -c ', num2str(CC(maxdex)), ' -q'];
% %                 model = train(double(TrainY), sparse(TrainX'), option);
% %                 [predicted_label, accuracy, prob_estimates] = predict(double(TestY), sparse(TestX'), model, ' -b 1');
% %                 ClassificationAccuracyShow{Tcurrent}(trialii,kii)    =  accuracy(1)   ;
% %                 nameddd     =   ['MultiClass_K_',num2str(K0),'_T_',num2str(T),'_Tcurrent_',num2str(Tcurrent),'_Trial_',num2str(trial),'_Data_After_Classification.mat']     ;
% %                 save(nameddd)   ;
% %                 %%==========================  CLASSIFICATION  END   ====================================================
% %             end
% %         end
% %     end
% % end
% % % save('NEWSGROUP20_ClassificationAccuracy_Data.mat','ClassificationAccuracyShow')   ;  %   ClassificationAccuracyShow{Tcurrent}(trialii,kii)    =  accuracy(1)   ;
% % %%================================================  Save Data  End  ============================================================
% %
% % %%================================================  Reshow Figure Begin  ===========================================================
% % % load('NEWSGROUP20_ClassificationAccuracy_Data.mat','ClassificationAccuracyShow')   ;
% % H = figure;set(H,'color','w'),
% % for Tcurrent    =   1:T
% %     subplot(1,2,1)  ;   hold on     ;
% %     boxplot(ClassificationAccuracyShow{Tcurrent})   ;  %    ylim([77,85]);
% % end
% % legend('Layer 1','Layer 2','Layer 3','Layer 4','Layer 5','Layer 6','Layer 7','Layer 8')
% % set(gca,'fontsize',20);set(findobj(gca,'type','text'),'FontSize',20);
% % txt = findobj(gca,'type','text');set(txt(1:end),'VerticalAlignment','Middle')
% %
% % for Tcurrent    =   1:T
% %     subplot(1,2,2);    hold on     ;
% %     errorbar(mean(ClassificationAccuracyShow{Tcurrent},1),std(ClassificationAccuracyShow{Tcurrent},1));     %   ylim([77,85]);
% % end
% % legend('Layer 1','Layer 2','Layer 3','Layer 4','Layer 5','Layer 6','Layer 7','Layer 8')
% % set(gca,'fontsize',20);
% % %%================================================  Reshow Figure Begin  ===========================================================
% 
% if i1==5
% mean(ClassificationAccuracyShow(:,:,K0Set111==128),1)
% std(ClassificationAccuracyShow(:,:,K0Set111==128),1)
% mean(ClassificationAccuracyShow(:,:,K0Set111==256),1)
% std(ClassificationAccuracyShow(:,:,K0Set111==256),1)
% mean(ClassificationAccuracyShow(:,:,K0Set111==512),1)
% std(ClassificationAccuracyShow(:,:,K0Set111==512),1)
% 
% else
%     mean(ClassificationAccuracyShow(:,:,K0Set111==200),1)
% std(ClassificationAccuracyShow(:,:,K0Set111==200),1)
% mean(ClassificationAccuracyShow(:,:,K0Set111==400),1)
% std(ClassificationAccuracyShow(:,:,K0Set111==400),1)
% mean(ClassificationAccuracyShow(:,:,K0Set111==800),1)
% std(ClassificationAccuracyShow(:,:,K0Set111==800),1)
% end