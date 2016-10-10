%%================================================  Reshow Figure Begin  ===========================================================

path = '/Users/zhoum/Box Sync/GBN_results/results2/'
T=5
K0_all = [25,50,100,200,400,600,800]%,800]%,800] %,200]; %,400,600,800]
eta_all=[0.01,0.05,0.1];

i1=2;


eta0 = 0.05;
Settings.TrimStrategy  =   5 %4,5
PerplexityShowTest = cell(T,1);
Settings.AllForPhi=1;
for i2 = 1:length(K0_all)
    K0 = K0_all(i2);
    for trial =1:5
        
        switch i1
            case 1
                dataname = '20newsTop2000'
            case 2
                dataname = 'NIPSTop2000'
        end
        
        K0 = K0_all(i2) %[50,100,200,400,600]
        
        Settings.PreTrainYes = 0
        
        Settings.percentage=0.3;
        
        dataname = [dataname,'_',num2str(Settings.percentage*100)]
        
        nameccc     =   [path,'Perplexity_',dataname,'_K0_',num2str(K0),'_T_',num2str(T),'_eta0',num2str(round(eta0*1000)),...
            '_PreTrain',num2str(Settings.PreTrainYes),'_AllForPhi',num2str(Settings.AllForPhi),'_TrimStrategy',num2str(Settings.TrimStrategy),'_Trial',num2str(trial),'_Data.mat']     ;
        
        load(nameccc,'Perp');
        
        for t=1:T
            PerplexityShowTest{t}(trial,i2) = Perp(t);
        end
        
    end
end

T   =   5
Settings.TrainBurnin  =   [1000,500*ones(1,T-1)]     ;       Settings.TrainCollection  =   [500,500*ones(1,T-1)];            Settings.TrainSampleSpace   =   1   ;


%load('NEWS20_Perplextiy_Data.mat')    ;
%   PerplexityShowTest{Tcurrent}(trial,K0)    =  ave.loglike(end)   ;
T   =   5   ; 
K0Set  =   K0_all; % 300 400]  ;
TrialSet    =     0:4     ;

H = figure;set(H,'color','w'),hold on     ;
 for Tcurrent    =   1:T
     aaaaa   =  (abs(PerplexityShowTest{Tcurrent}))  ; 
     errorbar(K0Set,mean(aaaaa,1),std(aaaaa,1),'LineWidth',1.5,'MarkerSize',10);     %   ylim([77,85]);
end
xlabel('K_{1max}','FontSize',20),ylabel('Perplexity','FontSize',20),legend('Layer 1','Layer 2','Layer 3','Layer 4','Layer 5')
%title('20NEWS'),xlim([0 420])
%set(gca,'fontsize',16);
set(gca,'xscale','log')
figure

%colors={':b*','--ko','-.ms','-rd',':>'}
colors={'--*','-.o',':s','-d','-->','-.h',':x','-p'}
%{':*','--o','-.s','-d',':>'}

layer1   =  mean(((PerplexityShowTest{1})),1)  ; 

%H = figure;set(H,'color','w'),hold on     ;
 for Tcurrent    =   1:T
     subplot(1,2,1)
     aaaaa   =  (abs(PerplexityShowTest{Tcurrent}))  ; 
     
     errorbar(K0Set,mean(aaaaa,1),std(aaaaa,1),colors{Tcurrent},'LineWidth',1.5,'MarkerSize',10);     %   ylim([77,85]);
     %plot(K0Set,mean(aaaaa,1),colors{Tcurrent})
     %set(gca,'xscale','log')
     hold on;
     %set(gca,'xscale','log')
     subplot(1,2,2)
     aaaaa   =  -bsxfun(@plus,-((PerplexityShowTest{Tcurrent})),layer1)  ; 
     %semilogx(K0Set,mean(aaaaa,1),colors{Tcurrent}); hold on   %   ylim([77,85]);
     errorbar(K0Set,mean(aaaaa,1),std(aaaaa,1),colors{Tcurrent},'LineWidth',1.5,'MarkerSize',10);    %   ylim([77,85]);
     %plot(K0Set,mean(aaaaa,1),colors{Tcurrent})
     %set(gca,'xscale','log')
     hold on;
 end
 
 K0Set(2)=[];
subplot(1,2,1);
legend('T = 1','T = 2','T = 3','T = 4','T = 5')
xlabel('K_{1max}'), ylabel('Perplexity');
%xlim([10,410])
ax=gca;
xlim([25,800])
set(ax,'XTick',K0Set)
set(ax,'XTickLabel',K0Set)
title('(a)')
subplot(1,2,2);
legend('T = 1','T = 2','T = 3','T = 4','T = 5')
xlabel('K_{1max}'), ylabel('Relative perplexity');
ax=gca;
set(ax,'XTick',K0Set)
set(ax,'XTickLabel',K0Set)
title('(b)')
xlim([25,800])
set(gcf,'color','w')

box on

%xlabel('K_{1max}','FontSize',20),ylabel('Perplexity','FontSize',20),legend('Layer 1','Layer 2','Layer 3','Layer 4','Layer 5')
%xlabel('K_{1max}'),ylabel('Perplexity'),
%legend('Layer 1','Layer 2','Layer 3','Layer 4','Layer 5')
%legend('T = 1','T = 2','T = 3','T = 4','T = 5')

%title('20NEWS'),%xlim([0 420])
%set(gca,'fontsize',16);


path = '/Users/zhoum/Box Sync/GBN_results/results1_global/'
for i2 = 1:length(K0_all)
    K0 = K0_all(i2);
    trial = 1;
    
    switch i1
        case 1
            dataname = '20newsTop2000';
        case 2
            dataname = 'NIPSTop2000';
    end
    
    K0 = K0_all(i2); %[50,100,200,400,600]
    
    Settings.PreTrainYes = 0;
    
    Settings.percentage=0.3;
    
    dataname = [dataname,'_',num2str(Settings.percentage*100)];
    
    nameccc     =   [path,'ZGlobalPerplexity_',dataname,'_K0_',num2str(K0),'_T_',num2str(T),'_eta0',num2str(round(eta0*1000)),...
        '_PreTrain',num2str(Settings.PreTrainYes),'_AllForPhi',num2str(Settings.AllForPhi),'_TrimStrategy',num2str(Settings.TrimStrategy),'_Trial',num2str(trial),'_Data.mat']     ;
    
    load(nameccc);
    ParaGlobal{5}.Phi;
    ['[',num2str(size(ParaGlobal{5}.Phi{1},2)),',' ...
        num2str(size(ParaGlobal{5}.Phi{2},2)),',' ...
        num2str(size(ParaGlobal{5}.Phi{3},2)),',' ...
        num2str(size(ParaGlobal{5}.Phi{4},2)),',' ...
        num2str(size(ParaGlobal{5}.Phi{5},2)),']']
    for t=1:T
        PerplexityShowTest{t}(trial,i2) = Perp(t);
    end
    
    
end

%%================================================  Reshow Figure Begin  ===========================================================



%%================================================  Save Data  Begin  ===========================================================
% T   =   8   ;
% K0Set  =   [25 50 100 200 300 400]  ;
% TrialSet    =     0:4     ;
% 
% PerplexityShowDataSet  =   cell(T,1)   ;
% PerplexityShowTrain       =   cell(T,1)   ;
% PerplexityShowTest        =   cell(T,1)   ;
% for i   =   1:T
%     PerplexityShowTrain{i} 	=   zeros(length(TrialSet),length(K0Set))   ;   
%     PerplexityShowTest{i} 	=   zeros(length(TrialSet),length(K0Set))   ;   
% end
% for Tcurrent    =   1:T
%     for kii     =   1:length(K0Set)
%         K0      =   K0Set(kii)  ;
%         for trialii     =   1:length(TrialSet)
%             trial       =   TrialSet(trialii)   ;
%             nameccc     =   ['NIPS_collapsed_trim_K0_',num2str(K0),'_T_',num2str(T),'_Tcurrent_',num2str(Tcurrent),'_Trial_',num2str(trial),'_Data_Saving.mat']     ;
%             if exist(nameccc,'file')
%                 load(nameccc)   ;
%                 avecun      =   ave     ;
%                 avecun.PhiTheta     =   []  ;
%                 PerplexityShowTrain{Tcurrent}(trialii,kii)    =  avecun.loglikeTrain(end)   ;
%                 PerplexityShowTest{Tcurrent}(trialii,kii)    =  avecun.loglike(end)   ;
%                 PerplexityShowDataSet{Tcurrent}{trialii,kii}    =   avecun ;
%             else
%                 PerplexityShowDataSet{Tcurrent}{trialii,kii}    =   [] ;
%             end
%         end
%     end
% end
% % save('NIPS_Perplextiy_Data.mat','PerplexityShowDataSet','PerplexityShowTest','PerplexityShowTrain','-v7.3')       % PerplexityShowDataSet{Tcurrent}{trialii,kii}    =   ave ;
% save('NEWS20_Perplextiy_Data.mat','PerplexityShowDataSet','-v7.3')       % PerplexityShowDataSet{Tcurrent}{trialii,kii}    =   ave ;
%%================================================  Save Data  End  ============================================================

