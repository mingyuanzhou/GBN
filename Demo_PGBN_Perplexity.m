dataset = 2; %[1,2]
K0 = 400; %[25,50,100,200,400,600,800]
trial = 1; %[1,2,3,4,5]

%Run the Demo Code with datset = 2, K0 in [25,50,100,200,400,600,800], and
%trial in [1,2,3,4,5] to reproduce the results used to plot Figure 12 of
%"Agumentable Gamma Belief Networks"

%Run the Demo Code with datset = 1, K0 in [25,50,100,200,400,600,800], and
%trial in [1,2,3,4,5] to reproduce the results used to plot Figure 13 of
%"Agumentable Gamma Belief Networks"

switch dataset
    case 1
        addpath('data/20news-bydate')
        load train.data;
        load test.data;
        test(:,1)=max(train(:,1))+test(:,1);
        train_test = [train;test];
        X_all =sparse(train_test(:,2),train_test(:,1),train_test(:,3));
        load train.label;
        load test.label;
        Para.trindx = 1 : length(train);
        Para.teindx = length(train) + (1: length(test));
        Para.Y = [train;test];
        
        fid=fopen('data/stop-word-list.txt');
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
        tmp = (sum(X_all,2)>=5); %words appear at least 5 times
        %tmp = (sum(X_all>0,2)>=5); %words appear in at lest 5 documents
        WO = WO(tmp);
        X_all=X_all(tmp,:);
        
        [~,dex]=sort(sum(X_all,2),'descend');
        X_all = X_all(dex(1:2000),:);
        WO = WO(dex(1:2000));
        
        Para.WO = WO;
        Para.stopwords = stopwords;
        
        clear train test dex stopwords tmp train_test WO;
        
        dataname = '20newsTop2000'
    case 2
        load data/nips12raw_str602
        X_all = counts;
        dataname = 'NIPSTop2000'
        prepar=1:size(X_all,2);
        X_all = X_all(1:2000,:);
        Para.WO = wl(1:2000);
end

%% 
eta= 0.05; %[0.01,0.05,0.1];
T = 5;
K = ones(1,T)*K0;

Para.TrainBurnin = [1000,500*ones(1,T-1)];
Para.TrainCollection = [500,500*ones(1,T-1)];

%% For comparison, you may run the same total number of iteration using a single layer model by setting the parameters as follows:
% T = 1;
% K = K0;
% Para.TrainBurnin = sum(Para.TrainBurnin)+sum(Para.TrainCollection)-500;
% Para.TrainCollection=500;

%%
Para.percentage=0.3;
Para.dataname = [dataname,'_',num2str(Para.percentage*100)];
Para.CollectionStep=5;
Para.train_idx = 1:size(X_all,2);

Para

%% Data
Para.DataType='Count';
rng(trial,'twister');
tic
[ParaGlobal,ParaLocal,Perplexity] = PGBN_Layerwise_perplexity(X_all,K,T,eta,Para);
TimeAll = toc;
name_save = ['Perplexity_',Para.dataname,'_K0_',num2str(K0),'_T_',num2str(T),'_eta',num2str(round(eta*1000)),...
   '_Trial',num2str(trial),'.mat'];
if trial==1
    save(['results2/Global',name_save],'ParaGlobal','K','T','eta','Para');
    %save(['results2/ZLocal',nameccc],'ParaLocal');
end
save(['results2/',name_save],'Perplexity');
