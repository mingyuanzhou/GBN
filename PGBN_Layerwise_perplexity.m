function [ParaGlobal,ParaLocal,Perplexity] = PGBN_Layerwise_perplexity(X_all,K,T,eta,Para)
%*************************************************************************
%Matlab code for
%[1] Mingyuan Zhou, Yulai Cong, and Bo Chen, "Augmentable gamma belief
%networks," Journal of Machine Learning Research, vol. 17, pp. 1-44, Sept.
%2016.

%[2] Mingyuan Zhou, Yulai Cong, and Bo Chen, "The Poisson gamma belief
%network," Neural Information Processing Systems (NIPS2015), Montreal,
%Canada, Dec. 2015.
%
% First version: March 2015 Second version: Sept 2015 Current version: July
% 2016
%
% Written by Mingyuan Zhou, http://mingyuanzhou.github.io/ Contributed by
% Yulai Cong, yulai_cong@163.com
%
%*************************************************************************

if ~exist('K','var')
    K = [256,128,64,32];
end
if ~exist('T','var')
    T = length(K);
end
if ~exist('eta','var')
    eta = 0.05;
end
if ~exist('Para','var')
    Para.TrainBurnin = [1000,500*ones(1,T-1)];
    Para.TrainCollection = [500,500*ones(1,T-1)];
    Para.percentage=0.3;
    Para.dataname = ['unknown_',num2str(Para.percentage*100)];
    Para.CollectionStep = 5;
    Para.train_idx = 1:size(X_all,2);
end

Para.DataType = 'Count';


%% Initialization
a0 = 0.01;   
b0 = 0.01;
e0 = 1;   
f0 = 1;
ac = 1;
bc = 1;

X = X_all(:,Para.train_idx);
[V,N] = size(X);

%Representing a word-document freqeuncy count matrix as word tokens,
%use (100*percentage)% of the tokens to form the training count matrix, 
%and use the others to form the testing count matrix.
[Xtrain,Xtest,WS,DS,WordTrainS,DocTrainS] = PartitionX_v1(X,Para.percentage);
WS = WS(WordTrainS); %word indices for training tokens
DS = DS(WordTrainS); %document indices for training tokens
ZS = DS-DS; %topic indices for training tokens

Yflagtrain = Xtrain>0;
Yflagtest = Xtest>0;
loglikeTrain = []; loglike=[];
ave.loglike=[];
ave.K=[];  
ave.gamma0=[];
Xmask=sparse(X);

Phi = cell(T,1);    
Eta = cell(T,1);
Theta = cell(T,1); 
c_j = cell(T+1,1);
for t=1:T+1
    c_j{t}=ones(1,N);
end
Xt_to_t1 = cell(T,1);   WSZS = cell(T,1);

ParaGlobal = cell(T,1);      
ParaLocal = cell(T,1);  
Perplexity = zeros(1,T);

%% GBN Gibbs sampling
for Tcurrent = 1:T
    
    ave.PhiTheta = sparse(V,N); ave.Count = 0;
    ave.PhiThetaSum = zeros(1,N);
    
    ParaGlobal{Tcurrent}.TrainTime =cputime;
    %%   Initial the Parameters of layer Tcurrent
    if Tcurrent == 1
        %Initilize the first hidden layer
        ZS = randi(K(1),length(DS),1); %Initilize the topic indices at random
        ZSDS = full(sparse(ZS,DS,1,K(1),N));
        ZSWS = full(sparse(ZS,WS,1,K(1),V));
        WSZS{1}=ZSWS';
        Xt_to_t1{1}=ZSDS;
        n_dot_k = sum(ZSDS,2);
        Eta{1} = eta(1);
        p_j = Calculate_pj(c_j,Tcurrent);
        r_k = 1/K(Tcurrent)*ones(K(Tcurrent),1);
        gamma0 = 1;       c0 = 1;
    else
        %Initilize hidden layer Tcurrent for Tcurrent>1
        K(Tcurrent) = K(Tcurrent-1); %Intilize the maximum width of a new layer using the width of the previous layer
        if K(Tcurrent)  <= 4
            break;
        end
        Eta{Tcurrent} = eta(min(length(eta),Tcurrent));
        Phi{Tcurrent} = rand(K(Tcurrent-1),K(Tcurrent));            
        Phi{Tcurrent} = bsxfun(@rdivide, Phi{Tcurrent}, max(realmin,sum(Phi{Tcurrent},1)));
        Theta{Tcurrent} = ones(K(Tcurrent),N)/K(Tcurrent);
        p_j = Calculate_pj(c_j,Tcurrent);
        r_k = 1/K(Tcurrent)*ones(K(Tcurrent),1);
        gamma0 = K(Tcurrent)/K(1);       c0 = 1;   c = 1;
    end
    
    %% Joint Learning
    for iter = 1 : (Para.TrainBurnin(Tcurrent) + Para.TrainCollection(Tcurrent) )
        tic
        %[Tcurrent,iter]
        if iter == Para.TrainBurnin(Tcurrent)
            TrimTcurrent; %Prune the inactive factors of Layer Tcurrent
        end
        
        %% ==================================== Upward Pass ===================================
        for t = 1:Tcurrent
            % ==================== Upward propogate latent counts ==================== 
            if t == 1
                % Collapsed Gibbs sampling for the first layer 
                dex111=randperm(length(ZS));
                ZS=ZS(dex111); DS=DS(dex111); WS=WS(dex111); %random permutation of the ordering of the word tokens
                if Tcurrent==1
                    shape = r_k*ones(1,N);
                else
                    shape = Phi{2}*Theta{2};
                end
                [ZSDS,ZSWS,n_dot_k,ZS] = GNBP_mex_collapsed_deep(ZSDS,ZSWS,n_dot_k,ZS,WS,DS,shape,Eta{1}(1));
                WSZS{t}=ZSWS';
                Xt_to_t1{t}=ZSDS;
            else
                % Blocked Gibbs sampling for layer t>1
                [Xt_to_t1{t},WSZS{t}] = CRT_Multrnd_Matrix(sparse(Xt_to_t1{t-1}),Phi{t},Theta{t});
            end
            % ==================== Sample Phi ========================
            if t > 1
                Phi{t} = SamplePhi(WSZS{t},Eta{t}); 
                if nnz(isnan(Phi{t}))
                    warning('Phi Nan');
                    Phi{t}(isnan(Phi{t})) = 0;
                end
            end
        end
        
        %% Sample r_k, gamma_0, and c_0
        Xt = CRT_sum_mex_matrix_v1(sparse(Xt_to_t1{Tcurrent}'),r_k')';
        [r_k,gamma0,c0]=Sample_rk(full(Xt),r_k,p_j{Tcurrent+1},gamma0,c0);
        
        %% ==================== Sample c_j and/or p_j ========================
        if iter>10
            if Tcurrent > 1
                p_j{2} = betarnd(  sum(Xt_to_t1{1},1)+a0   ,   sum(Theta{2},1)+b0  );
            else
                p_j{2} = betarnd(  sum(Xt_to_t1{1},1)+a0   ,   sum(r_k,1)+b0  );
            end
            p_j{2} = min( max(p_j{2},eps) , 1-eps);
            c_j{2} = (1-p_j{2})./p_j{2};
            for t = 3:(Tcurrent+1)
                if t == Tcurrent+1
                    c_j{t} = randg(sum(r_k)*ones(1,N)+e0) ./ (sum(Theta{t-1},1)+f0);
                else
                    c_j{t} = randg(sum(Theta{t},1)+e0) ./ (sum(Theta{t-1},1)+f0);
                end
            end
            p_j_temp = Calculate_pj(c_j,Tcurrent);
            p_j(3:end)=p_j_temp(3:end);
        end
        
        %% ==================== Sample Theta ========================
        for t = Tcurrent:-1:1
            if t == Tcurrent
                shape = r_k;
            else
                shape = Phi{t+1}*Theta{t+1};
            end
            if t > 1 || (~strcmp(Para.DataType, 'Count'))
                Theta{t} = bsxfun(@times, randg(bsxfun(@plus,shape,Xt_to_t1{t})),  1 ./ (c_j{t+1}-log(max(1-p_j{t},realmin))) );
                if nnz(isnan(Theta{t}))
                    warning('Theta Nan');
                    Theta{t}(isnan(Theta{t}))=0;
                end
            end
        end
        
        Timetmp = toc;
        
        ave.K(end+1) = nnz(sum(WSZS{1},1));
        %             ave.gamma0(iter) = gamma0;
        %             ave.eta(iter) = eta;
       
        if mod(iter,Para.CollectionStep)==0 && Para.percentage<1
            %Calculate Perplexity
            Phi{1} = SamplePhi(WSZS{1},Eta{1});
            if Tcurrent==1
                shape = r_k*ones(1,N);
            else
                shape = Phi{2}*Theta{2};
            end
            Theta{1} = bsxfun(@times, randg(bsxfun(@plus,shape,Xt_to_t1{1})),  p_j{2});
            
            X1 = Mult_Sparse(Xmask,Phi{1},Theta{1});
            X1sum = sum(Theta{1},1);
            
            X2 = bsxfun(@rdivide, X1,X1sum);
            loglike(end+1)=sum(Xtest(Yflagtest).*log(X2(Yflagtest)))/sum(Xtest(:));
            loglikeTrain(end+1)=sum(Xtrain(Yflagtrain).*log(X2(Yflagtrain)))/sum(Xtrain(:));
            
            if iter>Para.TrainBurnin(Tcurrent)
                ave.PhiTheta = ave.PhiTheta + X1;
                ave.PhiThetaSum = ave.PhiThetaSum + X1sum;
                ave.Count = ave.Count+1;
                X1 = ave.PhiTheta/ave.Count;
                X1sum = ave.PhiThetaSum/ave.Count;
                X1= bsxfun(@rdivide, X1,X1sum);
                ave.loglike(end+1) = sum(Xtest(Yflagtest).*log(X1(Yflagtest)))/sum(Xtest(:));
            else
                ave.loglike(end+1)  = NaN;
            end
                        
            clear X1 X2;
        end
        if mod(iter,100)==0
            fprintf('JointTrain Layer: %d, iter: %d, K: %d, TimePerIter: %d seconds. \n',Tcurrent,iter,nnz(sum(Xt,2)),Timetmp); 
            plot(exp(-loglikeTrain),'r');hold on
            plot(exp(-loglike),'b');hold on
            plot(exp(-ave.loglike),'k');hold off;
            drawnow
        end
        

        %%========================== One Training iteration End ==============================

    end
    ParaGlobal{Tcurrent}.Perplexitylexity=exp(-ave.loglike(end));
    Perplexity(Tcurrent) = exp(-ave.loglike(end))
    
    
    %% Save Global
    
    for t = 1:Tcurrent
        Phi{t} = SamplePhi(WSZS{t},Eta{t},true);
    end
    
    %  [r_k,gamma0,c0] = Sample_rk(full(sum(Xt,2)),r_k,p_j{Tcurrent+1},gamma0,c0,true);
    
    ParaGlobal{Tcurrent}.Phi = Phi;
    ParaGlobal{Tcurrent}.r_k = r_k;
    ParaGlobal{Tcurrent}.gamma0 = gamma0;
    ParaGlobal{Tcurrent}.c0 = c0;
    ParaGlobal{Tcurrent}.K = K(1:Tcurrent);
    ParaGlobal{Tcurrent}.Eta = Eta;
    ParaGlobal{Tcurrent}.TrainTime =cputime - ParaGlobal{Tcurrent}.TrainTime;
    
    for t = 1:Tcurrent
        tmp = sum(Theta{t},2);
        [~,ParaGlobal{Tcurrent}.Popularity{t}] = sort(tmp,'descend');
        tmp = sum(bsxfun(@rdivide,Theta{t},sum(Theta{t},1)),2);
        [~,ParaGlobal{Tcurrent}.PopularityFreq{t}] = sort(tmp,'descend');
    end
    for t = 1:Tcurrent+1
        ParaGlobal{Tcurrent}.cjmedian{t} = median(c_j{t});
    end
    
end

ave.PhiTheta = [];
ave.PhiThetaSum = [];
ParaGlobal{Tcurrent}.ave=ave;
ParaGlobal{Tcurrent}.loglike=loglike;
ParaGlobal{Tcurrent}.loglikeTrain=loglikeTrain;
