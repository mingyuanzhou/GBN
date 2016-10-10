function         [ParaGlobal,ParaLocal,Accuracy_all] = GBN_Layerwise(X_all,K,T,eta,Para)
%*************************************************************************
%Matlab code for
%[1] Mingyuan Zhou, Yulai Cong, and Bo Chen, "Augmentable gamma belief
%networks," Journal of Machine Learning Research, vol. 17, pp. 1-44, Sept.
%2016.

%[2] Mingyuan Zhou, Yulai Cong, and Bo Chen, "The Poisson gamma belief
%network," Neural Information Processing Systems (NIPS2015), Montreal,
%Canada, Dec. 2015.

% First version: March 2015
% Second version: Sept 2015
% Current version: July 2016
%
% Written by Mingyuan Zhou, http://mingyuanzhou.github.io/
% Contributed by Yulai Cong, yulai_cong@163.com
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
    Para.percentage=101;
    Para.dataname = 'unknown';
    Para.CollectionStep = 5;
    Para.train_idx = 1:size(X_all,2);
    Para.DataType = 'Count';
    
    Para.IsDisplay = 0;          Para.FigureGap = 1000;
    Para.SaveDataGap = 200;   Para.TestGap = 1e5;
    Para.IsTestInside = 0;
    Para.ParallelProcessing = 0;  Para.NumBlockParallel = 8;
    Para.TestBurnin = 500;       Para.TestCollection = 500;
end


%% Initialization
a0 = 0.01;   
b0 = 0.01;
e0 = 1;   
f0 = 1;
ac=1;
bc=1;

Epsilon = 0;


%% Initialization
X = X_all(:,Para.train_idx);
[V,N] = size(X);
switch Para.DataType
    case 'Positive'
        [ii,jj,M]=find(X>eps);    %Treat values smaller than eps as 0
        iijj=find(X>eps);
        Xmask = sparse(ii,jj,X(iijj),size(X,1),size(X,2));
    case 'Binary'
        [ii,jj,M] = find(X);
        iijj=find(X);
        Xmask = sparse(ii,jj,M,V,N);
    case 'Count'
        [~,~,WS,DS,~,~]= PartitionX_v1(X,101);
    otherwise
        error('Wrong Input Datatype!');
end
Phi = cell(T,1);    
Eta = cell(T,1);
Theta = cell(T,1); 
c_j = cell(T+1,1);
for t=1:T+1
    c_j{t}=ones(1,N);
end
Xt_to_t1 = cell(T,1);   WSZS = cell(T,1);

%% GBN Gibbs sampling
ParaGlobal = cell(T,1);      
ParaLocal = cell(T,1);  
Accuracy_all = cell(T,1);

for Tcurrent = 1:T
    ParaGlobal{Tcurrent}.TrainTime = cputime;
    %%   Initial Tcurrent Parameters
    if Tcurrent == 1
        %Initilize the first hidden layer
        if strcmp(Para.DataType, 'Count')
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
            Eta{1} = eta(1);
            Phi{Tcurrent} = rand(V,K(Tcurrent));            
            Phi{Tcurrent} = bsxfun(@rdivide, Phi{Tcurrent}, max(realmin,sum(Phi{Tcurrent},1)));
            Theta{Tcurrent} = ones(K(Tcurrent),N)/K(Tcurrent);
            p_j = Calculate_pj(c_j,Tcurrent);
            r_k = 1/K(Tcurrent)*ones(K(Tcurrent),1);
            gamma0 = 1;       c0 = 1;   
            c = Para.c;
        end
    else
        %Initilize hidden layer Tcurrent for Tcurrent>1
        K(Tcurrent) = K(Tcurrent-1); %Intilize the maximum width of a new layer using the width of the previous layer
        if K(Tcurrent)  <= 4
            break;
        end
        Eta{Tcurrent} = eta(min(length(eta),Tcurrent));
        % Phi{Tcurrent} = dirrnd(Eta{Tcurrent}*ones(1,K(Tcurrent)));
        Phi{Tcurrent} = rand(K(Tcurrent-1),K(Tcurrent));            
        Phi{Tcurrent} = bsxfun(@rdivide, Phi{Tcurrent}, max(realmin,sum(Phi{Tcurrent},1)));
        Theta{Tcurrent} = ones(K(Tcurrent),N)/K(Tcurrent);
        p_j = Calculate_pj(c_j,Tcurrent);
        r_k = 1/K(Tcurrent)*ones(K(Tcurrent),1);
        gamma0 = K(Tcurrent)/K(1);       c0 = 1;   
        %c = 1;
    end
    
    %% Joint Learning
    for iter = 1 : (Para.TrainBurnin(Tcurrent) + Para.TrainCollection(Tcurrent) )
        tic
        if iter == Para.TrainBurnin(Tcurrent)
            TrimTcurrent; %Prune the inactive factors of Layer Tcurrent
        end
        
        %% ==================================== Upward Pass ===================================
        for t = 1:Tcurrent
            if t == 1 
                if strcmp(Para.DataType, 'Count')
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
                    if strcmp(Para.DataType, 'Positive')
                        %%=============== Positive Data =================
                        if length(c)==1
                            %using the same scale for all data in the
                            %Poisson randomized gamma distributions
                            if Para.ParallelProcessing
                                Rate = mtimes_par(Phi{1},Theta{1},NumBlockParallel);
                            else
                                Rate = Phi{1}*Theta{1};
                            end
                            Rate = 2*sqrt(c*X(iijj).*Rate(iijj));
                            %                             Rate = Mult_Sparse(Xmask,Phi{1},Theta{1});
                            %                             Rate = nonzeros(2*sqrt(c*Xmask.*Rate));
                            if Para.ParallelProcessing
                                M = Truncated_bessel_rnd_par( Rate,NumBlockParallel);
                            else
                                %M = bessel_rnd( 2*sqrt(c*X(iijj).*Rate(iijj)) , Epsilon-1);
                                M = Truncated_bessel_rnd( Rate );
                            end
                        else
                            %using sample-specific scales in the
                            %Poisson randomized gamma distributions
                            if Para.ParallelProcessing
                                Rate = mtimes_par(Phi{1},Theta{1},NumBlockParallel);
                            else
                                Rate = Phi{1}*Theta{1};
                            end
                            Rate = 2*sqrt(c(jj)'.*X(iijj).*Rate(iijj));
                            if Para.ParallelProcessing
                                M = Truncated_bessel_rnd_par( Rate,NumBlockParallel);
                            else
                                %M = bessel_rnd( 2*sqrt(Rate(iijj)) , Epsilon-1);
                                M = Truncated_bessel_rnd( Rate );
                            end
                        end
                        if length(c)==1
                            c = gamrnd(sum(M(:))+numel(X)*Epsilon+ac , 1./(bc+sum(X(:))));
                        else
                            %c = randg(full(sparse(ii,1,M,V,1))+Epsilon*sum(X>0,2)+ac) ./ (bc+sum(X,2));
                            c = randg(full(sparse(1,jj,M,1,N))+ac) ./ (bc+sum(X,1));
                        end
                    else
                        %%=============== Binary Data =================
                        if Para.ParallelProcessing
                            %Rate = mtimes_par(Phi{1},Theta{1},NumBlockParallel);
                            Rate = Mult_Sparse_par(Xmask,Phi{1},Theta{1});
                        else
                            %Rate = Phi{1}*Theta{1};
                            Rate = Mult_Sparse(Xmask,Phi{1},Theta{1});
                        end
                        %Rate = Rate(iijj);
                        
                        if Para.ParallelProcessing
                            %M = truncated_Poisson_rnd_Parallel_cong(Rate(iijj), Para.NumBlockParallel);
                            M = Truncated_Poisson_rnd_par(full(Rate(iijj)), Para.NumBlockParallel);
                        else
                            M = truncated_Poisson_rnd(full(Rate(iijj)));
                        end
                    end
                    Xt = sparse(ii,jj,M,V,N);
                    if Para.ParallelProcessing
                        %[Xt_to_t1{t},WSZS{t}] = Multrnd_Matrix_Parallel_cong(Xt,Phi{t},Theta{t},NumBlockParallel);
                        %[Xt_to_t1{t},WSZS{t}] = parfun(@Multrnd_Matrix_mex_fast_v1,NumBlockParallel,[1,3],1,2,Xt,Phi{t},Theta{t});
                        [Xt_to_t1{t},WSZS{t}] = Multrnd_Matrix_mex_fast_v1_par(Xt,Phi{t},Theta{t},NumBlockParallel); 
                    else
                        [Xt_to_t1{t},WSZS{t}] = Multrnd_Matrix_mex_fast_v1(Xt,Phi{t},Theta{t});
                    end
                end
            else
                % Blocked Gibbs sampling for layer t>1
                if Para.ParallelProcessing
                    %[Xt_to_t1{t},WSZS{t}] = parfun(@CRT_Multrnd_Matrix,NumBlockParallel,[1,3],1,2,sparse(Xt_to_t1{t-1}),Phi{t},Theta{t});
                    [Xt_to_t1{t},WSZS{t}] = CRT_Multrnd_Matrix_par(sparse(Xt_to_t1{t-1}),Phi{t},Theta{t},NumBlockParallel);
                else
                    [Xt_to_t1{t},WSZS{t}] = CRT_Multrnd_Matrix(sparse(Xt_to_t1{t-1}),Phi{t},Theta{t});
                end
            end
            %% ==================== Sample Phi ========================
            if t > 1 || (~strcmp(Para.DataType, 'Count'))
                Phi{t} = SamplePhi(WSZS{t},Eta{t});            
                if nnz(isnan(Phi{t}))
                    warning('Phi Nan');
                    Phi{t}(isnan(Phi{t})) = 0;
                end
            end
        end

        %%==================================== Downward Pass ===================================
        
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
                %p_j{t} = 1 ./ (1 + c_j{t}./(-log(max(1-p_j{t-1},realmin))) );
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
                if Para.ParallelProcessing
                    temp = randg_par(bsxfun(@plus,shape,Xt_to_t1{t}),NumBlockParallel);
                    Theta{t} = bsxfun(@times, temp ,  1 ./ (c_j{t+1}-log(max(1-p_j{t},realmin))) );
                else
                    Theta{t} = bsxfun(@times, randg(bsxfun(@plus,shape,Xt_to_t1{t})),  1 ./ (c_j{t+1}-log(max(1-p_j{t},realmin))) );
                end        
                if nnz(isnan(Theta{t}))
                    warning('Theta Nan');
                    Theta{t}(isnan(Theta{t}))=0;
                end
            end
        end
        Timetmp = toc;
        if mod(iter,10)==0
            fprintf('JointTrain Layer: %d, iter: %d, K: %d, TimePerIter: %d seconds. \n',Tcurrent,iter,nnz(sum(Xt,2)),Timetmp);
        end
        %%==================== Figure ====================
        if Para.IsDisplay  &&  (mod(iter,Para.FigureGap)==0) && (~strcmp(Para.DataType, 'Count'))
            tmp=1;
            for t = 1:Tcurrent
                tmp = tmp*Phi{t};
                if strcmp(Para.DataType,'Binary')
                    if V==128
                        figure(80+t);DispDictionary_asym(1-exp(-tmp),16,8);drawnow;
                    else
                        figure(80+t);DispDictionary(1-exp(-tmp));drawnow;
                    end
                else
                    if V==128
                        figure(90+t);DispDictionary_asym(tmp,16,8);drawnow;
                    else
                        figure(90+t);DispDictionary(tmp);drawnow;
                    end
                end
            end
        end
        %%========================== One Training iteration End ==============================
    end
    
    %% Save Global
    if strcmp(Para.DataType, 'Positive')
        if length(c)==1
            c = (sum(M(:))+numel(X)*Epsilon+ac) ./ (bc+sum(X(:)));
        else
            c = (full(sparse(1,jj,M,1,N))+ac) ./ (bc+sum(X,1));
        end
    end
    
    for t = 1:Tcurrent
        Phi{t} = SamplePhi(WSZS{t},Eta{t},true);
    end
    
    %  [r_k,gamma0,c0] = Sample_rk(full(sum(Xt,2)),r_k,p_j{Tcurrent+1},gamma0,c0,true);
    
    if strcmp(Para.DataType, 'Positive')
        ParaGlobal{Tcurrent}.c = c;
    end
    ParaGlobal{Tcurrent}.Phi = Phi;
    ParaGlobal{Tcurrent}.r_k = r_k;
    ParaGlobal{Tcurrent}.gamma0 = gamma0;
    ParaGlobal{Tcurrent}.c0 = c0;
    ParaGlobal{Tcurrent}.K = K(1:Tcurrent);
    ParaGlobal{Tcurrent}.Eta = Eta;
    ParaGlobal{Tcurrent}.TrainTime =cputime - ParaGlobal{Tcurrent}.TrainTime;

    for t = 1:Tcurrent+1
        ParaGlobal{Tcurrent}.cjmedian{t} = median(c_j{t});
    end
    
    %% Testing
    ParaGlobal{Tcurrent}.TestTime =cputime;
    if Para.IsTestInside
        c_jmean = zeros(1,Tcurrent+1);
        for t = 1:(Tcurrent+1)
            c_jmean(t) = median(c_j{t});
        end
        [Accuracy_all{Tcurrent},ParaLocal{Tcurrent}] = GBN_Testing(X_all,ParaGlobal{Tcurrent},Tcurrent,Para,c_jmean);
    end
    ParaGlobal{Tcurrent}.TestTime =cputime - ParaGlobal{Tcurrent}.TestTime;
    
    save(['results2/',Para.name_save],'Accuracy_all');

    
end  %%=============== Tcurrent End =========================




