function         [ParaGlobal,ParaLocal,Accuracy_all]  =  GBN_Layerwise_IncreaseK(X_all,prepar,K,T,trial,DataType,dataname,SuPara,Settings)
%% 20150908
%
%*************************************************************************
%%%%%%%%%%%%%%%%%% Demo %%%%%%%%%%%%%%%%%%%%%%%%
% load mnist_gray.mat
% ndx = []; m = 500;
% for i = 0:1:9
%     % for i=[1,8]
%     tmp = find(train_label==i);
%     ndx = [ndx; tmp(1:m)];
% end
% X = train_mnist(:,ndx);
% Xlabel = train_label(ndx);
% label = Xlabel;
% rng(0,'twister');
% Y = sparse(label+1,1:size(X,2),1);
% Datatype = 'Positive';
%
% % Datatype = 'Binary';
% % X = double(X>rand(size(X)));

% Mingyuan Zhou,
% First version: March 2015
% Current version: May 29 2015
%*************************************************************************
%*************************************************************************
% if ~exist('eta','var')
%     eta=0.01;
% end
% if ~exist('K','var')
%     K = [256,128,64,32];
% end
% if ~exist('T','var')
%     T = length(K);
% end
% if ~exist('IsTest','var')
%     IsTest =false;
%     TestInput = [];
% end
% if ~exist('IsInferEta','var')
%     IsInferEta =true;
% end
% if ~exist('IsInfer_cj','var')
%     IsInfer_cj = true;
% end
% if ~exist('Burnin','var')
%     Burnin = 400;
% end
% if ~exist('Collection','var')
%     Collection = 100;
% end
% if ~exist('maxiter_preTrain','var')
%     maxiter_preTrain = 	200;
% end
% if min(X_all(:))<0 && strcmp(DataType,'Positive')
%     X_all(:) = (X_all(:)-min(X_all(:)))./(max(X_all(:))-min(X_all(:)));
% end

% if ~exist('Settings.nameccc','var')
%     Settings.nameccc = 'tempmat';
% end

%%   Settings
Epsilon     =   SuPara.Epsilon    ;
eta  =   SuPara.eta    ;
ac  =   SuPara.ac   ;   bc      =   SuPara.bc   ;
a0gamma  =   SuPara.a0gamma   ;   b0gamma      =   SuPara.b0gamma   ;
a0pj    =   SuPara.a0pj   ;   b0pj    =   SuPara.b0pj   ;
e0cj    =   SuPara.e0cj   ;   f0cj    =   SuPara.f0cj   ;
e0c0    =   SuPara.e0c0    ;   f0c0    =   SuPara.f0c0    ;
super1or2   =   SuPara.super1or2  ;

ParallelProcessing  =   Settings.ParallelProcessing     ;    NumBlockParallel   =   Settings.NumBlockParallel     ;
IsDisplay   =   Settings.IsDisplay    ; FigureGap  =   Settings.FigureGap      ;
SaveDataGap     =   Settings.SaveDataGap    ;       TestGap         =   Settings.TestGap     ;
IsInferEta  =   Settings.IsInferEta;     IsTestInside  =   Settings.IsTestInside     ;
PreTrainYes     =   Settings.PreTrainYes     ;      maxiter_preTrain    =   Settings.maxIter_preTrain       ;
TrainBurnin     =   Settings.TrainBurnin    ;      TrainCollection  =   Settings.TrainCollection    ;        TrainSampleSpace   =   Settings.TrainSampleSpace   ;
TestBurnin      =   Settings.TestBurnin         ;     TestCollection    =   Settings.TestCollection     ;        TestSampleSpace    =   Settings.TestSampleSpace    ;
ProcessMethod   =   Settings.ProcessMethod  ;

AllForPhi  =   Settings.AllForPhi      ;
TrimStrategy    =   Settings.TrimStrategy        ;
%  1=>TrimAllAfBurnin      % 2=>TrimTcurrentAfBurnin     %   3=>TrimAllAfEnd      %  4=>TrimTcurrentAfEnd

%% Initialization
if AllForPhi
    X = X_all   ;
else
    X = X_all(:,prepar.trindx);
end
[V,N] = size(X);
Kori    =   K(1)   ;    Tori    =   T   ;

switch DataType
    case 'Positive'
        [ii,jj,M]=find(X>eps);    %Treat values smaller than eps as 0
        iijj=find(X>eps);
        Xmask = sparse(ii,jj,X(iijj),size(X,1),size(X,2));
    case 'Binary'
        [ii,jj,M] = find(X);
        iijj=find(X);
    case 'Count'
        [~,~,WS,DS,~,~]= PartitionX_v1(X,101);
    otherwise
        error('Wrong Input Datatype!');
end

Phi = cell(T,1);    Eta = cell(T,1);
Theta = cell(T,1);   %   ThetaP = cell(T,1);    ThetaC = cell(T,1);
c_j = cell(T+1,1);
for t=1:T+1
    c_j{t}=ones(1,N);
end
Xt_to_t1 = cell(T,1);   WSZS    =   cell(T,1);
%c0t     =    ones(T,1)  ;   gamma0t     =   ones(T,1);

%% GBN Gibbs sampling
ParaGlobal  =   cell(T,1);      ParaLocal  =   cell(T,1);   Accuracy_all = cell(T,1);
MeanLogLikelihoodInCun   =   cell(T,1)   ;

K00=zeros(T,1);
for Tcurrent    =   1:T
    ParaGlobal{Tcurrent}.TrainTime =cputime;
    Phia    =   0   ;   Thetaa  =   0   ;   ca  =   0   ;
    %%   Initial Tcurrent Parameters
    if Tcurrent  ==  1
        if strcmp(DataType, 'Count')
            %Initilize the topic indices at random
            ZS = randi(K(1),length(DS),1);
            ZSDS = full(sparse(ZS,DS,1,K(1),N));
            ZSWS = full(sparse(ZS,WS,1,K(1),V));
            WSZS{1}=ZSWS';
            Xt_to_t1{1}=ZSDS;
            n_dot_k = sum(ZSDS,2);
            %Eta{Tcurrent}   =   eta(min(length(eta),Tcurrent))*ones(V,1)   ;
            Eta{1}   =   eta(1)   ;
            p_j = Calculate_pj(c_j,Tcurrent);
            r_k     =   1/K(Tcurrent)*ones(K(Tcurrent),1);
            gamma0  =   1   ;       c0  =   1   ;
        else
            %Eta{Tcurrent}   =   eta(min(length(eta),Tcurrent))*ones(V,1)   ;
            Eta{1}   =   eta(1)   ;
            Phi{Tcurrent}   =   rand(V,K(Tcurrent));            Phi{Tcurrent}   =   bsxfun(@rdivide, Phi{Tcurrent}, max(realmin,sum(Phi{Tcurrent},1)));
            Theta{Tcurrent} = ones(K(Tcurrent),N)/K(Tcurrent);
            p_j = Calculate_pj(c_j,Tcurrent);
            r_k     =   1/K(Tcurrent)*ones(K(Tcurrent),1);
            gamma0  =   1   ;       c0  =   1   ;   
            
            c = Settings.c;
            
            
            
        end
        K00(1)=K(1);
    else
        %Intilize the width of a new layer using the width of the previous layer
        if K(Tcurrent-1) >= K00(Tcurrent-1)*4/5 %&& TrimStrategy==6
            K(Tcurrent)     =   K(Tcurrent-1)*3;
        else
            K(Tcurrent)     =   K(Tcurrent-1);
        end
        K00(Tcurrent)=K(Tcurrent);
        if K(Tcurrent)  <=  4
            break   ;
        end
        %Eta{Tcurrent}   =   eta(min(length(eta),Tcurrent))*ones(K(Tcurrent-1),1)   ;
        Eta{Tcurrent}   =   eta(min(length(eta),Tcurrent)); %*ones(K(Tcurrent-1),1)   ;
        % Phi{Tcurrent}   =   dirrnd(Eta{Tcurrent}*ones(1,K(Tcurrent)));
        Phi{Tcurrent}   =   rand(K(Tcurrent-1),K(Tcurrent));            Phi{Tcurrent}   =   bsxfun(@rdivide, Phi{Tcurrent}, max(realmin,sum(Phi{Tcurrent},1)));
        Theta{Tcurrent} = ones(K(Tcurrent),N)/K(Tcurrent);
        % Phi{Tcurrent}   =   max(r_k,1e-2)/sum(max(r_k,1e-2))*ones(1,K(Tcurrent));           %Phi{Tcurrent} = bsxfun(@rdivide, Phi{Tcurrent}, sum(Phi{Tcurrent},1));
        % Theta{Tcurrent} = sum(r_k)*ones(K(Tcurrent),N)/K(Tcurrent);
        p_j = Calculate_pj(c_j,Tcurrent);
        r_k     =   1/K(Tcurrent)*ones(K(Tcurrent),1);
        gamma0  =   K(Tcurrent)/K(1)   ;       c0  =   1   ;   
        %c   =   1   ;
    end
    %% Pretrain Layer Tcurrent with the previous layers fixed.
    if PreTrainYes
        %%====================  Pretrain  ======================
        if Tcurrent   >=  2
            for iter    =   1:maxiter_preTrain
                %tic
                %%====================  Sample latent counts for the newly added layer  ======================
                %                     if ParallelProcessing
                %                         Xt  =   sparse( CRT_matrix_Parallel_cong( Xt_to_t1{Tcurrent-1}  ,  max(Phi{Tcurrent}*Theta{Tcurrent},eps)  ,  NumBlockParallel) ) ;
                %                     else
                % %                         if nnz(Xt_to_t1{Tcurrent-1})/numel(Xt_to_t1{Tcurrent-1})>0.5
                % %                             Xt = sparse(CRT_matrix( Xt_to_t1{Tcurrent-1} ,  max(Phi{Tcurrent}*Theta{Tcurrent},eps) ));
                % %                         else
                % %                             Xt = CRT_Mult_Sparse(sparse(Xt_to_t1{Tcurrent-1}),Phi,Theta);
                % %                         end
                %                         Xt = CRT_Mult_Sparse(sparse(Xt_to_t1{Tcurrent-1}),Phi{Tcurrent},Theta{Tcurrent});
                %
                %                     end
                %                     if ParallelProcessing
                %                         [Xt_to_t1{Tcurrent},WSZS{Tcurrent}]  =   Multrnd_Matrix_Parallel_cong(Xt,Phi{Tcurrent},Theta{Tcurrent},NumBlockParallel)     ;
                %                     else
                %                         [Xt_to_t1{Tcurrent},WSZS{Tcurrent}]   =   Multrnd_Matrix_mex_fast_v1(Xt,Phi{Tcurrent},Theta{Tcurrent});
                %                     end
                
                [Xt_to_t1{Tcurrent},WSZS{Tcurrent}] = CRT_Multrnd_Matrix(sparse(Xt_to_t1{Tcurrent-1}),Phi{Tcurrent},Theta{Tcurrent});
                
                %%====================  Sample Phi  ======================
                Phi{Tcurrent} = SamplePhi(WSZS{Tcurrent},Eta{Tcurrent});    %                figure(25),DispDictionaryImagesc(Phi{1});
                if nnz(isnan(Phi{Tcurrent})) || nnz(isinf(Phi{Tcurrent}))
                    Phia=1;
                end
                if nnz(isnan(Phi{Tcurrent}))
                    warning('Phi Nan');
                    Phi{Tcurrent}(isnan(Phi{Tcurrent}))   =   realmin   ;
                end
                %%====================  Sample r_k  ======================
                %[iii,jjj,sss]   =   find(Xt_to_t1{Tcurrent});
                %Xt = sparse(iii,jjj,CRT_vector(sss(:)',max(r_k(iii),eps)'),K(Tcurrent),N);
                Xt = CRT_sum_mex_matrix_v1(sparse(Xt_to_t1{Tcurrent}'),r_k')';
                [r_k,gamma0,c0]=Sample_rk(full(Xt),r_k,p_j{Tcurrent+1},gamma0,c0,false);
                %%====================  Sample Theta  ======================
                if Tcurrent==2
                    p_j{2} = betarnd(a0pj+sum(Xt_to_t1{1},1),b0pj+sum(Theta{2},1));
                    p_j{2} = min( max(p_j{2},eps) , 1-eps);
                    c_j{2} = (1-p_j{2})./p_j{2};
                    if nnz(isnan(c_j{2})) || nnz(isinf(c_j{2}))
                        ca=1;
                    end
                else
                    c_j{Tcurrent} = randg(e0cj+sum(Theta{Tcurrent},1))./(f0cj+sum(Theta{Tcurrent-1},1));
                    %                         p_j{Tcurrent} = 1./(1+c_j{Tcurrent}./(-log(max(1-p_j{Tcurrent-1},realmin))));
                    if nnz(isnan(c_j{Tcurrent})) || nnz(isinf(c_j{Tcurrent}))
                        ca=1;
                    end
                end
                c_j{Tcurrent+1} = randg(e0cj+sum(r_k)*ones(1,N))./(f0cj+sum(Theta{Tcurrent},1));
                if nnz(isnan(c_j{Tcurrent+1})) || nnz(isinf(c_j{Tcurrent+1}))
                    ca=1;
                end
                %                     p_j{Tcurrent+1}=1./(1+c_j{Tcurrent+1}./(-log(max(1-p_j{Tcurrent},realmin))));
                p_j_temp = Calculate_pj(c_j,Tcurrent);
                p_j(3:end)=p_j_temp(3:end);
                
                
                shape = r_k;
                if ParallelProcessing
                    temp     =  randg_Parallel_cong( bsxfun(@plus,shape,Xt_to_t1{Tcurrent}) , NumBlockParallel )   ;
                    Theta{Tcurrent} = bsxfun(@times, temp ,  1 ./ (c_j{Tcurrent+1}-log(max(1-p_j{Tcurrent},realmin))) );
                else
                    Theta{Tcurrent} = bsxfun(@times, randg(bsxfun(@plus,shape,Xt_to_t1{Tcurrent})),  1 ./ (c_j{Tcurrent+1}-log(max(1-p_j{Tcurrent},realmin))) );
                end     %                   figure(26),imagesc(Theta{1}),
                if nnz(isnan(Theta{Tcurrent})) || nnz(isinf(Theta{Tcurrent}))
                    Thetaa=1;
                end
                if nnz(isnan(Theta{Tcurrent}))
                    warning('Theta Nan');
                    Theta{Tcurrent}(isnan(Theta{Tcurrent}))     =   realmin ;
                end
                
                
                Timetmp     =   toc     ;
                if mod(iter,50)==0
                    fprintf('Pretain Layer: %d, iter: %d, TimePerIter: %d seconds. \n',Tcurrent,iter,Timetmp);
                end
                %%====================  One iteration End  ======================
            end
        end
    end
    %% Joint Learning
    
    
    
    for iter    =   1 : (TrainBurnin(Tcurrent) + TrainCollection(Tcurrent) * TrainSampleSpace )
        tic
        
        if TrimStrategy==5
            if iter == TrainBurnin(Tcurrent)
                TrimTcurrent     ;
            end
        end
        
        
        if TrimStrategy==6
            if (iter >=100) && (iter<=(TrainBurnin(Tcurrent) + TrainCollection(Tcurrent) -100) )
                TrimTcurrent     ;
            end
        end
        
        if TrimStrategy==1
            if iter == TrainBurnin(Tcurrent)
                TrimAll     ;
            end
        end
        
        
        %%==================================== Upward Pass ===================================
        for t   =   1:Tcurrent
            if t    ==  1 %&& Tcurrent==1
                if strcmp(DataType, 'Count')
                    dex111=randperm(length(ZS));
                    ZS=ZS(dex111);
                    DS=DS(dex111);
                    WS=WS(dex111);
                    if Tcurrent==1
                        shape = r_k*ones(1,N);
                        %[WSZS,DSZS,n_dot_k,r_k,ZS] = PFA_GNBP_partial_fixK(WSZS,DSZS,n_dot_k,r_k,ZS,WS,DS,c,gamma0,eta);
                    else
                        shape = Phi{2}*Theta{2};
                    end
                    [ZSDS,ZSWS,n_dot_k,ZS] = GNBP_mex_collapsed_deep(ZSDS,ZSWS,n_dot_k,ZS,WS,DS,shape,Eta{1}(1));
                    WSZS{t}=ZSWS';
                    Xt_to_t1{t}=ZSDS;
                else
                    if strcmp(DataType, 'Positive')
                        %%===============  Positive Data  =================
                        if length(c)==1
                            %Rate = Phi{1}*Theta{1};
                            if 1
                                if ParallelProcessing
                                    %Rate = parfun(@mtimes,NumBlockParallel,2,1,[],Phi{1},Theta{1});
                                    Rate = mtimes_par(Phi{1},Theta{1},NumBlockParallel);
                                else
                                    Rate = Phi{1}*Theta{1};
                                end
                                
                                %Rate = Phi{1}*Theta{1};
                                
                                Rate = 2*sqrt(c*X(iijj).*Rate(iijj));
                                
                                %                             Rate = Mult_Sparse(Xmask,Phi{1},Theta{1});
                                %                             Rate = nonzeros(2*sqrt(c*Xmask.*Rate));
                                %
                                %                             Rate = parfun(@Mult_Sparse,NumBlockParallel,[1,3],1,[],Xmask,Phi{1},Theta{1});
                                %                             Rate = nonzeros(2*sqrt(c*Xmask.*Rate));
                                
                                if ParallelProcessing
                                    %M  =   bessel_rnd_Parallel_cong(  2*sqrt(c*X(iijj).*Rate(iijj))  ,  Epsilon-1  ,NumBlockParallel)    ;
                                    %M = parfun(@Truncated_bessel_rnd,NumBlockParallel,reshape(Rate,1,[]));
                                    %M = parfun(@Truncated_bessel_rnd,NumBlockParallel,1,1,[],reshape(Rate,1,[]));
                                    M  = Truncated_bessel_rnd_par( Rate,NumBlockParallel);
                                else
                                    %M  =   bessel_rnd( 2*sqrt(c*X(iijj).*Rate(iijj)) , Epsilon-1) ;
                                    M  = Truncated_bessel_rnd( Rate );
                                end
                            else
                                
                                if ParallelProcessing
                                    M  = parfun(@Truncated_bessel_rnd,NumBlockParallel,[2,4],1,[], c,Xmask,Phi{1},Theta{1});
                                else
                                    M  = Truncated_bessel_rnd( c,Xmask,Phi{1},Theta{1});
                                end
                            end
                            
                            %                             if ParallelProcessing
                            %                                 %M  =   bessel_rnd_Parallel_cong(  2*sqrt(c*X(iijj).*Rate(iijj))  ,  Epsilon-1  ,NumBlockParallel)    ;
                            %                                 M = parfun(@Truncated_bessel_rnd,NumBlockParallel,reshape(2*sqrt(Rate(iijj)),1,[]));
                            %                             else
                            %                                 %M  =   bessel_rnd( 2*sqrt(c*X(iijj).*Rate(iijj)) , Epsilon-1) ;
                            %
                            %                                 M  = Truncated_bessel_rnd( 2*sqrt(c*X(iijj).*Rate(iijj)) );
                            %
                            %                             end
                        else
                            
                            if ParallelProcessing
                                %Rate = parfun(@mtimes,NumBlockParallel,2,1,[],Phi{1},Theta{1});
                                Rate = mtimes_par(Phi{1},Theta{1},NumBlockParallel);
                            else
                                Rate = Phi{1}*Theta{1};
                            end
                            
                            
                            %Rate = 2*sqrt(c(ii).*X(iijj).*Rate(iijj));
     
                            Rate = 2*sqrt(c(jj)'.*X(iijj).*Rate(iijj));
                            
                            if ParallelProcessing
                                M  = Truncated_bessel_rnd_par( Rate,NumBlockParallel);
                                % M  =   bessel_rnd_Parallel_cong( 2*sqrt(Rate(iijj))  ,  Epsilon-1  ,NumBlockParallel)    ;
                                %M = parfun(@Truncated_bessel_rnd,NumBlockParallel,reshape(2*sqrt(Rate(iijj)),1,[]));
                            else
                                %M   =   bessel_rnd( 2*sqrt(Rate(iijj)) , Epsilon-1);
                                M  = Truncated_bessel_rnd( Rate );
                                %M  = Truncated_bessel_rnd( 2*sqrt(Rate(iijj)) );
                                
                                
                                %                                                 BB = @Truncated_bessel_rnd;
                                %                                                  TotalLen = length(iijj);
                                %                                                  M = cell(0);
                                %                                                  NumBlock = 48;
                                %                                                  BlockSize = ceil(TotalLen/NumBlock);
                                %                                                  parfor block=1:NumBlock
                                %                                                      blockdex = (block-1)*BlockSize+1 : min(BlockSize,TotalLen);
                                %                                                      Rateiijj = iijj(blockdex);
                                %                                                      M{block} = feval(BB,Rateiijj);
                                %                                                  end
                                %                                                  M = cell2mat(M);
                            end
                        end
                        if length(c)==1
                            c = gamrnd(sum(M(:))+numel(X)*Epsilon+ac , 1./(bc+sum(X(:))));
                        else
                            %c = randg(full(sparse(ii,1,M,V,1))+Epsilon*sum(X>0,2)+ac) ./ (bc+sum(X,2));
                            c = randg(full(sparse(1,jj,M,1,N))+ac) ./ (bc+sum(X,1));
                        end
                    else
                        %%===============  Binary Data  =================
                        if ParallelProcessing
                            %Rate = parfun(@mtimes,NumBlockParallel,2,1,[],Phi{1},Theta{1});
                            Rate = mtimes_par(Phi{1},Theta{1},NumBlockParallel);
                        else
                            Rate = Phi{1}*Theta{1};
                        end
                        if ParallelProcessing
                            %M = truncated_Poisson_rnd_Parallel_cong(Rate(iijj), NumBlockParallel);
                            M = Truncated_Poisson_rnd_par(Rate(iijj), NumBlockParallel);
                        else
                            M = truncated_Poisson_rnd(Rate(iijj));
                        end
                    end
                    Xt = sparse(ii,jj,M,V,N);
                    if ParallelProcessing
                        %[Xt_to_t1{t},WSZS{t}]  =   Multrnd_Matrix_Parallel_cong(Xt,Phi{t},Theta{t},NumBlockParallel)     ;
                        %[Xt_to_t1{t},WSZS{t}]     =   parfun(@Multrnd_Matrix_mex_fast_v1,NumBlockParallel,[1,3],1,2,Xt,Phi{t},Theta{t});
                        [Xt_to_t1{t},WSZS{t}]     =   Multrnd_Matrix_mex_fast_v1_par(Xt,Phi{t},Theta{t},NumBlockParallel); 
                    else
                        [Xt_to_t1{t},WSZS{t}]     =   Multrnd_Matrix_mex_fast_v1(Xt,Phi{t},Theta{t});
                    end
                end
            else
                %                         if ParallelProcessing
                %                             Xt  =   sparse( CRT_matrix_Parallel_cong( Xt_to_t1{t-1}  ,  max(Phi{t}*Theta{t},eps)  ,  NumBlockParallel) ) ;
                %                         else
                %                            % Xt  =   sparse( CRT_matrix(Xt_to_t1{t-1} , max(Phi{t}*Theta{t},eps)) );
                %                             Xt = CRT_Mult_Sparse(sparse(Xt_to_t1{t-1}),Phi{t},Theta{t});
                %                         end
                %                         if ParallelProcessing
                %                             [Xt_to_t1{t},WSZS{t}]  =   Multrnd_Matrix_Parallel_cong(Xt,Phi{t},Theta{t},NumBlockParallel)     ;
                %                         else
                %                             [Xt_to_t1{t},WSZS{t}] = Multrnd_Matrix_mex_fast_v1(Xt,Phi{t},Theta{t});
                %                         end
                
                if ParallelProcessing
                    %[Xt_to_t1{t},WSZS{t}]     =   parfun(@CRT_Multrnd_Matrix,NumBlockParallel,[1,3],1,2,sparse(Xt_to_t1{t-1}),Phi{t},Theta{t});
                    [Xt_to_t1{t},WSZS{t}]     =  CRT_Multrnd_Matrix_par(sparse(Xt_to_t1{t-1}),Phi{t},Theta{t},NumBlockParallel);
                else
                    [Xt_to_t1{t},WSZS{t}] = CRT_Multrnd_Matrix(sparse(Xt_to_t1{t-1}),Phi{t},Theta{t});
                end
            end
            %%====================  Sample Eta  ========================
            if IsInferEta(t)
                %[Eta{t},gamma0t(t),c0t(t)] = SampleEta(WSZS{t},Eta{t},gamma0t(t),c0t(t));
                q_k = betarnd(sum(WSZS{t},1),Eta{t}*size(WSZS{t},1));
                Lv = CRT_sum_mex_v1(nonzeros(WSZS{t}),Eta{t});
                Eta{t} = gamrnd(0.01+Lv,1./(0.01-size(WSZS{t},1)*sum(log(max(1-q_k,realmin)))));
            end
            %%====================   Sample Phi  ========================
            if t > 1 || (~strcmp(DataType, 'Count'))
                Phi{t} = SamplePhi(WSZS{t},Eta{t}); %             figure(25),DispDictionaryImagesc(Phi{1});drawnow;
                if nnz(isnan(Phi{t})) || nnz(isinf(Phi{t}))
                    Phia=1;
                end
                if nnz(isnan(Phi{t}))
                    warning('Phi Nan');
                    Phi{t}(isnan(Phi{t}))   =   0   ;
                end
            end
        end
        %             [iii,jjj,sss]   =   find(Xt_to_t1{Tcurrent});
        %             Xt = sparse(iii,jjj,CRT_vector(sss(:)',max(r_k(iii),eps)'),K(Tcurrent),N);
        %%==================================== Downward Pass ===================================
        %[r_k,gamma0,c0]     =   Sample_rk(Xt,r_k,p_j{Tcurrent+1},gamma0,c0);
        
        Xt = CRT_sum_mex_matrix_v1(sparse(Xt_to_t1{Tcurrent}'),r_k')';
        
        
        if Tcurrent==1 || TrimStrategy~=7
            [r_k,gamma0,c0]=Sample_rk(full(Xt),r_k,p_j{Tcurrent+1},gamma0,c0);
        else
            
            if iter>=(TrainBurnin(Tcurrent) + TrainCollection(Tcurrent) * TrainSampleSpace )-50
                K_star=0;
            else
                K_star = 20;
            end
            
            c0 = randg(1+gamma0)/(1+sum(r_k));
            sumlogpi = sum(log(max(1-p_j{Tcurrent+1},realmin)));
            p_prime = -sumlogpi./(c0-sumlogpi);
            
            %gamma0 = randg(0 + CRT_sum_mex_v1(full(Xt),gamma0/KT))/(b0 - log(max(1-p_prime,realmin)));
            
            
            L_k=Xt;
            dexk = find(L_k>0);
            if K_star==0
                gamma0 = randg(0.01 + CRT_sum_mex_v1(full(Xt),gamma0/length(Xt)))/(0.01 - log(max(1-p_prime,realmin)));
               % L_k=L_k(dexk);
               % Xt = Xt(dexk);
               % Xt_to_t1{Tcurrent} = [Xt_to_t1{Tcurrent}(dexk,:);zeros(K_star,N)];
               % WSZS{Tcurrent}=[WSZS{Tcurrent}(:,dexk),zeros(K(Tcurrent-1),K_star)];
                r_k = randg(L_k+gamma0/length(Xt))/(-sumlogpi+ c0);
                if iter==(TrainBurnin(Tcurrent) + TrainCollection(Tcurrent) * TrainSampleSpace )-50
                    %TrimTcurrent;
                    TrimAll; 
                end
            else
                %gamma0 = gamrnd(0.01 + length(dexk),1/(0.01 - log(max(1-p_prime,realmin))));
                gamma0 = randg(0.01 + CRT_sum_mex_v1(full(Xt),gamma0/K(Tcurrent)))/(0.01 - log(max(1-p_prime,realmin)));

                L_k=L_k(dexk);
                Xt = Xt(dexk);
                Xt_to_t1{Tcurrent} = [Xt_to_t1{Tcurrent}(dexk,:);zeros(K_star,N)];
                WSZS{Tcurrent}=[WSZS{Tcurrent}(:,dexk),zeros(K(Tcurrent-1),K_star)];
               % r_k = [randg(L_k); randg(gamma0/K_star*ones(K_star,1))]/(-sumlogpi+ c0);
                
                r_k = [randg(L_k + gamma0/(K_star+length(L_k))); randg(gamma0/(K_star+length(L_k))*ones(K_star,1))]/(-sumlogpi+ c0);
                
                Phi{Tcurrent} = SamplePhi(WSZS{Tcurrent},Eta{Tcurrent}); %             figure(25),DispDictionaryImagesc(Phi{1});drawnow;
                if nnz(isnan(Phi{Tcurrent})) || nnz(isinf(Phi{Tcurrent}))
                    Phia=1;
                end
                if nnz(isnan(Phi{Tcurrent}))
                    warning('Phi Nan');
                    Phi{Tcurrent}(isnan(Phi{Tcurrent}))   =   0   ;
                end
                
            end
            
            %K = length(r_k);
            K(Tcurrent) = length(r_k);
            
            
            
            
        end
        
        
        
        %%====================   Sample c_j and/or p_j  ========================
        
        if iter>10
            if Tcurrent > 1
                p_j{2} = betarnd(  sum(Xt_to_t1{1},1)+a0pj   ,   sum(Theta{2},1)+b0pj  );
            else
                p_j{2} = betarnd(  sum(Xt_to_t1{1},1)+a0pj   ,   sum(r_k,1)+b0pj  );
            end
            p_j{2} = min( max(p_j{2},eps) , 1-eps);
            c_j{2} = (1-p_j{2})./p_j{2};
            if nnz(isnan(c_j{2})) || nnz(isinf(c_j{2}))
                ca=1;
            end
            for t   =   3:(Tcurrent+1)
                if t    ==  Tcurrent+1
                    c_j{t} = randg(sum(r_k)*ones(1,N)+e0cj) ./ (sum(Theta{t-1},1)+f0cj);
                else
                    c_j{t} = randg(sum(Theta{t},1)+e0cj) ./ (sum(Theta{t-1},1)+f0cj);
                end
                if nnz(isnan(c_j{t})) || nnz(isinf(c_j{t}))
                    ca=1;
                end
                %                         p_j{t} = 1 ./ (1 + c_j{t}./(-log(max(1-p_j{t-1},realmin))) );
            end
            p_j_temp = Calculate_pj(c_j,Tcurrent);
            p_j(3:end)=p_j_temp(3:end);
        end
        
        %%====================   Sample Theta  ========================
        for t  =   Tcurrent:-1:1
            if t    ==  Tcurrent
                shape = r_k;
            else
                shape = Phi{t+1}*Theta{t+1};
            end
            if t > 1 || (~strcmp(DataType, 'Count'))
                if ParallelProcessing
                    %temp     =  randg_Parallel_cong( bsxfun(@plus,shape,Xt_to_t1{t}) , NumBlockParallel )   ;
                    
                   % temp =  parfun(@randg,NumBlockParallel,1,1,[],bsxfun(@plus,shape,Xt_to_t1{t}));
                    temp =  randg_par(bsxfun(@plus,shape,Xt_to_t1{t}),NumBlockParallel);
                    
                    Theta{t} = bsxfun(@times, temp ,  1 ./ (c_j{t+1}-log(max(1-p_j{t},realmin))) );
                else
                    Theta{t} = bsxfun(@times, randg(bsxfun(@plus,shape,Xt_to_t1{t})),  1 ./ (c_j{t+1}-log(max(1-p_j{t},realmin))) );
                end         %            figure(26),imagesc(Theta{1}),drawnow
                %                         figure(10),hold on,plot(iter,max(max(Theta{t})),'*')
                if nnz(isnan(Theta{t})) || nnz(isinf(Theta{t}))
                    Thetaa=1;
                end
                if nnz(isnan(Theta{t}))
                    warning('Theta Nan');
                    Theta{t}(isnan(Theta{t}))=0;
                end
            end
        end
        
        Timetmp     =   toc     ;
         if mod(iter,10)==0
                fprintf('JointTrain Layer: %d, iter: %d, K: %d, TimePerIter: %d seconds. \n',Tcurrent,iter,nnz(sum(Xt,2)),Timetmp);
          end
        %%====================  Figure  ====================
        if IsDisplay  &&  (mod(iter,FigureGap)==0) && (~strcmp(DataType, 'Count'))
            tmp     =   1   ;
            for t   =   1:Tcurrent
                tmp     =   tmp * Phi{t}    ;
                if strcmp(DataType,'Binary')
                    %figure(90+t),DispDictionaryImagesc(1-exp(-tmp));drawnow;
                    if V==128
                        figure(80+t);DispDictionary_asym(1-exp(-tmp),16,8);drawnow;
                    else
                        figure(80+t);DispDictionary(1-exp(-tmp));drawnow;
                    end
                else
                    %figure(90+t),DispDictionaryImagesc(tmp);drawnow;
                    if V==128
                        figure(90+t);DispDictionary_asym(tmp,16,8);drawnow;
                    else
                        figure(90+t);DispDictionary(tmp);drawnow;
                    end
                end
            end
        end
        %             tmp     =   1   ;tmp2 = 1;
        %             for t   =   1:Tcurrent
        %                 tmp     =   tmp * Phi{t}    ;
        %                 tmp2    =   tmp2 .* c_j{t+1} ;
        %             end
        %             X_re = bsxfun(@rdivide, tmp * Theta{Tcurrent},tmp2) ;
        %             figure(11),hold on,plot(iter,max(sum((X-X_re).^2)))
        %%==========================  One Training iteration End  ==============================
    end
    
    %
    %         %%================   TrimStrategy   ======================
    %         %  1=>TrimAllAfBurnin      % 2=>TrimTcurrentAfBurnin     %   3=>TrimAllAfEnd      %  4=>TrimTcurrentAfEnd
    if TrimStrategy==2
        %                 if iter == TrainBurnin(Tcurrent)
        %                     TrimTcurrent     ;
        %                 end
    end
    %         switch TrimStrategy
    %             case 1  %% TrimStrategy == 1
    %
    %                     TrimAll  ;
    %
    %             case 2  %% TrimStrategy == 2
    %
    %                     TrimTcurrent     ;
    %
    %             case 3  %% TrimStrategy == 3
    %                 %if iter == (TrainBurnin(Tcurrent) + TrainCollection(Tcurrent) * TrainSampleSpace )-50
    %
    %                     TrimAll     ;
    %
    %             case 4  %% TrimStrategy == 4
    %                 %if iter == (TrainBurnin(Tcurrent) + TrainCollection(Tcurrent) * TrainSampleSpace )-50
    %
    %                     TrimTcurrent    ;
    %
    %         end
    
    if TrimStrategy==4
        TrimTcurrent    ;
    end
    
    
    %% Save Global
    if strcmp(DataType, 'Positive')
        if length(c)==1
            c = (sum(M(:))+numel(X)*Epsilon+ac) ./ (bc+sum(X(:))) ;
        else
            %c = (full(sparse(ii,1,M,V,1))+Epsilon*sum(X>0,2)+ac) ./ (bc+sum(X,2));
            c = (full(sparse(1,jj,M,1,N))+ac) ./ (bc+sum(X,1));
        end
    end
    
    %         for t = 1:Tcurrent
    %             if IsInferEta(t)
    %                 q_k = betarnd(sum(WSZS{t},1),Eta{t}*size(WSZS{t},1));
    %                 Lv = CRT_sum_mex(nonzeros(WSZS{t}),Eta{t});
    %                 Eta{t} = gamrnd(0.01+Lv,1./(0.01-size(WSZS{t},1)*sum(log(max(1-q_k,realmin)))));
    %                 % [Eta{t},gamma0t(t),c0t(t)] = SampleEta(WSZS{t},Eta{t},gamma0t(t),c0t(t),true);
    %             end
    %         end
    for t = 1:Tcurrent
        Phi{t} = SamplePhi(WSZS{t},Eta{t},true);
    end
    
    %  [r_k,gamma0,c0]     =   Sample_rk(full(sum(Xt,2)),r_k,p_j{Tcurrent+1},gamma0,c0,true);
    
    if strcmp(DataType, 'Positive')
        ParaGlobal{Tcurrent}.c    =   c   ;
    end
    ParaGlobal{Tcurrent}.Phi    =   Phi   ;
    ParaGlobal{Tcurrent}.r_k    =   r_k   ;
    ParaGlobal{Tcurrent}.gamma0    =   gamma0   ;
    ParaGlobal{Tcurrent}.c0    =   c0   ;
    ParaGlobal{Tcurrent}.K    =   K(1:Tcurrent)   ;
    ParaGlobal{Tcurrent}.Eta    =   Eta     ;
    %         ParaGlobal{Tcurrent}.gamma0t    =   gamma0t     ;
    %         ParaGlobal{Tcurrent}.c0t    =   c0t     ;
    ParaGlobal{Tcurrent}.TrainTime =cputime - ParaGlobal{Tcurrent}.TrainTime;
    
    for t   =   1:Tcurrent
        %tmp     =   sum(ThetaPAver{t},2)  ;
        tmp     =   sum(Theta{t},2)  ;
        [~,ParaGlobal{Tcurrent}.Popularity{t}]    =   sort(tmp,'descend')     ;
        tmp     =   sum(bsxfun(@rdivide,Theta{t},sum(Theta{t},1)),2)  ;
        [~,ParaGlobal{Tcurrent}.PopularityFreq{t}]    =   sort(tmp,'descend')     ;
    end
    
    for t   =   1:Tcurrent+1
        ParaGlobal{Tcurrent}.cjmedian{t}  = median(c_j{t});
    end
    
    %% Testing
    ParaGlobal{Tcurrent}.TestTime =cputime;
    if IsTestInside
        c_jmean = zeros(1,Tcurrent+1);
        for t   =   1:(Tcurrent+1)
            c_jmean(t)  =    median(c_j{t})   ;
        end
        [Accuracy_all{Tcurrent},ParaLocal{Tcurrent}]   =   GBN_Testing(X_all,prepar,ParaGlobal{Tcurrent},Tcurrent,DataType,c_jmean,SuPara,Settings)  ;
    end
    ParaGlobal{Tcurrent}.TestTime =cputime - ParaGlobal{Tcurrent}.TestTime;
    ParaLocal{Tcurrent}.Phia    =   Phia    ;
    ParaLocal{Tcurrent}.Thetaa    =   Thetaa    ;
    ParaLocal{Tcurrent}.ca    =   ca    ;
    
    
    save(['results2/',Settings.nameccc],'Accuracy_all');

    
end  %%=============== Tcurrent End =========================




