function [Accuracy_all,ParaLocal] = GBN_Testing(X,ParaGlobal,Tcurrent,Para,c_jmean)
a0 = 0.01;
b0 = 0.01;
e0 = 1;
f0 = 1;
ac=1;
bc=1;

%%
[V,N] = size(X);
switch Para.DataType
    case 'Positive'
        [ii,jj,M]=find(X>eps);    %Treat values smaller than eps as 0
        iijj=find(X>eps);
        Xmask = sparse(ii,jj,X(iijj),size(X,1),size(X,2));
    case 'Binary'
        [ii,jj,M] = find(X);
        iijj=find(X);
    case 'Count'
        %         [~,~,WS,DS,~,~]= PartitionX(X,101);
    otherwise
        error('Wrong Input Datatype!');
end
%% Initial ParaGlobal
if strcmp(Para.DataType, 'Positive')
    c = ParaGlobal.c;
    if length(c)>1
        c = mean(c)*ones(1,N);
    end
end
Phi = ParaGlobal.Phi;
r_k = ParaGlobal.r_k;
% gamma0 = ParaGlobal.gamma0;
% c0 = ParaGlobal.c0;
%% Initial ParaLocal
c_j = cell(Tcurrent+1,1);
for t=1:Tcurrent+1
    c_j{t}=ones(1,N)*c_jmean(t);
end
p_j = Calculate_pj(c_j,Tcurrent);
Theta = cell(Tcurrent,1);
for t=Tcurrent:-1:1
    if t==Tcurrent
        shape = r_k*ones(1,N);
    else
        if size(Phi{t+1},2)~=size(Theta{t+1},1)
            save('GBNTesting_wrong.mat','Phi','Theta');
        end
        shape = Phi{t+1}*Theta{t+1};
    end
    Theta{t} = max( bsxfun(@times, randg(shape) , 1./c_j{t+1}) , 1e-2 );
end
%% Initialization
%ThetaAver = cell(Tcurrent,1);
ThetaFreqAver = cell(Tcurrent,1);
c_jAver = cell(Tcurrent+1,1);
p_jAver = cell(Tcurrent+1,1);
for t = 1:(Tcurrent+1)
    if t    <= Tcurrent
        ThetaAver{t}=0;
        ThetaFreqAver{t}=0;
    end
    c_jAver{t}=0;
    p_jAver{t}=0;
end
Accuracy_all.Accuracy_Theta = [];
Xt_to_t1=cell(Tcurrent,1);
Accuracy_all.LogLikelihood =0;
%Accuracy_all.LogLikelihood_CSL=0;
%Accuracy_all.LogLikelihood_Harmonic=0;

Xmask=sparse(X);

%% Sample ParaLocal
for iter = 1 : (Para.TestBurnin + Para.TestCollection)
    tic
    %%==================================== Upward Pass ===================================
    for t = 1:Tcurrent
        if t == 1 %&& Tcurrent==1
            switch Para.DataType
                case 'Positive'
                    if length(c)==1
                        
                        if Para.ParallelProcessing
                            Rate = mtimes_par(Phi{1},Theta{1}, Para.NumBlockParallel);
                        else
                            Rate = Phi{1}*Theta{1};
                        end
                        Rate = 2*sqrt(c*X(iijj).*Rate(iijj));
                        if Para.ParallelProcessing
                            M = Truncated_bessel_rnd_par( Rate, Para.NumBlockParallel);
                        else
                            M = Truncated_bessel_rnd( Rate );
                        end
                    else
                        if Para.ParallelProcessing
                            Rate = mtimes_par(Phi{1},Theta{1}, Para.NumBlockParallel);
                        else
                            Rate = Phi{1}*Theta{1};
                        end
                        Rate = 2*sqrt(c(jj)'.*X(iijj).*Rate(iijj));
                        if Para.ParallelProcessing
                            M = Truncated_bessel_rnd_par( Rate, Para.NumBlockParallel);
                        else
                            M = Truncated_bessel_rnd( Rate );
                        end
                        c = randg(full(sparse(1,jj,M,1,N))+ac) ./ (bc+sum(X,1));
                    end
                    Xt = sparse(ii,jj,M,V,N);
                case 'Binary'
                    if Para.ParallelProcessing
                        Rate = Mult_Sparse_par(Xmask,Phi{1},Theta{1}, Para.NumBlockParallel);
                    else
                        Rate = Mult_Sparse(Xmask,Phi{1},Theta{1});
                    end
                    if Para.ParallelProcessing
                        M = Truncated_Poisson_rnd_par(full(Rate(iijj)), Para.NumBlockParallel);
                    else
                        M = truncated_Poisson_rnd(full(Rate(iijj)));
                    end
                    Xt = sparse(ii,jj,M,V,N);
                case 'Count'
                    Xt = sparse(X);
            end
            if Para.ParallelProcessing
                Xt_to_t1{t} = Multrnd_Matrix_mex_fast_v1_par(Xt,Phi{t},Theta{t}, Para.NumBlockParallel);
            else
                Xt_to_t1{t} = Multrnd_Matrix_mex_fast_v1(Xt,Phi{t},Theta{t});
            end
        else
            if Para.ParallelProcessing
                Xt_to_t1{t} = CRT_Multrnd_Matrix_par(sparse(Xt_to_t1{t-1}),Phi{t},Theta{t}, Para.NumBlockParallel);
            else
                Xt_to_t1{t} = CRT_Multrnd_Matrix(sparse(Xt_to_t1{t-1}),Phi{t},Theta{t});
            end
        end
    end
    %%==================================== Downward Pass ===================================
    %%==================== Sample Theta ========================
    if iter>10
        if Tcurrent > 1
            p_j{2} = betarnd(  sum(Xt_to_t1{1},1)+a0   ,   sum(Theta{2},1)+b0  );
        else
            p_j{2} = betarnd(  sum(Xt_to_t1{1},1)+a0   ,   sum(r_k,1)+b0  );
        end
        p_j{2} = min( max(p_j{2},realmin) , 1-realmin);
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
    
    for t = Tcurrent:-1:1
        if t == Tcurrent
            shape = r_k;
        else
            shape = Phi{t+1}*Theta{t+1};
        end
        if Para.ParallelProcessing
            temp = randg_par(bsxfun(@plus,shape,Xt_to_t1{t}), Para.NumBlockParallel);
            Theta{t} = bsxfun(@rdivide, temp ,  c_j{t+1}-log(max(1-p_j{t},realmin)) );
        else
            Theta{t} = bsxfun(@times, randg(bsxfun(@plus,shape,Xt_to_t1{t})), 1 ./ (c_j{t+1}-log(max(1-p_j{t},realmin))) );
        end
        if nnz(isnan(Theta{t}))
            warning('Theta Nan');
            Theta{t}(isnan(Theta{t}))=0;
        end
    end
    Timetmp = toc;
    if mod(iter,10)==0
        fprintf('Testing Layer: %d, iter: %d, TimePerIter: %d seconds. \n',Tcurrent,iter,Timetmp);
    end
    %%==================================== Average ===================================
    if iter > Para.TestBurnin %&& (mod(iter-Para.TestBurnin,TestSampleSpace)==0)
        for t = 1:(Tcurrent+1)
            if t <= Tcurrent
                %ThetaAver{t} 	= ThetaAver{t} + Theta{t} / Para.TestCollection;
                ThetaFreqAver{t} 	= ThetaFreqAver{t} + bsxfun(@rdivide,Theta{t},max(sum(Theta{t},1),realmin)) / Para.TestCollection;
            end
            c_jAver{t} = c_jAver{t} + c_j{t} / Para.TestCollection;
            p_jAver{t} = p_jAver{t} + p_j{t} / Para.TestCollection;
        end
        %         if  strcmp(DataType, 'Binary')
        %             Rate1=Phi{1}*Theta{1}(:,prepar.teindx);
        %             PBinary = 1-exp(-Rate1);
        %             tmp = X(:,prepar.teindx) .* log(max(realmin, PBinary) ) + (1-X(:,prepar.teindx)).*(-Rate1);
        %             LogLikelihood = sum(tmp,1);
        %             Accuracy_all.LogLikelihood_CSL = Accuracy_all.LogLikelihood_CSL + exp(LogLikelihood)/TestCollection;
        %             Accuracy_all.LogLikelihood_Harmonic = Accuracy_all.LogLikelihood_Harmonic + exp(-LogLikelihood)/TestCollection;
        %         end
    end
    
    
end     %%======================= One Testing iteration End ===========================


%%============== Calculate Accuracy ========================
if ~Para.ParallelProcessing
    %[AccuracyTmp] = TestingInside(ThetaFreqAver{1},Para);
    Accuracy_all = TestingInside(ThetaFreqAver{1},Para);
else
    %[AccuracyTmp] = TestingInside(ThetaFreqAver{1},Para,true,3);
    Accuracy_all = TestingInside(ThetaFreqAver{1},Para,true,3);
end
%Accuracy_all.Accuracy_ThetaFreqAver = AccuracyTmp.Accuracy_ThetaFreq_accuracy(1);
%Accuracy_all.Accuracy_ThetaFreqAver_AccAver = AccuracyTmp.Accuracy_ThetaFreq_AccAver;

% if  strcmp(DataType, 'Binary')
%     Accuracy_all.LogLikelihood_CSL = mean(log(max(Accuracy_all.LogLikelihood_CSL,realmin)));
%     Accuracy_all.LogLikelihood_Harmonic = mean(log(1./Accuracy_all.LogLikelihood_Harmonic));
% end

%%====== Save Local Parameters =============
%ParaLocal.ThetaAver = ThetaAver;
%ParaLocal.ThetaFreqAver = ThetaFreqAver;
ParaLocal.c_jAver = c_jAver;
ParaLocal.p_jAver = p_jAver;


