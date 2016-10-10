function ThetaAve=PGBN_test_new(X,T,TestInput,Burnin,Collection,e0,f0,e1,f1,IsPj2,pe0,pf0)

if nargin<10
    IsPj2=false;
end
Phi = TestInput.Phi;
r_k = TestInput.r_k;
N=size(X,2);
K = zeros(1,T);
for t=1:T
    K(t) = size(Phi{t},2);
end
Theta = cell(T,1);
c_j = cell(T+1,1);

for t=T:-1:1
    if t==T
        shape = r_k*ones(1,N);
    else
        shape = Phi{t+1}*Theta{t+1};
    end
    Theta{t} = max(randg(shape)./mean(TestInput.c_j{t+1}),1e-2);
end
for t=1:T+1
    c_j{t}=ones(1,N)*mean(TestInput.c_j{t});
end



Xt_to_t1 = cell(T,1);
WSZS=cell(T,1);

ThetaAve=cell(T,1);


text = [];
c0=1;
fprintf('\n Iteration: ');

a0=0.01;
b0=0.01;
%% Gibbs sampling
for iter=1:Burnin+Collection
    
    %Calculate p_j
    p_j = Calculate_pj(c_j,T);
    
    %% Upward Pass
    Xt = sparse(X);
    
    for t=1:T
        %Poisson factor analysis
        [Xt_to_t1{t},WSZS{t}] = Multrnd_Matrix_mex_fast(Xt,Phi{t},Theta{t});
        %Propogate counts to the next layer
        if t<T
            Xt = sparse(CRT_matrix(Xt_to_t1{t},max(Phi{t+1}*Theta{t+1},eps)));
        else
            [iii,jjj,sss]=find(Xt_to_t1{T});
            Xt = sparse(iii,jjj,CRT_vector(sss(:)',max(r_k(iii),eps)'),K(T),N);
        end
    end
    
    Tcurrent=T;
    for t=Tcurrent:-1:1
        
        if t==Tcurrent
            if Tcurrent>=2
                c_j{Tcurrent+1} = randg(e1+sum(r_k)*ones(1,N))./(f1+sum(Theta{Tcurrent},1));
                p_j{Tcurrent+1}=1./(1+c_j{Tcurrent+1}./(-log(max(1-p_j{Tcurrent},realmin))));
            else
                p_j{2} = betarnd(pe0+sum(Xt_to_t1{1},1),pf0+sum(r_k));
                c_j{2} = (1-p_j{2})./p_j{2};
            end

        elseif t>=2
            c_j{t+1} = randg(e1+sum(Theta{t+1},1))./(f1+sum(Theta{t},1));
            p_j{t+1}=1./(1+c_j{t+1}./(-log(max(1-p_j{t},realmin))));
        else
            p_j{2} = betarnd(pe0+sum(Xt_to_t1{1},1),pf0+sum(Theta{2},1));
            c_j{2} = (1-p_j{2})./p_j{2};
        end
        
        
        
        if t==Tcurrent
            shape = r_k;
        else
            shape = Phi{t+1}*Theta{t+1};
        end
        if t>1
            %Do not need to calculate the full Theta{1} for collapsed inference
            Theta{t} = bsxfun(@rdivide, randg(bsxfun(@plus,shape,Xt_to_t1{t})),  c_j{t+1}-log(max(1-p_j{t},realmin)));
            
        end
        
        %    end
        if nnz(isnan(Theta{t}))
            warning('Theta Nan');
            Theta{t}(isnan(Theta{t}))=1e-6;
        end
        
    end
    
    Theta{1} = bsxfun(@times, randg(bsxfun(@plus,shape,Xt_to_t1{1})),  p_j{2});
    
    
    
    
    if 0
        if iter>Burnin
            if iter==Burnin+1
                for t=1:T
                    ThetaAve{t}=Theta{t}/Collection;  %%%This is correct if Phi are fixed.
                end
            else
                for t=1:T
                    ThetaAve{t}=ThetaAve{t} + Theta{t}/Collection; %%%This is correct if Phi are fixed.
                end
            end
        end
    else
        if iter>Burnin
            if iter==Burnin+1
                for t=1:T
                    ThetaAve{t}=bsxfun(@rdivide,Theta{t},max(sum(Theta{t},1),realmin))/Collection;  %%%This is correct if Phi are fixed.
                end
            else
                for t=1:T
                    ThetaAve{t}=ThetaAve{t} + bsxfun(@rdivide,Theta{t},max(sum(Theta{t},1),realmin))/Collection; %%%This is correct if Phi are fixed.
                end
            end
        end
    end
    
    fprintf(repmat('\b',1,length(text)));
    text = sprintf('%d',iter);
    fprintf(text, iter);
    
end

