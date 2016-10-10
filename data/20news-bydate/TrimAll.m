if strcmp(DataType, 'Count') 
    %%=========    Count   ================
        [kk,kki,kkj] = unique(ZS); %,'stable');
        if Tcurrent==1
            gamma0=gamma0*length(kk)/K(1);
            r_k=r_k(kk);   
%             Xt  =   Xt(kk,:)    ;
        end
        K(1)=length(kk);
        ZS = kkj;
        ZSDS = full(sparse(ZS,DS,1,K(1),N));
        ZSWS = full(sparse(ZS,WS,1,K(1),V));
        n_dot_k = sum(ZSDS,2);
        WSZS{1}=ZSWS';
        Xt_to_t1{1}=ZSDS;    
        if Tcurrent > 1 
               % Eta{2}=Eta{2}(kk,:);
%                 WSZS{2}=WSZS{2}(kk,:);
                Phi{2}=Phi{2}(kk,:);
                Phi{2}=bsxfun(@rdivide,Phi{2},max(sum(Phi{2},1),realmin));
        end
        for t=2:Tcurrent
                dexK = find(sum(Xt_to_t1{t},2)==0);
                if ~isempty(dexK)
                    if t==Tcurrent
                        gamma0=gamma0*length(dexK)/K(t);
                        r_k(dexK)=[];  
%                         Xt(dexK)    =   []  ;
                    end
                    K(t)=K(t)-length(dexK);
                    Xt_to_t1{t}(dexK,:)=[];
                    Theta{t}(dexK,:)=[];
                    WSZS{t}(:,dexK)=[];
                    Phi{t}(:,dexK)=[];
                    if t<Tcurrent
                        %Eta{t+1}(dexK)=[];
%                         WSZS{t+1}(dexK,:)=[];
                        Phi{t+1}(dexK,:)=[];
                        Phi{t+1}=bsxfun(@rdivide,Phi{t+1},max(sum(Phi{t+1},1),realmin));
                    end
                end
        end
else
    %%=========    Positive  &  Binary   ================
        for t=1:Tcurrent
                dexK = find(sum(Xt_to_t1{t},2)==0);
                if ~isempty(dexK)
                    if t==Tcurrent
                        gamma0=gamma0*length(dexK)/K(t);
                        r_k(dexK)=[]; 
%                         Xt(dexK)    =   []  ;
                    end
                    K(t)=K(t)-length(dexK);
                    Xt_to_t1{t}(dexK,:)=[];
                    Theta{t}(dexK,:)=[];
                    WSZS{t}(:,dexK)=[];
                    Phi{t}(:,dexK)=[];
                    if t<Tcurrent
                        Eta{t+1}(dexK)=[];
%                         WSZS{t+1}(dexK,:)=[];
                        Phi{t+1}(dexK,:)=[];
                        Phi{t+1}=bsxfun(@rdivide,Phi{t+1},max(sum(Phi{t+1},1),realmin));
                    end
                end
        end
end