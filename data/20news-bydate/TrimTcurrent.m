if strcmp(DataType, 'Count')  && (Tcurrent==1)
        [kk,kki,kkj] = unique(ZS); %,'stable');
        gamma0=gamma0*length(kk)/K(1);
        r_k=r_k(kk);
%         Xt  =   Xt(kk,:)    ;
        K(1)=length(kk);
        ZS = kkj;
        ZSDS = full(sparse(ZS,DS,1,K(1),N));
        ZSWS = full(sparse(ZS,WS,1,K(1),V));
        n_dot_k = sum(ZSDS,2);
        WSZS{1}=ZSWS';
        Xt_to_t1{1}=ZSDS;
end
if (~strcmp(DataType, 'Count'))  || (Tcurrent > 1)
        for t=Tcurrent:Tcurrent
            dexK = find(sum(Xt_to_t1{t},2)==0);
            if ~isempty(dexK)
                gamma0=gamma0*length(dexK)/K(t);
                r_k(dexK)=[];
%                 Xt(dexK)    =   []  ;
                K(t)=K(t)-length(dexK);
                Xt_to_t1{t}(dexK,:)=[];
                Theta{t}(dexK,:)=[];        %  ThetaP{t}(dexK,:)    =   [];     ThetaC{t}(dexK,:)    =   [];     
                WSZS{t}(:,dexK)=[];
                Phi{t}(:,dexK)=[];
            end
        end
end