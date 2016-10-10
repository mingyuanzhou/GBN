%load('/Users/zhoum/Box Sync/GBN_results/results_PRG_GBN/ZGlobalMnist_K0_200_T_5_eta050_PreTrain0_AllForPhi0_TrimStrategy5_Trial1_Data.mat')
load('/Users/zhoum/Box Sync/GBN_results/results_PRG_GBN/ZGlobalMnist_CJ_K0_50_T_5_eta050_PreTrain0_AllForPhi0_TrimStrategy5_Trial1_Data')
Phi=ParaGlobal{5}.Phi;
Phit=1;
Phisort=Phi;
rsort = cell(5,1);
for layer=1:5
    Phit=Phit*Phi{layer};
    r_k = ParaGlobal{5}.r_k;
    for tt=4:-1:layer
        r_k = Phi{tt+1}*r_k;
    end
    
    [rsort{layer},dex]=sort(r_k,'descend');
    %dex = dex (rsort{layer}>0.0001);
    %[~,dex]=sort(sum(Theta{layer},2),'descend');
   % [rsort{layer},dex]=sort(sum(Theta{layer},2),'descend');
    
    
    
%     if layer==1
%         [~,dex]=sort(Phi{2}*Phi{3}*Phi{4}*Phi{5}*ParaGlobal{5}.r_k,'descend');
%     else
%         dex = ParaGlobal{5}.Popularity{layer};
%     end
    Phisort{layer} = Phisort{layer}(:,dex);
    if layer<5
        Phisort{layer+1} = Phisort{layer+1}(dex,:);
    end
figure;DispDictionary(Phit(:,dex));
end