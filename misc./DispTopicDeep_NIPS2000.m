
eta0 = 0.1;
%eta0 = 0.1;
if 1

    K=200; 
  KKK=K;  
    load('/Users/zhoum/Box Sync/GBN_results/results1_global/ZGlobalPerplexity_NIPSTop2000_30_K0_200_T_5_eta0100_PreTrain0_AllForPhi1_TrimStrategy5_Trial1_Data');
 Phi = ParaGlobal{5}.Phi;
   % wl=prepar.WO;
      
    i1=2
switch i1
    case 1
        SystemInRuningLinux=1;
        ToBeAnalized    =   9   ;
        IsBinaryClassificaiton  =   false   ;
        dataset=5;
        % IsBinaryClassificaiton  =  true  ;
        % dataset=1;
        LoadData_GBN    ;
        dataname = '20newsTop2000'
    case 2
        load /Users/zhoum/Dropbox/Zhou_Code/GammaBeliefNet_v1/data/nips12raw_str602
        X_all = counts;
        dataname = 'NIPSTop2000'
        prepar=1:size(X_all,2);
        X_all = X_all(1:2000,:);
        wl=wl(1:2000);
end

    
   % load(['/Users/zhoum/Box Sync/GBN_results/results1_local/ZLocalTrim_20newsBydate_K0_',num2str(K),'_T_5_eta0',num2str(eta0*1000),'_PreTrain0_AllForPhi0_TrimStrategy5_Trial1_Data.mat'])
else
    K=1024
    load(['/Users/zhoum/Box Sync/GBN_results/results1_global/ZGlobalTrim_20newsBydateTop2000_K0_',num2str(K),'_T_5_eta0',num2str(eta0*1000),'_PreTrain0_AllForPhi0_TrimStrategy5_Trial1_Data.mat'])
    Phi = ParaGlobal{5}.Phi;
    wl=prepar.WO;
    
   % load(['/Users/zhoum/Box Sync/GBN_results/results1_local/ZLocalTrim_20newsBydateTop2000_K0_',num2str(K),'_T_5_eta0',num2str(eta0*1000),'_PreTrain0_AllForPhi0_TrimStrategy5_Trial1_Data.mat'])
    
end

addpath('/Users/zhoum/Dropbox/Zhou_Code/Deep_PFA/topictoolbox')
%corpusname = '20newsBydate';
corpusname = 'NIPSTop2000';
corpusname =[corpusname,num2str(eta0*1000)]
%Theta = ParaLocal{5}.ThetaFreqAver;

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
    S{layer} = WriteTopics( Phit(:,dex) , 0.1 , wl,30);
end

for layer=1:5
    figure;plot(rsort{layer},'o')
end

topic1 = strsplit(S{1}{1});
StopWordDex = [];
for i=1:30 %length(topic1)
    StopWordDex=[StopWordDex,find(strcmp(wl,topic1{i}))];
end
wl(StopWordDex)=[];
Phi{1}(StopWordDex,:)=[];
% 
 topic2 = strsplit(S{1}{2});
StopWordDex = [];
for i=1:20 %length(topic2)
    StopWordDex=[StopWordDex,find(strcmp(wl,topic2{i}))];
end
wl(StopWordDex)=[];
Phi{1}(StopWordDex,:)=[];


topic3 = strsplit(S{1}{3});
StopWordDex = [];
for i=1:10 %length(topic3)
    StopWordDex=[StopWordDex,find(strcmp(wl,topic3{i}))];
end
wl(StopWordDex)=[];
Phi{1}(StopWordDex,:)=[];



Phit=1;
Phisort=Phi;
rsort = cell(5,1);
for layer=1:5
    Phit=Phit*Phi{layer};
    %[~,dex]=sort(sum(Theta{layer},2),'descend');
    % [rsort{layer},dex]=sort(sum(Theta{layer},2),'descend');
    %[rsort{layer},dex]=sort(ParaGlobal{layer}.r_k,'descend');
    %     if layer==1
    %         [rsort{layer},dex]=sort(Phi{2}*Phi{3}*Phi{4}*Phi{5}*ParaGlobal{5}.r_k,'descend');
    %     else
    %     dex = ParaGlobal{5}.Popularity{layer};
    %     end
    
    r_k = ParaGlobal{5}.r_k;
    for tt=4:-1:layer
        r_k = Phi{tt+1}*r_k;
    end
    
    [rsort{layer},dex]=sort(r_k,'descend');
    
    %[rsort{layer},dex]=sort(sum(Theta{layer},2),'descend');
    Phisort{layer} = Phisort{layer}(:,dex);
    if layer<5
        Phisort{layer+1} = Phisort{layer+1}(dex,:);
    end
    S{layer} = WriteTopics( Phit(:,dex) , 0.1 , wl,30);
end

%for nodeFrom=[2,4,6,13,15,22,24,26]% tech size(Phisort{startlayer},2)
%for nodeFrom=[2,4,6,15,22,24]% tech size(Phisort{startlayer},2)
%for nodeFrom=[2,4,13]; %,26]% softwares, windows, graphics
%for nodeFrom=[6,15,22,24]% hardwares, mac dirve
%for nodeFrom=13 % internet
%for nodeFrom=[8,14,23,9] %hokey baseball, med
%for nodeFrom=[10,12,16,19]% car and bike
%for nodeFrom=[1,11,21]% religion
%for nodeFrom=[27,30] %Science
%for nodeFrom=[18] %,28]% Israel
%for nodeFrom=[5,7,18,25]; %,25] %,13]% independent
%for nodeFrom = [3,17] % goverment
%      for nodeFrom = 3 %gun
%for nodeFrom=10 %for sale
%for nodeFrom=1: size(Phisort{startlayer},2)


startlayer=5;
endlayer=1;
Colors = {'black', 'blue','orange','green','red'};
for nodeFrom=1:min(40,size(Phisort{startlayer},2))
    filename=[corpusname,'_K',num2str(K),'_start',num2str(startlayer),'_node',num2str(nodeFrom),'_end',num2str(endlayer),'updown.dot'];
    fid = fopen(filename,'W');
    fprintf(fid,'digraph G {\n');
    fprintf(fid, 'size="18,18!";ranksep=4; ratio = auto; layers="1:2:3:4:5";\n');
    fprintf(fid,  'graph [fontname = "helvetica"];')
        fprintf(fid,  'node [fontname = "helvetica"];')
        fprintf(fid, 'edge [fontname = "helvetica"];')
        %fprintf(fid, 'size="15,15!";ranksep=4; ratio = auto; layers="1:2:3:4:5";\n');
        fprintf(fid, 'ranksep=2.8;\n');
    
    %fprintf(fid,'graph [rankdir = "LR"];')
    edgelist={};
    nodelist={};
    MaxSize=[20,10,10,10,10]*40;
    
    NetWorkOnly = false;
    minWeight=[10,20,5,5,5,5,5];
    %,10,15] %size(Phisort{startlayer},2)  %size(Phisort{end},2) %1:min(size(Phisort,2),5)
    [edgelist,nodelist]=dotwrite_updown_weight(fid,nodeFrom,startlayer,Phisort,S,edgelist,nodelist,endlayer,MaxSize,rsort,NetWorkOnly,minWeight,Colors);
    
     
    
    for i=1:5
        temp=[];
        for j=1:length(nodelist)
            if strcmp(nodelist{j}(2),num2str(i))
                ccc=nodelist{j};
                temp=[temp ' ' ccc];
            end
        end
        fprintf(fid, ['{rank=same;' temp '}']);
    end
    
    fprintf(fid,'}');
    fclose(fid);
    system(['/usr/local/bin/dot -Gepsilon=.0000001 -Tpdf ', filename, ' -o ',filename(1:end-4),'.pdf'])

end

startlayer=5;
endlayer=2;
Colors = {'black', 'blue','orange','green','red'};
for nodeFrom=1:min(40,size(Phisort{startlayer},2))
    filename=[corpusname,'_K',num2str(K),'_start',num2str(startlayer),'_node',num2str(nodeFrom),'_end',num2str(endlayer),'updown.dot'];
    fid = fopen(filename,'W');
    fprintf(fid,'digraph G {\n');
    fprintf(fid, 'size="18,18!";ranksep=4; ratio = auto; layers="1:2:3:4:5";\n');
    fprintf(fid,  'graph [fontname = "helvetica"];')
        fprintf(fid,  'node [fontname = "helvetica"];')
        fprintf(fid, 'edge [fontname = "helvetica"];')
        %fprintf(fid, 'size="15,15!";ranksep=4; ratio = auto; layers="1:2:3:4:5";\n');
        fprintf(fid, 'ranksep=2.8;\n');
    
    %fprintf(fid,'graph [rankdir = "LR"];')
    edgelist={};
    nodelist={};
    MaxSize=[20,10,10,10,10]*40;
    
    NetWorkOnly = false;
    minWeight=[10,10,1,1,1,1,1];
    %,10,15] %size(Phisort{startlayer},2)  %size(Phisort{end},2) %1:min(size(Phisort,2),5)
    [edgelist,nodelist]=dotwrite_updown_weight(fid,nodeFrom,startlayer,Phisort,S,edgelist,nodelist,endlayer,MaxSize,rsort,NetWorkOnly,minWeight,Colors);
    
     
    
    for i=1:5
        temp=[];
        for j=1:length(nodelist)
            if strcmp(nodelist{j}(2),num2str(i))
                ccc=nodelist{j};
                temp=[temp ' ' ccc];
            end
        end
        fprintf(fid, ['{rank=same;' temp '}']);
    end
    
    fprintf(fid,'}');
    fclose(fid);
    system(['/usr/local/bin/dot -Gepsilon=.0000001 -Tpdf ', filename, ' -o ',filename(1:end-4),'.pdf'])

end

startlayer=3;
endlayer=1;
for nodeFrom=1:min(50,size(Phisort{startlayer},2))
    filename=[corpusname,'_K',num2str(K),'_start',num2str(startlayer),'_node',num2str(nodeFrom),'_end',num2str(endlayer),'updown.dot'];
    fid = fopen(filename,'W');
    fprintf(fid,'digraph G {\n');
    fprintf(fid, 'size="18,18!";ranksep=4; ratio = auto; layers="1:2:3:4:5";\n');
    fprintf(fid,  'graph [fontname = "helvetica"];')
        fprintf(fid,  'node [fontname = "helvetica"];')
        fprintf(fid, 'edge [fontname = "helvetica"];')
        %fprintf(fid, 'size="15,15!";ranksep=4; ratio = auto; layers="1:2:3:4:5";\n');
        fprintf(fid, 'ranksep=2.8;\n');
    %fprintf(fid,'graph [rankdir = "LR"];')
    edgelist={};
    nodelist={};
    MaxSize=[20,10,10,10,10]*40;
    
    NetWorkOnly = false;
    minWeight=1;
    %,10,15] %size(Phisort{startlayer},2)  %size(Phisort{end},2) %1:min(size(Phisort,2),5)
    [edgelist,nodelist]=dotwrite_updown(fid,nodeFrom,startlayer,Phisort,S,edgelist,nodelist,endlayer,MaxSize,rsort,NetWorkOnly,minWeight,Colors);
    
    
    for i=1:5
        temp=[];
        for j=1:length(nodelist)
            if strcmp(nodelist{j}(2),num2str(i))
                ccc=nodelist{j};
                temp=[temp ' ' ccc];
            end
        end
        fprintf(fid, ['{rank=same;' temp '}']);
    end
    
    fprintf(fid,'}');
    fclose(fid);
     system(['/usr/local/bin/dot -Gepsilon=.0000001 -Tpdf ', filename, ' -o ',filename(1:end-4),'.pdf'])

end

if 0
    %TopicsAll={'god','independent','sports','mideast','gun','car_sales','hardware','software_graphics','nasa','security','med','isreal'};
    TopicsAll={'car','god','independent','sports','tax','gun','car_sales','hardware','software_graphics','nasa','security','med','isreal'};

    for ii=1:length(TopicsAll)
        MetaTopic = TopicsAll{ii};
        switch MetaTopic
            case 'windows'
                NodeList=1;
            case 'god'
                 NodeList=[2,8,28] %,20,26,28]
            case 'car'
                %NodeList = [1,11,19]; %2
                %NodeList = [3,11];
                NodeList = [11];
            case 'sports'
                %NodeList = [6,13,17,15]
                %NodeList = [8,9,14,23]
                NodeList = [7,13,26]
            case 'tax'
                %NodeList = [19,25]
                NodeList = [19]
            case 'gun'
                %NodeList = [8]
                NodeList = [16,18]
            case 'car_sales'
                %NodeList = [5,9,16,20]
                NodeList = [3,11]
            case 'hardware'
                % NodeList = [10,18,23,24]
                %NodeList = [5,10,12,17,21]
                NodeList = [10,17]
            case 'software_graphics'
                %NodeList = [4,12,14,26]
                NodeList =[4];
            case 'nasa'
                NodeList = 6;
            case 'security'
                NodeList = 9;
            case 'med'
                NodeList = 14;
            case 'isreal'
                NodeList = [15 23];
        end
        startlayer=5;
        endlayer=1;
        fid = fopen(['20newsupdown.dot'],'W');
        fprintf(fid,'digraph G {\n');
        fprintf(fid,  'graph [fontname = "helvetica"];')
        fprintf(fid,  'node [fontname = "helvetica"];')
        fprintf(fid, 'edge [fontname = "helvetica"];')
        %fprintf(fid, 'size="15,15!";ranksep=4; ratio = auto; layers="1:2:3:4:5";\n');
        fprintf(fid, 'ranksep=2.8;\n');
        
        %fprintf(fid, 'ranksep=3.5;\n');
        %fprintf(fid,'graph [rankdir = "LR"];')
        
        edgelist={};
        nodelist={};
        MaxSize=[10,10,10,10,10]*500;
        
        NetWorkOnly = false;
        minWeight=3;
        for nodeFrom=NodeList
            %,10,15] %size(Phisort{startlayer},2)  %size(Phisort{end},2) %1:min(size(Phisort,2),5)
            [edgelist,nodelist]=dotwrite_updown(fid,nodeFrom,startlayer,Phisort,S,edgelist,nodelist,endlayer,MaxSize,rsort,NetWorkOnly,minWeight,Colors);
            
     
            % [edgelist,nodelist]=dotwrite_updown_LR(fid,nodeFrom,startlayer,Phisort,S,edgelist,nodelist,endlayer,MaxSize,rsort,NetWorkOnly);
            
        end
        for i=1:5
            temp=[];
            for j=1:length(nodelist)
                if strcmp(nodelist{j}(2),num2str(i))
                    ccc=nodelist{j};
                    temp=[temp ' ' ccc];
                end
            end
            fprintf(fid, ['{rank=same;' temp '}']);
        end
        
        fprintf(fid,'}');
        fclose(fid);
        system(['/usr/local/bin/dot -Gepsilon=.0000001 -Tpdf 20newsupdown.dot -o ',corpusname, '5to1_',MetaTopic ,'.pdf'])
    end
end



if 0
    %TopicsAll={'god','independent','sports','mideast','gun','car_sales','hardware','software_graphics','nasa','security','med','isreal'};
    TopicsAll={'god','god1','god2','god3','independent','sports','tax','gun','car_sales','hardware','software_graphics','nasa','security','med','isreal','islamic','windows_apple','mideast'};

    for ii=1:length(TopicsAll)
        MetaTopic = TopicsAll{ii};
        switch MetaTopic
            case 'windows'
                NodeList=[5,10];
            case 'god'
                 %NodeList=[3,8,18,23,26,33]
                 NodeList=[3]
            case 'god1'
                NodeList=[8]
            case 'god2'
                 %NodeList=[3,8,18,23,26,33]
                 NodeList=[3,24]
            case 'god3'
                 %NodeList=[3,8,18,23,26,33]
                 NodeList=[8,24]
            case 'car'
                %NodeList = [1,11,19]; %2
                NodeList = [4,9];
            case 'sports'
                %NodeList = [6,13,17,15]
                %NodeList = [8,9,14,23]
                NodeList = [6,7]
            case 'tax'
                %NodeList = [19,25]
                NodeList = [25]
            case 'gun'
                %NodeList = [8]
                NodeList = 13
            case 'car_sales'
                %NodeList = [5,9,16,20]
                NodeList = [22]
            case 'hardware'
                % NodeList = [10,18,23,24]
                NodeList = [2,15,16]
            case 'software_graphics'
                %NodeList = [4,12,14,26]
                NodeList =[1,5];
            case 'nasa'
                NodeList = [17,20,28];
            case 'security'
                NodeList = [12,31];
            case 'med'
                NodeList = 14;
            case 'isreal'
                NodeList = [18,30];
            case 'islamic'
                NodeList = 33;
            case 'mideast'
                NodeList=[18,30,33];
            case 'windows_apple'
                NodeList = [2,5];
        end
        startlayer=3;
        endlayer=1;
        fid = fopen(['20newsupdown.dot'],'W');
        fprintf(fid,'digraph G {\n');
        fprintf(fid,  'graph [fontname = "helvetica"];')
        fprintf(fid,  'node [fontname = "helvetica"];')
        fprintf(fid, 'edge [fontname = "helvetica"];')
        %fprintf(fid, 'size="15,15!";ranksep=4; ratio = auto; layers="1:2:3:4:5";\n');
        fprintf(fid, 'ranksep=2.8;\n');
        
        %fprintf(fid, 'ranksep=3.5;\n');
        %fprintf(fid,'graph [rankdir = "LR"];')
        
        edgelist={};
        nodelist={};
        MaxSize=[10,10,10,10,10]*500;
        
        NetWorkOnly = false;
        minWeight=1;
        for nodeFrom=NodeList
            %,10,15] %size(Phisort{startlayer},2)  %size(Phisort{end},2) %1:min(size(Phisort,2),5)
            [edgelist,nodelist]=dotwrite_updown(fid,nodeFrom,startlayer,Phisort,S,edgelist,nodelist,endlayer,MaxSize,rsort,NetWorkOnly,minWeight,Colors);
            % [edgelist,nodelist]=dotwrite_updown_LR(fid,nodeFrom,startlayer,Phisort,S,edgelist,nodelist,endlayer,MaxSize,rsort,NetWorkOnly);
            
        end
        for i=1:5
            temp=[];
            for j=1:length(nodelist)
                if strcmp(nodelist{j}(2),num2str(i))
                    ccc=nodelist{j};
                    temp=[temp ' ' ccc];
                end
            end
            fprintf(fid, ['{rank=same;' temp '}']);
        end
        
        fprintf(fid,'}');
        fclose(fid);
        system(['/usr/local/bin/dot -Gepsilon=.0000001 -Tpdf 20newsupdown.dot -o ',corpusname, '3to1_',MetaTopic ,'.pdf'])
    end
end


if 0
    %TopicsAll={'god','independent','sports','mideast','gun','car_sales','hardware','software_graphics','nasa','security','med','isreal'};
    TopicsAll={'god','god1','god2','independent','sports','tax','gun','car_sales','hardware','software_graphics','nasa','security','med','isreal'};

    for ii=1:length(TopicsAll)
        MetaTopic = TopicsAll{ii};
        switch MetaTopic
            case 'windows'
                NodeList=1;
            case 'god'
                 NodeList=[2,7] %,20,26,28]
            case 'god1'
                 NodeList=[2] %,20,26,28]
            case 'god2'
                 NodeList=[7] %,20,26,28]
            case 'car'
                %NodeList = [1,11,19]; %2
                NodeList = [3,11];
            case 'sports'
                %NodeList = [6,13,17,15]
                %NodeList = [8,9,14,23]
                NodeList = [8,14,25]
            case 'tax'
                %NodeList = [19,25]
                NodeList = [19]
            case 'gun'
                %NodeList = [8]
                NodeList = [18]
            case 'car_sales'
                %NodeList = [5,9,16,20]
                NodeList = [24]
            case 'hardware'
                % NodeList = [10,18,23,24]
                %NodeList = [5,10,12,17,21]
                NodeList = [10,17]
            case 'software_graphics'
                %NodeList = [4,12,14,26]
                NodeList =[4];
            case 'nasa'
                NodeList = 6;
            case 'security'
                NodeList = 9;
            case 'med'
                NodeList = 13;
            case 'isreal'
                NodeList = [16 23];
        end
        startlayer=5;
        endlayer=1;
        fid = fopen(['20newsupdown.dot'],'W');
        fprintf(fid,'digraph G {\n');
        fprintf(fid,  'graph [fontname = "helvetica"];')
        fprintf(fid,  'node [fontname = "helvetica"];')
        fprintf(fid, 'edge [fontname = "helvetica"];')
        %fprintf(fid, 'size="15,15!";ranksep=4; ratio = auto; layers="1:2:3:4:5";\n');
        %fprintf(fid, 'ranksep=2.8;\n');
        
        fprintf(fid, 'ranksep=3.5;\n');
        fprintf(fid,'graph [rankdir = "LR"];')
        
        edgelist={};
        nodelist={};
        MaxSize=[10,10,10,10,10]*500;
        
        NetWorkOnly = false;
        minWeight=0;
        for nodeFrom=NodeList
            %,10,15] %size(Phisort{startlayer},2)  %size(Phisort{end},2) %1:min(size(Phisort,2),5)
           % [edgelist,nodelist]=dotwrite_updown(fid,nodeFrom,startlayer,Phisort,S,edgelist,nodelist,endlayer,MaxSize,rsort,NetWorkOnly,minWeight,Colors);
            
     
            [edgelist,nodelist]=dotwrite_updown_LR(fid,nodeFrom,startlayer,Phisort,S,edgelist,nodelist,endlayer,MaxSize,rsort,NetWorkOnly,minWeight,Colors);
            
        end
        for i=1:5
            temp=[];
            for j=1:length(nodelist)
                if strcmp(nodelist{j}(2),num2str(i))
                    ccc=nodelist{j};
                    temp=[temp ' ' ccc];
                end
            end
            fprintf(fid, ['{rank=same;' temp '}']);
        end
        
        fprintf(fid,'}');
        fclose(fid);
        system(['/usr/local/bin/dot -Gepsilon=.0000001 -Tpdf 20newsupdown.dot -o ',corpusname, '5to1_',MetaTopic ,'LR.pdf'])
    end
    
    
    %TopicsAll={'god','independent','sports','mideast','gun','car_sales','hardware','software_graphics','nasa','security','med','isreal'};
    TopicsAll={'god','independent','sports','tax','gun','car_sales','hardware','software_graphics','nasa','security','med','isreal'};

    for ii=1:length(TopicsAll)
        MetaTopic = TopicsAll{ii};
        switch MetaTopic
            case 'windows'
                NodeList=1;
            case 'god'
                 NodeList=[2,7] %,20,26,28]
            case 'car'
                %NodeList = [1,11,19]; %2
                NodeList = [3,11];
            case 'sports'
                %NodeList = [6,13,17,15]
                %NodeList = [8,9,14,23]
                NodeList = [8,14,25]
            case 'tax'
                %NodeList = [19,25]
                NodeList = [19]
            case 'gun'
                %NodeList = [8]
                NodeList = [18]
            case 'car_sales'
                %NodeList = [5,9,16,20]
                NodeList = [24]
            case 'hardware'
                % NodeList = [10,18,23,24]
                %NodeList = [5,10,12,17,21]
                NodeList = [10,17]
            case 'software_graphics'
                %NodeList = [4,12,14,26]
                NodeList =[4];
            case 'nasa'
                NodeList = 6;
            case 'security'
                NodeList = 9;
            case 'med'
                NodeList = 13;
            case 'isreal'
                NodeList = [16 23];
        end
        startlayer=5;
        endlayer=2;
        fid = fopen(['20newsupdown.dot'],'W');
        fprintf(fid,'digraph G {\n');
        fprintf(fid,  'graph [fontname = "helvetica"];')
        fprintf(fid,  'node [fontname = "helvetica"];')
        fprintf(fid, 'edge [fontname = "helvetica"];')
        %fprintf(fid, 'size="15,15!";ranksep=4; ratio = auto; layers="1:2:3:4:5";\n');
        %fprintf(fid, 'ranksep=2.8;\n');
        
        fprintf(fid, 'ranksep=3.5;\n');
        fprintf(fid,'graph [rankdir = "LR"];')
        
        edgelist={};
        nodelist={};
        MaxSize=[10,10,10,10,10]*500;
        
        NetWorkOnly = false;
        minWeight=0.01;
        for nodeFrom=NodeList
            %,10,15] %size(Phisort{startlayer},2)  %size(Phisort{end},2) %1:min(size(Phisort,2),5)
           % [edgelist,nodelist]=dotwrite_updown(fid,nodeFrom,startlayer,Phisort,S,edgelist,nodelist,endlayer,MaxSize,rsort,NetWorkOnly,minWeight,Colors);
            
     
            [edgelist,nodelist]=dotwrite_updown_LR(fid,nodeFrom,startlayer,Phisort,S,edgelist,nodelist,endlayer,MaxSize,rsort,NetWorkOnly,minWeight,Colors);
            
        end
        for i=1:5
            temp=[];
            for j=1:length(nodelist)
                if strcmp(nodelist{j}(2),num2str(i))
                    ccc=nodelist{j};
                    temp=[temp ' ' ccc];
                end
            end
            fprintf(fid, ['{rank=same;' temp '}']);
        end
        
        fprintf(fid,'}');
        fclose(fid);
        system(['/usr/local/bin/dot -Gepsilon=.0000001 -Tpdf 20newsupdown.dot -o ',corpusname, '5to2_',MetaTopic ,'LR.pdf'])
    end
    
    
    %TopicsAll={'god','independent','sports','mideast','gun','car_sales','hardware','software_graphics','nasa','security','med','isreal'};
    TopicsAll={'god','independent','sports','tax','gun','car_sales','hardware','software_graphics','nasa','security','med','isreal'};

    for ii=1:length(TopicsAll)
        MetaTopic = TopicsAll{ii};
        switch MetaTopic
            case 'windows'
                NodeList=1;
            case 'god'
                 NodeList=[2,7] %,20,26,28]
            case 'car'
                %NodeList = [1,11,19]; %2
                NodeList = [3,11];
            case 'sports'
                %NodeList = [6,13,17,15]
                %NodeList = [8,9,14,23]
                NodeList = [8,14,25]
            case 'tax'
                %NodeList = [19,25]
                NodeList = [19]
            case 'gun'
                %NodeList = [8]
                NodeList = [18]
            case 'car_sales'
                %NodeList = [5,9,16,20]
                NodeList = [24]
            case 'hardware'
                % NodeList = [10,18,23,24]
                %NodeList = [5,10,12,17,21]
                NodeList = [10,17]
            case 'software_graphics'
                %NodeList = [4,12,14,26]
                NodeList =[4];
            case 'nasa'
                NodeList = 6;
            case 'security'
                NodeList = 9;
            case 'med'
                NodeList = 13;
            case 'isreal'
                NodeList = [16 23];
        end
        startlayer=5;
        endlayer=3;
        fid = fopen(['20newsupdown.dot'],'W');
        fprintf(fid,'digraph G {\n');
        fprintf(fid,  'graph [fontname = "helvetica"];')
        fprintf(fid,  'node [fontname = "helvetica"];')
        fprintf(fid, 'edge [fontname = "helvetica"];')
        %fprintf(fid, 'size="15,15!";ranksep=4; ratio = auto; layers="1:2:3:4:5";\n');
        fprintf(fid, 'ranksep=2.8;\n');
        
        %fprintf(fid, 'ranksep=3.5;\n');
        %fprintf(fid,'graph [rankdir = "LR"];')
        
        edgelist={};
        nodelist={};
        MaxSize=[10,10,10,10,10]*500;
        
        NetWorkOnly = false;
        minWeight=0.01;
        for nodeFrom=NodeList
            %,10,15] %size(Phisort{startlayer},2)  %size(Phisort{end},2) %1:min(size(Phisort,2),5)
           % [edgelist,nodelist]=dotwrite_updown(fid,nodeFrom,startlayer,Phisort,S,edgelist,nodelist,endlayer,MaxSize,rsort,NetWorkOnly,minWeight,Colors);
            
     
            [edgelist,nodelist]=dotwrite_updown(fid,nodeFrom,startlayer,Phisort,S,edgelist,nodelist,endlayer,MaxSize,rsort,NetWorkOnly,minWeight,Colors);
            
        end
        for i=1:5
            temp=[];
            for j=1:length(nodelist)
                if strcmp(nodelist{j}(2),num2str(i))
                    ccc=nodelist{j};
                    temp=[temp ' ' ccc];
                end
            end
            fprintf(fid, ['{rank=same;' temp '}']);
        end
        
        fprintf(fid,'}');
        fclose(fid);
        system(['/usr/local/bin/dot -Gepsilon=.0000001 -Tpdf 20newsupdown.dot -o ',corpusname, '5to3_',MetaTopic ,'.pdf'])
    end
    
    
    %TopicsAll={'god','independent','sports','mideast','gun','car_sales','hardware','software_graphics','nasa','security','med','isreal'};
    TopicsAll={'god','independent','sports','tax','gun','car_sales','hardware','software_graphics','nasa','security','med','isreal'};

    for ii=1:length(TopicsAll)
        MetaTopic = TopicsAll{ii};
        switch MetaTopic
            case 'windows'
                NodeList=1;
            case 'god'
                 NodeList=[2,7] %,20,26,28]
            case 'car'
                %NodeList = [1,11,19]; %2
                NodeList = [3,11];
            case 'sports'
                %NodeList = [6,13,17,15]
                %NodeList = [8,9,14,23]
                NodeList = [8,14,25]
            case 'tax'
                %NodeList = [19,25]
                NodeList = [19]
            case 'gun'
                %NodeList = [8]
                NodeList = [18]
            case 'car_sales'
                %NodeList = [5,9,16,20]
                NodeList = [24]
            case 'hardware'
                % NodeList = [10,18,23,24]
                %NodeList = [5,10,12,17,21]
                NodeList = [10,17]
            case 'software_graphics'
                %NodeList = [4,12,14,26]
                NodeList =[4];
            case 'nasa'
                NodeList = 6;
            case 'security'
                NodeList = 9;
            case 'med'
                NodeList = 13;
            case 'isreal'
                NodeList = [16 23];
        end
        startlayer=5;
        endlayer=2;
        fid = fopen(['20newsupdown.dot'],'W');
        fprintf(fid,'digraph G {\n');
        fprintf(fid,  'graph [fontname = "helvetica"];')
        fprintf(fid,  'node [fontname = "helvetica"];')
        fprintf(fid, 'edge [fontname = "helvetica"];')
        %fprintf(fid, 'size="15,15!";ranksep=4; ratio = auto; layers="1:2:3:4:5";\n');
        fprintf(fid, 'ranksep=2.8;\n');
        
        %fprintf(fid, 'ranksep=3.5;\n');
        %fprintf(fid,'graph [rankdir = "LR"];')
        
        edgelist={};
        nodelist={};
        MaxSize=[10,10,10,10,10]*500;
        
        NetWorkOnly = false;
        minWeight=0.01;
        for nodeFrom=NodeList
            %,10,15] %size(Phisort{startlayer},2)  %size(Phisort{end},2) %1:min(size(Phisort,2),5)
           % [edgelist,nodelist]=dotwrite_updown(fid,nodeFrom,startlayer,Phisort,S,edgelist,nodelist,endlayer,MaxSize,rsort,NetWorkOnly,minWeight,Colors);
            
     
            [edgelist,nodelist]=dotwrite_updown(fid,nodeFrom,startlayer,Phisort,S,edgelist,nodelist,endlayer,MaxSize,rsort,NetWorkOnly,minWeight,Colors);
            
        end
        for i=1:5
            temp=[];
            for j=1:length(nodelist)
                if strcmp(nodelist{j}(2),num2str(i))
                    ccc=nodelist{j};
                    temp=[temp ' ' ccc];
                end
            end
            fprintf(fid, ['{rank=same;' temp '}']);
        end
        
        fprintf(fid,'}');
        fclose(fid);
        system(['/usr/local/bin/dot -Gepsilon=.0000001 -Tpdf 20newsupdown.dot -o ',corpusname, '5to2_',MetaTopic ,'.pdf'])
    end
end

if 0
    %% draw disconnected networks
    endlayer=1;
    for startlayer=5:-1:2
        dex=find(sum(Phisort{startlayer}>0.01,2)==0);
        for nodeFrom=dex(:)'
            fid = fopen([corpusname,'_',num2str(KKK),'_rootlayer',num2str(startlayer-1),'_',num2str(nodeFrom),'updown.dot'],'W');
            fprintf(fid,'digraph G {\n');
            fprintf(fid,  'graph [fontname = "helvetica"];')
            fprintf(fid,  'node [fontname = "helvetica"];')
            fprintf(fid, 'edge [fontname = "helvetica"];')
            %fprintf(fid, 'size="15,15!";ranksep=4; ratio = auto; layers="1:2:3:4:5";\n');
            fprintf(fid, 'ranksep=2.8;\n');
            
            %fprintf(fid,'graph [rankdir = "LR"];')
            edgelist={};
            nodelist={};
            MaxSize=[10,10,10,10,10]*10;
            
            NetWorkOnly = false;
            %,10,15] %size(Phisort{startlayer},2)  %size(Phisort{end},2) %1:min(size(Phisort,2),5)
            if startlayer-1==1
                [edgelist,nodelist]=dotwrite_updown_bottomelayer(fid,nodeFrom,startlayer-1,Phisort,S,edgelist,nodelist,endlayer,MaxSize,rsort,NetWorkOnly);
                
            else
                [edgelist,nodelist]=dotwrite_updown(fid,nodeFrom,startlayer-1,Phisort,S,edgelist,nodelist,endlayer,MaxSize,rsort,NetWorkOnly);
            end
            
            for i=1:5
                temp=[];
                for j=1:length(nodelist)
                    if strcmp(nodelist{j}(2),num2str(i))
                        ccc=nodelist{j};
                        temp=[temp ' ' ccc];
                    end
                end
                fprintf(fid, ['{rank=same;' temp '}']);
            end
            
            fprintf(fid,'}');
            fclose(fid);
        end
    end
    
end
%%
if 0
    
    
    fid = fopen([corpusname,'downup.dot'],'W');
    fprintf(fid,'digraph G {\n');
    fprintf(fid, 'ranksep=6; size ="10,10"; ratio = "compress"; pack=true; layers="1:2:3:4:5";\n');
    %fprintf(fid,'graph [rankdir = "LR"];')
    edgelist={};
    nodelist={};
    endlayer=3;
    startlayer=length(Phi);
    MaxSize=[10,10,10,10,10]*3;
    startlayer=5;
    endlayer=1;
    NetWorkOnly = false;
    for nodeFrom=1:5 %size(Phisort{endlayer},2)  %size(Phisort{end},2) %1:min(size(Phisort,2),5)
        [edgelist,nodelist]=dotwrite_downup(fid,nodeFrom,startlayer,Phisort,S,edgelist,nodelist,endlayer,MaxSize,rsort,NetWorkOnly);
    end
    
    for i=1:5
        temp=[];
        for j=1:length(nodelist)
            if strcmp(nodelist{j}(2),num2str(i))
                ccc=nodelist{j};
                temp=[temp ' ' ccc];
            end
        end
        fprintf(fid, ['{rank=same;' temp '}']);
    end
    
    fprintf(fid,'}');
    fclose(fid);
    
    
    
    fid = fopen([corpusname,num2str(KKK),'downup_network.dot'],'W');
    fprintf(fid,'digraph G {\n');
    fprintf(fid, 'ranksep=4; ratio = "auto"; layers="1:2:3:4:5";\n');
    %fprintf(fid,'graph [rankdir = "LR"];')
    edgelist={};
    nodelist={};
    endlayer=3;
    startlayer=length(Phi);
    MaxSize=[10,10,10,10,10]*400;
    startlayer=5;
    endlayer=1;
    NetWorkOnly = true;
    for nodeFrom=1:size(Phisort{endlayer},2)  %size(Phisort{end},2) %1:min(size(Phisort,2),5)
        [edgelist,nodelist]=dotwrite_downup(fid,nodeFrom,startlayer,Phisort,S,edgelist,nodelist,endlayer,MaxSize,rsort,NetWorkOnly);
    end
    
    for i=1:5
        temp=[];
        for j=1:length(nodelist)
            if strcmp(nodelist{j}(2),num2str(i))
                ccc=nodelist{j};
                temp=[temp ' ' ccc];
            end
        end
        fprintf(fid, ['{rank=same;' temp '}']);
    end
    
    fprintf(fid,'}');
    fclose(fid);
    
    
    fid = fopen([corpusname,num2str(KKK),'updown_network.dot'],'W');
    fprintf(fid,'digraph G {\n');
    fprintf(fid, 'ranksep=4; ratio = auto; layers="1:2:3:4:5";\n');
    %fprintf(fid,'graph [rankdir = "LR"];')
    edgelist={};
    nodelist={};
    MaxSize=[10,10,10,10,10]*400;
    startlayer=5;
    endlayer=1;
    NetWorkOnly = true;
    for nodeFrom=1:size(Phisort{startlayer},2)  %size(Phisort{end},2) %1:min(size(Phisort,2),5)
        [edgelist,nodelist]=dotwrite_updown(fid,nodeFrom,startlayer,Phisort,S,edgelist,nodelist,endlayer,MaxSize,rsort,NetWorkOnly);
    end
    
    for i=1:5
        temp=[];
        for j=1:length(nodelist)
            if strcmp(nodelist{j}(2),num2str(i))
                ccc=nodelist{j};
                temp=[temp ' ' ccc];
            end
        end
        fprintf(fid, ['{rank=same;' temp '}']);
    end
    
    fprintf(fid,'}');
    fclose(fid);
    
end
