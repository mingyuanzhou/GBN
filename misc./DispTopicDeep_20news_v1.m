
%load nips12raw_str602
KKK=400;%25;50;100

corpusname = '20news';
addpath('/Users/zhoum/Box Sync/PGBN_results_lastLayer/PGNB_V5_Experiment_1_2_20150607/Eta0_1')
load(['MultiClass_collapsed_trim_20news_K0_',num2str(KKK),'_T_5_Trial_1_Burnin_1000_Collection_500.mat'])
Phi=PhiAve{5};
addpath('/Users/zhoum/Dropbox/Zhou_Code/Deep_PFA/topictoolbox')
S=cell(5,1);

wl=WO;
Theta{1}=ThetaAve{1}(:,~dexTest);
for t=1:5
    Theta{t}=ThetaAve{t}(:,~dexTest);
end

Phit=1;
Phisort=Phi;
rsort = cell(5,1);
for layer=1:5
    Phit=Phit*Phi{layer};
    %[~,dex]=sort(sum(Theta{layer},2),'descend');
    [rsort{layer},dex]=sort(sum(Theta{layer},2),'descend');
    Phisort{layer} = Phisort{layer}(:,dex);
    if layer<5
        Phisort{layer+1} = Phisort{layer+1}(dex,:);
    end
    S{layer} = WriteTopics( Phit(:,dex) , 0.1 , wl,30);
end
topic1 = strsplit(S{1}{1});
StopWordDex = [];
for i=1:length(topic1)
    StopWordDex=[StopWordDex,find(strcmp(wl,topic1{i}))];
end
wl(StopWordDex)=[];
Phi{1}(StopWordDex,:)=[];


Phit=1;
Phisort=Phi;
rsort = cell(5,1);
for layer=1:5
    Phit=Phit*Phi{layer};
    %[~,dex]=sort(sum(Theta{layer},2),'descend');
    [rsort{layer},dex]=sort(sum(Theta{layer},2),'descend');
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
TopicsAll={'god','independent','sports','mideast','gun','car_sales','hardware','software_graphics','nasa','security','med','isreal'};
for ii=1:length(TopicsAll)
    MetaTopic = TopicsAll{ii};
    switch MetaTopic
        case 'god'
            %NodeList = [3,7,21]
            NodeList=[1,11,21]
        case 'independent'
            %NodeList = [1,11,19]; %2
            NodeList = [5,7,18];
        case 'sports'
            %NodeList = [6,13,17,15]
            %NodeList = [8,9,14,23]
            NodeList = [8,14]
        case 'mideast'
            %NodeList = [19,25]
            NodeList = [18,28]
        case 'gun'
            %NodeList = [8]
            NodeList = [3]
        case 'car_sales'
            %NodeList = [5,9,16,20]
            NodeList = [10,12,16,19]
        case 'hardware'
           % NodeList = [10,18,23,24]
             NodeList = [6,15,22,24]
        case 'software_graphics'
            %NodeList = [4,12,14,26]
            NodeList =[2,4,13];
        case 'nasa'
            NodeList = 5;
        case 'security'
            NodeList = 7;
        case 'med'
            NodeList = 9;
        case 'isreal'
            NodeList = 18;
    end
    startlayer=5;
    endlayer=1;
    fid = fopen([corpusname,'updown.dot'],'W');
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
    
    for nodeFrom=NodeList
        %,10,15] %size(Phisort{startlayer},2)  %size(Phisort{end},2) %1:min(size(Phisort,2),5)
        [edgelist,nodelist]=dotwrite_updown(fid,nodeFrom,startlayer,Phisort,S,edgelist,nodelist,endlayer,MaxSize,rsort,NetWorkOnly);
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
    system(['/usr/local/bin/dot -Gepsilon=.0000001 -Tpdf 20newsupdown.dot -o ',corpusname, 'all_',MetaTopic ,'.pdf'])
end

startlayer=5;
endlayer=1;
for nodeFrom=11:20 %size(Phisort{startlayer},2)
    fid = fopen([corpusname,'_',num2str(KKK),'_',num2str(nodeFrom),'updown.dot'],'W');
    fprintf(fid,'digraph G {\n');
    fprintf(fid, 'size="18,18!";ranksep=4; ratio = auto; layers="1:2:3:4:5";\n');
    
    %fprintf(fid,'graph [rankdir = "LR"];')
    edgelist={};
    nodelist={};
    MaxSize=[10,10,10,10,10]*50;
    
    NetWorkOnly = false;
    %,10,15] %size(Phisort{startlayer},2)  %size(Phisort{end},2) %1:min(size(Phisort,2),5)
    [edgelist,nodelist]=dotwrite_updown(fid,nodeFrom,startlayer,Phisort,S,edgelist,nodelist,endlayer,MaxSize,rsort,NetWorkOnly);
    
    
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
