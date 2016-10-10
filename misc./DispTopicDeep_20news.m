
%load nips12raw_str602
KKK=400;%25;50;100
%load('/Users/zhoum/Box Sync/PGBN_results_lastLayer/NIPS100_0_30_25.mat')
%load(['/Users/zhoum/Box Sync/PGBN_results_lastLayer/NIPS',num2str(KKK),'_0_30_25.mat'])
%corpusname = 'NIPS';

% 
% topic1 = {'network' 'innput' 'neural' 'networks' 'output' 'function' 'data' 'algorithm' 'set' 'model' 'learning'};
% StopWordDex = [];
% for i=1:length(topic1)
%     StopWordDex=[StopWordDex,find(strcmp(wl,topic1{i}))];
% end
%  wl(StopWordDex)=[];
%  Phi{1}(StopWordDex,:)=[];
if 1
    corpusname = '20news';
addpath('/Users/zhoum/Box Sync/PGBN_results_lastLayer/PGNB_V5_Experiment_1_2_20150607')
IsBinaryClassificaiton =false
read20news

wl=WO;


load(['MultiClass_collapsed_trim_20news_K0_',num2str(KKK),'_T_5_Trial_1_Burnin_1000_Collection_500.mat'])
Phi=PhiAve{5};
Theta=ThetaAve;


%  dex=[find(strcmp(wl,'edu')),find(strcmp(wl,'com')), find(strcmp(wl,'writes')),find(strcmp(wl,'article')),find(strcmp(wl,'don'))];
%  dex=[dex,find(strcmp(wl,'just')),find(strcmp(wl,'like')),find(strcmp(wl,'people')),find(strcmp(wl,'think')),find(strcmp(wl,'does'))];
%  dex=[dex,find(strcmp(wl,'know')),find(strcmp(wl,'time')),find(strcmp(wl,'good')),find(strcmp(wl,'apr'))];
% % 
%  wl(dex)=[];
%  Phi{1}(dex,:)=[];

topic1 = {'don' 'just' 'edu' 'writes' 'like' 'think' 'article' 'people' 'know' 'com' 'time' 'good' 'way' 'make' 'does' 'say' 'want' 've' 'really' 'right' ...
    'did' 'use' 'going' 'apr' 'point' 'things' 'thing' 'll' 'better' 'sure'}
%topic1 = {'don' 'just' 'edu' 'writes' 'like' 'think' 'article' 'people' 'know' 'com' 'does' 'did' 'apr'}
StopWordDex = [];
for i=1:length(topic1)
    StopWordDex=[StopWordDex,find(strcmp(wl,topic1{i}))];
end
 wl(StopWordDex)=[];
 Phi{1}(StopWordDex,:)=[];



for t=1:5
Theta{t}=ThetaAve{t}(:,~dexTest);
end
end

Phit=1;
for t=1:5
    Phit=Phit*Phi{t};
end

[~,dex]=sort(sum(Theta{5},2),'descend');
for t=1:3
[~,dext]=sort(Phit(:,dex(t)),'descend');
wl(dext(1:6))'
end


[~,dex]=sort(sum(Theta{1},2),'descend');
for t=1:3
[~,dext]=sort(Phi{1}(:,dex(t)),'descend');
wl(dext(1:6))'
end


addpath('/Users/zhoum/Dropbox/Zhou_Code/Deep_PFA/topictoolbox')

S=cell(5,1);

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
    S{layer} = WriteTopics( Phit(:,dex) , 0.1 , wl,20 );
end



for i=1:length(S{1})
i
S{1}(i)
end

for i=1:length(S{5})
i
S{5}(i)
end

% Phisort=Phi;
% rsort = cell(5,1);
% for layer=length(Phi):-1:1
%     [rsort{layer},dex]=sort(sum(Theta{layer},2),'descend');
%     Phisort{layer} = Phi{layer}(:,dex);
% end


% 
% fid = fopen([corpusname,'.dot'],'W');
% fprintf(fid,'digraph G {\n ranksep=5; \n ratio = auto;\n');
% %fprintf(fid,'graph [rankdir = "LR"];')
% edgelist={};
% endlayer=3;
% startlayer=length(Phi);
% 
% startlayer=5;
% endlayer=1;
% for nodeFrom=1:5 %5:size(Phisort{startlayer},2)  %size(Phisort{end},2) %1:min(size(Phisort,2),5)
%     edgelist=dotwrite(fid,nodeFrom,startlayer,Phisort,S,edgelist,endlayer);
% end
% fprintf(fid,'}');
% fclose(fid);



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



startlayer=5;
endlayer=1;
fid = fopen([corpusname,'updown.dot'],'W');
fprintf(fid,'digraph G {\n');
fprintf(fid,  'graph [fontname = "helvetica"];')
fprintf(fid,  'node [fontname = "helvetica"];')
fprintf(fid, 'edge [fontname = "helvetica"];')
%fprintf(fid, 'size="15,15!";ranksep=4; ratio = auto; layers="1:2:3:4:5";\n');
fprintf(fid, 'ranksep=2.8;\n');

%fprintf(fid,'graph [rankdir = "LR"];')
edgelist={};
nodelist={};
MaxSize=[10,10,10,10,10]*5;

NetWorkOnly = false;
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
for nodeFrom=[5,7,18,25]; %,25] %,13]% independent
  %for nodeFrom = [3,17] % goverment
%      for nodeFrom = 3 %gun
%for nodeFrom=10 %for sale
%for nodeFrom=1: size(Phisort{startlayer},2)
%for nodeFrom=1:size(Phisort{startlayer},2)
    %,10,15] %size(Phisort{startlayer},2)  %size(Phisort{end},2) %1:min(size(Phisort,2),5)
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

system(['/usr/local/bin/dot -Gepsilon=.0000001 -Tpdf 20newsupdown.dot -o ',corpusname,num2str(startlayer) ,'.pdf'])


startlayer=5;
endlayer=1;
for nodeFrom=1:size(Phisort{startlayer},2)
    fid = fopen([corpusname,'_',num2str(KKK),'_',num2str(nodeFrom),'updown.dot'],'W');
    fprintf(fid,'digraph G {\n');
    fprintf(fid, 'size="18,18!";ranksep=4; ratio = auto; layers="1:2:3:4:5";\n');
    
    %fprintf(fid,'graph [rankdir = "LR"];')
    edgelist={};
    nodelist={};
    MaxSize=[10,10,10,10,10]*5;
    
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



% 
% for nodeFrom=1:min(size(Phisort,2),5)
%     for layer=length(Phi):-1:2
%         dex=find(Phisort{layer}(:,nodeFrom)>1/size(Phisort,1));
%         for nodeTo = dex(:)'
%             fprintf(fid,['"' num2str(layer) '_' num2str(nodeFrom) '"' ' -> ' '"' num2str(layer-1) '_' num2str(nodeTo) '"' '\n']);
%         end
%     end
% end
% 
% 
% 
% fprintf(fid,'}');
% fclose(fid);
% 
% for layer=length(Phi):-1:2
%     for nodeFrom=1:min(size(Phisort{layer},2),10)
%         dex=find(Phisort{layer}(:,nodeFrom)>1/size(Phisort,1));
%         for nodeTo = dex(:)'
%             %fprintf(fid,[num2str(layer) '_' num2str(nodeFrom) ' -> ' num2str(layer-1) '_' num2str(nodeTo) '\n']);
%             fprintf(fid,[num2str(layer*500+nodeFrom) ' -> ' num2str((layer-1)*500 + nodeTo) '\n']);
%         end
%     end
% end
% fprintf(fid,'}');
% fclose(fid);
% 
% 
% 
%  digraph G {
%  "main" -> "parse dfad" -> execute;
% main -> init;
%  main -> cleanup;
%  execute -> make_string;
% execute -> printf
%  init -> make_string;
%  main -> printf;
%  execute -> compare;
%  }