load('GBN_MNIST_Results.mat')   ;
load('Main_PretrainTrim_MNIST_T_5_SuPara2_Trial_1_Popularity.mat')   ;

ToBeAnalysed    =   3  ;   % 1:K0=50       2:K0=100        3:K0=200        4:K0=400        5:K0=600
% SortResults     =   1   ;
% ParaGlobal  =   ParaGlobal_Set{ToBeAnalysed}   ;
% for Tcurrent    =   T
%         for tnow    =   1:Tcurrent
%             if tnow == 1
%                 Phitmp  =   ParaGlobal{Tcurrent}.Phi{tnow} ;
%             else
%                 Phitmp  =   Phitmp * ParaGlobal{Tcurrent}.Phi{tnow} ;
%             end
%             if SortResults
%                 H=figure(Tcurrent*10+tnow);set(H,'color','w'),DispDictionaryImshowNegative(Phitmp(:, Popularity{ToBeAnalysed}{tnow} ));
%             else
%                 H=figure(Tcurrent*10+tnow);set(H,'color','w'),DispDictionaryImshowNegative(Phitmp);
%             end
%         end
% end

ParaGlobal  =   ParaGlobal_Set{ToBeAnalysed} ;
%Accuracy_all  =  Accuracy_all_Set{ToBeAnalysed} ;
Phi = cell(5,1);
for t=1:5
    Phi{t}=ParaGlobal{5}.Phi{t} ;
end
K0=K0Set(ToBeAnalysed);

S=cell(5,1);
Phit=1;
Phisort=Phi;
rsort = cell(5,1);
filepath = 'digits_images';
for layer=1:5
    Phit=Phit*Phi{layer};
    %[~,dex]=sort(sum(Theta{layer},2),'descend');
    dex = Popularity{ToBeAnalysed}{layer};
    rsort{layer} = (1./(1:length(dex)));
    %[rsort{layer},dex]=sort(sum(Theta{layer},2),'descend');
    Phisort{layer} = Phisort{layer}(:,dex);
    if layer<5
        Phisort{layer+1} = Phisort{layer+1}(dex,:);
    end
end
if 0
    Phit=1;
    for layer=1:5
        Phit=Phit*Phi{layer};
        %[~,dex]=sort(sum(Theta{layer},2),'descend');
        dex = Popularity{ToBeAnalysed}{layer};
        for i=1:length(dex)
            I=DispDictionary(Phit(:,dex(i)));
            II = zeros(size(I,1)+2,size(I,2)+2,3);
            switch layer
                case 5
                    II(:,:,1)=1;
                case 4
                    II(:,:,2)=1;
                case 3
                    II(:,:,1)=1;
                    II(:,:,2)=1/2;
                case 2
                    II(:,:,3)=1;
                case 1
                    II(:,:,1)=1/2;
                    II(:,:,2)=1/2;
                    II(:,:,3)=1/2;
            end
            
            II(3:end-2,3:end-2,:)=1-I(2:end-1,2:end-1,:);
            %II(2:end-1,2:end-1,:)=1-I(2:end-1,2:end-1,:);
            %I = I(2:end-1,2:end-1,:);
            imwrite(II,[filepath '/K0' num2str(K0) '_layer' num2str(layer) '_' num2str(i) '.png'])
        end
        %S{layer} = WriteTopics( Phit(:,dex) , 0.1 , wl,20 );
    end
end
for startlayer=2:5
    endlayer=startlayer-1;
    switch endlayer
        case 1
            %NodeList = 1:10:100;
            %NodeList = [1,7,27,30,41,58,68,69,74,75]
            NodeList = [1,7,27,30,41,68,74,75,15]
        case 2
            NodeList = [7,8,12,10,16,34,35,15,1:4:50];
            %NodeList = [7,8,12,10,16,34,35,1:4:80];
        case 3
            NodeList = [10,1:4:100];
        case 4
            NodeList = [2,3,8,11,12,1:4:100];
    end
    NodeList(NodeList>K0)=randi(K0,1,nnz(NodeList>K0));
    fid = fopen([dataname,'updown.dot'],'W');
    fprintf(fid,'digraph G {\n');
    %fprintf(fid,  'graph [fontname = "helvetica"];')
    %fprintf(fid,  'node [fontname = "helvetica"];')
    fprintf(fid,  'node[label=""];')
    fprintf(fid, 'edge [color=gray19 ];'); %[fontname = "helvetica"];')
    %fprintf(fid, 'size="15,15!";ranksep=4; ratio = auto; layers="1:2:3:4:5";\n');
    fprintf(fid, 'ranksep=3;ratio = auto;\n');
    
    %fprintf(fid,'graph [rankdir = "LR"];')
    edgelist={};
    nodelist={};
    MaxSize=[10,10,10,10,10]*100;
    NetWorkOnly = false;
    for nodeFrom=NodeList; %,25] %,13]% independent
        %for nodeFrom = [3,17] % goverment
        %      for nodeFrom = 3 %gun
        %   for nodeFrom=1:size(Phisort{startlayer},2)
        %,10,15] %size(Phisort{startlayer},2)  %size(Phisort{end},2) %1:min(size(Phisort,2),5)
        [edgelist,nodelist]=dotwrite_updown_images(fid,nodeFrom,startlayer,Phisort,S,edgelist,nodelist,endlayer,MaxSize,rsort,NetWorkOnly,K0,filepath);
        
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
    system(['/usr/local/bin/neato -Gepsilon=.0000001 -Tpdf MNISTupdown.dot -o ',dataname,'_K0',num2str(K0),'_',num2str(startlayer) ,'.pdf'])
    system(['/usr/local/bin/neato -Gepsilon=.0000001 -Tpng MNISTupdown.dot -o ',dataname,'_K0',num2str(K0),'_',num2str(startlayer) ,'.png'])
    
    %system(['/usr/local/bin/neato -Gepsilon=.0000001 -Tpdf MNISTupdown.dot -o ',dataname,num2str(startlayer) ,'a.pdf'])
end
if 0
    startlayer=5;
    endlayer=1;
    for nodeFrom=1:size(Phisort{startlayer},2)
        fid = fopen([dataname,'_',num2str(K0),'_',num2str(nodeFrom),'updown.dot'],'W');
        fprintf(fid,'digraph G {\n');
        fprintf(fid,  'node[label=""];')
        fprintf(fid, 'edge [color=orange];'); %[fontname = "helvetica"];')
        %fprintf(fid, 'size="15,15!";ranksep=4; ratio = auto; layers="1:2:3:4:5";\n');
        fprintf(fid, 'ranksep=1.5;\n');
        %fprintf(fid, 'size="18,18!";ranksep=4; ratio = auto; layers="1:2:3:4:5";\n');
        
        %fprintf(fid,'graph [rankdir = "LR"];')
        edgelist={};
        nodelist={};
        MaxSize=[10,10,10,10,10]*300;
        
        NetWorkOnly = false;
        %,10,15] %size(Phisort{startlayer},2)  %size(Phisort{end},2) %1:min(size(Phisort,2),5)
        [edgelist,nodelist]=dotwrite_updown_images(fid,nodeFrom,startlayer,Phisort,S,edgelist,nodelist,endlayer,MaxSize,rsort,NetWorkOnly,K0,filepath);
        
        
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
if 0
    fid = fopen([dataname,num2str(K0),'downup_network.dot'],'W');
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
    
    
    fid = fopen([dataname,num2str(K0),'updown_network.dot'],'W');
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
        [edgelist,nodelist]=dotwrite_updown_images(fid,nodeFrom,startlayer,Phisort,S,edgelist,nodelist,endlayer,MaxSize,rsort,NetWorkOnly,K0,filepath);
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
