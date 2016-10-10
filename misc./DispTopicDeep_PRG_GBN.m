path = '/Users/zhoum/Box Sync/GBN_results/results_PRG_GBN_results2/'

for K0=[50,100,200,400]
    load([path,'ZGlobalMnist_CJ_K0_',num2str(K0),'_T_5_eta050_PreTrain0_AllForPhi0_TrimStrategy5_Trial1_Data'])
    
    filepath = 'digits_images';
    
    
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
        figure;I = DispDictionary(-Phit(:,dex));
        imwrite(I,[filepath '/K1max_' num2str(K0) '_SortedPhi' num2str(layer) '_5.png'])
    end
end

K0=400;
load([path,'ZGlobalMnist_CJ_K0_',num2str(K0),'_T_5_eta050_PreTrain0_AllForPhi0_TrimStrategy5_Trial1_Data'])
T = 5;
theta = randg(ParaGlobal{T}.r_k*ones(1,800))./ParaGlobal{T}.cjmedian{T+1};
tave=3;
for t=T:-1:tave
    
     theta = gamrnd(ParaGlobal{T}.Phi{t}*theta,1./ParaGlobal{T}.cjmedian{t});
    
end

for t=tave-1:-1:1
        theta = (ParaGlobal{T}.Phi{t}*theta)./ParaGlobal{T}.cjmedian{t};
end

theta(:,sum(theta,1)<5)=[];
%Phit  = ParaGlobal{5}.Phi{1}*theta;
%theta = bsxfun(@rdivide,theta,sum(theta,1));
%Phit = 1-exp(-ParaGlobal{4}.Phi{1}*ParaGlobal{4}.Phi{2}*ParaGlobal{4}.Phi{3}*theta);
%Phit = 1-exp(-ParaGlobal{T}.Phi{1}*theta);
%Phit = rand(size(Phit))<Phit;
Phit=theta;
figure
DispDictionary(Phit);

if 0
    for  K0=[50,100,200,400]
        load([path,'ZGlobalMnist_CJ_K0_',num2str(K0),'_T_5_eta050_PreTrain0_AllForPhi0_TrimStrategy5_Trial1_Data'])
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
            Phisort{layer} = Phisort{layer}(:,dex);
            if layer<5
                Phisort{layer+1} = Phisort{layer+1}(dex,:);
            end
            %[~,dex]=sort(sum(Theta{layer},2),'descend');
            % dex = Popularity{ToBeAnalysed}{layer};
            if 0
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
            end
            %S{layer} = WriteTopics( Phit(:,dex) , 0.1 , wl,20 );
        end
        
        dataname = 'MNIST_2015';
        S=[];
        for startlayer=2:5
            endlayer=startlayer-1;
            switch endlayer
                case 1
                    NodeList = 1:10:120;
                    %NodeList = [1,7,27,30,41,58,68,69,74,75]
                    %NodeList = [1,7,27,30,41,68,74,75,15,8,50]
                case 2
                    %NodeList = [7,8,12,10,16,34,35,15,48,50,208,56,62,74,1,5,11,17,30]; %1:4:50,
                    %NodeList = [7,8,12,10,16,34,35,1:4:80];
                    NodeList = [1:5:110]; %1:4:50,
                case 3
                    %NodeList = [10,1:4:100];
                    NodeList = [1:5:150]; 
                case 4
                    %NodeList = [2,3,8,11,12,1:4:100,30];
                    NodeList = [1:5:150]; 
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
            system(['/usr/local/bin/neato -Gepsilon=.0000001 -Tpdf ' , dataname, 'updown.dot -o ',dataname,'_K0',num2str(K0),'_',num2str(startlayer) ,'.pdf'])
            system(['/usr/local/bin/neato -Gepsilon=.0000001 -Tpng ' , dataname, 'updown.dot -o ',dataname,'_K0',num2str(K0),'_',num2str(startlayer) ,'.png'])
            
            %system(['/usr/local/bin/neato -Gepsilon=.0000001 -Tpdf MNISTupdown.dot -o ',dataname,num2str(startlayer) ,'a.pdf'])
        end
    end
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
