%Reproduce Figures 3-5 of Mingyuan Zhou, Yulai Cong, and Bo Chen, ?Augmentable gamma belief networks," Journal of Machine Learning Research, vol. 17, pp. 1-44, Sept. 2016.
%please first run DispTopicDeep_20news_v2.m
for iter=1:6
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
    startlayer=5;
    endlayer=4;
    Phisort1=Phisort;
    rsort1=rsort;
    Phisort1{5}=Phisort1{5}*0;
    
    %for nodeFrom=(iter-1)*1+(1:5:30)
    for nodeFrom=(iter-1)*5+(1:5)
        %,10,15] %size(Phisort{startlayer},2)  %size(Phisort{end},2) %1:min(size(Phisort,2),5)
        [edgelist,nodelist]=dotwrite_updown(fid,nodeFrom,startlayer,Phisort1,S,edgelist,nodelist,endlayer,MaxSize,rsort1,NetWorkOnly,minWeight,Colors);
        
        
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
    system(['/usr/local/bin/dot -Gepsilon=.0000001 -Tpdf 20newsupdown.dot -o ',corpusname,num2str(startlayer),'_',num2str(iter), '.pdf'])
end


for iter=1:6
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
    startlayer=3;
    endlayer=2;
    Phisort1=Phisort;
    rsort1=rsort;
    Phisort1{startlayer}=Phisort1{startlayer}*0;
    
    %for nodeFrom=(iter-1)*1+(1:5:30)
    for nodeFrom=(iter-1)*5+(1:5)
        %,10,15] %size(Phisort{startlayer},2)  %size(Phisort{end},2) %1:min(size(Phisort,2),5)
        [edgelist,nodelist]=dotwrite_updown(fid,nodeFrom,startlayer,Phisort1,S,edgelist,nodelist,endlayer,MaxSize,rsort1,NetWorkOnly,minWeight,Colors);
        
        
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
    system(['/usr/local/bin/dot -Gepsilon=.0000001 -Tpdf 20newsupdown.dot -o ',corpusname,num2str(startlayer),'_',num2str(iter), '.pdf'])
end



for iter=1:9
    fid = fopen('20newsdownup.dot','W');
    fprintf(fid,'digraph G {\n');
    fprintf(fid, 'ranksep=4; ratio = "auto"; layers="1:2:3:4:5";\n');
    %fprintf(fid,'graph [rankdir = "LR"];')
    edgelist={};
    nodelist={};
    endlayer=3;
    startlayer=length(Phi);
    MaxSize=[10,10,10,10,10]*400;
    startlayer=2;
    endlayer=1;
    NetWorkOnly = false;
    Phisort1=Phisort;
    rsort1=rsort;
    Phisort1{2}=Phisort1{2}*0;
    for nodeFrom=(iter-1)*1+(1:10:90) %size(Phisort{end},2) %1:min(size(Phisort,2),5)
        [edgelist,nodelist]=dotwrite_downup(fid,nodeFrom,startlayer,Phisort1,S,edgelist,nodelist,endlayer,MaxSize,rsort1,NetWorkOnly);
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
    system(['/usr/local/bin/dot -Gepsilon=.0000001 -Tpdf 20newsdownup.dot -o ',corpusname,num2str(endlayer),'_',num2str(iter), '.pdf'])

end