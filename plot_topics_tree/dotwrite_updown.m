function [edgelist,nodelist] = dotwrite_updown(fid,nodeFrom,layer,Phisort,S,edgelist,nodelist,endlayer,MaxSize,rsort,NetWorkOnly,minWeight,Colors)
%Write topic tree from top to bottom
%Mingyuan Zhou
%August, 2015
if ~exist('minWeight','var')
    minWeight = 1;
end
if ~exist('Colors','var')
    Colors = cell(10,1);
    for i=1:10
        Colors{i}='black';
    end
end

if layer>endlayer
        %dex=find(Phisort{layer}(:,nodeFrom)>max(1/max(size(Phisort{layer},1),20),0.01));
        %dex=find(Phisort{layer}(:,nodeFrom)>max(2/size(Phisort{layer},1),minWeight));
        dex=find(Phisort{layer}(:,nodeFrom)>minWeight/size(Phisort{layer},1));
%         if layer>2
%             dex=find(Phisort{layer}(:,nodeFrom)>1/max(size(Phisort{layer},1)));
%         else
%             dex=find(Phisort{layer}(:,nodeFrom)>0.01);
%         end
        dex=dex(dex<=MaxSize(layer));
        SS = strsplit(S{layer}{nodeFrom});
        FontSize = [',fontsize =' num2str((50*(rsort{layer}(nodeFrom)/rsort{layer}(1))^(1/10)))];
       if layer>endlayer
        Str=['"' num2str(layer) '_' num2str(nodeFrom) '"' ' [ rank =' num2str(layer) FontSize ', shape=box, style="rounded", color=' Colors{layer}  ', label=' '"' num2str(nodeFrom) ' ' strjoin(SS(1:4)) '\n' strjoin(SS(5:8)) '\n' strjoin(SS(9:12)) '"' ']\n' ];
       else
        Str=['"' num2str(layer) '_' num2str(nodeFrom) '"' ' [ rank =' num2str(layer) FontSize ', shape=box, style="rounded", color=' Colors{layer}  ', label=' '"' num2str(nodeFrom) ' ' strjoin(SS(1)) '\n' strjoin(SS(2:3)) '\n' strjoin(SS(4:5)) '\n' strjoin(SS(6:7)) '"' ']\n' ];
       end
        if NetWorkOnly
            FontSize = [',fontsize =' num2str(max(50*(rsort{layer}(nodeFrom)/rsort{layer}(1))^(1/10),1))];
            Str=['"' num2str(layer) '_' num2str(nodeFrom) '"' ' [ rank =' num2str(layer) FontSize ',width=.01, height=.01, shape=ellipse, label=' '"' num2str(nodeFrom) '"' ']\n' ];

        end
        if nnz(strcmp(nodelist,Str))==0
            nodelist{end+1}=['"' num2str(layer) '_' num2str(nodeFrom) '"'];
            fprintf(fid,Str);
        end

        %fprintf(fid,['"' num2str(layer) '_' num2str(nodeFrom) '"' ' [ label=' '"' S{layer}{nodeFrom} '"' ']\n' ]);
        for nodeTo = dex(:)'
             if NetWorkOnly
            Str = [ '"' num2str(layer) '_' num2str(nodeFrom) '"' ' -> '  '"' num2str(layer-1) '_' num2str(nodeTo) '"' ' [penwidth=' num2str(max((5*(Phisort{layer}(nodeTo,nodeFrom)).^(1/2)),0)) ',dir=forward]' '\n'];
             else
                 Str = [ '"' num2str(layer) '_' num2str(nodeFrom) '"' ' -> '  '"' num2str(layer-1) '_' num2str(nodeTo) '"' ' [penwidth=' num2str(max((15*(Phisort{layer}(nodeTo,nodeFrom)).^(1/2)),0)) ',dir=forward]' '\n'];
             end
            if nnz(strcmp(edgelist,Str))==0
                edgelist{end+1}=Str;
                fprintf(fid,Str);
            end
            
            
            SS = strsplit(S{layer-1}{nodeTo});
            FontSize = [',fontsize =' num2str((50*(rsort{layer-1}(nodeTo)/rsort{layer-1}(1)).^(1/10)))];
            if layer-1>endlayer
            Str = ['"' num2str(layer-1) '_' num2str(nodeTo) '"' ' [rank =' num2str(layer+1)  FontSize ', shape=box, style="rounded", color=' Colors{layer-1}  ',label=' '"' num2str(nodeTo) ' ' strjoin(SS(1:4)) '\n' strjoin(SS(5:8)) '\n' strjoin(SS(9:12)) '"' ']\n' ];
            else
            Str = ['"' num2str(layer-1) '_' num2str(nodeTo) '"' ' [rank =' num2str(layer+1)  FontSize ', shape=box, style="rounded", color=' Colors{layer-1}  ',label=' '"' num2str(nodeTo) ' ' strjoin(SS(1)) '\n' strjoin(SS(2)) '\n' strjoin(SS(3)) '\n' strjoin(SS(4))...
                 '\n' strjoin(SS(5))  '\n' strjoin(SS(6))  '\n' strjoin(SS(7))  '\n' strjoin(SS(8)) '\n' strjoin(SS(9)) '\n' strjoin(SS(10)) '\n' strjoin(SS(11)) '\n' strjoin(SS(12)) '"' ']\n' ];
            end
            if NetWorkOnly
                FontSize = [',fontsize =' num2str(max(50*(rsort{layer-1}(nodeTo)/rsort{layer-1}(1)).^(1/10),1))];
                Str = ['"' num2str(layer-1) '_' num2str(nodeTo) '"' ' [rank =' num2str(layer+1)  FontSize ' ,width=.01, height=.01, shape=ellipse, label=' '"' num2str(nodeTo) '"' ']\n' ];
            end
            if nnz(strcmp(nodelist,Str))==0
                nodelist{end+1}=['"' num2str(layer-1) '_' num2str(nodeTo) '"'];
                fprintf(fid,Str);
            end
            %fprintf(fid,['"' num2str(layer-1) '_' num2str(nodeTo) '"' ' [ label=' '"' S{layer-1}{nodeTo} '"' ']\n' ]);
            %fprintf(fid,['"' S{layer}{nodeFrom} '"' ' -> ' '"' S{layer-1}{nodeFrom} '"' '\n']);
        end
        for j_in=dex(:)'
            [edgelist,nodelist]=dotwrite_updown(fid,j_in,layer-1,Phisort,S,edgelist,nodelist,endlayer,MaxSize,rsort,NetWorkOnly,minWeight,Colors);
        end
else
%     nodeTo=nodeFrom;
%     SS = strsplit(S{layer}{nodeTo});
%     FontSize = [',fontsize =' num2str((50*(rsort{layer}(nodeTo)/rsort{layer}(1)).^(1/10)))];
%     if layer-1>endlayer
%         Str = ['"' num2str(layer) '_' num2str(nodeTo) '"' ' [rank =' num2str(layer+1)  FontSize ', shape=box, style="rounded", label=' '"' num2str(nodeTo) ' ' strjoin(SS(1:4)) '\n' strjoin(SS(5:8)) '\n' strjoin(SS(9:12)) '"' ']\n' ];
%     else
%         Str = ['"' num2str(layer) '_' num2str(nodeTo) '"' ' [rank =' num2str(layer+1)  FontSize ', shape=box, style="rounded", label=' '"' num2str(nodeTo) ' ' strjoin(SS(1)) '\n' strjoin(SS(2)) '\n' strjoin(SS(3)) '\n' strjoin(SS(4))...
%             '\n' strjoin(SS(5))  '\n' strjoin(SS(6))  '\n' strjoin(SS(7))  '\n' strjoin(SS(8)) '\n' strjoin(SS(9)) '\n' strjoin(SS(10)) '\n' strjoin(SS(11)) '\n' strjoin(SS(12)) '"' ']\n' ];
%     end
%     if NetWorkOnly
%         FontSize = [',fontsize =' num2str(max(50*(rsort{layer}(nodeTo)/rsort{layer}(1)).^(1/10),1))];
%         Str = ['"' num2str(layer) '_' num2str(nodeTo) '"' ' [rank =' num2str(layer+1)  FontSize ' ,width=.01, height=.01, shape=ellipse, label=' '"' num2str(nodeTo) '"' ']\n' ];
%     end
%     if nnz(strcmp(nodelist,Str))==0
%         nodelist{end+1}=['"' num2str(layer) '_' num2str(nodeTo) '"'];
%         fprintf(fid,Str);
%     end
    return;
end
