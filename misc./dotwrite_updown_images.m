function [edgelist,nodelist] = dotwrite_updown_images(fid,nodeFrom,layer,Phisort,S,edgelist,nodelist,endlayer,MaxSize,rsort,NetWorkOnly,K0,filepath)
%Write topic tree from top to bottom
%Mingyuan Zhou
%August, 2015

if layer>endlayer
        dex=find(Phisort{layer}(:,nodeFrom)>max(1/max(size(Phisort{layer},1),20),0.01));
        dex=dex(dex<=MaxSize(layer));
        i=nodeFrom;
        SS = [filepath '/K0' num2str(K0) '_layer' num2str(layer) '_' num2str(i) '.png'];
        FontSize = [',fontsize =' num2str((50*(rsort{layer}(nodeFrom)/rsort{layer}(1))^(1/10)))];
       if layer>endlayer
        Str=['"' num2str(layer) '_' num2str(nodeFrom) '"' ' [ rank =' num2str(layer) FontSize ', shape=none, image="' SS  '"]\n' ];
       else
        Str=['"' num2str(layer) '_' num2str(nodeFrom) '"' ' [ rank =' num2str(layer) FontSize ', shape=none, image="' SS  '"]\n' ];
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
                Str = [ '"' num2str(layer) '_' num2str(nodeFrom) '"' ' -> '  '"' num2str(layer-1) '_' num2str(nodeTo) '"' ' [penwidth=' num2str(max((8*(Phisort{layer}(nodeTo,nodeFrom)).^(1/2)),0)) ',dir=forward]' '\n'];
             else
                 Str = [ '"' num2str(layer) '_' num2str(nodeFrom) '"' ' -> '  '"' num2str(layer-1) '_' num2str(nodeTo) '"' ' [penwidth=' num2str(max((8*(Phisort{layer}(nodeTo,nodeFrom)).^(1/2)),0)) ',dir=forward]' '\n'];
             end
            if nnz(strcmp(edgelist,Str))==0
                edgelist{end+1}=Str;
                fprintf(fid,Str);
            end
            
            i=nodeTo;
            SS = [filepath '/K0' num2str(K0) '_layer' num2str(layer-1) '_' num2str(i) '.png'];
            FontSize = [',fontsize =' num2str((50*(rsort{layer-1}(nodeTo)/rsort{layer-1}(1)).^(1/10)))];
            if layer-1>endlayer
            Str = ['"' num2str(layer-1) '_' num2str(nodeTo) '"' ' [rank =' num2str(layer+1)  FontSize ', shape=none, image="' SS  '"]\n' ];
            else
            Str = ['"' num2str(layer-1) '_' num2str(nodeTo) '"' ' [rank =' num2str(layer+1)  FontSize ', shape=none, image="' SS  '"]\n' ];
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
            [edgelist,nodelist]=dotwrite_updown_images(fid,j_in,layer-1,Phisort,S,edgelist,nodelist,endlayer,MaxSize,rsort,NetWorkOnly,K0,filepath);
        end
else
    return;
end
