function [edgelist,nodelist] = dotwrite_downup(fid,nodeFrom,startlayer,Phisort,S,edgelist,nodelist,layer,MaxSize,rsort,NetWorkOnly)
%Write topic tree from bottom to top
%Mingyuan Zhou
%August, 2015

if layer<startlayer
        %dex=find(Phisort{layer+1}(nodeFrom,:)>max(1/max(size(Phisort{layer+1},1),20),0.001));
        dex=find(Phisort{layer+1}(nodeFrom,:)>0.01);
        dex=dex(dex<=MaxSize(layer));
        SS = strsplit(S{layer}{nodeFrom});
        FontSize = [',fontsize =' num2str((50*(rsort{layer}(nodeFrom)/rsort{layer}(1))^(1/10)))];
      %  Str=['"' num2str(layer) '_' num2str(nodeFrom) '"' ' [ rank =' num2str(layer) FontSize  ', shape=box, style="rounded", label='  '"' num2str(nodeFrom) ' ' strjoin(SS(1)) '\n' strjoin(SS(2:3)) '\n' strjoin(SS(4:5)) '"' ']\n' ];
        
       % if layer>endlayer
       % Str=['"' num2str(layer) '_' num2str(nodeFrom) '"' ' [ rank =' num2str(layer) FontSize ', shape=box, style="rounded",  label=' '"' num2str(nodeFrom) ' ' strjoin(SS(1:4)) '\n' strjoin(SS(5:8)) '\n' strjoin(SS(9:12)) '"' ']\n' ];
      
        Str = ['"' num2str(layer) '_' num2str(nodeFrom) '"' ' [rank =' num2str(layer+1)  FontSize ', shape=box, style="rounded", label=' '"' num2str(nodeFrom) ' ' strjoin(SS(1)) '\n' strjoin(SS(2)) '\n' strjoin(SS(3)) '\n' strjoin(SS(4))...
                 '\n' strjoin(SS(5))  '\n' strjoin(SS(6))  '\n' strjoin(SS(7))  '\n' strjoin(SS(8)) '\n' strjoin(SS(9)) '\n' strjoin(SS(10)) '\n' strjoin(SS(11)) '\n' strjoin(SS(12)) '"' ']\n' ];

        %else
       % Str=['"' num2str(layer) '_' num2str(nodeFrom) '"' ' [ rank =' num2str(layer) FontSize ', shape=box, style="rounded", color=' Colors{layer}  ', label=' '"' num2str(nodeFrom) ' ' strjoin(SS(1)) '\n' strjoin(SS(2:3)) '\n' strjoin(SS(4:5)) '\n' strjoin(SS(6:7)) '"' ']\n' ];
       %end
        
        
        if NetWorkOnly
            Str=['"' num2str(layer) '_' num2str(nodeFrom) '"' ' [ rank =' num2str(layer) FontSize  ',width=.01, height=.01, shape=ellipse, label='  '"' num2str(nodeFrom) '"' ']\n' ];
        end
        
        if nnz(strcmp(nodelist,Str))==0
            nodelist{end+1}=['"' num2str(layer) '_' num2str(nodeFrom) '"'];
            fprintf(fid,Str);
        end

        %fprintf(fid,['"' num2str(layer) '_' num2str(nodeFrom) '"' ' [ label=' '"' S{layer}{nodeFrom} '"' ']\n' ]);
        for nodeTo = dex(:)'
            if NetWorkOnly
                Str = [ '"' num2str(layer+1) '_' num2str(nodeTo) '"' ' -> '  '"' num2str(layer) '_' num2str(nodeFrom) '"' ' [penwidth=' num2str(max(5*(Phisort{layer+1}(nodeFrom,nodeTo)).^(1/2),0)) ']' '\n'];
            else
                Str = [ '"' num2str(layer+1) '_' num2str(nodeTo) '"' ' -> '  '"' num2str(layer) '_' num2str(nodeFrom) '"' ' [penwidth=' num2str(max(15*(Phisort{layer+1}(nodeFrom,nodeTo)).^(1/2),0)) ']' '\n'];
                
            end
            if nnz(strcmp(edgelist,Str))==0
                edgelist{end+1}=Str;
                fprintf(fid,Str);
            end
            
            
            SS = strsplit(S{layer+1}{nodeTo});
            FontSize = [',fontsize =' num2str((50*(rsort{layer+1}(nodeTo)/rsort{layer+1}(1)).^(1/10)))];
         %   Str = ['"' num2str(layer+1) '_' num2str(nodeTo) '"' ' [rank =' num2str(layer+1) FontSize  ', shape=box, style="rounded", label=' '"' num2str(nodeTo) ' ' strjoin(SS(1)) '\n' strjoin(SS(2:3)) '\n' strjoin(SS(4:5)) '"' ']\n' ];
            
            Str = ['"' num2str(layer+1) '_' num2str(nodeTo) '"' ' [rank =' num2str(layer+1)  FontSize ', shape=box, style="rounded", label=' '"' num2str(nodeTo) ' ' strjoin(SS(1)) '\n' strjoin(SS(2)) '\n' strjoin(SS(3)) '\n' strjoin(SS(4))...
                 '\n' strjoin(SS(5))  '\n' strjoin(SS(6))  '\n' strjoin(SS(7))  '\n' strjoin(SS(8)) '\n' strjoin(SS(9)) '\n' strjoin(SS(10)) '\n' strjoin(SS(11)) '\n' strjoin(SS(12)) '"' ']\n' ];

            if NetWorkOnly
                Str = ['"' num2str(layer+1) '_' num2str(nodeTo) '"' ' [rank =' num2str(layer+1) FontSize  ',width=.01, height=.01, shape=ellipse, label=' '"' num2str(nodeTo)  '"' ']\n' ];
            end
            if nnz(strcmp(nodelist,Str))==0
                nodelist{end+1}=['"' num2str(layer+1) '_' num2str(nodeTo) '"'];
                fprintf(fid,Str);
            end
            %fprintf(fid,['"' num2str(layer-1) '_' num2str(nodeTo) '"' ' [ label=' '"' S{layer-1}{nodeTo} '"' ']\n' ]);
            %fprintf(fid,['"' S{layer}{nodeFrom} '"' ' -> ' '"' S{layer-1}{nodeFrom} '"' '\n']);
        end
        for j_in=dex(:)'
            [edgelist,nodelist]=dotwrite_downup(fid,j_in,startlayer,Phisort,S,edgelist,nodelist,layer+1,MaxSize,rsort,NetWorkOnly);
        end
else
    return;
end
