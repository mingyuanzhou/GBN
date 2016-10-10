function edgelist = dotwrite(fid,nodeFrom,layer,Phisort,S,edgelist,endlayer)
if layer>endlayer
        dex=find(Phisort{layer}(:,nodeFrom)>max(5/max(size(Phisort{layer},1),20),0.05));
        SS = strsplit(S{layer}{nodeFrom});
        fprintf(fid,['"' num2str(layer) '_' num2str(nodeFrom) '"' ' [ label=' '"' num2str(nodeFrom) ' ' strjoin(SS(1)) '\n' strjoin(SS(2:3)) '\n' strjoin(SS(4:5)) '"' ']\n' ]);
        %fprintf(fid,['"' num2str(layer) '_' num2str(nodeFrom) '"' ' [ label=' '"' S{layer}{nodeFrom} '"' ']\n' ]);
        for nodeTo = dex(:)'
            
            Str = ['"' num2str(layer) '_' num2str(nodeFrom) '"' ' -> ' '"' num2str(layer-1) '_' num2str(nodeTo) '"' '\n'];
            if nnz(strcmp(edgelist,Str))==0
                edgelist{end+1}=Str;
                fprintf(fid,Str);
            end
            SS = strsplit(S{layer-1}{nodeTo});
            fprintf(fid,['"' num2str(layer-1) '_' num2str(nodeTo) '"' ' [ label=' '"' num2str(nodeTo) ' ' strjoin(SS(1)) '\n' strjoin(SS(2:3)) '\n' strjoin(SS(4:5)) '"' ']\n' ]);
            %fprintf(fid,['"' num2str(layer-1) '_' num2str(nodeTo) '"' ' [ label=' '"' S{layer-1}{nodeTo} '"' ']\n' ]);
            %fprintf(fid,['"' S{layer}{nodeFrom} '"' ' -> ' '"' S{layer-1}{nodeFrom} '"' '\n']);
        end
        for j_in=dex(:)'
            edgelist=dotwrite(fid,j_in,layer-1,Phisort,S,edgelist,endlayer);
        end
else
    return;
end
