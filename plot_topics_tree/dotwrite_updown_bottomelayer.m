function [edgelist,nodelist] = dotwrite_updown_bottomelayer(fid,nodeFrom,layer,Phisort,S,edgelist,nodelist,endlayer,MaxSize,rsort,NetWorkOnly)
%Write topic tree from top to bottom
%Mingyuan Zhou
%August, 2015

nodeTo=nodeFrom;
SS = strsplit(S{layer}{nodeTo});
FontSize = [',fontsize =' num2str((50*(rsort{layer}(nodeTo)/rsort{layer}(1)).^(1/10)))];
if layer-1>endlayer
    Str = ['"' num2str(layer) '_' num2str(nodeTo) '"' ' [rank =' num2str(layer+1)  FontSize ', shape=box, style="rounded", label=' '"' num2str(nodeTo) ' ' strjoin(SS(1:4)) '\n' strjoin(SS(5:8)) '\n' strjoin(SS(9:12)) '"' ']\n' ];
else
    Str = ['"' num2str(layer) '_' num2str(nodeTo) '"' ' [rank =' num2str(layer+1)  FontSize ', shape=box, style="rounded", label=' '"' num2str(nodeTo) ' ' strjoin(SS(1)) '\n' strjoin(SS(2)) '\n' strjoin(SS(3)) '\n' strjoin(SS(4))...
        '\n' strjoin(SS(5))  '\n' strjoin(SS(6))  '\n' strjoin(SS(7))  '\n' strjoin(SS(8)) '\n' strjoin(SS(9)) '\n' strjoin(SS(10)) '\n' strjoin(SS(11)) '\n' strjoin(SS(12)) '"' ']\n' ];
end
if NetWorkOnly
    FontSize = [',fontsize =' num2str(max(50*(rsort{layer}(nodeTo)/rsort{layer}(1)).^(1/10),1))];
    Str = ['"' num2str(layer) '_' num2str(nodeTo) '"' ' [rank =' num2str(layer+1)  FontSize ' ,width=.01, height=.01, shape=ellipse, label=' '"' num2str(nodeTo) '"' ']\n' ];
end
if nnz(strcmp(nodelist,Str))==0
    nodelist{end+1}=['"' num2str(layer) '_' num2str(nodeTo) '"'];
    fprintf(fid,Str);
end
