function X_tf=Term_Freq(X)
X_tf = bsxfun(@rdivide,X,max(sum(X,1),realmin));
