function X_idf = idf(X)
X_idf = log(size(X,2)./(sum(X>0,2)+1));
