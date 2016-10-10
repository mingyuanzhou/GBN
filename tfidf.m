function X_tfidf=tfidf(X)

X_tfidf = bsxfun(@rdivide,X,max(sum(X,1),realmin));
X_tfidf = bsxfun(@times,X_tfidf,log(size(X,2))./(sum(X>0,2)+1));
