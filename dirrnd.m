function G = dirrnd(A,option)
if nargin<2
    option = 'Gibbs';
end
if strcmp(option,'Gibbs')
    G = gamrnd(A,1);
    G = bsxfun(@rdivide, G, sum(G,1));
else
    G = bsxfun(@rdivide, A, sum(A,1));
end