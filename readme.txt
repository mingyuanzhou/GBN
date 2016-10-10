Matlab code for the following two papers:

[1] Mingyuan Zhou, Yulai Cong, and Bo Chen, “Augmentable gamma belief networks," Journal of Machine Learning Research, vol. 17, pp. 1-44, Sept. 2016. 

[2] Mingyuan Zhou, Yulai Cong, and Bo Chen, "The Poisson gamma belief network," Neural Information Processing Systems (NIPS2015), Montreal, Canada, Dec. 2015.

The code is written by Mingyuan Zhou and contributed by Yulai Cong.


(1) Compile c mex files if the provided mex files do not work on your computers:
make.m

(2) Deep topic modeling with Poisson gamma belief network (PGBN), using perplexity on heldout word tokens to evaluate the performance:

Demo_PGBN_Perplexity.m	%produce the results used to plot Figures 12-13 of [1]		
PGBN_Layerwise_perplexity.m
	
(3) Unsupervised multilayer feature learning on a word-document count matrix with the Poisson gamma belief network (PGBN). Perform logistic regression with the extracted features at the first layer to evaluate the performance:

Demo_PGBN_FeatureExtraction.m	%produce the results used to plot Figures 3-7, 9-11, and 17-21 and Table 3 of [1]	
GBN_Layerwise.m
GNB_Testing.m

(4) Unsupervised multilayer feature learning on a binary matrix with the Bernoulli-Poisson gamma belief network (BerPo-GBN). Perform logistic regression with the extracted features at the first layer to evaluate the performance:

Demo_BerPo_GBN_FeatureExtraction.m	%produce the results used to plot Figures 14 of [1]	
GBN_Layerwise.m
GNB_Testing.m
truncated_Poisson_rnd.m

(5) Unsupervised multilayer feature learning on a nonnegative real matrix with the Poisson randomized gamma gamma belief network (PRG-GBN). Perform logistic regression with the extracted features at the first layer to evaluate the performance:

Demo_PRG_GBN_FeatureExtraction.m	%produce the results used to plot Figures 15 and 22-24 of [1]	
GBN_Layerwise.m
GNB_Testing.m
Truncated_bessel_rnd.m

(6) Logistic regression to produce the results for Table 1 of [1]
Classify_Direct.m



Folders:
(i) plot_images_tree: the folder contains the files and code to reproduce Figures 22-24 of [1]

(ii) plot_topics_tree: the folder contains the files and code to reproduce Figures 3-7 and 17-21 of [1]

(iii) misc.: files that have not been cleaned



Supporting functions:

(a) Sample from the truncated Bessel distribution using its probability mass function, starting from its mode and sampling along both sides of the mode:
Truncated_bessel_rnd.m
Rand_Truncated_PMF_bessel.c

(b) Prune the inactive factors of the current top hidden layer:
TrimTcurrent.m

(c) Assign the word tokens to latent topics using the multinomial distribution and/or upward propagate latent counts using the CRT distribution:
GNBP_mex_collapsed_deep.c (collapsed Gibbs sampling for the first layer)
Multrnd_Matrix_mex_fast_v1.c (blocked Gibbs sampling for the first layer)
CRT_Multrnd_Matrix.c 
CRT_sum_mex_matrix_v1.c
CRT_sum_mex_v1.c,

(d) Efficiently calculate and store the matrix Xmask.*(Phi*Theta) if Xmask is a sparse indicator matrix:
Mult_Sparse.c

(e) Plot the probability distribution functions of the Poisson randomized gamma and truncated Bessel distributions
PRG_Bessel_Figure.m (used to produce Figure 16)		

(f) Display the inferred basis vectors at various hidden layers as images:
DispDictionary.m	
DispDictionary_asym.m	(specifying image sizes)

(g)*_par.m (a parallel version of the function *)

Copyright (c), 2016, Mingyuan Zhou
http://mingyuanzhou.github.io/

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.