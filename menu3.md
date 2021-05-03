@def title = "Co-expression measures"
@def hascode = true
@def rss = "Running through various types of coexpression measures that we use to analyse genetic data."
@def rss_title = "Coexpression measures"
@def rss_pubdate = Date(2021, 4, 26)

@def tags = ["syntax", "code", "image"]

# Co-expression measures

\toc

## Gene expression data

One obtainable output of RNA sequencing experiments is gene expression count data that represent how present the gene is in a particular sample. 
This can be represented as an $n \times m$ count matrix, where $n$ is the number of genes, and $m$ is the number of samples in the experiment:

|Gene | Sample 1 | Sample 2 | ... | Sample $m$
|-----|----------|----------|-----|----------
|Gene 1 | 12 | 100 |... | 23
|Gene 2 | 0 | 11 |... | 341
|Gene 3 | 223 | 1 |... | 1045
|  .   |   .   |   .   |   . 
|  .   |   .   |   .   |   . 
|  .   |   .   |   .   |   . 
|Gene $n$ | 0 | 0 |... | 19

For the purposes of our analysis, we will view each of these genes as a random variable $X_g$ with $m$ observations:

$$
X_g = \{x_1,x_2,x_3,...,x_m\} 
$$


## Calculating gene variance

Often the first thing we want to identify in the data is a subset of genes which are of interest, i.e. the genes that have the most change in their expression across samples.
One way to do this is to look at the variance of each of our random variables $X_g$. Recall that the variance can be estimated using the formula

$$
\sigma_{X_g}^2 =\frac{1}{m-1}\sum_{i=1}^m (x_i - \bar{x})^2  
$$

If we then identify those genes that have the highest variance, we have a good way to restrict our analsyis to ony those genes that may be responding to changes in experimental conditions.

## Identifying interactions between genes

We then want to see to what degree these genes of interest are interacting with one another. 
The most simple approach here is to look at the (Pearson) correlation between each of our random variables: 

\begin{align} 
\rho_{XY} &= \frac{\sum_{i=1}^m (x_i - \bar{x})(y_i - \bar{y})}{\sqrt{\sum_{i=1}^m (x_i - \bar{x})^2}\sqrt{\sum_{i=1}^m (y_i - \bar{y})^2}}  \\
	&= \frac{cov(X,Y)}{\sigma_X \sigma{y}}
\end{align}

i.e. the correlation is the covariance of $X$ and $Y$ normalised by their standard deviations.
This means that the correlation gives a value that is normalised between $-1$ and $1$, giving us a relative measure of how closely each random variable $X$ and $Y$ correspond to each other across the sample set.

## Forming a co-expression network 

The correlation coefficient for every pair of genes can be arranged as an $n\times n$ similarity matrix $S$.   
Note that $S$ is quite similar to an adjacency matrix representation of a network, and indeed if we were to take the correlation values as the weights for edges then it could be viewed as one.
In our case, we want an unweighted network, so the simplest approach is to apply a hard threshold on the (absolute) correlation values to decide whether or not an edge exists between each gene pair.

## Indirect associations in the network
One problem here is that the correlation between two genes may be largely dependent on their shared relationship with another gene.
In our network, this will result in a large overrepresentation of triangle type structures.     
Partial correlation seeks to eliminate this dependency on a third variable:
$$
\rho_{XY|Z} = \frac{\rho_{XY} - \rho_{XZ}\rho_{YZ}}{\sqrt{1-\rho_{XZ}^2}\sqrt{1-\rho_{YZ}^2}}
$$

## Example

```julia:./code_pg1/ex1
using Statistics
## 4 gene vectors; B and C are both randomly perturbated from A
A = LinRange(1,10,100) + randn(100,1)
B = A + randn(100,1)
C = A + randn(100,1)

##generate correlation matrix for genes
S = cor(hcat(A,B,C))  
@show S
## partial correlation of B and C when we exclude the influence of A:
parBCA = (S[2,3]-S[2,1]*S[3,1])/(sqrt(1-S[2,1]^2)*sqrt(1-S[3,1]^2))
@show parBCA
## partial correlation of B and A when we exclude the influence of C:
parBAC = (S[1,2]-S[2,3]*S[1,3])/(sqrt(1-(S[2,3])^2)*sqrt(1-(S[1,3])^2))
@show parBAC
``` 
\output{./code_pg1/ex1}
