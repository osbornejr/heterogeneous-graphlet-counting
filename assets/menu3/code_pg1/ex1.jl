# This file was generated, do not modify it. # hide
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