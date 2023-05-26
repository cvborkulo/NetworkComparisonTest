# NetworkComparisonTest

Statistical comparison of two networks with respect to three invariance measures

Network approaches to psychometric constructs, in which constructs are modeled in terms of interactions between their constituent factors, have rapidly gained popularity in psychology. Applications of such network approaches to various psychological constructs have moved from a descriptive stance, in which the goal is to estimate the network structure that pertains to a construct, to a more comparative stance, in which the goal is to compare network structures across populations. 

This package facilitates this recent movements by providing the necessary methodological tools. The Network Comparison Test (NCT) uses resampling-based permutation testing to compare network structures from two independent, cross-sectional data sets on invariance of 1) network structure, 2) edge (connection) strength, and 3) global strength. 

# Example

```r
library("IsingSampler")
library("IsingFit")

### Simulate binary datasets under null hypothesis:
N <- 6 # Number of nodes
nSample <- 500 # Number of samples

# Ising parameters:
Graph <- matrix(sample(0:1,N^2,TRUE,prob = c(0.8, 0.2)),N,N) * runif(N^2,0.5,2)
Graph <- pmax(Graph,t(Graph))
diag(Graph) <- 0
Thresh <- -rowSums(Graph) / 2

# Simulate:
data1 <- IsingSampler(nSample, Graph, Thresh)
data2 <- IsingSampler(nSample, Graph, Thresh)
colnames(data1) <- colnames(data2) <- c('V1', 'V2', 'V3', 'V4', 'V5', 'V6')

# Testing the three aspects that are validated (network invariance, global strength, edge weight)
# 2 edges are tested here: between variable 1 and 2, 
# and between 3 and 6 (can be list(c(2,1),c(6,3)) as well)
Res_1 <- NCT(data1, data2, gamma=0, it=1000, binary.data = TRUE, 
test.edges=TRUE, edges=list(c(1,2),c(3,6)))

# Plot results of the network structure invariance test:
plot(Res_1, what="network")

# Plot results of global strength invariance test:
plot(Res_1, what="strength")

# Plot results of the edge invariance test:

plot(Res_1, what="edge")
```

# Background Information
For more information on the Network Comparison Test, take a look at:

Van Borkulo, C. D., van Bork, R., Boschloo, L., Kossakowski, J. J., Tio, P., Schoevers, R. A., Borsboom, D., & Waldorp, L. J. (2022). Comparing network structures on three aspects: A permutation test. Psychological Methods. DOI: 10.1037/met0000476

# Bug Reports, Feature Request, or Contributing
If you encounter any bugs or have ideas for new features, you can submit them by creating an issue on Github. Additionally, if you want to contribute to the development of NCT, you can initiate a branch with a pull request; we can review and discuss the proposed changes.
