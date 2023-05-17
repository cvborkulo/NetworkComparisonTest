suppressPackageStartupMessages({
  library(IsingSampler)
  library(IsingFit)
  library(bootnet)
})

### Simulate binary datasets under null hypothesis:
### underlying network structures are similar
# Input:
N <- 6 # Number of nodes
nSample <- 500 # Number of samples

# Ising parameters:
set.seed(123)
Graph <- matrix(sample(0:1,N^2,TRUE,prob = c(0.8, 0.2)),N,N) * runif(N^2,0.5,2)
Graph <- pmax(Graph,t(Graph))
Graph[4,1] <- Graph[4,1]*-1
Graph[1,4] <- Graph[1,4]*-1
Graph[5,1] <- Graph[5,1]*-1
Graph[1,5] <- Graph[1,5]*-1
Graph[6,1] <- Graph[6,1]*-1
Graph[1,6] <- Graph[1,6]*-1
diag(Graph) <- 0
Thresh <- -rowSums(Graph) / 2

# Simulate:
data1 <- IsingSampler(nSample, Graph, Thresh)
data2 <- IsingSampler(nSample, Graph, Thresh)
colnames(data1) <- colnames(data2) <- c('V1', 'V2', 'V3', 'V4', 'V5', 'V6')


set.seed(123)
NCT_a <- NCT(data1, data2, gamma=0, it=10, binary.data = TRUE, 
             test.edges=TRUE, edges=list(c(1,2),c(3,6)), progressbar = FALSE)

test_that("NCT works when directly inputting data", {
  expect_snapshot(summary(NCT_a))
})

## Plot results of global strength invariance test (not reliable with only 10 
# permutations!)
# plot(NCT_a, what="strength") # TODO: test this!

## (2) Feeding the estimateNetwork() output into NCT
est_1 <- estimateNetwork(data1, default = "IsingFit", tuning = 0, verbose = FALSE)
est_2 <- estimateNetwork(data2, default = "IsingFit", tuning = 0, verbose = FALSE)
## When using estimateNetwork() output, there is no need to specify gamma and binary.data 
## This yields similar output as NCT_a
set.seed(123)
NCT_b <- NCT(est_1, est_2, it=10, test.edges=TRUE, 
             edges=list(c(1,2),c(3,6)), progressbar = FALSE, verbose = FALSE)
test_that("NCT works with output from estimateNetwork", {
  expect_snapshot(summary(NCT_b))
})




## Next, an example of testing whether there are differences in node strength 
# when data is paired (e.g., a group which is measured pre- and post-treatement). 
# Also, here you can see how to specify that you want to take the sign of node strength 
# into account (by default, the absolute value is taken and, therefore, the sign is 
# ignored).

## abs = FALSE
set.seed(123)
NCT_c = NCT(est_1, est_2, paired = TRUE, abs = FALSE, test.edges = TRUE, 
            edges = list(c(1,2),c(3,6)), test.centrality = TRUE, 
            centrality = c("strength"), nodes = "all", it=10, progressbar = FALSE)

test_that("NCT works when testing node strength", {
  expect_snapshot(summary(NCT_c))
})

## Finally, an example how to test for differences in centrality (e.g., expectedInfluence)

set.seed(123)
NCT_d = NCT(est_1, est_2, paired = TRUE, abs = FALSE, test.edges = TRUE, 
            edges = list(c(1,2),c(3,6)), test.centrality = TRUE, 
            centrality = c("expectedInfluence"), nodes = "all", it=10, progressbar = FALSE)

test_that("NCT works when testing expected influence", {
  expect_snapshot(summary(NCT_d))
})
