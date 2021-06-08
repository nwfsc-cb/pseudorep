### this code simulates data like that shown in Figure S8.
Ne = 20
NGens1 = 6
NReps = 1000
NSets = 100
NLoci = 10
NPairs = NLoci*(NLoci-1)/2
PairPairs = NPairs*(NPairs-1)/2
bottom = NLoci*(NLoci-1)/2 - 1
top = (NLoci-2)*(NLoci-3)/2

FMiss = top/bottom # fraction of PairPairs that have no overlap
NMiss = FMiss * PairPairs # number of PairPairs that have no overlap
NHits = PairPairs - NMiss # number of PairPairs that have overlap of one locus
NPairs
PairPairs
NHits
NMiss

BigRSQ = array(NA, dim = c(NLoci,NLoci,NReps))
CorHits = rep(NA,NHits)
CorMiss = rep(NA,NMiss)
PairLoci = matrix(NA,2,NPairs)
meanHits = 1:NSets
meanMiss = 1:NSets

counter = 0
for (j in 1:(NLoci-1))  {
for (k in (j+1):NLoci)  {
counter = counter + 1
PairLoci[1,counter] = j
PairLoci[2,counter] = k
}}  ##2 end for j,k


###########                  
reproduce <- function(parents){
  # N = number of inds (rows)
  n_parents = dim(parents)[[1]]
  # L = number of loci (columns)
  L = dim(parents)[[2]]
  
  # construct matrices of the parental genotypes
  parents_geno1 = parents[parents1, ]
  parents_geno2 = parents[parents2, ]
  
  # offspring are a combination of the two parents
    offspring <- rbinom(n=length(parents_geno1), size = 1, p = parents_geno1 / 2) + 
      rbinom(n = length(parents_geno2), size = 1, p = parents_geno2 / 2)
    # convert offspring vector back to a matrix
  offspring <- matrix(data = offspring, nrow =Ne, ncol = L )
  
  return(offspring) 
                                        } # end function                                                                   
#############
####### #######
GetLD <- function(Geno)     { 

mat = data.matrix(Geno)
mat = t(mat)
# Center each variable
mat = mat - rowMeans(mat);
# Standardize each variable
mat = mat / sqrt(rowSums(mat^2));   
# Calculate correlations
r = tcrossprod(mat);
rsq = r^2 

 return(rsq)      } # end function  	        

##############

#### Main program                       
## We are simulating a monoecious population with random selfing = original Wright-Fisher ideal population

for (q in 1:NSets)  {
print(paste0("Set = ",q))
flush.console()    

for (jj in 1:NReps)  {

##Initialize population as 100% heterozygotes so starting P = 0.5 at each locus
Geno0 = matrix(1,Ne,NLoci)

##Burnin with NGens1 generations of random mating and genetic drift to equilibrate LD

for (kk in 1:NGens1)    {
 parents1 <- sample.int(n = Ne, size = Ne, replace = TRUE)
 parents2 <- sample.int(n = Ne, size = Ne, replace = TRUE)
## create Ne offspring in generation k, each with NLoci genotypes
Geno1 = reproduce(Geno0)

## offpsring become parents of next generation
Geno0 = Geno1
                      } # end for kk

BigRSQ[,,jj] = GetLD(Geno0)

}  ## end for jj

A = c(NA,NA)
B = A
counthits = 0
countmiss = 0
for (j in 1:(NPairs-1))  {
A = PairLoci[,j]
for (k in (j+1):NPairs)  {
B = PairLoci[,k]
## check for hits
   if (A[1] == B[1] | A[1] == B[2] | A[2] == B[1] | A[2] == B[2])  { 
     counthits = counthits+1
     CorHits[counthits] = cor(BigRSQ[A[1],A[2],],BigRSQ[B[1],B[2],],use="pairwise.complete.obs") }
     else {
     countmiss = countmiss+1
     CorMiss[countmiss] = cor(BigRSQ[A[1],A[2],],BigRSQ[B[1],B[2],],use="pairwise.complete.obs") }
  }}  ## end for j,k
  
meanHits[q] = mean(CorHits)
meanMiss[q] = mean(CorMiss)

}  ## end for q

##meanHits  # data for each set
##meanMiss

mean(meanHits) # correlations of correlations of pairs that share one locu
mean(meanMiss) # pairs that don't
   

