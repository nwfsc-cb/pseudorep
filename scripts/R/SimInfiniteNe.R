## simulating infinite Ne, with each sample using different loci

NLoci = c(2000,500,100)
Sampsize = c(100,50)
NSamps = 5
Multiplier = 2
BigLoci = NLoci[1]*Multiplier*NSamps
X = c(2,1,0)
p = c(0.05,0.15,0.25,0.35,0.45)
proportions = c(9,5,3,2,1)

NReps = 500

MRSQ = list()
for (LL in 1:length(NLoci))  {
  MRSQ[[LL]] = array(NA, dim = c(NReps,NSamps,length(Sampsize)))
  }
  
############                                                                             
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

 TotL = dim(Geno)[[2]]
 for (i in 1:TotL)  {  
 for (j in i:TotL)  {  rsq[j,i] = NA}}    

return(rsq)
 } # end function  	
###############################

for (k in 1:NReps)  {

Freqs = sample(p,BigLoci,replace=T,prob=proportions)  ## get L-shaped distribution of parametric allele frequencies
 
  print(paste0("Replicate = ",k))
  flush.console() 
  
for (q in 1:NSamps) {

  for (qq in 1:length(Sampsize))  {
  
   Genos = matrix(NA,Sampsize[qq],(2*NLoci[1]))
   start = (q-1)*Multiplier*NLoci[1] + 1
   stop = start + Multiplier*NLoci[1] - 1

   count = 0
    for (m in start:stop)  {  
    count=count+1
    ## pick each genotype at each locus randomly for Sampsize individuals 
    A = Freqs[m]^2
    B = 2*Freqs[m]*(1-Freqs[m])
    C = (1-Freqs[m])^2
    G = sample(X,Sampsize[qq],replace=T,prob=c(A,B,C))
    Genos[,count] = G
    }  ### end for m

  cut = 1/Sampsize[qq]
  P = colMeans(Genos)/2
  GenoS1 = Genos[, P > cut & P < (1-cut), drop = F]  
  GenoS1 = GenoS1[,1:NLoci[1]]
 
BigD = GetLD(GenoS1)

    for (LL in 1:length(NLoci))  {
    MRSQ[[LL]][k,q,qq] = mean(BigD[1:NLoci[LL],1:NLoci[LL]],na.rm = T)
    }   ## end for LL
  
  }   ## end for qq
   
 }   ## end for q
 
}  # end for k
  
 Phi = array(NA, dim = c(NReps,length(Sampsize),length(NLoci)))
 
  for (LL in 1:length(NLoci))  {
  for (qq in 1:length(Sampsize))  {
    for (k in 1:NReps)  {
    Phi[k,qq,LL] = var(MRSQ[[LL]][k,,qq])/mean(MRSQ[[LL]][k,,qq])^2
    }  ## end for k
    }  ## end for qq
    }  ## end for LL
  
nprime = 2/colMeans(Phi)  
colnames(nprime) = NLoci
rownames(nprime) = Sampsize
nprime

nprimen = nprime
for (j in 1:length(Sampsize))  {
nprimen[j,]=2*nprime[j,]/(NLoci*(NLoci-1))
}
nprimen

