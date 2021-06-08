#########This code simulates unlinked loci and allows multiple samples of size S <= Ne from the final population
Ne = 200
NGens = 7
NReps = 100
NSamps = 10
phi1 = matrix(NA,NReps,length(Sampsize))

Sampsize = c(100,50,25)
NLociStart =20000  ## should be at least 1.5 times NSamps*max(NLoci)
NLoci = c(1000,500,250,100)
##NLoci = c(3000,1000,500,250,100,50)
TargetLoci = NSamps*NLoci[1]

tempphi = 1:length(NLoci)
FinalPhi = matrix(NA,length(Sampsize),length(NLoci))
colnames(FinalPhi) = NLoci
rownames(FinalPhi) = Sampsize
FinalMean = FinalPhi

###########                  
reproduce <- function(parents){
  noffspring = length(parents1)
  # L = number of loci (columns)
  L = dim(parents)[[2]]
  
  # construct matrices of the parental genotypes
  parents_geno1 = parents[parents1, ]
  parents_geno2 = parents[parents2, ]
  
  # offspring are a combination of the two parents
    offspring <- rbinom(n=length(parents_geno1), size = 1, p = parents_geno1 / 2) + 
      rbinom(n = length(parents_geno2), size = 1, p = parents_geno2 / 2)
    # convert offspring vector back to a matrix
  offspring <- matrix(data = offspring, nrow =noffspring, ncol = L )
  
  return(offspring) 
                                        } # end function                                                                   
##############
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

##return(mean(rsq,na.rm=T)) 
 return(rsq)
 } # end function  	        

#####################

## set up Pop matrices

PopPhi = list()
for (P in 1:NReps)  {
PopPhi[[P]] = matrix(NA,length(Sampsize),length(NLoci)) 
rownames(PopPhi[[P]]) = Sampsize
colnames(PopPhi[[P]]) = NLoci
}

## We are simulating a monoecious population with random selfing = original Wright-Fisher ideal population

for (jj in 1:NReps) {

##Initialize population as 100% heterozygotes so starting P = 0.5 at each locus
Geno1 = matrix(1,Ne,NLociStart)

 print(paste0("Population = ",jj))
  flush.console() 

for (k in 1:NGens)    {

##For each replicate, burnin with NGens generations of random mating and genetic drift to equilibrate LD
parents1 <- sample.int(n = Ne, size = Ne, replace = TRUE)
parents2 <- sample.int(n = Ne, size = Ne, replace = TRUE)

## create Ne offspring in generation k, each with NLociStart genotypes
Geno2 = reproduce(Geno1)

## offpsring become parents of next generation
Geno1 = Geno2 
                      } # end for k
                                            
                      
  ### Set up sample size matrices

SampRsq = list()
for (SS in 1:length(Sampsize))  {
SampRsq[[SS]] = matrix(NA,length(NLoci),NSamps) 
}

for (SS in 1:length(Sampsize))  {

         print(paste0("Working on samplesize",Sampsize[SS]))
	  flush.console() 
     
     BigR = matrix(NA,TargetLoci,TargetLoci)	  
     
 for (q in 1:NSamps)  {	  
      ## Draw a random sample from the Ne individuals in final generation
         GenoS = Geno1[sample(nrow(Geno1),Sampsize[SS],replace=F),]   
         P = colMeans(GenoS)/2 
      
      ##eliminate loci monomorphic or nearly so  
      cut = 1/Sampsize[SS] ## eliminates singletons and doubletons
      GenoS= GenoS[, P > cut & P < (1-cut), drop = F] 
      ## reduce to target loci
      GenoS = GenoS[,1:TargetLoci]     
 
     ## get r^2 for a chunk of loci
           startR = 1+(q-1)*NLoci[1]  
           stopR = startR + NLoci[1] - 1
           BigR[startR:stopR,startR:stopR] = GetLD(GenoS[,startR:stopR])
             print(paste0("Got Rsq in subset ",q))
             flush.console()    
   }  # end for q	        
	          
     ## Get mean rsq for subsets of loci
        for (L in 1:length(NLoci))  {
          MeanR = 1:NSamps
              
          for (q in 1:NSamps)  {
               startL = 1+(q-1)*NLoci[L]  
               stopL = startL + NLoci[L] - 1
               MeanR[q] = mean(BigR[startL:stopL,startL:stopL],na.rm=T)  
           } ## end for q
           
          SampRsq[[SS]][L,] = MeanR        
        }   ## end for L
        
      for(L in 1:length(NLoci))  {
      tempphi[L] = var(SampRsq[[SS]][L,])/mean(SampRsq[[SS]][L,])^2
      }  # end for L

    PopPhi[[jj]][SS,] = tempphi
    
 }  ## end for SS   

}  ## end for jj                   
                      
  ##average phi across Pop matrices
  dummy = matrix(0,length(Sampsize),length(NLoci)) 
  for (PD in 1:NReps)  {
  dummy = dummy + PopPhi[[PD]]  
 }
  FinalPhi = dummy/NReps

X=2/FinalPhi  ## nprime
X
Y = NLoci*(NLoci-1)/2 ## n = total number of pairs of loci
t(t(X)/Y)
