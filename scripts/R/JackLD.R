## This file reads simulated genotype files, generates new samples, and estimates nprime for each
## sample using the Jones et al. jackknife method

beginning = "Jack-N"
path1 = "mypath1"
path2 = "mypath2"
pathout="mypath3"

NPeds = 2
NPops = 2
NMut = 4

mapend = ".chrom_counts"
end = ".geno"

Ne = 200
Chr = 4  ## number of chromosomes
NLoci = c(300,100)
TargetLoci = NLoci[1]
Sampsize = c(40,25)
NSamples = 2
NSets = NPeds*NPops*NMut*NSamples
RepData = matrix(NA,NSets,4)
colnames(RepData) = c("Ped","Pop","Mut","Sample")

JackVar = array(NA, dim = c(NSets,length(Sampsize),length(NLoci)))
colnames(JackVar) = Sampsize
JackPhi = JackVar
JackMean = JackVar
JackVarSeparate = JackVar
JackPhiSeparate = JackVar
JackMeanSeparate = JackVar
FullR = JackVar
FullRSeparate = JackVar

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

return(rsq)      } # end function  	        

#####################

sets = 0

for (ped in 1:NPeds)  {
for (pop in 1:NPops)   {
for (mut in 1:NMut)     {

MeanRsq = list()
  for (j in 1:length(Sampsize))  {
  MeanRsq[[j]] = rep(NA,Sampsize[j])
  } # end j 

MeanRsqSeparate = MeanRsq
 
 print(paste0("Population = ",ped,pop,mut))
  flush.console() 

filetoread = paste("Ne_",Ne,".PED_",(ped-1),".Chr_",Chr,".NPOP_4",".MUT_",(mut-1),".SINGLEPOP_",pop,end,sep="")
filetoread2 = paste("Ne_",Ne,".PED_",(ped-1),".Chr_",Chr,".NPOP_4",".MUT_",(mut-1),".SINGLEPOP_",pop,mapend,sep="")
maps = read.table(file=paste(path2,filetoread2,sep=""))

Geno1 = read.table(file=paste(path1,filetoread,sep=""))
Geno1 = t(Geno1)

cvec = 999
for (j in 1:Chr)  {
a = rep(j,maps[j,1])
cvec = c(cvec,a)
}
cvec = cvec[-1]   ## vector giving chromosome for each locus

z = rbind(cvec,Geno1)

z = z[,sample(1:dim(z)[[2]])]  ## randomize order of loci

z = as.matrix(z)
dimnames(z) = NULL

cvec2 = z[1,]   ### sorted vector of chromosome data   
z = z[-1,]  ## first remove first row of chromosome data

 for (q in 1:NSamples)  {  
     sets = sets + 1
     RepData[sets,1] = ped
     RepData[sets,2] = pop
     RepData[sets,3] = mut
     RepData[sets,4] = q
          
  for (qq in 1:length(Sampsize))  {  
     print(paste0("Sample ",q,"Sample size = ",  Sampsize[qq]))
     flush.console()     
    
           GenoS = z[sample(nrow(z),Sampsize[qq],replace=F),]  ### Get a random sample of individuals
           P = colMeans(GenoS)/2     
	       ##eliminate loci monomorphic or nearly so  
	       GenoS = rbind(cvec2,GenoS)  ## replace chr vector for screening loci
	       cut = 1/Sampsize[qq]
	       z2 = GenoS[, P > cut & P < (1-cut), drop = F]
	       z2 = z2[,sample(1:dim(z2)[[2]])]  ## randomize order of loci
	       ## reduce to target loci
	       z2 = z2[,1:TargetLoci] 	     
	       c3 = z2[1,]   ### save sorted and trimmed vector of loci that will be used
	       z2 = z2[-1,]
	       
	       TempR = GetLD(z2)  ## Get data for full sample size of individuals and loci
	       SepR = TempR
	         for (L1 in 1:(NLoci[1]-1))   {
	         for (L2 in (L1+1):NLoci[1])  {
	         if(c3[L1]==c3[L2]) { SepR[L1,L2] = NA  }
	          }}  ## end for L1, L2
		       
       for (LL in 1:length(NLoci))  {
	       FullR[sets,qq,LL] = mean(TempR[1:NLoci[LL],1:NLoci[LL]],na.rm=T)
               FullRSeparate[sets,qq,LL] = mean(SepR[1:NLoci[LL],1:NLoci[LL]],na.rm=T)
               
         for (j in 1:Sampsize[qq])  {
          ##    print(paste0("Jacksample = ",j))
	  ##    flush.console() 

            BigR = GetLD(z2[-j,1:NLoci[LL]])
            MeanRsq[[qq]][j] = mean(BigR,na.rm=T)  
         
            BigRSeparate = BigR
                 for (L1 in 1:(NLoci[LL]-1))   {
	          for (L2 in (L1+1):NLoci[LL])  {
	          if(c3[L1]==c3[L2]) { BigRSeparate[L1,L2] = NA  }
	          }}  ## end for L1, L2
             MeanRsqSeparate[[qq]][j] = mean(BigRSeparate,na.rm=T)  ## mean rsq excluding within-chromosome comparisons
         }  ## end for j
         
         JackVar[sets,qq,LL] = (0.84^2)*var(MeanRsq[[qq]])*(Sampsize[qq]-1)^2/(Sampsize[qq])  ### 0.84 is Jones et al adjustment
         JackMean[sets,qq,LL] = mean(MeanRsq[[qq]])
         JackPhi[sets,qq,LL] = JackVar[sets,qq,LL]/JackMean[sets,qq,LL]^2
         JackVarSeparate[sets,qq,LL] = (0.84^2)*var(MeanRsqSeparate[[qq]])*(Sampsize[qq]-1)^2/(Sampsize[qq])  ### 0.84 is Jones et al adjustment
         JackMeanSeparate[sets,qq,LL] = mean(MeanRsqSeparate[[qq]])
         JackPhiSeparate[sets,qq,LL] = JackVarSeparate[sets,qq,LL]/JackMeanSeparate[sets,qq,LL]^2
       
       }  # end for LL

     } ## end for qq        

  } ## end for q
     
}  # end for mut
       
}  ## end for pop
}  ## end for ped


for (LL in 1:length(NLoci))  {

y=cbind(RepData,NLoci[LL],2/JackPhi[,,LL], 2/JackPhiSeparate[,,LL],FullR[,,LL],FullRSeparate[,,LL])

write.csv(y,file=paste(pathout,beginning,Ne,"Chr",Chr,"L",NLoci[LL],".csv",sep=""))

}  ## end for LL

## starting in column G, output file has 4 sets of columns, with each set having a column for each sample size
## set 1: nprime for all pairwise comparisons
## set 2: nprime for comparisons on different chromosomes
## set 3: mean r^2 for all pairwise comparisons
## set 4: mean r^2 for comparisons on different chromosomes