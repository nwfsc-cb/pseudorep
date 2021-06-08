### This code reads in simulated genotypes and computes nprime for LD analyses
### it expects genotype files with the format and naming convention of the file 
### 'Ne_200.PED_0.CHR_64.NPOP_4.MUT_1.SINGLEPOP_1.ts.geno'
### for each genotype file, it also reads a second file with the same name and extension 'ts.chrom_counts'

path = "mypath1"

NPeds = 4
NPops = 4
NMut = 8
NSamps = 4
Ne = 200
Chr = 64  ## number of chromosomes

beginning = "Z-SLiM"
mapend = ".ts.chrom_counts"
end = ".ts.geno"

Sampsize = c(50,25)
NLoci = c(5000,3000,2000,1000,500,250,100)

BigPhi = list()
for (SS in 1:length(Sampsize))  {
BigPhi[[SS]] = matrix(NA,length(NLoci),NPops*NPeds) 
rownames(BigPhi[[SS]]) = NLoci
}
BigPhiSeparate = BigPhi


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

counter = 0
for (ped in 1:NPeds)  {
for (pop in 1:NPops)   {

counter = counter + 1

### Set up r^2 arrays 

PRSQ = list()
for (SS in 1:length(Sampsize))  {
PRSQ[[SS]] = array(NA, dim = c(NSamps,NMut,length(NLoci)))
} 
PRSQSeparate = PRSQ

for (mut in 1:NMut)   {

 print(paste0("                  Population = ",ped,pop,mut))
  flush.console() 

filetoread = paste("Ne_",Ne,".PED_",(ped-1),".CHR_",Chr,".NPOP_4.MUT_",(mut-1),".SINGLEPOP_",pop,end,sep="")
filetoread2 = paste("Ne_",Ne,".PED_",(ped-1),".CHR_",Chr,".NPOP_4.MUT_",(mut-1),".SINGLEPOP_",pop,mapend,sep="")

z = read.table(file=paste(path,filetoread,sep=""))

maps = read.table(file=paste(path,filetoread2,sep=""))

cvec = 999
for (j in 1:Chr)  {
a = rep(j,maps[j,1])
cvec = c(cvec,a)
}
cvec = cvec[-1]   ## vector giving chromosome for each locus

z = rbind(cvec,z)

z = z[,sample(1:dim(z)[[2]])]  ## randomize order of loci

z = as.matrix(z)
dimnames(z) = NULL

cvec2 = z[1,]   ### sorted vector of chromosome data   


for (q in 1:NSamps)  {	  

for (SS in 1:length(Sampsize))  {

         print(paste0("Working on samplesize",Sampsize[SS]))
	  flush.console()      

      ## Draw a random sample from the Ne individuals in final generation
      z2 = z[-1,]  ## first remove first row of chromosome data
      z2 = z2[sample(nrow(z2),Sampsize[SS],replace=F),] 
      P = colMeans(z2)/2 
            
      ##eliminate loci monomorphic or nearly so  
      z2 = rbind(cvec2,z2)  ## replace chr vector for screening loci
      cut = 1/Sampsize[SS]
      z2 = z2[, P > cut & P < (1-cut), drop = F] 
        
      cvec3 = z2[1,]   ### save sorted and trimmed vector of loci that will be used
      z2 = z2[-1,]    
        
           start = (q-1)*NLoci[1] + 1
           stop = start + NLoci[1] - 1
           BigRsq = GetLD(z2[,start:stop])
            ### For Separate matrix, put NA for any comparisons of loci on the same chromosome
            BigRSeparate = BigRsq         
            c3 = cvec3[start:stop]
              for (L1 in 1:(NLoci[1]-1))   {
              for (L2 in (L1+1):NLoci[1])  {
              if(c3[L1]==c3[L2]) { BigRSeparate[L1,L2] = NA  }
              }}  ## end for L1, L2   
                 
            for (LL in 1:length(NLoci))  {
              PRSQ[[SS]][q,mut,LL] = mean(BigRsq[1:NLoci[LL],1:NLoci[LL]],na.rm=T)
              PRSQSeparate[[SS]][q,mut,LL] = mean(BigRSeparate[1:NLoci[LL],1:NLoci[LL]],na.rm=T)
              }   ## end for LL        
            
 }  ## end for SS                
                    
}  # end for q	        

}  ## end for mut


  for (SS in 1:length(Sampsize))  {
  for (LL in 1:length(NLoci))  {
       A = PRSQ[[SS]][,,LL]
       B = PRSQSeparate[[SS]][,,LL]
       A1 = as.vector(A)
       B1 = as.vector(B)
       BigPhi[[SS]][LL,counter] = var(A1)/mean(A1)^2
       BigPhiSeparate[[SS]][LL,counter] = var(B1)/mean(B1)^2
  }  ## end for LL
  }  ## end for SS     
   
      for (LL in 1:length(NLoci))  {
        C = matrix(NA,NMut*NSamps,length(Sampsize))
        C1 = C
          for (SS in 1:length(Sampsize))  {  
           A = PRSQ[[SS]][,,LL]
           B = PRSQSeparate[[SS]][,,LL]
          A1 = as.vector(A)
          B1 = as.vector(B)
          C[,SS] = A1
          C1[,SS] = B1
          }  # end for SS
         C2 = cbind(C,C1) 
       outfile = paste(pathout,beginning,"-Ne",Ne,"-Chr",Chr,"-Loci",NLoci[LL],".out",sep="") 
       write.table(C2, file = outfile,append = TRUE,quote=FALSE,col.names=FALSE)
     }  ## end for LL        
          
}   ## end for pop

}  ## end for ped

FinalPhi = matrix(NA,length(NLoci),length(Sampsize))
rownames(FinalPhi) = NLoci
colnames(FinalPhi) = Sampsize
FinalPhiSeparate = FinalPhi

for (SS in 1:length(Sampsize))  {
 FinalPhi[,SS] = rowMeans(BigPhi[[SS]])
 FinalPhiSeparate[,SS] = rowMeans(BigPhiSeparate[[SS]])
  }

2/FinalPhi
2/FinalPhiSeparate

