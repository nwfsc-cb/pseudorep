## bigld <- read.csv(file="c:/users/robin.waples/dropbox/pseudo-replication/SLiM/SLiMnprimeForModeling.csv", stringsAsFactors=FALSE)
bigld <- read.csv(file="SLiMnprimeForModeling2.csv", stringsAsFactors=FALSE)

 ### drop data for Chr = 1
 bigld = subset(bigld,bigld$Chr > 1)

 ### drop data for all comparisons
 bigld2 = cbind(bigld[,1:4],bigld[,8:10])

  ## restrict to target sample sizes
   smallld = subset(bigld2, bigld2$S ==25 | bigld2$S ==50 | bigld2$S ==100)
 
## restrict to S < Ne
   smallld2 = subset(smallld, smallld$S < smallld$Ne)

 ldtrim = as.data.frame(smallld2)
 colnames(ldtrim) = c("Chr", "Ne","S", "Loci","n","nprime", "nprimen")
 
 #### ensure nprime does not exceed n
# for (i in 1:dim(ldtrim)[[1]])  {
# ldtrim$nprime[i] = min(ldtrim$nprime[i],ldtrim$n[i])
# }
     
ld = ldtrim

## remove NA for nprime
ld = ld[!is.na(ld$nprime),]

 
 
 ##########################################

 ## Martin function using Y = npairs rather than Loci
 
 NeX = c(50,200,800,3200)
 ChrX= c(4,16,64)
 SX = c(25,50,100)
 nsets = (length(NeX)-1)*length(ChrX)*length(SX) + length(ChrX)
  
 ###########
 ff <- function(x,p){
   x/(1/p[1]^p[3]+(x/p[2])^p[3])^(1/p[3])
}

ffFit <- function(p,x,y){
  sum((ff(x,p)-y)^2)
}
############

 nlines = dim(ld)[[1]]
 
 fits = matrix(NA,nsets,6)
 colnames(fits) = c("Chr","Ne","S","p1","p2","p3")
 
 predict4 =seq(1:nlines)
 
 ### fit 3 parameters to each set of Ne*Chr*S nprime values
 
 counter = 0
 
 for (a in 1:length(NeX))   {
 for (b in 1:length(ChrX))  {
 for (c in 1:length(SX))    {
    temp = subset(ld,ld$Ne==NeX[a] & ld$Chr == ChrX[b] & ld$S == SX[c])
    if(nrow(temp) > 0)  {
      counter = counter + 1
      L = temp$Loci
      N = temp$nprime
      m1 <- nlm(ffFit,p=c(log10(N[1]-L[1]),max(log10(N)),2),x=log10(L),y=log10(N))
      p = m1$estimate
      fits[counter,1] = temp[1,1]
      fits[counter,2] = temp[1,2]
      fits[counter,3] = temp[1,3]
      fits[counter,4] = p[1]
      fits[counter,5] = p[2]
      fits[counter,6] = p[3]   }  ## end if
   }}}  ## end for a,b,c
   

## find best overall covariates for each of 3 parameters
IChr = 1/fits[,1]
INe = 1/fits[,2]
IS = 1/fits[,3]
SNe = fits[,3]/fits[,2]
z1 = lm(fits[,4] ~ log(fits[,1]) + log(fits[,2]) + log(fits[,3]) )
z1a = lm(fits[,4] ~ IChr + INe + IS)
summary(z1)
summary(z1a)
z1b = lm(fits[,4] ~ log(fits[,1]) * log(fits[,2]) + log(fits[,3]) )
summary(z1b)
z1c = lm(fits[,4] ~ log(fits[,1]) * log(fits[,2]) + log(fits[,3]) + SNe )
summary(z1c)
z1d = lm(fits[,4] ~ log(fits[,1]) * log(fits[,2]) + log(fits[,3])* log(fits[,2]) )
summary(z1d)
z1e = lm(fits[,4] ~ log(fits[,1]) * log(fits[,2]) + log(fits[,3])* log(fits[,2]) + log(fits[,3])* log(fits[,1]))
summary(z1e)

## Get AIC for various models
aicmat <- AIC(z1,z1a,z1b,z1c,z1d,z1e)
aicmat <- aicmat[order(aicmat[,2]), ]  #sort by AIC
aicmat   

z2 = lm(fits[,5] ~ log(fits[,1]) * log(fits[,2]) + log(fits[,3])* log(fits[,2]) )
summary(z2)
z2a = lm(fits[,5] ~ log(fits[,1]) * log(fits[,2]) + log(fits[,3])* log(fits[,2]) + log(fits[,3])* log(fits[,1]))
summary(z2a)
z2b = lm(fits[,5] ~ log(fits[,1]) * log(fits[,2]) )
summary(z2b)
z2c = lm(fits[,5] ~ log(fits[,1]) * log(fits[,2]) +  log(fits[,3]) )
summary(z2c)

aicmat <- AIC(z2,z2a,z2b,z2c)
aicmat <- aicmat[order(aicmat[,2]), ]  #sort by AIC
aicmat   


z3 = lm(fits[,6] ~ log(fits[,1]) * log(fits[,2]) + log(fits[,3])* log(fits[,2]) )
summary(z3)
z3a = lm(fits[,6] ~ log(fits[,1]) * log(fits[,2]) + log(fits[,3])* log(fits[,2]) + log(fits[,3])* log(fits[,1]))
summary(z3a)
z3b = lm(fits[,6] ~ log(fits[,1]) * log(fits[,2]) )
summary(z3b)
z3c = lm(fits[,6] ~ log(fits[,1]) * log(fits[,2]) +  log(fits[,3]) )
summary(z3c)
z3d = lm(fits[,6] ~ IChr*INe + IChr*IS + SNe)
summary(z3d)

aicmat <- AIC(z3,z3a,z3b,z3c,z3d)
aicmat <- aicmat[order(aicmat[,2]), ]  #sort by AIC
aicmat   

q1 = z1b$coefficients
q2 = z2b$coefficients
q5 = z3$coefficients

######check modeling prediction to full fits
#Efits = matrix(NA,nrow(fits),3)
#for (j in 1:nrow(fits))  {
#N=fits[j,2]
#C=fits[j,1]
#SS=fits[j,3]
#p[1] = q1[1] + q1[2]*log(C) + q1[3]*log(N) + q1[4]*log(SS) + q1[5]*log(C)*log(N)  
#p[2] = q2[1] + q2[2]*log(C) + q2[3]*log(N) + q2[4]*log(C)*log(N) 
#p[3] = q5[1] + q5[2]*log(C) + q5[3]*log(N) + q5[4]*log(SS) + q5[5]*log(C)*log(N) + q5[6]*log(SS)*log(N)  
#Efits[j,] = c(p[1],p[2],p[3])
#} 

#ObsExp = fits[,4:6]/Efits
#median(ObsExp)

### predict p parameters and nprime, given Ne , Chr, and S
#for (j in 1:nlines)  {
#N=ld[j,2]
#C=ld[j,1]
#SS=ld[j,3]
#p[1] = q1[1] + q1[2]*log(C) + q1[3]*log(N) + q1[4]*log(SS) + q1[5]*log(C)*log(N)  
#p[2] = q2[1] + q2[2]*log(C) + q2[3]*log(N) + q2[4]*log(C)*log(N) 
#p[3] = q5[1] + q5[2]*log(C) + q5[3]*log(N) + q5[4]*log(SS) + q5[5]*log(C)*log(N) + q5[6]*log(SS)*log(N)  
# A = ff(log10(ld[j,4]),p)
#predict4[j] = 10^A
#} 

#### ensure predicted nprime does not exceed n
#for (i in 1:dim(ld)[[1]])  {
#predict4[i] = min(predict4[i],ld$n[i])
# }
 
#OE4 = predict4/ld$nprime
#c(max(OE4),min(OE4),median(OE4),var(OE4))

