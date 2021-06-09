###for all analyses

bigFst <- read.csv(file="Fst/FstnprimeFinal2.csv", stringsAsFactors=FALSE)

Fsttrim = as.data.frame(bigFst)
colnames(Fsttrim) = c("Chr", "Ne","Asc","S", "Loci","nprime")
 
    
Fst = Fsttrim

## remove NA for nprime
Fst = Fst[!is.na(Fst$nprime),]


############### restrict to Loci > minLoci

minLoci = 100

Fstshort = subset(Fst,Fst$Loci>=minLoci)
 
##########################################


#########################################################################

Fst = Fstshort

 ## Martin functions
 
 NeX = c(50,200,800,3200)
 ChrX= c(1,4,16,64)
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

 nlines = dim(Fst)[[1]]
 
 fits = matrix(NA,nsets,6)
 colnames(fits) = c("Chr","Ne","S","p1","p2","p3")
 
 predictF =seq(1:nlines)
 
  ### fit 3 parameters to each set of Ne*Chr*S nprime values
  
  counter = 0
  
  for (a in 1:length(NeX))   {
  for (b in 1:length(ChrX))  {
  for (c in 1:length(SX))    {
     temp = subset(Fst,Fst$Ne==NeX[a] & Fst$Chr == ChrX[b] & Fst$S == SX[c])
     if(nrow(temp) > 0)  {
       counter = counter + 1
       L = temp$Loci
       N = temp$nprime
       m1 <- nlm(ffFit,p=c(log10(20),max(log10(N)),2),x=log10(L),y=log10(N))
       p = m1$estimate
       fits[counter,1] = temp[1,1]
       fits[counter,2] = temp[1,2]
       fits[counter,3] = temp[1,4]
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
z1f = lm(fits[,4] ~ log(fits[,1]) * log(fits[,2]))
summary(z1f)

## Get AIC for various models
aicmat <- AIC(z1,z1a,z1b,z1c,z1d,z1e,z1f)
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
z2d = lm(fits[,5] ~ log(fits[,1]) * log(fits[,2]) + log(fits[,3])* log(fits[,2]) + SNe )
summary(z2d)
z2e = lm(fits[,5] ~ log(fits[,1]) * log(fits[,2]) + log(fits[,3])* log(fits[,2]) + SNe + IChr + INe + IS)
summary(z2e)
z2f = lm(fits[,5] ~ log(fits[,1]) * log(fits[,2]) + log(fits[,3])* log(fits[,2]) + SNe + IChr )
summary(z2f)
z2g = lm(fits[,5] ~ log(fits[,1]) * log(fits[,2]) + log(fits[,3])* log(fits[,2]) + log(fits[,3])* log(fits[,1])+ SNe )
summary(z2g)

aicmat <- AIC(z2,z2a,z2b,z2c,z2d,z2e,z2f,z2g)
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
z3e = lm(fits[,6] ~ log(fits[,1]) * log(fits[,2]) + log(fits[,3])* log(fits[,2]) + log(fits[,3])* log(fits[,1]) + SNe + IChr)
summary(z3e)
z3f = lm(fits[,6] ~ log(fits[,1]) + log(fits[,2]) + log(fits[,3]))
summary(z3f)

aicmat <- AIC(z3,z3a,z3b,z3c,z3d,z3e)
aicmat <- aicmat[order(aicmat[,2]), ]  #sort by AIC
aicmat   

# get coefficients from each of the best models
q1 = z1f$coefficients
q2 = z2g$coefficients
q5 = z3a$coefficients

######check modeling prediction to full fits
Efits = matrix(NA,nrow(fits),3)
for (j in 1:nrow(fits))  {
N=fits[j,2]
C=fits[j,1]
SS=fits[j,3]
p[1] = q1[1] + q1[2]*log(C) + q1[3]*log(N) + q1[4]*log(C)*log(N)  
p[2] = q2[1] + q2[2]*log(C) + q2[3]*log(N) + q2[4]*log(SS) + q2[5]*SS/N + q2[6]*log(C)*log(N) + q2[7]*log(SS)*log(N) + q2[8]*log(SS)*log(C)
p[3] = q5[1] + q5[2]*log(C) + q5[3]*log(N) + q5[4]*log(SS) + q5[5]*log(C)*log(N) + q5[6]*log(SS)*log(N) + q5[7]*log(C)*log(SS)  
Efits[j,] = c(p[1],p[2],p[3])
} 

ObsExp = fits[,4:6]/Efits
median(ObsExp)

predict4 = rep(NA,nrow(Fst))
### predict p parameters and nprime, given Ne , Chr, and S
for (j in 1:nlines)  {
N=Fst[j,2]
C=Fst[j,1]
SS=Fst[j,4]
p[1] = q1[1] + q1[2]*log(C) + q1[3]*log(N) + q1[4]*log(C)*log(N)  
p[2] = q2[1] + q2[2]*log(C) + q2[3]*log(N) + q2[4]*log(SS) + q2[5]*SS/N + q2[6]*log(C)*log(N) + q2[7]*log(SS)*log(N) + q2[8]*log(SS)*log(C)
p[3] = q5[1] + q5[2]*log(C) + q5[3]*log(N) + q5[4]*log(SS) + q5[5]*log(C)*log(N) + q5[6]*log(SS)*log(N) + q5[7]*log(C)*log(SS)  
 A = ff(log10(Fst[j,5]),p)
predict4[j] = 10^A
} 

# #### ensure predicted nprime does not exceed number of loci
# for (i in 1:dim(Fst)[[1]])  {
# predict4[i] = min(predict4[i],Fst$Loci[i])
#  }
#  
# OE4 = predict4/Fst$nprime
# c(max(OE4),min(OE4),median(OE4),var(OE4))
# 	
# ### now delete Chr = 1 from comparisons of predicted and true
# Ratio = predict4/Fst$nprime
# B = cbind(Fst[,1:6],predict4,Ratio)
# B = subset(B,B$Chr>1)
# c(max(B$Ratio),min(B$Ratio),median(B$Ratio),var(B$Ratio))
# plot(log10(B$predict4),log10(B$nprime))
# cor(log10(B$predict4),log10(B$nprime))
# 
# jpeg("Fst/Results/ObsExp-Fstnprime.jpg",width = 500, height = 400)
# plot(log10(B$nprime),log10(B$predict4),xlab = "log10(true nprime)",ylab = "log10(fitted nprime)" )
# dev.off()

