#### This code predicts nprime for Fst as a function of covariates [Ne,Chr,S,NLoci]
#### User can supply vectors for Ne and NLoci to capture uncertainty or to explore experimental design options
### nprime is the effective number of diallelic loci

######### User-supplied input data
Ne = c(100,200,400,800,1600,3200,10000)  ### should span the plausible range of true Ne
   ### Note: Ne > 3200 is outside the range of values we evaluated; use results with caution
Chr = 18.5 ## equivalent number of chromosomes of length 50 Mb
S = 50  ## diploid number of individuals; use harmonic mean if it varies across loci; should be <Ne
NLoci = c(100,250,500,1000,2500,5000,10000,25000,50000,200000)  ### a range of interest to the user
##############

Inputs = 1:4
OutPut = matrix(NA,length(NLoci),length(Ne))
rownames(OutPut) = NLoci
colnames(OutPut) = Ne

###########
GetNprime <- function(Inputs)     { 

######
 ff <- function(x,p){
   x/(1/p[1]^p[3]+(x/p[2])^p[3])^(1/p[3])}
######
N = Inputs[1]  
C = Inputs[2] 
SS= Inputs[3]  
Loci= Inputs[4] 
p = 1:3

q1 = c(0.958633285 ,0.013469143 , 0.007135748, -0.001863156) 
q2 = c(2.08304896, -0.24271278, 0.00923448, -0.57920739, 1.21567379, 0.05929167 , 0.12727216, 0.09241921) 
q5 = c(-8.8963881, 2.7271474, 1.8244296, 3.8598313, -0.1480143, -0.4126190, -0.3886260 )

p[1] = q1[1] + q1[2]*log(C) + q1[3]*log(N) + q1[4]*log(C)*log(N)  
p[2] = q2[1] + q2[2]*log(C) + q2[3]*log(N) + q2[4]*log(SS) + q2[5]*SS/N + q2[6]*log(C)*log(N) + q2[7]*log(SS)*log(N) + q2[8]*log(SS)*log(C)
p[3] = q5[1] + q5[2]*log(C) + q5[3]*log(N) + q5[4]*log(SS) + q5[5]*log(C)*log(N) + q5[6]*log(SS)*log(N) + q5[7]*log(C)*log(SS)  
A = ff(log10(Loci),p)
Enprime = 10^A

#### ensure predicted nprime does not exceed number of loci
Enprime = min(Enprime, Loci)

return(Enprime)
 } # end function  
############

for (j in 1:length(NLoci))  {
for (k in 1:length(Ne))    {

Inputs[1] = Ne[k]
Inputs[2] = Chr
Inputs[3] = S
Inputs[4] = NLoci[j]

OutPut[j,k] = GetNprime(Inputs)

}}  ## end j,k

## rows are numbers of loci; columns are Ne
round(OutPut)
