#### This code predicts nprime for LD as a function of covariates [Ne,Chr,S,NLoci]
#### User can supply vectors for Ne and NLoci to capture uncertainty or to explore experimental design options
### nprime is the effective number of pairwise comparisons of diallelic loci

######### User-supplied input data
Ne = c(100,200,400,800,1600,3200,10000)  ### should span the plausible range of true Ne
Chr = 18.5  ## equivalent number of chromosomes of length 50 Mb
S = 50  ## diploid number of individuals; use harmonic mean if it varies across loci; should be <Ne
NLoci = c(100,500,1000,5000,10000,50000)  ### a range of interest to the user
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

q1 = c(1.53450230,  0.12518487,  0.06692614, -0.06710636,-0.24470293, -0.01474650) 
q2 = c(0.52949689,0.12559264,0.61411490,0.24493530,0.50730040,-1.38178677,-0.01697811,-0.03484065)
q5 =  c(1.2317425, 0.5675758,-0.3671743, 1.2211834, 0.3815077, 0.1881963,-0.6756503 )
p[1] = q1[1] + q1[2]*log(C) + q1[3]*log(N) + q1[4]*log(SS) + q1[5]*SS/N + q1[6]*log(C)*log(N)  
p[2] = q2[1] + q2[2]*log(C) + q2[3]*log(N) + q2[4]*log(SS) + q2[5]*SS/N + q2[6]/C + q2[7]*log(C)*log(N) + q2[8]*log(SS)*log(N)
p[3] = q5[1] + q5[2]*log(C) + q5[3]*log(N) + q5[4]*log(SS) + q5[5]*log(C)*log(N) + q5[6]*log(SS)*log(N) + q5[7]*log(C)*log(SS)  
A = ff(log10(Loci),p)
Enprime = 10^A

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
