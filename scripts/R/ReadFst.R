path = "c:/users/lukes/dropbox/pseudo-replication/SLiM/Final Files/"
filetoread = paste(path,"FstNe200.csv",sep="")
BigFstData = read.csv(filetoread)

P1 = subset(BigFstData,BigFstData$POP == "1_2")
P1[,4] = 1
P2 = subset(BigFstData,BigFstData$POP == "1_3")
P2[,4] = 2
P3 = subset(BigFstData,BigFstData$POP == "1_4")
P3[,4] = 3
P4 = subset(BigFstData,BigFstData$POP == "2_3")
P4[,4] = 4
P5 = subset(BigFstData,BigFstData$POP == "2_4")
P5[,4] = 5
P6 = subset(BigFstData,BigFstData$POP == "3_4")
P6[,4] = 6

PP = rbind(P1,P2,P3,P4,P5,P6)

FstData=PP

Ne = 200
ChrX = c(1,4,16,64)
Sampsize = c(25,50,100)
nped = 4
npop = 6
nmut = 6
nsamps = 8
nreps = 10  ## number of random permutations of samps x muts
R = min(nmut,nsamps)
NLoci = c(100,250,500,1000,2500,5000,10000,25000,50000,100000,200000)
Ascer = c(0,1,2)

nsets = nped*npop*nreps

numrows = length(Ascer)*length(ChrX)*length(Sampsize)*length(NLoci)
OutPutData = matrix(NA,numrows,9)
colnames(OutPutData) = c("Chr","Ne","Asc","Sampsize","NLoci","Fst-U","nprime-U","Fst-W","nprime-W")

bigcount = 0

for (A in 1:3)  {
for (Chrom in 1:length(ChrX))  {
for (SS in 1:length(Sampsize))  {
for (LL in 1:length(NLoci))  {

  print(paste0("Working on ",A,Chrom,SS,LL))
  flush.console() 
  
bigcount = bigcount + 1
mydata = subset(FstData,FstData$ASC == Ascer[A] & FstData$CHR == ChrX[Chrom] & FstData$S == Sampsize[SS] & FstData$LOCI == NLoci[LL])

tempU = matrix(NA,nsets,2)
colnames(tempU) = c("mean","phi")
tempW = tempU

 count = 0
  for (ped in 1:nped)  {
  for (pop in 1:npop)  {
  littledata = subset(mydata,mydata$PED == (ped-1) & mydata$POP == pop)

  for (REP in 1:nreps)  {
    count = count+1
   pick1 = sample(nsamps,R,replace=F)
   pick2 = sample(nmut,R,replace=F)
    AA = 1:nmut
    BB = AA
    for (j in 1:R)  {
    X = subset(littledata, littledata$SAMP == (pick1[j]-1) & littledata$MUT == (pick2[j]-1))
    AA[j] = X$GSTU
    BB[j] = X$GSTW
    }  ## end for j
     tempU[count,1] = mean(AA)
     tempU[count,2] = var(AA)/mean(AA)^2
     tempW[count,1] = mean(BB)
     tempW[count,2] = var(BB)/mean(BB)^2
  } ## end for REP
  
    }}  ## end for Ped,Pop
  
OutPutData[bigcount,1] = ChrX[Chrom]
OutPutData[bigcount,2] = Ne
OutPutData[bigcount,3] = Ascer[A]
OutPutData[bigcount,4] = Sampsize[SS]
OutPutData[bigcount,5] = NLoci[LL]
OutPutData[bigcount,6] = mean(tempU[,1])
OutPutData[bigcount,7] = 2/mean(tempU[,2])
OutPutData[bigcount,8] = mean(tempW[,1])
OutPutData[bigcount,9] = 2/mean(tempW[,2])


}}}}  # end for A,Chr,S,L

##OutPutData

write.csv(OutPutData, file=paste(path,"FstNe200Out.csv",sep=""))
