##This code estimates Qbar for use in Equation S13
## It assumes Poisson variance in CO number and random placement of COs
NReps = 500000
COs = 0.5  ## mean number of crossovers per chromosome
qbar = rep(1,NReps)
for (q in 1:NReps)  {
  n = rpois(1,COs)
  if (n > 0)  {
   j = runif(n)
   jj = c(0,sort(j),1)
   segments = diff(jj)
    qbar[q]  = sum(segments^2)
   }  # end if
 } # end for q
 mean(qbar)
 
