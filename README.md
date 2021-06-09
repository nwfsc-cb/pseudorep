# pseudorep
Repo for Waples et al. "Pseudoreplication in genomics-scale datasets", Molecular Ecology Resources


### Simulated data
Output from our simulations used as input to the regression modeling is found in the data/ foler, https://github.com/nwfsc-cb/pseudorep/tree/main/data

nprimeLD.csv, nprimeFst.csv – these files give the nprime values for LD and Fst, respectively.  These data were used to make many of the tables and figures, and these files were the input data used for developing the models to predict nprime as a function of covariates (Ne, Chr, L, S).

FstNe200.csv – This is an example output file for simulations for Fst with Ne = 200.  Each of the >460K lines gives mean Fst values, averaged over the specified number of loci, for one replicate sample for one pair of populations, for the specified pedigree and mutational replicates.  Both weighted and unweighted Fst values are shown for both Nei and Hudson estimators, under 3 ascertainment schemes.  Most of our analyses used results for weighted Nei Fst and ascertainment = 2 (MAF = 0.05).   

Ne_200.PED_0.CHR_64.NPOP_4.MUT_1.SINGLEPOP_1.ts.geno – This file is an example of simulated genotypes (50K loci for Ne = 200 individuals) that are read in by the file ReadLD.R and used to compute nprime for LD.

Ne_200.PED_0.CHR_64.NPOP_4.MUT_1.SINGLEPOP_1.ts.chrom_counts – This file is a companion to the genotype file above; it gives the number of loci on each of the 64 chromosomes, in order.  This information is used to restrict comparisons to pairs of loci on different chromosomes.


### Scripts 
The /scripts folder contains several folders with code to replicate our analysis. We have included Readme.md files for each folder, which are organized by analysis and programming language, https://github.com/nwfsc-cb/pseudorep/tree/main/scripts

PredictLD.R, PredictFst.R – These files use results of our modeling exercise to predict nprime for any new combinations of the four covariates (Ne, Chr, L, S) specified by the user.

SimUnlinkedLD.R – This code simulates unlinked loci and computes nprime for LD analyses.

SimInfiniteNe.R – This code simulates sampling from a population with infinite Ne and computes nprime for LD analyses.  Although we did not present results for Fst with infinite Ne, that could easily be mimicked with this code by comparing two replicate samples from the same parametric allele frequencies.

TriangleNe.R – This code simulates the kind of data used for Figure S8.

ReadFst.R – This code reads in files with the same format as FstNe200.csv and computes nprime for each scenario.

ReadLD.R – This code reads in genotype files in the format of ‘Ne_200.PED_0.CHR_64.NPOP_4.MUT_1.SINGLEPOP_1.ts.geno’ and computes nprime.  This code also reads in the associated ‘ts.chrom_counts’ files, which allow it to compute nprime when restricting comparisons to pairs of loci on different chromosomes.

JackLD.R – This code reads in the same pairs of files as does ‘ReadLD.R’, generates new samples, and for each sample computes mean r2 and estimated effective degrees of freedom, using the Jones et al. jackknife method implemented in NeEstimator.    

FstFitting.R contains an example of iterating over multiple model forms to identify the most parsimonious model

LDFitting.R contains an example of iterating over multiple model forms to identify the most parsimonious model for the LD analyses
