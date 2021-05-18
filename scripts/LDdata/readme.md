## `recap_and_mutate.py` 
produces *.SINGLEPOP_{POPid}.ts files used for the LD analyses. These are made by loading SLiM output, recapitating, mutating and saving the output.  Next each daughter population is selected in turn. 

## `make_LD_data .py` 
produces the *.geno files used in the LD analyses. These files are made by loading the output of `recap_and_mutate.py`, and filtering by minor allele count.

