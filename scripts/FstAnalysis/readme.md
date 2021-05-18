## `analyze_FST.py` 
estimates Fst on subsets of individuals and loci from pairs of popualtions, used for the Fst analyses. 
These are made by loading the recapitated SLiM output produced by `recap_and_mutate.py`.  Then pairs of populations, and individuals are selected, mutations are placed on the tree-sequences and an ascertainment procedure is conducted. The output are summary statistics including Fst. 

## `analyze_FST_blockjackknife.py` 
The same basic procedure as in `analyze_FST.py`, but in addition two block jackknife procedures are performed to calculated std errors, with chromosome and 5mb-sized blocks. 

## `analyze_FST.example.sh`
This shell script gives an example how to call a`nalyze_FST.py` and `analyze_FST_blockjackknife.py`.
