This directory contains three types of files.

`{CHR}_chr.slim` this contains SLiM code for simualtions with {CHR} number of chromosomes.  This is a forward simualtion of one ancestral popualtion splitting into four daughter populations, and produces and tree-seqeunce file as output. 

`sim.slim.Ne_{Ne}.sh` this is a shell script invoking the SLiM scripts for different values of {Ne}.

`sim.slim.Ne_{Ne}.log` this is a log of the SLiM output from the forward simualtions, and contains the random seeds used for simulation.
