#!/usr/bin/env python
import os
import sys
import glob
import scipy as sp
import random
import numpy as np
import pandas as pd
import itertools
import msprime
import pyslim
import allel
import tskit
sys.path.append(os.path.abspath('../python'))
import simNe_v5_functions

import matplotlib.pyplot as plt
import seaborn as sns

# # Analysis
NPED = 4
NMUT = 8
NCHROMS = [1,4,16,64]
NES = [50, 200, 800]


# combines the haploid genotypes into diploid genotypes
def convert_gt(gt):
    return(gt[:,::2] + gt[:,1::2])

# ### Just spit out the genotypes
for NCHROM in NCHROMS:
    # make genome map here
    genome_map = simNe_v5_functions.make_genome_map(NCHROM * 50e6, NCHROM)
    for NE in NES:    
        for PED in range(NPED):
            for MUT in range(NMUT):
                for POPid in [1,2,3,4]:
                    single_ts_path = f"/home/ryan/simNe/paper/sim/slim/singlepop/Ne_{NE}.PED_{PED}.CHR_{NCHROM}.NPOP_4.MUT_{MUT}.SINGLEPOP_{POPid}.ts"
                    ts = tskit.load(single_ts_path)
                    ascertained = simNe_v5_functions.strip_mac(ts, mac=3)
                    if ascertained.num_sites > 75000:
                        ascertained = simNe_v5_functions.downsample_snps(ascertained, 75000)
                    site_pos = np.array([site.position for site in ascertained.sites()])[:,None]
                    chrom_counts = np.bincount(simNe_v5_functions.find_chrom(site_pos, genome_map))[1:]
                    gt = convert_gt(ascertained.genotype_matrix())

                    geno_path = f"/home/ryan/simNe/paper/sim/slim/geno/Ne_{NE}.PED_{PED}.CHR_{NCHROM}.NPOP_4.MUT_{MUT}.SINGLEPOP_{POPid}"
                    np.savetxt(geno_path+'.geno', gt, 
                           delimiter = '\t', fmt = '%i')
                    np.savetxt(geno_path+'.chrom_counts', chrom_counts, 
                           delimiter = '\t', fmt = '%i')
                    ts_path = f"/home/ryan/simNe/paper/sim/slim/geno/Ne_{NE}.PED_{PED}.CHR_{NCHROM}.NPOP_4.MUT_{MUT}.SINGLEPOP_{POPid}.ascertained"
                    ascertained.dump(geno_path+'.ts')


