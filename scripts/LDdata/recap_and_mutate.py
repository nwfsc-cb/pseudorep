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
import tszip

sys.path.append(os.path.abspath('../python'))
import simNe_v5_functions


NPED = 12
NMUT = 12
NCHROMS = [1,4,16,64]
NES = [50, 200, 800, 3200]


rec_rate = 1e-8


# # Recapitate
for NCHROM in NCHROMS:
    #genome_map = simNe_v5_functions.make_genome_map(NCHROM * 50e6, NCHROM)
    for NE in NES:
        for PED in range(NPED):
            ped_ts_path = f"/home/ryan/simNe/paper/sim/slim/pedigree/Ne_{NE}.PED_{PED}.CHR_{NCHROM}.NPOP_4.pedigree.ts"
            slim_ts = tskit.load(ped_ts_path)
            recap = msprime.simulate(from_ts = slim_ts, recombination_rate=rec_rate, Ne = NE, 
                    population_configurations = [msprime.PopulationConfiguration(),
                                msprime.PopulationConfiguration(),
                                msprime.PopulationConfiguration(),
                                msprime.PopulationConfiguration(),
                                msprime.PopulationConfiguration(),]
                    )
            recap = recap.simplify()
            recap_ts_path = f"/home/ryan/simNe/paper/sim/slim/recap/Ne_{NE}.PED_{PED}.CHR_{NCHROM}.NPOP_4.recap.ts"
            recap.dump(recap_ts_path)



NPED * len(NCHROMS) * len(NES)


# # Mutate

min_sites = 500000
mut_rate = 1e-8

for NCHROM in NCHROMS:
    for NE in NES:
        for PED in range(NPED):
            recap_ts_path = f"/home/ryan/simNe/paper/sim/slim/recap/Ne_{NE}.PED_{PED}.CHR_{NCHROM}.NPOP_4.recap.ts"
            recap = tskit.load(recap_ts_path)
            for MUT in range(NMUT):
                ts = recap
                nsites = 0
                # keep mutating until we have at least min_sites
                while nsites < min_sites:
                    ts = msprime.mutate(ts, rate=mut_rate, random_seed=None, keep=True)
                    nsites = ts.num_sites
                # retain exactly min_sites
                ts = simNe_v5_functions.downsample_snps(ts, min_sites, fail = True)
                # compress as small as possible
                ts = ts.simplify(reduce_to_site_topology=True)
                mut_ts_path = f"/home/ryan/simNe/paper/sim/slim/mut/Ne_{NE}.PED_{PED}.CHR_{NCHROM}.NPOP_4.MUT_{MUT}.tsz"
                tszip.compress(ts, mut_ts_path)


NPED * len(NCHROMS) * len(NES) * NMUT


# ## Split out each separate pop

for NCHROM in NCHROMS:
    for NE in NES:    
        for PED in range(NPED):
            for MUT in range(NMUT):
                mut_ts_path = f"/home/ryan/simNe/paper/sim/slim/mut/Ne_{NE}.PED_{PED}.CHR_{NCHROM}.NPOP_4.MUT_{MUT}.ts"
                ascertained_1 = tskit.load(mut_ts_path)
                for POPid in [1,2,3,4]:
                    single_pop_ts = ascertained_1.simplify(samples = ascertained_1.samples(population_id = POPid), 
                                                   filter_populations = True, # remove other pops (default)
                                                   filter_individuals = True, # remove other inds (default)
                                                   filter_sites = False, # will keep all sites in each pop (not default)
                                                  )
                    single_ts_path = f"/home/ryan/simNe/paper/sim/slim/singlepop/Ne_{NE}.PED_{PED}.CHR_{NCHROM}.NPOP_4.MUT_{MUT}.SINGLEPOP_{POPid}.ts"
                    single_pop_ts.dump(single_ts_path)

