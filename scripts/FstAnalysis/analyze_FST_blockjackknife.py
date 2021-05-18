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
import tskit
import allel


os.environ['MKL_NUM_THREADS'] = '8'


sys.path.append(os.path.abspath('../python'))
import simNe_v5_functions


# Parsing and validating  cmd line arguements
narg = len(sys.argv) - 1
assert narg == 3

NE = int(sys.argv[1])
NCHROM = int(sys.argv[2])
PED = int(sys.argv[3])

assert (NE in [200, 800])
assert (NCHROM in [4, 16])
assert (PED in [0,1,2,3,4,5,6,7,8,9,10,11])

def get_genotype_array(ts):
    l = ts.num_sites
    n = ts.num_samples
    ga = allel.GenotypeArray(ts.genotype_matrix().reshape(l, int(n/2), 2), dtype='i1')
    return(ga)

# combines the haploid genotypes into diploid genotypes
def convert_gt(gt):
    return(gt[:,::2] + gt[:,1::2])

# fast method to calculate r2
def matmul_r2(snp_genos):
    norm = (snp_genos - snp_genos.mean(1, keepdims=True))/snp_genos.std(1, ddof=1, keepdims=True)
    return(((norm @ norm.T) /(norm.shape[1]-1))**2)

def non_diag_mean(mat):
    dim = mat.shape[0]
    return((mat.sum()-dim)/((dim*dim)-dim))

def get_fst(ts):

    # shared by later Fst calculations
    ga = get_genotype_array(ts)
    pop1_indices = (ts.samples(population=1)[::2]/2).astype(int)
    pop2_indices = (ts.samples(population=2)[::2]/2).astype(int)
    subpops = [pop1_indices, pop2_indices]
    ac1 = ga[:, pop1_indices].count_alleles()
    ac2 = ga[:, pop2_indices].count_alleles()

    # Nei Gst
    a,b,u = get_nei_gst(ts, ac1, ac2)
    nei_gst_unweighted = u
    nei_gst_weighted = np.sum(a)/np.sum(b)
    nei_gst_variance = np.nanvar(a/b)

    # Nei Fst (from Bhatia)
    a,b = get_nei_fst(ts, ac1, ac2)
    nei_fst_weighted = np.sum(a)/np.sum(b)
    nei_fst_unweighted = np.nanmean(a/b)

    # Hudson Fst
    a, b = get_hudson_fst(ac1, ac2)
    hudson_fst_weighted = np.sum(a)/np.sum(b)
    hudson_fst_unweighted = np.nanmean(a/b)
    hudson_fst_variance = np.nanvar(a/b)

    # WC Fst
    a,b,c = get_wc_fst(ga, subpops)
    wc_fst_weighted = np.sum(a)/(np.sum(a)+np.sum(b)+np.sum(c))
    wc_fst_unweighted = np.nanmean(a/(a+b+c))

    n = int(ts.num_samples/2)
    sample_sets = [range(0,n),range(n,n+n)]
    fst_branch_norm = np.float(ts.Fst(sample_sets = sample_sets, mode = 'branch', windows = None, span_normalise = True))
    fst_site_norm = np.float(ts.Fst(sample_sets = sample_sets, mode = 'site', windows = None, span_normalise = True))

    fst_branch_nonnorm = np.nanmean(ts.Fst(sample_sets = sample_sets, mode = 'branch', windows = 'trees', span_normalise = False))
    fst_site_nonnorm = np.nanmean(ts.Fst(sample_sets = sample_sets, mode = 'site', windows = 'sites', span_normalise = False))

    res = {
        'nei_gst_unweighted': nei_gst_unweighted,
        'nei_gst_weighted': nei_gst_weighted,

        'nei_fst_weighted': nei_fst_weighted,
        'nei_fst_unweighted': nei_fst_unweighted,

        'hudson_fst_weighted': hudson_fst_weighted,
        'hudson_fst_unweighted': hudson_fst_unweighted,

        'wc_fst_weighted': wc_fst_weighted,
        'wc_fst_unweighted': wc_fst_unweighted,

        'fst_branch_norm': fst_branch_norm,
        'fst_branch_nonnorm': fst_branch_nonnorm,
        'fst_site_norm': fst_site_norm,
        'fst_site_nonnorm': fst_site_nonnorm,
        'hudson_fst_variance':hudson_fst_variance,
        'nei_gst_variance': nei_gst_variance,
        }
    return(res)


def get_wc_fst(ga, subpops):
    a,b,c = allel.weir_cockerham_fst(g = ga, subpops = subpops)
    return(a,b,c)

def get_hudson_fst(ac1, ac2):
    num, denom = allel.hudson_fst(ac1, ac2)
    return(num, denom)

def get_patterson_f2(ac1, ac2):
    num, denom = allel.patterson_fst(ac1, ac2)
    return(num, denom)

def get_nei_gst(ts, ac1, ac2):
    # from robins code
    n1 = len(ts.samples(population=1))
    n2 = len(ts.samples(population=2))
    P1 = ac1[:,1]/n1
    P2 = ac2[:,1]/n2
    sP1 = 1.0 - (P1**2 + (1-P1)**2) # expected heterozygosity in P1
    sP2 = 1.0 - (P2**2 + (1-P2)**2) # expected heterozygosity in P2
    Hexp = (sP1+sP2)/2.0 # mean heterozygosity
    Pbar = (P1 + P2)/2.0 # mean freq
    Htot = 1 - (Pbar**2 + (1-Pbar)**2)
    F = 1.0 - Hexp/Htot
    Fst_u = np.nanmean(F) # unweighted mean over loci
    G_num = Htot - Hexp
    G_denom = Htot
    return(G_num, G_denom, Fst_u)

def get_nei_fst(ts, ac1, ac2):
    # from Bhatia et al, citing Nei 1986
    l = ts.num_sites
    n = ts.num_samples
    n1 = len(ts.samples(population=1))
    n2 = len(ts.samples(population=2))
    P1 = ac1[:,1]/n1
    P2 = ac2[:,1]/n2
    Pavg = (P1+P2)/2
    num = (P1-P2)**2 # squared dif in alele freqs
    denom = 2 * Pavg*(1-Pavg) # exp het for average allele freq
    return(num, denom)


def sample_inds_fst(Ne, nsamp):
    """selected the samples relating to nsamp diploid individuals from each population"""
    inds = np.array(range(int(Ne)))
    inds_chosen1 = np.random.choice(inds, nsamp, replace = False)
    inds_chosen2 = np.random.choice(inds, nsamp, replace = False) + int(Ne)
    inds_chosen = np.concatenate([inds_chosen1, inds_chosen2])
    a, b = inds_chosen * 2, inds_chosen * 2+1
    c = np.empty((a.size + b.size), dtype=a.dtype)
    c[0::2] = a
    c[1::2] = b
    samples = sorted(c)
    return(samples)



### Parameters
NMUT = 6
mut_rate = 2e-8
NREP_SAMP = 8
NSAMPLES = [25, 50, 100]
NSNP = [5000,20000,50000]
ASCERTAINED =[2]
fwd_of_Ne = {50:10, 100:20, 200:40, 800:160, 3200:640}


stat_order = ['nei_gst_unweighted','nei_gst_weighted','nei_fst_weighted','nei_fst_unweighted','hudson_fst_weighted','hudson_fst_unweighted','wc_fst_weighted','wc_fst_unweighted','fst_branch_norm','fst_branch_nonnorm','fst_site_norm','fst_site_nonnorm', 'hudson_fst_variance', 'nei_gst_variance']

# Results file
path = f'/home/ryan/simNe/paper/results/fst_jackknife/fst.jackknife2.NE_{NE}.NCHROM_{NCHROM}.PED_{PED}.results.txt'
# write header
header = 'NE\tNCHROM\tPED\tPOPid\tS\tSAMP_REP\tLOCI\tMUT\tASCERTAIN\tJK_BINSIZE\tJK_BIN_EXCLUDE\tBIN_LOCI\tR2' + '\t' +'\t'.join(stat_order) + '\n'
with open(path, 'w') as OUTFILE:
    OUTFILE.write(header)



# Analysis
genome_map = simNe_v5_functions.make_genome_map(NCHROM * 50e6, NCHROM)
recap_ts_path = f"/home/ryan/simNe/paper/sim/slim/recap/Ne_{NE}.PED_{PED}.CHR_{NCHROM}.NPOP_4.recap.ts"
ts = tskit.load(recap_ts_path)


# now sampling a population pair
for poppair in itertools.combinations([1,2,3,4], 2):
    pop_ts = ts.simplify(samples = list(ts.samples(population_id = poppair[0])) +
                         list(ts.samples(population_id = poppair[1])))
    for S in NSAMPLES:
    # take S individuals from each population
        if S > NE:
                # can't sample if S is greater than Ne
                break
        for SAMP_REP in range(NREP_SAMP):
            samples = sample_inds_fst(NE, S)
            ts_samp = pop_ts.simplify(samples = samples)
            for ASCERTAIN in ASCERTAINED:
                for LOCI in NSNP:
                    # mutate
                    nsites = 0
                    ts_analyze = ts_samp.tables.tree_sequence()
                    # keep mutating until we have at least min_sites
                    mut_tries = 0
                    while nsites < (LOCI*NMUT):
                        if ASCERTAIN == 0:
                            ts_analyze = msprime.mutate(ts_analyze, rate=mut_rate,
                                random_seed=None, keep=True, start_time = 0)
                        elif ASCERTAIN == 1:
                            ts_analyze = msprime.mutate(ts_analyze, rate=mut_rate,
                                random_seed=None, keep=True, start_time = fwd_of_Ne[NE])
                        elif ASCERTAIN == 2:
                            ts_analyze = msprime.mutate(ts_analyze, rate=mut_rate,
                                random_seed=None, keep=True, start_time = 0)
                            # only keep sites above 5% overall
                            five_percent = int(np.floor(S*4)*0.05)
                            ts_analyze = simNe_v5_functions.strip_mac(ts_analyze, mac=five_percent)
                        else:
                            assert False
                        nsites = ts_analyze.num_sites
                        #mut_tries += 1
                        #print(f'S: {S} | nmuts: {mut_tries} | nsites: {nsites}')

                    # retain exactly min_sites
                    ts_analyze = simNe_v5_functions.downsample_snps(ts_analyze, (LOCI*NMUT), fail = True)
                    all_sites = np.array([s.id for s in ts_analyze.sites()])
                    locus_sets = np.random.choice(a = all_sites, size = (NMUT, LOCI), replace = False)
                    for MUT in range(NMUT):
                        # remove the rest of the sites
                        to_delete = np.delete(locus_sets, MUT, axis = 0).flatten()
                        ts_analyze2 = ts_analyze.delete_sites(to_delete)

                        # get the position of the remaining sites
                        site_pos = np.array([site.position for site in ts_analyze2.sites()])
                        jk_site_ids = np.array([s.id for s in ts_analyze2.sites()])
                        # assign sites to 5mb/1chr groups

                        for JK_BINS in ['5mb', '1chr']:
                            mean_r2 = 0
                            if JK_BINS == '5mb':
                                bin_ends = np.arange(0, int(NCHROM * 50e6 + 5e6), int(5e6))
                                #jk_groups = pd.cut(site_pos, bins = bins_5mb, labels = False)
                            elif JK_BINS == '1chr':
                                bin_ends = np.arange(0, int(NCHROM * 50e6 + 5e6), int(50e6))
                                #jk_groups = simNe_v5_functions.find_chrom(site_pos[:,None], genome_map)
                            else:
                                assert False, 'incorrect specification of JK_BINS'
                            for ibin, (start, stop) in enumerate(zip(bin_ends, bin_ends[1:])):
                                interval = np.array([start,stop]).reshape(1, -1)
                                ts_analyze_jk = ts_analyze2.delete_intervals(interval)
                                fst_res = get_fst(ts_analyze_jk)
                                bin_sites = ts_analyze_jk.num_sites

                                # write out the results
                                with open(path, 'a') as OUTFILE:
                                    pair_str = '_'.join([str(x) for x in poppair])
                                    OUTFILE.write(f'{NE}\t{NCHROM}\t{PED}\t{pair_str}\t{S}\t{SAMP_REP}\t{LOCI}\t{MUT}\t{ASCERTAIN}')
                                    OUTFILE.write(f'\t{JK_BINS}\t{ibin}\t{bin_sites}\t{mean_r2}')
                                    for stat in stat_order:
                                        OUTFILE.write(f'\t{fst_res[stat]}')
                                    OUTFILE.write('\n')
                        #without jackknife
                        JK_BINS = 'None'
                        ibin = 0
                        fst_res = get_fst(ts_analyze2)
                        bin_sites = ts_analyze2.num_sites

                        # get mean r2 within each pop
                        mean_r2 = 0
                        for popid in [1,2]:
                            r2_ts = ts_analyze2.simplify(samples = ts_analyze2.samples(population_id = popid))
                            r2_ts = simNe_v5_functions.strip_mac(r2_ts, mac=2)
                            r2_ts = simNe_v5_functions.downsample_snps(r2_ts, 5000, fail = False)
                            gt = convert_gt(r2_ts.genotype_matrix())
                            ldmat = matmul_r2(gt)
                            mean_r2 += non_diag_mean(ldmat)
                        # report the average of the two pops
                        mean_r2 = mean_r2 / 2

                        # write out the results
                        with open(path, 'a') as OUTFILE:
                            pair_str = '_'.join([str(x) for x in poppair])
                            OUTFILE.write(f'{NE}\t{NCHROM}\t{PED}\t{pair_str}\t{S}\t{SAMP_REP}\t{LOCI}\t{MUT}\t{ASCERTAIN}')
                            OUTFILE.write(f'\t{JK_BINS}\t{ibin}\t{bin_sites}\t{mean_r2}')
                            for stat in stat_order:
                                OUTFILE.write(f'\t{fst_res[stat]}')
                            OUTFILE.write('\n')
