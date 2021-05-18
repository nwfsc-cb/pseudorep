import numpy as np
import scipy as sp
from numba import jit


def strip_mac(ts, mac):
    """removes sites with <= mac from the tree-sequence"""
    # also get the other threshold
    mac2 = ts.num_samples - mac
    tables = ts.dump_tables()
    tables.sites.clear()
    tables.mutations.clear()
    for tree in ts.trees():
        for site in tree.sites():
            nmut = len(site.mutations)
            assert nmut <= 1  # Only supports infinite sites muts.
            if nmut == 1:
                mut = site.mutations[0]
                if (tree.num_samples(mut.node) > mac) & (tree.num_samples(mut.node) < mac2):
                    site_id = tables.sites.add_row(
                        position=site.position,
                        ancestral_state=site.ancestral_state)
                    tables.mutations.add_row(
                        site=site_id, node=mut.node, derived_state=mut.derived_state)
    return(tables.tree_sequence())
    
    
def downsample_snps(ts, nsnps, fail = True):
    """downsample a tree-sequence to a maximum number of snps"""
    if ts.num_mutations== nsnps:
        return(ts)
    elif ts.num_mutations < nsnps:
        if fail:
            assert False
        else:
            return(ts)
    else:
        keep = frozenset(np.random.choice(a = ts.num_mutations, size = nsnps, replace = False))
        tables = ts.dump_tables()
        tables.sites.clear()
        tables.mutations.clear()
        mut_ix = 0
        for tree in ts.trees():
            for site in tree.sites():
                assert len(site.mutations) == 1  # Only supports infinite sites muts.
                if mut_ix in keep:
                    mut = site.mutations[0]
                    site_id = tables.sites.add_row(
                        position=site.position,
                        ancestral_state=site.ancestral_state)
                    tables.mutations.add_row(
                        site=site_id, node=mut.node, derived_state=mut.derived_state)
                mut_ix +=1 
        return tables.tree_sequence()
    
def downsample_individuals(ts, max_inds, fail = True):
    ind_array = ts.samples().reshape(-1,2) # which nodes correspond to which individuals
    keep_inds = np.random.choice(a = int(ts.num_samples/2), size = max_inds, replace = False)
    keep_samples = ind_array[sorted(keep_inds)].ravel()
    return(ts.simplify(samples = keep_samples))




def annotate_chromosomes(ts, genome_map):
    """uses the genome map to annotate each site in the ts with its chromosome
    chromosome is saved in the metadata column
    use int.from_bytes(site.metadata, byteorder = 'little') to recover it
    """
    tables = ts.dump_tables()
    tables.sites.clear()
    site_pos = np.array([site.position for site in ts.sites()])[:,None]
    chroms = find_chrom(site_pos, genome_map)
    
    for chrom, site in zip(chroms, ts.sites()):
        site_id = tables.sites.add_row(
                        position=site.position,
                        ancestral_state=site.ancestral_state, 
                        metadata = chrom.tobytes())
    return(tables.tree_sequence())




@jit(nopython=True, nogil=True, cache=False)
def find_chrom(sites, genome_map):
    _, chroms = np.where((sites >= genome_map[0] ) & (sites <= genome_map[1]))
    chroms += 1
    return(chroms)

@jit(nopython=True, nogil=True, cache=True)
def make_genome_map(length, nchrom):
    starts = np.arange(0,length, length/nchrom)
    ends = np.arange(length/nchrom,length+1, length/nchrom)
    return(starts, ends)


@jit(nopython=True, nogil=True, cache=False)
def make_mask(chroms):
    # relies on chroms being in sorted order
    n = len(chroms)
    nmax = chroms.max()
    mask = np.ones(n*n)
    mask = mask.reshape(n, n).astype(np.bool_)
    offset = 0
    for x in np.bincount(chroms):
        if x>0:
            mask[offset:offset+x, offset:offset+x] = False
        offset +=x
    return(mask)

def convert_gt(gt):
    gt = (gt[::2] + gt[1::2])
    return(gt)

def get_overall_mean(ld_mat):
    return(ld_mat[np.triu_indices_from(ld_mat, k=1)].mean())

def get_r2_fast(geno_mat):
    norm_snps = (geno_mat - sp.mean(geno_mat, 0)) / sp.std(geno_mat, 0, ddof=1)
    norm_snps = norm_snps.T
    num_snps, num_indivs = norm_snps.shape
    ld_mat = (sp.dot(norm_snps, norm_snps.T) / float(num_indivs-1))**2
    return(ld_mat)