import re
import numpy as np
import os 
import glob
import subprocess as sp
import p4
from rpy2.robjects.packages import importr
import rpy2.robjects as ro

def get_seq(aln):
	seq = {}
	id = ''
	f = open(aln, 'r')
	for l in f:
		if re.search('>', l):
			id = re.search('>(\S+)', l).group(1)
			id = re.sub('^_R_', '', id)
			seq[id] = ''
		else:
			seq[id] += l.rstrip().upper()
	f.close()

	return seq

def get_miss(seq, length):
	vals = [s.count('-') + s.count('N') for s in seq.values()]
	miss = np.mean(vals) / float(length)
	return miss

def get_het(seq, length):
	hets = [len(re.findall('[r|y|s|w|k|m|R|Y|S|W|K|M]', s)) for s in seq.values()]
	return np.mean(hets) / float(length)

def get_gc(seq, length):
        vals = [s.count('G') + s.count('C') for s in seq.values()]
        gc = np.mean(vals) / float(length)
        return gc

def get_pic(seq, length):
	pic = 0
	redund = 0

	for i in range(0, length):
		bases = [s[i] for s in seq.values()]
		bases = [bp for bp in bases if bp in ['A', 'T', 'C', 'G']]
		bases = list(set(bases))
		if len(bases) > 1:
			pic += 1
		if len(bases) > 2:
			redund += 1

	return [pic / float(length), redund / float(length)]

def get_tree_support(treefile):
	ro.r('t1 = read.tree("%s")' % treefile)

	support = ro.r('mean(as.numeric(t1$node.label), na.rm=T)')[0]
	ro.r('root_tip = diag(vcv.phylo(t1))')
	root_tip_var = ro.r('var(root_tip, na.rm=T)')[0]
	tree_len = ro.r('sum(t1$edge.length)')[0]
	ro.r('blen = t1$edge.length[1:Ntip(t1)]')
        ro.r('weird = t1$edge[which(blen > mean(blen) * %s), 2]' % 5)
	branch_out = ro.r('length(weird)')[0]

	return [support, root_tip_var, tree_len, branch_out]

	
importr('ape')
alndir = '/scratch/drabosky_flux/sosi/SqCL_July2017/family_phylogeny/trim_300_0.05_0.3/'
alns = glob.glob(alndir + '*aln')
treedir = '/scratch/drabosky_flux/sosi/SqCL_July2017/family_phylogeny/rooted_300_0.05_0.3/'
outdir = '/scratch/drabosky_flux/sosi/SqCL_July2017/family_phylogeny/'

res = {}
for aln in alns:
	locus = re.sub('^.*/', '', aln)
	locus = re.sub('.fasta.aln', '', locus)
	res[locus] = {}

	# get seq
	seq = get_seq(aln)

	# get length
	res[locus]['length'] = len(seq.values()[0])
	
	# get completeness
	res[locus]['occupancy'] = len(seq)

	# get missing
	res[locus]['missing'] = get_miss(seq, res[locus]['length'])

	# get heterozygosity
	res[locus]['heterozygosity'] = get_het(seq, res[locus]['length'])

	# GC content
	res[locus]['GC'] = get_gc(seq, res[locus]['length'])

	# pics / saturation
	vals = get_pic(seq, res[locus]['length'])
	res[locus]['PICs'] = vals[0]
	res[locus]['saturation'] = vals[1]

	# get tree
	# get average SH, root-tip var, tree len, tree height
	treefile = os.path.join(treedir, '%s.SH.tre' % locus)
	if os.path.isfile(treefile):
		vals = get_tree_support(treefile)
	else:
		vals = ['NA', 'NA', 'NA', 'NA']

	res[locus]['SH'] = vals[0]
	res[locus]['root_tip_var'] = vals[1]
	res[locus]['tree_len'] = vals[2]
	res[locus]['branch_outliers'] = vals[3]


loc1 = res.keys()[0]
types = sorted([type1 for type1 in res[loc1]])
out = open(os.path.join(outdir, 'loci_profiling1.csv'), 'w')
out.write('locus,%s\n' % (','.join(types)))
for locus in res:
	vals = [str(res[locus][type1]) for type1 in types] 
	out.write('%s,%s\n' % (locus, ','.join(vals)))
out.close()
