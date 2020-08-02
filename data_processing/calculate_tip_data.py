import re
import pandas as pd
import glob
import gzip
import os
import numpy as np
import argparse
from rpy2.robjects.packages import importr
import rpy2.robjects as ro

keys = ['heterozygosity', 'num_AHE_contigs',
        'num_UCE_contigs', 'num_total',
        'avg_cov', 'avg_length',
  	'num_sites', 'per_GC',
	'avg_tip_len', 'avg_tip_sd']

def get_seq(f):
        seq = {}
        id = ''
        
        f = open(f, 'r')
        for l in f:
                if re.search('>', l):
                        id = re.search('>(\S+)', l).group(1)
			id = re.sub('^_R_', '', id)
                        seq[id] = ''
                else:
                        seq[id] += l.rstrip().upper()
        f.close()

        return seq

d = pd.read_csv("/scratch/drabosky_flux/sosi/phylogeny/squamate_phylogenomics_v11.csv")
d = d.ix[d.family_level == True]
qual = {}

for s in d['sample'].tolist():
	print(s)
	prg = '/nfs/turbo/lsa-rabosky/Lab/SqCL_July2017/PRG/%s.fasta' % s
	vcf = '/nfs/turbo/lsa-rabosky/Lab/SqCL_July2017/variants/%s.qual_filtered20.cov_filtered2.vcf.gz' % s

	qual[s] = {}
	for key in keys:
		qual[s][key] = None	

	if os.path.isfile(prg):
		seq = get_seq(prg)
		qual[s]['num_total'] = len(seq)
		lens = [len(y) for x, y in seq.items()]
		if len(lens) > 0:
			qual[s]['avg_length'] = np.mean(lens)
		qual[s]['num_AHE_contigs'] = len([x for x in seq if re.search('AHE', x)])
		qual[s]['num_UCE_contigs'] = len([x for x in seq if re.search('uce', x)])

	if os.path.isfile(vcf):
		het = {'diff': 0, 'denom': 0}
		cov = 0
		f = gzip.open(vcf, 'r')
		for l in f:
			if not re.search("^#", l):
				d = re.split('\t', l.rstrip())
				geno = re.search('^(\S\S\S)', d[9]).group(1)
				geno = re.split('/', geno)
				dp = re.search('DP=(\d+)', l).group(1)
				cov += int(dp)
				if int(dp) > 9:
					het['denom'] += 1
					if geno[0] != geno[1]:
                                        	het['diff'] += 1
		
		if het['denom'] > 0:
			qual[s]['heterozygosity'] = het['diff'] / float(het['denom'])
			qual[s]['avg_cov'] = cov / float(het['denom'])

# alignment
aln = {}
for s in qual:
	aln[s] = {'GC_count': 0, 'tot_seq': 0, 'miss_count': 0}
alns = glob.glob('/scratch/drabosky_flux/sosi/family_phylogeny/trim_300_0.05_0.3/*aln')
for a in alns:
	print(a)
	seq = get_seq(a)
	seq2 = {}
	for ind, s in seq.items():
		ind = re.sub('_R_', '', ind)
		seq2[ind]  = s

	for ind, s in seq2.items():
		if ind in aln:
			aln[ind]['tot_seq'] += len(s)
			aln[ind]['miss_count'] += (s.count('N') + s.count('-'))
			aln[ind]['GC_count'] += (s.count('G') + s.count('C'))

for ind in aln:
	if aln[ind]['tot_seq'] > 0:
		qual[ind]['num_sites'] = aln[ind]['miss_count'] / float(aln[ind]['tot_seq'])
		qual[ind]['per_GC'] = aln[ind]['GC_count'] / float(aln[ind]['tot_seq'])

# tree length
trees = glob.glob('/scratch/drabosky_flux/sosi/family_phylogeny/rooted_300_0.05_0.3/*tre')
importr('phangorn')
importr('ape')

tre = {}
for s in qual:
        tre[s] = {'tip_len': 0, 'tot_trees': 0, 'tip_sd': 0}
for tree in trees:
	print(tree)
	ro.r('t1 = read.tree("%s")' % tree)
	ro.r('root_tip = diag(vcv.phylo(t1))')
	ro.r('mrt = median(root_tip, na.rm=T)')
	ro.r('sd_tip = root_tip / mrt')
	root_tip = ro.r('root_tip')
	sd_tip = ro.r('sd_tip')
	names = ro.r('names(root_tip)')

	for a, b, c in zip(names, root_tip, sd_tip):
		if a in tre:
			tre[a]['tot_trees'] += 1
			tre[a]['tip_len'] += b
			tre[a]['tip_sd'] += c
for ind in tre:
	if tre[ind]['tot_trees'] > 0:
		qual[ind]['avg_tip_len'] = tre[ind]['tip_len'] / float(tre[ind]['tot_trees'])
		qual[ind]['avg_tip_sd'] = tre[ind]['tip_sd'] / float(tre[ind]['tot_trees'])

o = open('/scratch/drabosky_flux/sosi/family_phylogeny/tip_data.csv', 'w')
o.write('%s,%s\n' % ('ind', ','.join(keys)))
for ind in qual:
	o.write('%s,%s\n' % (ind, ','.join([str(qual[ind][x]) for x in keys])))
o.close()
