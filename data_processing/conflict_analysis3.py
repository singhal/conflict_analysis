import re
import glob
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-n", "--name", help="name to search")
args = parser.parse_args()

outdir = '/scratch/drabosky_flux/sosi/family_phylogeny/conflicts_v2/'
bpdir = '/scratch/drabosky_flux/sosi/family_phylogeny/conflicts_v2/bp/'
trees = '/scratch/drabosky_flux/sosi/family_phylogeny/rooted_21June18.csv'
lnldir = os.path.join('/scratch/drabosky_flux/sosi/family_phylogeny/conflicts_v2/lnl_tests/', args.name)

def get_bp(bpf, trees, o):
	stem = re.sub('^.*/', '', bpf)
	d = re.split('\.', stem)
	conf = d[0]
	edge = d[1]
	
	f = open(bpf, 'r')
	head = f.next()
	yay = re.search(': (.*)$', f.next().rstrip()).group(1)
	no = re.search(': (.*)$', f.next().rstrip()).group(1)
	yay = [trees[int(x) - 1] for x in re.split('\s+', yay)]
	no = [trees[int(x) - 1] for x in re.split('\s+', no)]

	for tree in trees:
		if tree in yay:
			value = 'TRUE'
		elif tree in no:
			value = 'FALSE'
		else:
			value = 'NA'

		o.write('%s,%s,%s,CONFLICT,%s\n' % (conf, edge, tree, value))
	f.close()

def get_trees(trees):
	f = open(trees, 'r')
	t = []
	header = f.next()

	for l in f:
		d = re.sub("\"", "", l.strip())
		t.append(d)
	f.close()

	return t

def get_lnl(subdir, trees, o):
	files = glob.glob(subdir + '/*info*')
	lnls = {}

	d = re.split('\.', args.name)
        conf = d[0]
        edge = d[1]

	for file in files:
		f = open(file, 'r')
		lnl = 'NA'
		for l in f:
			if re.search('Final GAMMA-based Score of best tree', l):
				lnl = re.search('(\S+)$', l.rstrip()).group(1)
		f.close()

		locus = re.search('info\.(\S+)_%s' % args.name, file).group(1)
		lnls[locus] = lnl

	for tree in trees:
		if tree in lnls:
			o.write('%s,%s,%s,LNL,%s\n' % (conf, edge, tree, lnls[tree]))
		else:
			o.write('%s,%s,%s,LNL,%s\n' % (conf, edge, tree, 'NA'))

	
trees = get_trees(trees)
bpf = glob.glob(bpdir + '*%s*txt' % args.name)
o = open('%s%s.results.csv' % (outdir, args.name), 'w')
for bp in bpf:
 	get_bp(bp, trees, o)
get_lnl(lnldir, trees, o)
