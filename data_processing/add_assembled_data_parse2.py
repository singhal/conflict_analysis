import re
import glob
import pandas as pd

def parse_fasta(file):
	f = open(file, 'r')
	seq = {}
	id = ''
	for l in f:
		if re.search('>', l):
			id = re.search('>(\S+)', l).group(1)
			seq[id] = ''
		else:
			seq[id] += l.strip()
	return seq

d = pd.read_csv("~/Desktop/samples/added_inds.csv")
counts = []
for ix, row in d.iterrows():
	file = '/Users/sonal/Desktop/new/%s.fasta' % row['sample']
	x = parse_fasta(file)
	counts.append(len(x))
d['num_contigs'] = counts
d.to_csv("~/Desktop/samples/added_inds2.csv", index=False)