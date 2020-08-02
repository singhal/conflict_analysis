import re
import glob
from Bio import SeqIO
import pandas as pd

def parse_phy(file):
	f = open(file, 'r')
	seq = {}
	head = f.readline()
	for l in f:
		l = l.rstrip()
		if len(l) > 5:
			d = re.split('\s+', l)
			seq[d[0]] = d[1]
	return seq

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

d = pd.read_csv("~/Dropbox/Oz_Crown_Ages/samples/previous/added_inds3.csv", encoding="latin-1")
inds = {}
for x, y in d.iterrows():
	inds[y['original']] = y['sample']

'''
# 10.1186/s12862-016-0611-6
dir = '/Users/sonal/Dropbox/Oz_Crown_Ages/prev_published_data/s12862-016-0611-6/'
files = glob.glob(dir + '*nexus')
ids = {}
seqs = {}

for ix, file in enumerate(files):
	locname = 'locus%s' % ix
	x = SeqIO.parse(file, "nexus")
	for record in x:
		id = record.id
		seq = str(record.seq)
		orig = len(seq)

		seq = re.sub('N', '', seq)
		seq = re.sub('\?', '', seq)
		seq = re.sub('\-', '', seq)

		if len(seq) / float(orig) > 0.2:
			if id not in ids:
				ids[id] = 1
			if id not in seqs:
				seqs[id] = {}
			seqs[id][locname] = seq

for id in seqs:
	out = '/Users/sonal/Desktop/new/%s.fasta' % inds[id]
	o = open(out, 'w')
	for name, s in seqs[id].items():
		o.write('>%s\n%s\n' % (name, s))
	o.close()
'''

dir = '/Users/sonal/Dropbox/Oz_Crown_Ages/prev_published_data/CH-15-248/'
files = glob.glob(dir + '*nex')
seqs = {}

for ix, file in enumerate(files):
	locname = 'locus%s' % ix
	x = SeqIO.parse(file, "nexus")
	for record in x:
		id = record.id
		seq = str(record.seq)
		orig = len(seq)

		seq = re.sub('N', '', seq)
		seq = re.sub('\?', '', seq)
		seq = re.sub('\-', '', seq)
		seq = seq.upper()

		if len(seq) / float(orig) > 0.2:
			if id not in seqs:
				seqs[id] = {}
			seqs[id][locname] = seq

for id in seqs:
	out = '/Users/sonal/Desktop/new/%s.fasta' % inds[id]
	o = open(out, 'w')
	for name, s in seqs[id].items():
		o.write('>%s\n%s\n' % (name, s))
	o.close()

dir = '/Users/sonal/Dropbox/Oz_Crown_Ages/prev_published_data/evv026/'
files = glob.glob(dir + '*nex')
seqs = {}

for ix, file in enumerate(files):
	locname = 'locus%s' % ix
	x = SeqIO.parse(file, "nexus")
	for record in x:
		id = record.id
		seq = str(record.seq)
		orig = len(seq)

		seq = re.sub('N', '', seq)
		seq = re.sub('\?', '', seq)
		seq = re.sub('\-', '', seq)
		seq = seq.upper()

		if len(seq) / float(orig) > 0.2:
			if id not in seqs:
				seqs[id] = {}
			seqs[id][locname] = seq

for id in seqs:
	out = '/Users/sonal/Desktop/new/%s.fasta' % inds[id]
	o = open(out, 'w')
	for name, s in seqs[id].items():
		o.write('>%s\n%s\n' % (name, s))
	o.close()

'''
dir = '/Users/sonal/Dropbox/Oz_Crown_Ages/prev_published_data/s12862-015-0503-1/'
files = glob.glob(dir + '*phylip')
seqs = {}

for ix, file in enumerate(files):
	locname = 'locus%s' % ix
	x = parse_phy(file)
	for id, seq in x.items():
		orig = len(seq)
		seq = re.sub('N', '', seq)
		seq = re.sub('\?', '', seq)
		seq = re.sub('\-', '', seq)
		seq = seq.upper()

		if len(seq) / float(orig) > 0.2:
			if id not in seqs:
				seqs[id] = {}
			seqs[id][locname] = seq


for id in seqs:
	out = '/Users/sonal/Desktop/new/%s.fasta' % inds[id]
	o = open(out, 'w')
	for name, s in seqs[id].items():
		o.write('>%s\n%s\n' % (name, s))
	o.close()

dir = '/Users/sonal/Dropbox/Oz_Crown_Ages/prev_published_data/syw001/'
files = glob.glob(dir + '*phy')
seqs = {}

for ix, file in enumerate(files):
	locname = 'locus%s' % ix
	x = parse_phy(file)
	for id, seq in x.items():
		orig = len(seq)
		seq = re.sub('N', '', seq)
		seq = re.sub('\?', '', seq)
		seq = re.sub('\-', '', seq)
		seq = seq.upper()

		if len(seq) / float(orig) > 0.2:
			if id not in seqs:
				seqs[id] = {}
			seqs[id][locname] = seq

for id in seqs:
	out = '/Users/sonal/Desktop/new/%s.fasta' % inds[id]
	o = open(out, 'w')
	for name, s in seqs[id].items():
		o.write('>%s\n%s\n' % (name, s))
	o.close()

dir = '/Users/sonal/Dropbox/Oz_Crown_Ages/prev_published_data/j.ympev.2017.03.017/'
files = glob.glob(dir + '*fasta')
seqs = {}

for ix, file in enumerate(files):
	locname = 'locus%s' % ix
	x = parse_fasta(file)
	for id, seq in x.items():
		orig = len(seq)
		seq = re.sub('N', '', seq)
		seq = re.sub('\?', '', seq)
		seq = re.sub('\-', '', seq)
		seq = seq.upper()

		if len(seq) / float(orig) > 0.2:
			if id not in seqs:
				seqs[id] = {}
			seqs[id][locname] = seq

for id in seqs:
	out = '/Users/sonal/Desktop/new/%s.fasta' % inds[id]
	o = open(out, 'w')
	for name, s in seqs[id].items():
		o.write('>%s\n%s\n' % (name, s))
	o.close()


dir = '/Users/sonal/Dropbox/Oz_Crown_Ages/prev_published_data/rsbl.2012.0331/'
files = glob.glob(dir + '*nex')
seqs = {}

for ix, file in enumerate(files):
	locname = 'locus%s' % ix
	x = SeqIO.parse(file, "nexus")
	for record in x:
		id = record.id
		seq = str(record.seq)
		orig = len(seq)

		seq = re.sub('N', '', seq)
		seq = re.sub('\?', '', seq)
		seq = re.sub('\-', '', seq)
		seq = seq.upper()

		if len(seq) / float(orig) > 0.2:
			if id not in seqs:
				seqs[id] = {}
			seqs[id][locname] = seq

for id in seqs:
	out = '/Users/sonal/Desktop/new/%s.fasta' % inds[id]
	o = open(out, 'w')
	for name, s in seqs[id].items():
		o.write('>%s\n%s\n' % (name, s))
	o.close()


dir = '/Users/sonal/Dropbox/Oz_Crown_Ages/prev_published_data/rsbl20170393/'
files = glob.glob(dir + '*phy')
seqs = {}

for ix, file in enumerate(files):
	locname = 'locus%s' % ix
	x = parse_phy(file)
	for id, seq in x.items():
		orig = len(seq)
		seq = re.sub('N', '', seq)
		seq = re.sub('\?', '', seq)
		seq = re.sub('\-', '', seq)
		seq = seq.upper()

		if len(seq) / float(orig) > 0.2:
			if id not in seqs:
				seqs[id] = {}
			seqs[id][locname] = seq

for id in seqs:
	out = '/Users/sonal/Desktop/new/%s.fasta' % inds[id]
	o = open(out, 'w')
	for name, s in seqs[id].items():
		o.write('>%s\n%s\n' % (name, s))
	o.close()


file = "/Users/sonal/Dropbox/Oz_Crown_Ages/prev_published_data/j.ympev.2014.08.023.nex"
f = open(file, 'r')
locs = {}
locnum = 1
for l in f:
	if re.search('CHARSET', l):
		start = re.search('(\d+)-', l).group(1)
		start = int(start) - 1
		end = re.search('-(\d+);', l).group(1)
		end = int(end)
		locs['locus%s' % locnum] = [start, end]
		locnum += 1
f.close()	

x = SeqIO.parse(file, "nexus")
for record in x:
	id = record.id
	seq = str(record.seq)

	out = '/Users/sonal/Desktop/new/%s.fasta' % inds[id]
	o = open(out, 'w')
	for loc in locs:
		o.write('>%s\n%s\n' % (loc, seq[locs[loc][0]:locs[loc][1]]))
	o.close()

file = "/Users/sonal/Dropbox/Oz_Crown_Ages/prev_published_data/j.ympev.2014.06.013.nex.txt"
f = open(file, 'r')
locs = {}
locnum = 1
for l in f:
	if re.search('CHARSET', l):
		start = re.search('(\d+)-', l).group(1)
		start = int(start) - 1
		end = re.search('-(\d+);', l).group(1)
		end = int(end)
		locs['locus%s' % locnum] = [start, end]
		locnum += 1
f.close()	

x = SeqIO.parse(file, "nexus")
for record in x:
	id = record.id
	seq = str(record.seq)

	out = '/Users/sonal/Desktop/new/%s.fasta' % inds[id]
	o = open(out, 'w')
	for loc in locs:
		o.write('>%s\n%s\n' % (loc, seq[locs[loc][0]:locs[loc][1]]))
	o.close()


file1 = "/Users/sonal/Dropbox/Oz_Crown_Ages/prev_published_data/jbi.12989.txt"
f = open(file1, 'r')
locs = {}
locnum = 1
for l in f:
	if re.search('DNA', l):
		start = re.search('(\d+)-', l).group(1)
		start = int(start) - 1
		end = re.search('-(\d+)', l).group(1)
		end = int(end)
		locs['locus%s' % locnum] = [start, end]
		locnum += 1
f.close()	

file = "/Users/sonal/Dropbox/Oz_Crown_Ages/prev_published_data/jbi.12989.nex"
x = SeqIO.parse(file, "nexus")
for record in x:
	id = record.id
	seq = str(record.seq)
	seq = seq.upper()

	out = '/Users/sonal/Desktop/new/%s.fasta' % inds[id]
	o = open(out, 'w')
	for loc in locs:
		o.write('>%s\n%s\n' % (loc, seq[locs[loc][0]:locs[loc][1]]))
	o.close()

file1 = "/Users/sonal/Dropbox/Oz_Crown_Ages/prev_published_data/j.ympev.2x016.07.002.txt"
f = open(file1, 'r')
locs = {}
locnum = 1
for l in f:
	if re.search('DNA', l):
		start = re.search('(\d+)-', l).group(1)
		start = int(start) - 1
		end = re.search('-(\d+)', l).group(1)
		end = int(end)
		locs['locus%s' % locnum] = [start, end]
		locnum += 1
f.close()	

file = "/Users/sonal/Dropbox/Oz_Crown_Ages/prev_published_data/j.ympev.2x016.07.002.phy"
x = parse_phy(file)

for id, seq in x.items():
	out = '/Users/sonal/Desktop/new/%s.fasta' % inds[id]
	o = open(out, 'w')
	for loc in locs:
		o.write('>%s\n%s\n' % (loc, seq[locs[loc][0]:locs[loc][1]]))
	o.close()

file1 = "/Users/sonal/Dropbox/Oz_Crown_Ages/prev_published_data/zoj.12392.charSets"
f = open(file1, 'r')
locs = {}
locnum = 1
for l in f:
	if re.search('DNA', l):
		start = re.search('(\d+)-', l).group(1)
		start = int(start) - 1
		end = re.search('-(\d+)', l).group(1)
		end = int(end)
		locs['locus%s' % locnum] = [start, end]
		locnum += 1
f.close()	

file = "/Users/sonal/Dropbox/Oz_Crown_Ages/prev_published_data/zoj.12392.phylip"
x = SeqIO.parse(file, "phylip-relaxed")
for record in x:
	id = record.id
	seq = str(record.seq)
	seq = seq.upper()

	out = '/Users/sonal/Desktop/new/%s.fasta' % inds[id]
	o = open(out, 'w')
	for loc in locs:
		o.write('>%s\n%s\n' % (loc, seq[locs[loc][0]:locs[loc][1]]))
	o.close()
'''