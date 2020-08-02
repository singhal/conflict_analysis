import re
import os
import glob

def compare_inds(t1, t2, c):
	unique = True

	# print(set(t1.keys()))
	# print(set(c.keys()))

	if set(t1.keys()) != set(c.keys()):
		if len(t1.keys()) <= len(set(c.keys())):
			# need to make sure they aren't equal b/c of missing data
			clades = []
			for ind in c:
				if ind in t1:
					clades.append(1)
				elif ind in t2:
					clades.append(2)
			# print(clades)
			clades = list(set(clades))
			if len(clades) == 1 and clades[0] == 1:
				unique = False
	else:
		unique = False

	# print(unique)
	# print('***')

	return unique

def parse_file(nodename, file):
	f = open(file, 'r')
	fline = f.readline().rstrip()

	c1 = re.search('\(0\) \(\(([^)]+)', fline).group(1)
	c2 = re.search(',([^)]+)\);', fline).group(1)

	c1 = dict({(x, 1) for x in re.split(',', c1)})
	c2 = dict({(x, 1) for x in re.split(',', c2)})

	# number of trees support
	support1 = int(re.search('conctrees \[(\d+)\]', f.readline()).group(1))
	support2 = int(re.search('conftrees \[(\d+)\]', f.readline()).group(1))
	support = support1  + support2 

	conf = []


	for l in f:
		if re.search('^\s+\(\d+\)', l):
			t1 = re.search('^.*\(\(([^)]+)', l).group(1)
			t2 = re.search(',([^)]+)\);', l).group(1)

			t1 = dict({(x, 1) for x in re.split(',', t1)})
			t2 = dict({(x, 1) for x in re.split(',', t2)})

			trees = int(re.search('^\s+\(\d+\)\s+\d+\s+(\d+)', l).group(1))

			unique = True
			for ix, c in enumerate(conf):
				u = compare_inds(t1, t2, c[0])
				if not u:
					unique = False
					if trees > c[1]:
						c[1] = trees
						
			if unique:
				conf.append([t1, trees, t2])

	max_con = max([con[1] for con in conf])

	out2 = os.path.join(outdir, '%s.csv' % nodename)
	o2 = open(out2, 'w')
	o2.write('nodename,percent,max_support_percent,species1,species2\n')
	prop_s = support1 / float(support)
	for ix, con in enumerate(conf):
		sps1 = ':'.join(sorted(con[0].keys()))
		sps2 = ':'.join(sorted(con[2].keys()))
		prop = con[1] / float(support)
		if prop >= 0.5 * prop_s or prop > 0.05:
			o2.write('%s,%s,%s,%s,%s\n' % (nodename, prop, prop_s, sps1, sps2))
	# print(sps, round(prop, 3))
	o2.close()

	return ':'.join(c1.keys()), support1, support2, support, max_con
	
outdir = '/Users/sonal/Desktop/activeWork/conflict_analysis/concat_edges/'
files = glob.glob("%s*txt" % outdir)
out = '/Users/sonal/Desktop/activeWork/conflict_analysis/concat_edges80.csv'

o = open(out, 'w')
o.write('nodename,tips,support,conflict,total_inf,max_conflict\n')
for file in files:
	print(file)
	nodename = re.search('(node\d+).txt', file).group(1)
	d = parse_file(nodename, file)
	o.write('%s,%s\n' % (nodename, ','.join([str(x) for x in d])))
o.close()
# ideal output would generate tree with original clade
# tree with most conflict clade
