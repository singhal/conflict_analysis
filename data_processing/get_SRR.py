from Bio import Entrez
import urllib
import re
import pandas as pd
import lxml.etree

Entrez.email = 'sosi@umich.edu'

# https://gist.github.com/martijnvermaat/4619109

# ascp -k1 -Tr -l1000M -i ~/.aspera/cli/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-
# instant/reads/ByRun/sra/SRR/SRR304/SRR304976/SRR304976.sra SRR304976.sra

def get_runs(sample):
	handle = Entrez.esearch(db='sra', term=sample)
	record = Entrez.read(handle)

	if not len(record['IdList']) == 1:
		print('Found %d entries in SRA for "%s" instead of the expected 1'
			  % (len(record['IdList']), sample))
	if len(record['IdList']) > 0:
		result = record['IdList'][0]

		handle = Entrez.efetch(db='sra', id=result)
		entry = lxml.etree.parse(handle)
		return entry.xpath('//EXPERIMENT_PACKAGE_SET/EXPERIMENT_PACKAGE/RUN_SET/RUN/@accession')
	else:
		return []

def run_file(file):
	d = pd.read_csv(file)
	outfile = re.sub('.csv', '_v2.csv', file)
	srrs = []
	for ix, row in d.iterrows():
		x = get_runs(row['SRA accession'])
		if len(x) > 0:
			srrs.append(x[0])
		else:
			srrs.append('NA')
	d['SRR'] = srrs
	d.to_csv(outfile)
	return outfile

file1 = "/Users/sonal/Dropbox/Oz_Crown_Ages/prev_published_data/syv058.csv"
file2 = "/Users/sonal/Dropbox/Oz_Crown_Ages/prev_published_data/j.ympev.2016.04.015.csv"
#out1 = run_file(file1)
#out2 = run_file(file2)
out1 = "/Users/sonal/Dropbox/Oz_Crown_Ages/prev_published_data/syv058_v2.csv"
out2 = "/Users/sonal/Dropbox/Oz_Crown_Ages/prev_published_data/j.ympev.2016.04.015_v2.csv"

d1 = pd.read_csv(out1, index_col=0)
d1.columns = ['lineage', 'sample', 'SRA', 'SRR']
d2 = pd.read_csv(out2, index_col=0)
d2.columns = ['lineage', 'sample', 'SRA', 'SRR']

d = d1.append(d2)

d['lineage'] = [re.sub(' ', '_', x) for x in d.lineage]
d['sample'] = [re.sub(' ', '_', x) for x in d['sample']]

d = d[d.SRR.notnull()]
d.to_csv("~/Dropbox/Oz_Crown_Ages/prev_published_data/published_samples.csv", index=False)