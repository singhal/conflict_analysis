import re
import numpy as np
import random
import os 
import glob
import subprocess as sp
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
import argparse

parser = argparse.ArgumentParser()
# tree 1
parser.add_argument("-t", "--tree", required=True, help='tree to evaluate')
# alndir
parser.add_argument("-a", "--alndir", required=True, help='aln dir')
# outdir
parser.add_argument("-o", "--outdir", required=True, help='out dir')
args = parser.parse_args()

outs = ['taeGut2', 'galGal5', 'chrPic1', 'allMis1', 'hg38']

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

    newseq = {}
    for id, s in seq.items():
        if id not in outs:
            newseq[id] = s
    for out in outs:
        if out in seq:
            newseq['outgroup'] = seq[out]
            break

    return newseq

def print_reduced(outdir, seq, locus, treefile, stem):
    # get tips
    tips = [tip for tip in seq]
    ro.r.assign('tips', tips)
    ro.r('tips = unlist(tips)')
    
    # get union
    ro.r('g = read.tree("%s")' % treefile)
    ro.r('keep = intersect(g$tip.label, tips)')

    ro.r('g = drop.tip(g, setdiff(g$tip.label, keep))')
    gtree = os.path.join(outdir, '%s.%s.tre' % (locus, stem))
    class1 = ro.r('class(g)')[0]
    if class1 == "phylo":
        nodenum = int(ro.r('Nnode(g)')[0])
        if nodenum > 1:
            ro.r('write.tree(g, "%s")' % gtree)
        else:
            gtree = ''
    else:
        gtree = ''

    len_keep = int(ro.r('length(keep)')[0])
    keep = ro.r('keep')
    alnout = os.path.join(outdir, '%s.%s.aln' % (locus, stem))
    if len_keep > 0 and len(gtree) > 0:
        o = open(alnout, 'w')
        for ind in keep:
            o.write('>%s\n%s\n' % (ind, seq[ind]))
        o.close()
    
    return gtree, alnout
    
importr('ape')
alndir = args.alndir
alns = glob.glob(alndir + '/*aln')

stem = re.search('([^/]+).tre', args.tree).group(1)

if not os.path.isdir(args.outdir):
    os.mkdir(args.outdir)
outdir = os.path.join(args.outdir, stem)
if not os.path.isdir(outdir):
    os.mkdir(outdir)

out = open(os.path.join(args.outdir, '%s.commands.txt' % stem), 'w')

print(len(alns))
for aln in alns:
    locus = re.sub('^.*/', '', aln)
    locus = re.sub('.fasta.aln', '', locus)

    # get seq
    seq = get_seq(aln)

    # print out reduced alignment
    # print out reduced tree
    (tree, aln) = print_reduced(outdir, seq, locus, args.tree, stem)
    if (len(tree) > 0) :
        call1 = 'raxmlHPC -m GTRGAMMA -n %s_%s -s %s -g %s -p %s' % (locus, stem, aln, tree, random.randint(1, 1000))   
        out.write('cd %s\n' % outdir)
        out.write(call1 + '\n')
out.close()
