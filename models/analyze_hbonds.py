# read in all the hbonds files in a directory
# for each, use rex to extract the chain and residues.
#           remove duplicates
#           identify inter/intra
# find the contacts that appear in at least X% of the files
# sort them and output table S2!

import sys,os,re,glob
from collections import defaultdict
from Bio import SeqIO
def interface(ch1,ch2):
    d = {'CD':"intra",
         'DI':"inter"}
    return d[''.join(sorted(ch1+ch2))]

def get_hbonds(hdir):
    inter = defaultdict(int)
    intra = defaultdict(int)
    nfiles = 0
    for nfn,fn in enumerate(glob.glob(hdir+'*')):
        inf = open(fn,'r')
        nfiles += 1
        found_inter = set()
        found_intra = set()
        for nl,l in enumerate(inf):
            if nl<6:
                continue
            fields = l.split()
            m = re.search(r':(\d+)\.(\w)@[\w\s]+:(\d+)\.(\w)',l)
            res1 = int(m.group(1))
            chain1 = m.group(2)
            res2 = int(m.group(3))
            chain2 = m.group(4)
            i = interface(chain1,chain2)
            key = (min(res1,res2),max(res1,res2))
            if i=='inter':
                found_inter.add(key)
            else:
                found_intra.add(key)
        for f in found_inter:
            inter[f]+=1
        for f in found_intra:
            intra[f]+=1
        inf.close()
    for i in inter:
        inter[i] = float(inter[i])/nfiles*100
    for i in intra:
        intra[i] = float(intra[i])/nfiles*100
    return inter,intra


min_frac = 20
closed_dir = 'results/hbonds_closed/'
open_dir = 'results/hbonds_open/'
seq_fn = '../data/gtub_yeast.fasta'
closed_inter,closed_intra = get_hbonds(closed_dir)
open_inter,open_intra = get_hbonds(open_dir)

# get sequence for labeling
with open(seq_fn) as handle:
    record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
seq = record_dict[record_dict.keys()[0]]


# get hbonds
all_keys = set(closed_inter.keys()+closed_intra.keys()+open_inter.keys()+open_intra.keys())
all_keys = sorted(list(all_keys))

# print table
print '            cInter cIntra oInter oIntra'
for k in all_keys:
    ret = '%i%s-%i%s\t'%(k[0],seq[k[0]-1],k[1],seq[k[1]-1])
    found = False
    if closed_inter[k]>min_frac:
        ret+='X\t'
        found = True
    else:
        ret+='\t'
    if closed_intra[k]>min_frac:
        ret+='X\t'
        found = True
    else:
        ret+='\t'
    if open_inter[k]>min_frac:
        ret+='X\t'
        found = True
    else:
        ret+='\t'
    if open_intra[k]>min_frac:
        ret+='X\t'
        found = True
    else:
        ret+='\t'
    if found:
        print ret
