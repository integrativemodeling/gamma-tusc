from modeller import *
import sys
sys.path.append('../src/')
from sequence_tools import *

env = environ()
gcp4_pdblim = alignment(env,file='gcp4_pdblim.pir')

### make PDBlim alignments

aln2 = alignment(env,file='TUBGCP_promals.fasta',alignment_format="FASTA",
                 align_codes=['GCP4_HUMAN','GCP2_YEAST'])
ir = InsertionRemover(aln2,1e6)
ir.limit_to_pdb(gcp4_pdblim)
aln2_fix = ir.get_alignment()
write_align(aln2_fix,'gcp2/pdblim')
pir_to_pap(env,'gcp2/pdblim')

aln3 = alignment(env,file='TUBGCP_promals.fasta',alignment_format="FASTA",
                 align_codes=['GCP4_HUMAN','GCP3_YEAST'])
ir3 = InsertionRemover(aln3,1e6)
ir3.limit_to_pdb(gcp4_pdblim)
aln3_fix = ir3.get_alignment()
write_align(aln3_fix,'gcp3/pdblim')
pir_to_pap(env,'gcp3/pdblim')

### create PDBlim combos
# read the SSEs
s2 = get_seqs_from_pir('gcp2/pdblim.pir')
mo2 = ModelOptions('gcp2/GCP2.options')
s3 = get_seqs_from_pir('gcp3/pdblim.pir')
mo3 = ModelOptions('gcp3/GCP3.options')
len2 = len([s for s in s2['GCP2_YEAST'] if s not in '/-'])
len3 = len([s for s in s3['GCP3_YEAST'] if s not in '/-'])

# create new combo alignment and write it
aln_combo = alignment(env)
aln_combo.append_sequence(s2['GCP4_HUMAN']+'/'+s3['GCP4_HUMAN'])
aln_combo.append_sequence(s2['GCP2_YEAST']+'/'+s3['GCP3_YEAST'])
aln_combo[0].code = 'GCP4_HUMAN'
aln_combo[0].prottyp = 'structureX'
aln_combo[0].atom_file = 'gcp4_2copies.pdb'
aln_combo[0].source = 'HUMAN'
aln_combo[0].range = ('FIRST:A','LAST:B')
aln_combo[1].code='GTUSC'
write_align(aln_combo,'combo/pdblim23')

# adjust the starting position of GCP3 SSEs and combine it with GCP2's SSEs
sse2 = mo2.get_data()
sse3 = mo3.get_data()

s2s = sse2['sses']
s3s = sse3['sses']
for h in s3s['helices']:
    s2s['helices'].append([h[0]+len2,h[1]+len2])
for h in s3s['loops']:
    s2s['loops'].append([h[0]+len2,h[1]+len2])
for sheet in s3s['sheets']:
    tmp = []
    for h in sheet:
        tmp.append([h[0]+len2,h[1]+len2])
    s2s['sheets'].append(tmp)
stream = file('combo/pdblim.options','w')
yaml.dump(sse2,stream)
