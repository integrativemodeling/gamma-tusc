# this script creates the combos from individual alignments. kind of difficult to generalize so I'm leaving it here.

from modeller import *
import sys
sys.path.append('../src/')
from sequence_tools import *

def get_seq_len(seq):
    return len([s for s in seq if s not in '-/.'])

### config
template_fn = 'closed_template.pdb'

### read sequences etc
env = environ()
s2 = get_seqs_from_pir('gcp2/ins5_sse.pir')
s3 = get_seqs_from_pir('gcp3/ins5_sse.pir')
sG = get_seqs_from_pir('gtub/gtub_rp15_promals_ins6.pir')
s110 = get_seqs_from_fasta('../data/spc110_fragment.fasta')
mo2 = ModelOptions('gcp2/GCP2.options')
mo3 = ModelOptions('gcp3/GCP3.options')
aln2o = alignment(env,file='gcp2/ins5_to_orig.pir')
aln3o = alignment(env,file='gcp3/ins5_to_orig.pir')
out_dir = 'combo/'

### create alignment
aln_combo = alignment(env)
temp_seq = '/'.join([s2['GCP4_HUMAN'],s3['GCP4_HUMAN'],sG['GTUB_HUMAN'],sG['GTUB_HUMAN'],s110['Spc110'],s110['Spc110']])
temp_seq += '/'+temp_seq
targ_seq = '/'.join([s2['GCP2_YEAST'],s3['GCP3_YEAST'],sG['GTUB_YEAST'],sG['GTUB_YEAST'],s110['Spc110'],s110['Spc110']])
targ_seq += '/'+targ_seq
aln_combo.append_sequence(temp_seq)
aln_combo.append_sequence(targ_seq)
aln_combo[0].name = "TEMPLATE"
aln_combo[0].code = "TEMPLATE"
aln_combo[0].atom_file = template_fn
aln_combo[0].prottyp = 'structureX'
aln_combo[0].source = "HUMAN"
aln_combo[1].name = "TARGET"
aln_combo[1].code = "TARGET"
write_align(aln_combo,out_dir+'ins5')

### Combine ModelOptions files
# create shifters
shift2 = SequenceShifter(aln2o)
shift3 = SequenceShifter(aln3o)

# add insertions to MOs, then shift
print 'shift 2'
mo2s = shift2.shift_model_options(mo2) # this should ALSO add insertions
# also, should not return insertions on N or C termini! b/c won't be able to translate.
print 'shift 3'
mo3s = shift3.shift_model_options(mo3)


# combine all the MO files together into one
len2 = get_seq_len(s2['GCP2_YEAST'])
len3 = get_seq_len(s3['GCP3_YEAST'])
lenG = get_seq_len(sG['GTUB_YEAST'])
lenS = get_seq_len(s110['Spc110'])
print 'len2',len2,'len3',len3,'lenG',lenG,'lenS',lenS
moF = ModelOptions()
moF.append(mo2s)
moF.append(mo3s,len2)
moF.append(mo2s,len2+len3+2*lenG+2*lenS)
moF.append(mo3s,2*len2+len3+2*lenG+2*lenS)

# add symmetries to MO and write
lenGTUSC = len2+len3+2*lenG+2*lenS
moF.symmetries = [[[1,lenGTUSC],[lenGTUSC+1,lenGTUSC*2]]]
moF.write(out_dir+'ins5.options')
