from modeller import *
import sys

def run(aln_fn,pdb_fn,code,out_prefix):

    # read seq
    env = environ()
    aln = alignment(env,file=aln_fn,align_codes=code,alignment_format="FASTA")
    full_seq = ''.join([p.code for p in aln[0].residues])

    # read PDB (todo: add segment as option)
    mdl = model(env,file=pdb_fn,model_segment=('FIRST:A','LAST:A'))

    # align
    new_aln = alignment(env)
    new_aln.append_sequence(full_seq)
    new_aln[0].code = "ORIG_SEQ"
    new_aln.append_model(mdl,'PDB_SEQ')
    new_aln.salign(overhang=50, gap_penalties_1d=(-450, -50),
                   alignment_type='tree', output='ALIGNMENT')
    new_aln.write(out_prefix+'.pap',alignment_format='PAP')
    new_aln.write(out_prefix+'.pir',alignment_format='PIR')

if __name__=="__main__":
    if len(sys.argv)!=5:
        print "This script extracts a sequence from an alignment and aligns it to a PDB file"
        print "Usage: <align_fn> <pdb_fn> <code> <out_prefix>"
    else:
        run(*sys.argv[1:])
