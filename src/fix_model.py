import argparse
import IMP
import IMP.atom
from modeller import *
import string
import os
import sequence_tools

def parse_args():
    usage = """This script allows you to fix a model based on alignment files
    you input the model file, all the align_to_orig files, and the order that they appear (allowing repeats).
    output is fixed PDB file with the correct chain IDs and sequence numbering.
    """
    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument("-i",action="store",dest="model_fn",help="input model to fix",required=True)
    parser.add_argument("-o",action="store",dest="out_fn",help="output PDB file",required=True)
    parser.add_argument("-f",action="store",dest="aln_fns",required=True,nargs="+",
                        help="alignment files mapping molecule sequences to their originals."
                        "Can provide any number of them, and can just use a fasta file if you don't need to align.")
    parser.add_argument("-n",action="store",dest="order",required=True,
                        help="order that the (-f) sequences appear in the model (numbers start with 1)"
                        "e.g. for two repeated sequences you could write: -f aln1.pir aln2.pir -n 1212")
    args = parser.parse_args()
    return args

def run():
    args = parse_args()
    mdl = IMP.Model()
    mh = IMP.atom.read_pdb(args.model_fn,mdl)
    all_res = IMP.atom.get_by_type(mh,IMP.atom.RESIDUE_TYPE)
    mh_new = IMP.atom.Hierarchy.setup_particle(IMP.Particle(mdl))

    env = environ()
    chs_new = []
    cnames = string.ascii_uppercase
    cnames = list(cnames)
    cnames.reverse()
    n_final = 0
    for naln in args.order:
        ch = IMP.atom.Chain.setup_particle(IMP.Particle(mdl),cnames.pop())
        mh_new.add_child(ch)
        num_aln = int(naln)
        ext = os.path.splitext(args.aln_fns[num_aln-1])
        if ext[1]=='.pir':
            aln = alignment(env,file=args.aln_fns[num_aln-1])
            s_orig = sequence_tools.get_seq_from_aln(aln,0)
            s_final = sequence_tools.get_seq_from_aln(aln,1)
        elif ext[1]=='.fasta':
            aln = alignment(env,file=args.aln_fns[num_aln-1],alignment_format="FASTA")
            s_orig = sequence_tools.get_seq_from_aln(aln,0)
            s_final = sequence_tools.get_seq_from_aln(aln,0)
        else:
            print "Alignment file must be .pir or .fasta. You have",ext[1]
            exit()
        print 'reading',args.aln_fns[num_aln-1]


        for n0,(s0,s1) in enumerate(zip(s_orig,s_final)):
            if s1!='-':
                res = IMP.atom.create_clone(all_res[n_final])
                IMP.atom.Residue(res).set_index(n0+1)
                ch.add_child(res)
                n_final+=1
    IMP.atom.write_pdb(mh_new,args.out_fn)
    print 'wrote to',args.out_fn

if __name__=="__main__":
    run()
