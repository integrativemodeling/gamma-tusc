# read in the original alignment, specifying template and target
# identify insertions and remove them, writing final template-target alignment
# ALSO create an alignment of the final target sequence to the original one so you can fix indexes later
# modifiers:
#   file aligning template to a PDB sequence. this creates more insertions
#   ModelOptions file contains SSEs. these "rescue" insertion columns

from modeller import *
import sys
from optparse import OptionParser
import sequence_tools



def parse_args():
    usage = """usage %prog [options] <alignment_fn> <template_code> <target_code> <max_insertion_size> <out_prefix>
    Removes insertions, with options to
    """
    parser=OptionParser(usage)
    parser.add_option("-s","--sse_fn",dest="sse_fn",default='',
                      help="an input file for making sse restraints. use if you want to ensure"
                      " that sse-restrainted residues aren't removed with other insertions.")
    parser.add_option("-p","--pdb_aln",dest="pdb_aln",default='',
                      help="alignment file for the template. changes template column to gap "
                      "where there is no structure")
    parser.add_option("-d","--data_dir",dest="data_dir",default='',
                      help="directory with atom files")
    parser.add_option("-o","--out_aln_to_orig_prefix",dest="out_aln_to_orig_prefix",default='',
                      help="Output an alignment file for the final sequence to the original. "
                      "Use if you need to recover original sequence numbering.")
    (options, args) = parser.parse_args()
    if len(args) != 5:
        parser.error("incorrect number of arguments "+str(len(args))+" " +str(args))
    return [options,args]

def run():
    ### read data
    env = environ()
    options, args = parse_args()
    aln_fn, template_code, target_code, max_ins_size, out_prefix=args
    max_ins_size = int(max_ins_size)

    aln = alignment(env,aln_fn,align_codes=[template_code,target_code])
    iremove = sequence_tools.InsertionRemover(aln,max_ins_size)
    if options.pdb_aln!='':
        iremove.limit_to_pdb(options.pdb_aln)
    if options.sse_fn!='':
        mo = ModelOptions(options.sse_fn)
        iremove.enable_sses(mo)
    fix_aln = iremove.get_alignment()
    fix_aln.write(out_prefix+'.pir',alignment_format="PIR")
    fix_aln.write(out_prefix+'.pap',alignment_format="PAP")
    print 'wrote alignment',out_prefix+'.pir'
    if options.out_aln_to_orig_prefix!='':
        orig_aln = iremove.get_current_target_alignment()
        orig_aln.write(options.out_aln_to_orig_prefix+'.pir',alignment_format="PIR")
        orig_aln.write(options.out_aln_to_orig_prefix+'.pap',alignment_format="PAP")
        print 'wrote aln to orig file',options.out_aln_to_orig_prefix+'.pir'

if __name__=="__main__":
    run()
