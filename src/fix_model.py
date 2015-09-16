import argparse

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
    print args.aln_fns,args.order
if __name__=="__main__":
    run()
