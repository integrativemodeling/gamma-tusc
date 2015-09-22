from modeller import *
import sys

def run(pdb_fn,out_pdb_fn):
    env=environ()
    io_data.hydrogen=False
    new_mdl=model(env,file=pdb_fn)
    new_mdl.write(out_pdb_fn)

if __name__=="__main__":
    if len(sys.argv)!=3:
        print "Usage: <in_fn> <out_fn>"
    else:
        run(*sys.argv[1:])
