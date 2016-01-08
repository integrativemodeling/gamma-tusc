from chimera import runCommand as rc
import glob
import sys,os

wd = sys.argv[3]       # provide starting directory ($PWD). sorry.
in_match = sys.argv[4] # pdbs to analyzie
out_dir = sys.argv[5]  # directory to dump results

os.chdir(wd)

try:
    os.stat(out_dir)
except:
    os.mkdir(out_dir)

for nfn,fn in enumerate(glob.glob(in_match)):
    print 'reading',fn
    out_fn = os.path.join(out_dir,str(nfn)+'.hbonds')
    rc('open %s'%fn)
    rc('sel #0:.C #0:.D #0:.I')
    rc('hbonds intermodel false intramol false distSlop 2.0 angleSlop 90.0 selRestrict both saveFile "%s"'%out_fn)
    rc('close all')
rc('stop now')
