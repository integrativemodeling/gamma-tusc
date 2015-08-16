# does loop modeling, also adds chain breaks where the insertion file requests it

from modeller import *
from modeller.automodel import *
import os,sys
import string
from optparse import OptionParser
import sequence_tools
import random
import IMP
import IMP.multifit

def parse_args():
    usage = """usage %prog [options] <align_fn> <nmodels> <out_dir>
    """
    parser=OptionParser(usage)
    parser.add_option("-l","--nloops",dest="nloops",default=0,type='int',
                      help="number of loop models to make.")
    parser.add_option("-d","--data_dir",dest="data_dir",default='',
                      help="directory with atom files")
    parser.add_option("-s","--options_fn",dest="options_fn",default='',
                      help="file")
    parser.add_option("-b","--no_restrain_beta",dest="no_restrain_beta",action="store_true",
                      default=False,help="Flag if you want to ignore betas in the SSE file")
    parser.add_option("-r","--random_seed",dest="random_seed",default=0,type='int',
                      help="random seed. if 0 will pick one at random")
    (options, args) = parser.parse_args()
    if len(args) != 3:
        parser.error("incorrect number of arguments "+str(len(args))+" " + str(args))
    return [options,args]

class MyModel(automodel): #dope_loopmodel):
    breaks=[]
    sym_copies=[]
    sses=[]

    def set_data(self,inserts=None,sses=None,sym_pairs=None):
        self.inserts=inserts
        self.sses=sses
        self.sym_pairs=sym_pairs

    def special_restraints(self,aln):
        rsr = self.restraints
        at = self.atoms
        res = self.residues

        ### adding sses
        if not self.sses is None:
            for h in self.sses['helices']:
                print 'helix',h,''.join(res[r].code for r in xrange(h[0]-1,h[1]))
                rsr.add(secondary_structure.alpha(self.residue_range(h[0]-1,h[1]-1)))

            for sts in self.sses['sheets']:
                ends=[]
                for s in sts:
                    start=s[0]-1
                    stop=s[1]-1
                    sstrand=''.join([res[sss].code for sss in range(start,stop+1)])
                    print 'adding strand',start,sstrand
                    rsr.add(secondary_structure.strand(self.residue_range(start,stop)))
                    ends+=[start,stop]
                sslen=ends[1]-ends[0]+1
                print 'adding sheet of len',sslen
                print 'start',res[ends[0]].code,res[ends[0]].atoms['N']
                print 'stop',res[ends[-1]].code,res[ends[-1]].atoms['O']
                rsr.add(secondary_structure.sheet(res[ends[0]].atoms['N'],
                                                  res[ends[-1]].atoms['O'],
                                                  sheet_h_bonds=-1*sslen))

        ### adding distance restraint over inserts
        if not self.inserts is None:
            for ins in self.inserts:
                start,ilen=ins
                dist=ilen*3.5 #IMP.multifit.get_approximated_radius(ilen)
                print 'creating insert restraint from',res[start-1],'to',res[start],'dist',dist
                rsr.add(forms.upper_bound(group=physical.xy_distance,
                                          feature=features.distance(res[start-1].atoms['CA'],
                                                                    res[start].atoms['CA']),
                                          mean=dist, stdev=0.1))

        ### adding symmetry restraint
        if not self.sym_pairs is None:
            s1=selection(self.residue_range(self.sym_pairs[0][0],
                                            self.sym_pairs[0][1])).only_atom_types('CA')
            s2=selection(self.residue_range(self.sym_pairs[1][0],
                                            self.sym_pairs[1][1])).only_atom_types('CA')
            print ''.join([r.code for r in self.residue_range(self.sym_pairs[0][0],self.sym_pairs[0][1])])
            print ''.join([r.code for r in self.residue_range(self.sym_pairs[1][0],self.sym_pairs[1][1])])
            self.restraints.symmetry.append(symmetry(s1,s2,1.0))
        else:
            self.restraints.symmetry.report(1.0)


def run(env=None):
    options,args=parse_args()
    aln_fn,nmodels,out_dir=args

    ### Setup
    if options.random_seed==0:
        random_seed=random.randint(-50000,-2)
    print 'using random seed',random_seed
    env = environ(rand_seed=random_seed)
    aln_fn = os.path.abspath(aln_fn)
    out_dir = os.path.abspath(out_dir)
    nmodels = int(nmodels)
    env.libs.topology.read(file='$(LIB)/top_heav.lib')
    env.libs.parameters.read(file='$(LIB)/par.lib')
    if options.data_dir!='':
        atom_dir=os.path.abspath(options.data_dir)
        env.io.atom_files_directory.append(atom_dir)

    ### read in extra data
    sses=None
    symmetry_pairs=None
    inserts=None
    if options.options_fn!='':
        mo = sequence_tools.ModelOptions(options.options_fn)
        if options.no_restrain_beta:
            mo._data['sses']['sheets'] = []
        sses = mo.get_sses()
    #if options.symmetry_pairs!='':
    #    symmetry_pairs=[map(int,sp.split(':')) for sp in options.symmetry_pairs.split(',')]
    #if options.ins_fn!='':
    #    inserts=tools.parse_insert_file(options.ins_fn)

    ### create model, passing it various data
    a0 = alignment(env,file=aln_fn)
    a = MyModel(env, alnfile=aln_fn,
                knowns=a0[0].code, sequence=a0[1].code,
                assess_methods=(assess.DOPE,assess.GA341))
    a.starting_model = 1
    a.ending_model = nmodels
    a.set_data(inserts,sses,symmetry_pairs)

    ### perform modeling
    owd = os.getcwd()
    print 'starting in',owd
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    os.chdir(out_dir)
    print 'working from',os.getcwd()
    if options.nloops>0:
        #a.md_level=None
        a.loop.starting_model=1
        a.loop.ending_model=options.nloops
        #a.loop.md_level=refine.very_fast
    a.make()
    os.chdir(owd)

if __name__=="__main__":
    run()
