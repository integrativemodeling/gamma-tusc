# does loop modeling, also adds chain breaks where the insertion file requests it

from modeller import *
from modeller.automodel import *
import os,sys
import string
from optparse import OptionParser
import sequence_tools
import random

def parse_args():
    usage = """usage %prog [options] <align_fn> <nmodels>
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
    parser.add_option("-o","--out_dir",dest="out_dir",default='',
                      help="write files to this directory instead of current working one")
    (options, args) = parser.parse_args()
    if len(args) != 2:
        parser.error("incorrect number of arguments "+str(len(args))+" " + str(args))
    return [options,args]

class MyModel(automodel): #dope_loopmodel):
    breaks=[]
    sym_copies=[]
    sses=[]

    def set_data(self,inserts=[],sses=[],sym_pairs=[]):
        self.inserts = inserts
        self.sses = sses
        self.sym_pairs = sym_pairs

    #def select_loop_atoms(self):
    #    sel0 = super(MyModel,self).select_loop_atoms()
    #    sel1 = selection(self.residue_range(1894,1897))
    #    print '>>>refining additional loops',''.join(self.residues[r].code for r in xrange(1893,1897))
    #    return sel0 | sel1

    def special_restraints(self,aln):
        rsr = self.restraints
        at = self.atoms
        res = self.residues

        ### adding sses
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
        for ins in self.inserts:
            start,ilen=ins
            dist=ilen*3.5 #IMP.multifit.get_approximated_radius(ilen)
            print 'creating insert restraint from',res[start-1],'to',res[start],'dist',dist
            rsr.add(forms.upper_bound(group=physical.xy_distance,
                                      feature=features.distance(res[start-1].atoms['CA'],
                                                                res[start].atoms['CA']),
                                      mean=dist, stdev=0.1))

        ### adding symmetry restraint
        for spair in self.sym_pairs:
            s1=selection(self.residue_range(spair[0][0],
                                            spair[0][1])).only_atom_types('CA')
            s2=selection(self.residue_range(spair[1][0],
                                            spair[1][1])).only_atom_types('CA')
            print ''.join([r.code for r in self.residue_range(spair[0][0],spair[0][1])])
            print ''.join([r.code for r in self.residue_range(spair[1][0],spair[1][1])])
            self.restraints.symmetry.append(symmetry(s1,s2,1.0))
        if self.sym_pairs==[]:
            self.restraints.symmetry.report(1.0)


def run(env=None):
    options,args=parse_args()
    aln_fn,nmodels=args

    ### Setup
    if options.random_seed==0:
        random_seed=random.randint(-50000,-2)
    print 'using random seed',random_seed
    env = environ(rand_seed=random_seed)
    aln_fn = os.path.abspath(aln_fn)

    nmodels = int(nmodels)
    env.libs.topology.read(file='$(LIB)/top_heav.lib')
    env.libs.parameters.read(file='$(LIB)/par.lib')
    if options.data_dir!='':
        atom_dir=os.path.abspath(options.data_dir)
        env.io.atom_files_directory.append(atom_dir)

    ### read in extra data
    if options.options_fn!='':
        sses = []
        if options.options_fn!='':
            mo = sequence_tools.ModelOptions(options.options_fn)
            if options.no_restrain_beta:
                mo._data['sses']['sheets'] = []
            sses = mo.get_sses()
        symmetry_pairs = [] #mo.get_symmetries() #disabling for now, probably not needed
        inserts = mo.get_insertions()

    ### create model, passing it various data
    a0 = alignment(env,file=aln_fn)
    a = MyModel(env, alnfile=aln_fn,
                knowns=a0[0].code, sequence=a0[1].code,
                assess_methods=(assess.DOPE,assess.GA341))
    a.starting_model = 1
    a.ending_model = nmodels
    a.set_data(inserts,sses,symmetry_pairs)

    ### perform modeling
    if options.out_dir!='':
        out_dir = os.path.abspath(options.out_dir)
        owd = os.getcwd()
        print 'starting in',owd
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        os.chdir(out_dir)
        print 'working from',os.getcwd()
    if options.nloops>0:
        a.loop.starting_model=1
        a.loop.ending_model=options.nloops
        a.loop.md_level = refine.very_fast
    a.make()

    if options.out_dir!='':
        os.chdir(owd)

if __name__=="__main__":
    run()
