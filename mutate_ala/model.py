# for each in top ensemble build mutated model!

from modeller import *
from modeller.automodel import *
import os,sys

class MyModel(automodel):
    def select_atoms(self):
        sel = selection(self.residue_range('2075:Q', '2118:Q'),
                        self.residue_range('2031:P', '2074:P'),
                        self.residue_range('4193:7', '4236:7'),
                        self.residue_range('4149:6', '4192:6')).only_sidechain()
        return sel


aln_fn = os.path.abspath('aln.pir')
env = environ()
env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')
#env.edat.nonbonded_sel_atoms = 2
a=MyModel(env, alnfile=aln_fn,
            knowns='TEMPLATE', sequence='TARGET')
a.starting_model = 1
a.ending_model = 1
a.make()
