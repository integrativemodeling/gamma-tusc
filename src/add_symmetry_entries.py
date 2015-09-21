import IMP
import IMP.atom
import sys
from optparse import OptionParser

def parse_args():
    usage = """usage %prog [options] <input.pdb> <output.pdb> <num_chains_per_group>
    Add symmetry identifiers. Must be consecutive in PDB file.
    """
    parser=OptionParser(usage)
    (options, args) = parser.parse_args()
    if len(args) != 3:
        parser.error("incorrect number of arguments "+str(len(args))+" " +str(args))
    return [options,args]

def run():
    ### read data
    mdl = IMP.Model()
    options,args = parse_args()
    in_fn,out_fn,n_per_group=args
    n_per_group = int(n_per_group)

    mh = IMP.atom.read_pdb(in_fn,mdl)
    chs = mh.get_children()
    if len(chs)%n_per_group!=0:
        print "WARNING: Number per group does not divide evenly into num chains"
    #num_groups= len(chs)/n_per_group
    #print num_groups

    for p in IMP.core.get_leaves(mh):
        IMP.atom.Atom(p).set_occupancy(0)
        IMP.atom.Atom(p).set_temperature_factor(0)

    for nchain,c in enumerate(chs):
        sel = IMP.atom.Selection(c,atom_type=IMP.atom.AtomType("CA"))
        ngroup = nchain/n_per_group
        n_in_group = nchain%n_per_group
        print IMP.atom.Chain(c).get_id(),len(IMP.core.get_leaves(c)),ngroup+1,n_in_group+1
        for p in sel.get_selected_particles():
            IMP.atom.Atom(p).set_occupancy(1) #ngroup+1)
            IMP.atom.Atom(p).set_temperature_factor(ngroup+1) #n_in_group+1)

    IMP.atom.write_pdb(mh,out_fn)
if __name__=="__main__":
    run()
