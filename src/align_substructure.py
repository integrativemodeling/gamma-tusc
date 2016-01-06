# align two chains from two pdb files over a PDB range. write them both

import IMP
import IMP.atom
import IMP.algebra
import argparse
import re

def parse_sel(dat):
    m1 = re.search(r'^(\w+)-(\w+).(\w)$',dat)
    m2 = re.search(r'^(\w+)-.(\w)$',dat)
    m3 = re.search(r'^.(\w)$',dat)
    start=1
    stop=-1

    if m1:
        start,stop,chain=m1.group(1,2,3)
    elif m2:
        start,chain=m2.group(1,2)
    elif m3:
        chain=m3.group(1)
    else:
        print 'error with selection argument',dat
        exit()
    start=int(start)
    stop=int(stop)
    return (start,stop,chain)

def write_some_chains(hier,cnames,mdl,out_fn):
    chains=IMP.atom.get_by_type(hier,IMP.atom.CHAIN_TYPE)
    mh_out=IMP.atom.Hierarchy(IMP.Particle(mdl))
    for c in chains:
        if IMP.atom.Chain(c).get_id() in cnames:
            mh_out.add_child(IMP.atom.Chain(c))
    IMP.atom.write_pdb(mh_out,out_fn)

def parse_args():
    usage = """Align two structures. Can provide list of selections (see help with -h)
    """
    parser=argparse.ArgumentParser(usage)
    parser.add_argument("--f1",action="store",dest="f1",help="File 1")
    parser.add_argument("--f2",action="store",dest="f2",help="File 2")
    parser.add_argument("-w",action="store_true",dest="write_only_selected_chains",
                        help="Only write chains that were selected",
                        default=False)
    parser.add_argument("--out_f1",action="store",dest="out_f1",help="Write file 1 "
                        "same as orig unless you use flag -w, then it extracts that chain")
    parser.add_argument("--out_f2",action="store",dest="out_f2",help="Transformed file 2")
    parser.add_argument("-s",action="store",dest="selections",
                        nargs="+",
                        help="List of selections."
                        " Format is chimera style: start-stop.chain, start-.chain, .chain")
    result = parser.parse_args()
    if result.f1 is None or result.f2 is None:
        print "Must supply f1 and f2"
        exit()

    sels=[]
    if result.selections is None:
        pass
    else:
        for s in result.selections:
            sels.append(parse_sel(s))

    return result.f1,result.f2,sels,result.out_f1,result.out_f2,result.write_only_selected_chains

def run():
    f1,f2,selections,out_f1,out_f2,write_only_selected_chains = parse_args()
    print 'got selections',selections
    mdl=IMP.Model()
    mh1=IMP.atom.read_pdb(f1,mdl)
    mh2=IMP.atom.read_pdb(f2,mdl)

    ps1=[]
    ps2=[]
    if len(selections)==0:
        ps1=IMP.atom.get_leaves(mh1)
        ps2=IMP.atom.get_leaves(mh2)
    else:
        for sel in selections:
            start,stop,chain=sel
            if stop==-1:
                sss=IMP.atom.Selection(mh1,chain=chain,atom_type=IMP.atom.AtomType("CA"))
                stop=max([IMP.atom.get_residue(IMP.atom.Atom(a)).get_index() \
                          for a in sss.get_selected_particles()])
            s1=IMP.atom.Selection(mh1,chain=chain,residue_indexes=range(start,stop+1),
                                  atom_type=IMP.atom.AtomType("CA"))
            s2=IMP.atom.Selection(mh2,chain=chain,residue_indexes=range(start,stop+1),
                                  atom_type=IMP.atom.AtomType("CA"))
            ps1+=s1.get_selected_particles()
            ps2+=s2.get_selected_particles()
    coords1=[IMP.core.XYZ(p).get_coordinates() for p in ps1]
    coords2=[IMP.core.XYZ(p).get_coordinates() for p in ps2]
    trans=IMP.algebra.get_transformation_aligning_first_to_second(coords2,coords1)
    IMP.atom.transform(mh2,trans)

    if write_only_selected_chains:
        chs=[sel[2] for sel in selections]
        if out_f1:
            write_some_chains(mh1,chs,mdl,out_f1)
        if out_f2:
            write_some_chains(mh2,chs,mdl,out_f2)
    else:
        if out_f1:
            IMP.atom.write_pdb(mh1,out_f1)
        if out_f2:
            IMP.atom.write_pdb(mh2,out_f2)

if __name__=="__main__":
    run()
