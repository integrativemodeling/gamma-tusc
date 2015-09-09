# based on two states, calculate the local neighborhood overlap.
# for each residue
#  get all triangles in stateA, stateB
#  for each top triangle match:
#   align by that triangle
#   expand by greedily finding nearest neighbors
#   calculate metric: frac overlap, or size/DRMS?


import IMP
import IMP.atom
import IMP.multifit2
from optparse import OptionParser

def parse_args():
    usage = """usage %prog [options] <f1.pdb> <f2.pdb> <chains_for_analysis> <outf>
    Compare two states of the same protein, f1 and f2.
    For each residue, tries to overlap neighborhoods.
    Reports overlap, DRMS and size of the matching nearby particles.
    Chains for analysis: look only at these residues
    Chains for universe: all atoms for picking neighborhoods (will just use whole structure)
    """
    parser=OptionParser(usage)
    parser.add_option("-d","--dist",dest="dist",default=6.0,type='float',
                      help="neighborhood size for search")
    parser.add_option("-o","--out2",dest="out2",default='',
                      help="Paint neighborhood changes on file 2")
    parser.add_option("-s","--start_frac",dest="start_frac",default=0,
                      type=float,
                      help="Starting fraction of residues to analyze")
    parser.add_option("-t","--stop_frac",dest="stop_frac",default=1.0,
                      type=float,
                      help="Stop fraction of residues to analyze")

    (options, args) = parser.parse_args()
    if len(args) != 4:
        parser.error("incorrect number of arguments "+str(len(args))+" " + str(args))
    return [options,args]

def get_ps_for_chains(mh,chs):
    ret = []
    for cname in chs:
        print 'gathering',cname
        sel = IMP.atom.Selection(mh,chain=cname)
        ps = sel.get_selected_particles()
        ret += ps
    return ret

def run():
    options,args = parse_args()
    fn1,fn2,chs_analysis,out_fn = args
    mdl = IMP.Model()
    mh1 = IMP.atom.read_pdb(fn1,mdl)
    mh2 = IMP.atom.read_pdb(fn2,mdl)
    rb1 = IMP.atom.create_rigid_body(mh1)

    # PS for analysis
    ps1a = IMP.atom.Selection(mh1,chain_ids=list(chs_analysis),
                              atom_type=IMP.atom.AtomType("CA")).get_selected_particles()
    ps2a = IMP.atom.Selection(mh2,chain_ids=list(chs_analysis),
                              atom_type=IMP.atom.AtomType("CA")).get_selected_particles()

    # PS for universe
    ps1u = IMP.atom.Selection(mh1).get_selected_particles()
    ps2u = IMP.atom.Selection(mh2).get_selected_particles()
    vs1u = [IMP.core.XYZ(p).get_coordinates() for p in ps1u]
    vs2u = [IMP.core.XYZ(p).get_coordinates() for p in ps2u]
    nn1 = IMP.algebra.NearestNeighbor3D(vs1u) #use this to find particles that need to match
    nn2 = IMP.algebra.NearestNeighbor3D(vs2u)

    # config
    outf = open(out_fn,'w')
    outf.write('# <chain> <residue> <DRMS> <match_len> <len/DRMS>\n')
    num_analysis = len(ps1a)
    start = int(num_analysis * options.start_frac)
    stop = int(num_analysis * options.stop_frac)
    print 'start',start,'stop',stop
    native_dist = 2.0
    for np,(p1,p2) in enumerate(zip(ps1a,ps2a)):
        if np>=start and np<stop:
            # get info
            res = IMP.atom.get_residue(IMP.atom.Atom(p1))
            this_res = res.get_index()
            this_chain = IMP.atom.Chain(res.get_parent()).get_id()

            v1 = IMP.core.XYZ(p1).get_coordinates()
            v2 = IMP.core.XYZ(p2).get_coordinates()
            #print '\n>>>> analyzing',this_chain,this_res

            # get particles for neighbor search
            ball1 = nn1.get_in_ball(v1,options.dist)

            # create alignment along backbone and transform
            bb1 = IMP.atom.Selection(mh1,chain_id=this_chain,residue_indexes=[this_res-1,this_res,this_res+1]).get_selected_particles()
            bb2 = IMP.atom.Selection(mh2,chain_id=this_chain,residue_indexes=[this_res-1,this_res,this_res+1]).get_selected_particles()
            trans = IMP.core.get_transformation_aligning_first_to_second(bb1,bb2)
            IMP.core.transform(rb1,trans)

            # get native overlap
            ngood = 0
            for idx in ball1:
                pos = IMP.core.XYZ(ps1u[idx]).get_coordinates()
                ball2 = nn2.get_in_ball(pos,native_dist)
                if len(ball2)>0:
                    ngood+=1
            frac_good = float(ngood)/len(ball1)
            IMP.core.transform(rb1,trans.get_inverse())
            #print 'frac good',frac_good
            outf.write('%s %i %.5f\n'%(this_chain,this_res,1.0-frac_good))
    outf.close()
if __name__=="__main__":
    run()
2
