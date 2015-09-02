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
    usage = """usage %prog [options] <f1.pdb> <f2.pdb> <chains_for_universe> <chains_for_analysis> <outf>
    Compare two states of the same protein, f1 and f2.
    For each residue, tries to overlap neighborhoods.
    Reports overlap, DRMS and size of the matching nearby particles.
    Chains for universe: all atoms for picking neighborhoods
    Chains for analysis: look only at these residues
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
    if len(args) != 5:
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
    fn1,fn2,chs_univ,chs_analysis,out_fn = args
    mdl = IMP.Model()
    mh1u = IMP.atom.read_pdb(fn1,mdl)
    mh2u = IMP.atom.read_pdb(fn2,mdl)
    mh1a = IMP.atom.read_pdb(fn1,mdl,IMP.atom.CAlphaPDBSelector())
    mh2a = IMP.atom.read_pdb(fn2,mdl,IMP.atom.CAlphaPDBSelector())


    # PS for analysis
    ps1a = IMP.atom.Selection(mh1a,chain_ids=list(chs_analysis)).get_selected_particles()
    ps2a = IMP.atom.Selection(mh2a,chain_ids=list(chs_analysis)).get_selected_particles()

    # PS for universe
    ps1u = IMP.atom.Selection(mh1u,chain_ids=list(chs_univ)).get_selected_particles()
    ps2u = IMP.atom.Selection(mh2u,chain_ids=list(chs_univ)).get_selected_particles()
    vs1u = [IMP.core.XYZ(p).get_coordinates() for p in ps1u]
    vs2u = [IMP.core.XYZ(p).get_coordinates() for p in ps2u]
    nn1 = IMP.algebra.NearestNeighbor3D(vs1u)
    nn2 = IMP.algebra.NearestNeighbor3D(vs2u)

    # config
    param = IMP.multifit2.MultiFitParams()
    query_match_radius = 6.0 #options.dist #5.0
    dscores = []
    outf = open(out_fn,'w')
    outf.write('# <chain> <residue> <DRMS> <match_len> <len/DRMS>\n')
    num_analysis = len(ps1a)
    start = int(num_analysis * options.start_frac)
    stop = int(num_analysis * options.stop_frac)
    print 'start',start,'stop',stop
    for np,(p1,p2) in enumerate(zip(ps1a,ps2a)):
        if np>=start and np<stop:
            # get info
            res = IMP.atom.get_residue(IMP.atom.Atom(p1))
            this_res = res.get_index()
            this_chain = IMP.atom.Chain(res.get_parent()).get_id()
            v1 = IMP.core.XYZ(p1).get_coordinates()
            v2 = IMP.core.XYZ(p2).get_coordinates()
            print '\n>>>> analyzing',this_chain,this_res

            # get balls - NOT including the residue itself
            s1 = set(IMP.atom.Selection(mh1u,chain_id=this_chain,residue_index=this_res).get_selected_particles())
            s2 = set(IMP.atom.Selection(mh2u,chain_id=this_chain,residue_index=this_res).get_selected_particles())

            # try expanding radius til you have enough to match
            test = 0
            radius = options.dist
            while test == 0:
                ball1 = [ps1u[i] for i in nn1.get_in_ball(v1,radius) if ps1u[i] not in s1]
                ball2 = [ps2u[j] for j in nn2.get_in_ball(v2,radius) if ps2u[j] not in s2]
                ap = IMP.multifit2.PointAlignment(ball2,param) # TO these ps
                results = ap.align(ball1,query_match_radius)
                test = len(results)
                radius += 0.5


            trans,pairs = results[0]
            id2,id1 = zip(*pairs)
            drms = IMP.atom.get_drms([vs1u[i] for i in id1],
                                     [vs2u[j] for j in id2])
            print 'radius',radius,'drms',drms,'match size',len(id1),'size/drms',len(id1)/drms
            del ap
            #print 'trans',trans
            #print 'got',len(matches),'matches'
            #trans,pairs = IMP.multifit2.get_alignment_first_to_second(ball2,ball1)
            #print 'res',n1,'num pairs',len(pairs),'trans',trans
            outf.write('%s %i %.5f %i %.5f\n'%(this_chain,this_res,drms,len(id1),float(len(id1))/drms))
    outf.close()
if __name__=="__main__":
    run()
2
