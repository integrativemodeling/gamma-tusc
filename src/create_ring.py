# create a ring for some files
# this is useful for doing correct CC scoring!

import IMP
import IMP.em
import sys
import os
import string
import glob
from optparse import OptionParser


def add_transformed_clones(mh_global,mh,transform,ntrans):
    cur_trans = transform
    for nt in xrange(ntrans):
        mhc = IMP.atom.create_clone(mh)
        rbc = IMP.atom.create_rigid_body([mhc])
        IMP.core.transform(rbc,cur_trans)
        mh_global.add_child(mhc)
        cur_trans=transform*cur_trans


#prefix = sys.argv[1]
file_match = sys.argv[1]
#nstart = int(sys.argv[2])
#nstop = int(sys.argv[3])
state = 'closed'
data_dir = 'data/'

if state == 'closed':
    map_fn = os.path.join(data_dir,'ringmasked-crop.mrc')
    resolution  = 6.9
    apix = 1.88
    origin = [95,100,104]
    vt = [-0.0136204, -0.0500845, 18.8688]
    rot = [0.885275, 0.000900034, 3.91722e-05, -0.465068]
elif state == 'open':
    maps_fn = os.path.join(data_dir,'data/openmap-crop.mrc')
    resolution = 8.0
    apix = 1.19
    origin = [160,160,160]
    vt = [-0.524138, -0.870234, 22.2567]
    rot = [0.885414, 0.0164601, -0.0249601, -0.46384]

num_trans = 2
num_itrans = 4
trans = IMP.algebra.Transformation3D(IMP.algebra.Rotation3D(rot),
                                     IMP.algebra.Vector3D(vt))

itrans = trans.get_inverse()

### read proteins and add transforms
mdl = IMP.Model()
for fname in glob.glob(file_match):
    path,b = os.path.split(fname)
    if 'ring' in b or 'stat' in b:
        continue
    #fname = os.path.join(prefix,str(nfile)+'.pdb')
    print 'making ring from',fname
    mh = IMP.atom.read_pdb(fname,mdl)
    mh_ring = IMP.atom.Hierarchy(IMP.Particle(mdl))
    mh_ring.add_child(mh)
    add_transformed_clones(mh_ring,mh,trans,num_trans)
    add_transformed_clones(mh_ring,mh,itrans,num_itrans)
    IMP.atom.write_pdb(mh_ring,os.path.join(path,b.split('.')[0]+'_ring.pdb'))
    IMP.atom.destroy(mh_ring)
    del mh
    del mh_ring






'''
### read protein
mdl=IMP.Model()
sel=IMP.atom.NonWaterNonHydrogenPDBSelector()
if options.ca_only:
    sel=IMP.atom.CAlphaPDBSelector()
mh0=IMP.atom.read_pdb(prot_fn,mdl,sel)
if options.chains!='':
    mh1=IMP.atom.Hierarchy.setup_particle(IMP.Particle(mdl))
    for c in IMP.atom.get_by_type(mh0,IMP.atom.CHAIN_TYPE):
        if IMP.atom.Chain(c).get_id() in options.chains:
            mh1.add_child(IMP.atom.Chain(c))
else:
    mh1=mh0
rb1=IMP.atom.create_rigid_body([mh1],'rb0')

### if transforms were reqested, add to mh:
if options.num_trans!=0 or options.num_invtrans!=0:
    rot=IMP.algebra.get_rotation_from_matrix(0.568131,-0.822938,0.0,0.822938,0.56813,0.0,0.0,0.0,1.0)
    if options.open:
        trans=IMP.algebra.Transformation3D(rot,[0,0,-22.2])
    else:
        trans=IMP.algebra.Transformation3D(rot,[0,0,-18.8])
    invtrans=trans.get_inverse()
    mh_global=IMP.atom.Hierarchy(IMP.Particle(mdl))
    mh_global.add_child(mh1)
    add_transformed_clones(mh_global,mh1,invtrans,options.num_invtrans,1)
    add_transformed_clones(mh_global,mh1,trans,options.num_trans,options.num_invtrans)
    print 'added',options.num_trans+options.num_invtrans,'copies'
else:
    mh_global=mh1
if options.out_pdb!='':
    cnames=string.uppercase+string.lowercase
    for nc,c in enumerate(IMP.atom.get_by_type(mh_global,IMP.atom.CHAIN_TYPE)):
        IMP.atom.Chain(c).set_id(cnames[nc])
    IMP.atom.write_pdb(mh_global,options.out_pdb)

### read map and compute score. remember the returned scores are 1-CC
dmap=IMP.em.read_map(in_map_fn,IMP.em.MRCReaderWriter())
dmap.get_header().set_resolution(resolution)
dmap.update_voxel_size(voxel_size)

if options.open:
    o=-160*1.19
    dmap.set_origin(o,o,o)
else:
    pass
    #dmap.set_origin(-95*1.88,-100*1.88,-104*1.88)

### try sampling density map
ps=IMP.core.get_leaves(mh_global)


print 'creating map'

#smap=IMP.em.SampledDensityMap(ps,resolution,voxel_size)
smap=IMP.em.SampledDensityMap(dmap.get_header())
smap.set_particles(ps)
smap.resample()
smap.calcRMS()
IMP.em.write_map(smap,options.out_map,IMP.em.MRCReaderWriter())

#smap=IMP.em.read_map('tst_map.mrc',IMP.em.MRCReaderWriter())
print 'wrote map'
threshold = smap.get_header().dmin+0.0000001
print 'threshold',threshold
#localcc=IMP.em.CoarseCC.cross_correlation_coefficient(dmap,smap,threshold)
localcc=IMP.em.CoarseCC.cross_correlation_coefficient(smap,dmap,threshold)
print 'score',localcc
'''
'''


print 'num particles',len(ps)
score=IMP.em.compute_fitting_scores(ps,
                                    dmap,
                                    [IMP.algebra.get_identity_transformation_3d()],
                                    False,
                                    options.local)

print 'score',1-score.get_score(0)
'''
