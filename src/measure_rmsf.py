import IMP
import IMP.atom
import sys,os,glob
import numpy as np

file_match = sys.argv[1]
model_fn = sys.argv[2]
out_prefix = sys.argv[3]

mdl = IMP.Model()
all_coords = []

for nfn,fn in enumerate(glob.glob(file_match)):
    print fn
    mh = IMP.atom.read_pdb(fn,mdl,IMP.atom.CAlphaPDBSelector())
    sel = IMP.atom.Selection(mh,chains = 'ABCDEF')
    for nps,p in enumerate(sel.get_selected_particles()):
        if nfn==0:
            all_coords.append([])
        all_coords[nps].append(IMP.core.XYZ(p).get_coordinates())
    IMP.atom.destroy(mh)
    del mh

print 'calculating RMSF'
outf = open(out_prefix+'.dat','w')
rmsfs = []
for nps,vlist in enumerate(all_coords):
    meanpos = IMP.algebra.get_centroid(vlist)
    rmsf = np.sqrt(sum(np.linalg.norm(v - meanpos)**2/len(vlist) for v in vlist))
    rmsfs.append(rmsf)
    outf.write('%i %.6f\n'%(nps,rmsf))
outf.close()

print 'got',len(rmsfs),'rmsfs'
print 'writing bfactors'
mh0 = IMP.atom.read_pdb(model_fn,mdl)
mh1 = IMP.atom.Hierarchy.setup_particle(IMP.Particle(mdl))
nres = 0
for chain in mh0.get_children()[:6]:
    for res in chain.get_children():
        for atom in res.get_children():
            IMP.atom.Atom(atom).set_temperature_factor(rmsfs[nres])
        nres += 1
    mh1.add_child(chain)
print 'total res',nres
IMP.atom.write_pdb(mh1,out_prefix+'.pdb')
