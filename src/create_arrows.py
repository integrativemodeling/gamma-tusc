import IMP
import IMP.atom
import IMP.algebra
import sys,os

fn1 = sys.argv[1]    #stop
fn2 = sys.argv[2]    #start
out_fn = sys.argv[3]
chains = sys.argv[4]


mdl = IMP.Model()
mh1 = IMP.atom.read_pdb(fn1,mdl,IMP.atom.CAlphaPDBSelector())
mh2 = IMP.atom.read_pdb(fn2,mdl,IMP.atom.CAlphaPDBSelector())


ca1 = IMP.atom.Selection(mh1,chains=chains).get_selected_particles()
ca2 = IMP.atom.Selection(mh2,chains=chains).get_selected_particles()
#ca1 = IMP.core.get_leaves(mh1)
#ca2 = IMP.core.get_leaves(mh2)
outf = open(out_fn,'w')

#max_dist = 0.0
#for c1,c2 in zip(ca1,ca2):
#    x1 = IMP.core.XYZ(c1).get_coordinates()
#    x2 = IMP.core.XYZ(c2).get_coordinates()
#    dist = IMP.algebra.get_distance(x1,x2)
#    if dist>max_dist:
#        max_dist=dist

#cyl_rad = 0.4 #0.6
#cone_base = 0.55 #0.9
#frac_cyl = 0.7 #0.6

#cyl_rad = 0.6
#cone_base = 0.9
#frac_cyl = 0.6

### thicker, fewer
cyl_rad = 0.8
cone_base = 1.5
frac_cyl = 0.6
n_per_arrow = 10
n_skip = 2 #1 means don't skip
# new plan: get mean vector for every 10 CAs?
all_pairs = []
ct = 0
cur_x1 = IMP.algebra.Vector3D(0,0,0)
cur_x2 = IMP.algebra.Vector3D(0,0,0)

for c1,c2 in zip(ca1,ca2):
    x1 = IMP.core.XYZ(c1).get_coordinates()
    x2 = IMP.core.XYZ(c2).get_coordinates()
    cur_x1+=x1
    cur_x2+=x2
    ct+=1
    if ct==n_per_arrow:
        all_pairs.append([cur_x1/n_per_arrow,cur_x2/n_per_arrow])
        cur_x1 = IMP.algebra.Vector3D(0,0,0)
        cur_x2 = IMP.algebra.Vector3D(0,0,0)
        ct=0

for x1,x2 in all_pairs[0::n_skip]:
    #x2 = (x1+x2)/2
    #x2 = (x1+x2)/2
    #fdist = 4*IMP.algebra.get_distance(x1,x2)/max_dist #fidst1 is the biggest, fdist 0 is smallest dist
    # red goes from 0.75 to 1.0
    #r = fdist*0.45+0.55
    #g = 0.55*(1.0-fdist)
    #b = 0.55*(1.0-fdist)
    r = 1.0
    g = 0.0
    b = 0.0
    outf.write('.color %.2f %.2f %.2f\n'%(r,g,b))
    outf.write('.arrow %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f\n'%(x1[0],x1[1],x1[2],x2[0],x2[1],x2[2],cyl_rad,cone_base,frac_cyl))
    #outf.write('.arrow %i %i %i %i %i %i\n'%(x1[0],x1[1],x1[2],x2[0],x2[1],x2[2]))
outf.close()
