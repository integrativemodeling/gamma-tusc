import IMP
import IMP.atom
import sys,os


pdb_fn = sys.argv[1]
score_fn = sys.argv[2]
out_fn = sys.argv[3]
copy_map={'B':'H',  #assuming gtusc here...
          'G':'A',
          'D':'J',
          'I':'C'}

mdl = IMP.Model()
mh1 = IMP.atom.read_pdb(pdb_fn,mdl)
inf = open(score_fn,'r')
for l in inf:
    if l[0]=='#':
        continue
    #chain,res,drms,match_len,score = l.split()
    chain,res,score = l.split()
    res = int(res)
    #drms = float(drms)
    #match_len = int(match_len)
    score = float(score)
    sel = IMP.atom.Selection(mh1,chain_ids=[chain,copy_map[chain]],residue_index=res)
    for p in sel.get_selected_particles():
        IMP.atom.Atom(p).set_temperature_factor(score)
IMP.atom.write_pdb(mh1,out_fn)
