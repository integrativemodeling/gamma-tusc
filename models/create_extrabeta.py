# this is an unfortunate script. MDFF doesn't pick up on all the beta sheets because of a STUPID bug.
# run this script to add extra betas. you will have to modify it if the structure changes at all. sorry.

# numbering system is in raw residue number (output from modeller works)
extra_sheets=(((20,24),(27,31)),
              ((200,202),(208,210)),
              ((367,370),(382,385)),
              ((759,761),(768,770)))
tx = 2118

ss = ''
for sheet in extra_sheets:
    s = sheet[0][0],sheet[0][1],sheet[1][0],sheet[1][1]
    ss += '(residue %i to %i) or (residue %i to %i) or '%(s[0],s[1],s[2],s[3])
    ss += '(residue %i to %i) or (residue %i to %i) or '%(s[0]+tx,s[1]+tx,s[2]+tx,s[3]+tx)
ss = ss.strip(' or ')
o = 'ssrestraints -psf tusc_autopsf.psf -pdb tusc_autopsf.pdb -o tusc-morebeta.txt -k_prot 2000 -sel "%s" -hbdonorsel "name N and backbone and (%s)" -hbaccsel "name O and backbone and (%s)" -hbonds -hbacut 45.0 -hbdcut 3.9 -labels' %(ss,ss,ss)
print o
