package require autopsf
package require mdff
package require ssrestraints
package require cispeptide
package require chirality
mol new tusc.pdb
autopsf -mol 0
mol new tusc_autopsf.pdb
mol addfile tusc_autopsf.psf
mdff gridpdb -psf tusc_autopsf.psf -pdb tusc_autopsf.pdb -o tusc_autopsf-grid.pdb
ssrestraints -psf tusc_autopsf.psf -pdb tusc_autopsf.pdb -o tusc-extrabonds.txt -hbonds -k_prot 2000
cispeptide restrain -o tusc-extrabonds-cispeptide.txt
chirality restrain -o tusc-extrabonds-chirality.txt
ssrestraints -psf tusc_autopsf.psf -pdb tusc_autopsf.pdb -o tusc-morebeta.txt -k_prot 2000 -sel "(residue 20 to 24) or (residue 27 to 31) or (residue 2138 to 2142) or (residue 2145 to 2149) or (residue 200 to 202) or (residue 208 to 210) or (residue 2318 to 2320) or (residue 2326 to 2328) or (residue 367 to 370) or (residue 382 to 385) or (residue 2485 to 2488) or (residue 2500 to 2503) or (residue 759 to 761) or (residue 768 to 770) or (residue 2877 to 2879) or (residue 2886 to 2888)" -hbdonorsel "name N and backbone and ((residue 20 to 24) or (residue 27 to 31) or (residue 2138 to 2142) or (residue 2145 to 2149) or (residue 200 to 202) or (residue 208 to 210) or (residue 2318 to 2320) or (residue 2326 to 2328) or (residue 367 to 370) or (residue 382 to 385) or (residue 2485 to 2488) or (residue 2500 to 2503) or (residue 759 to 761) or (residue 768 to 770) or (residue 2877 to 2879) or (residue 2886 to 2888))" -hbaccsel "name O and backbone and ((residue 20 to 24) or (residue 27 to 31) or (residue 2138 to 2142) or (residue 2145 to 2149) or (residue 200 to 202) or (residue 208 to 210) or (residue 2318 to 2320) or (residue 2326 to 2328) or (residue 367 to 370) or (residue 382 to 385) or (residue 2485 to 2488) or (residue 2500 to 2503) or (residue 759 to 761) or (residue 768 to 770) or (residue 2877 to 2879) or (residue 2886 to 2888))" -hbonds -hbacut 45.0 -hbdcut 3.9
