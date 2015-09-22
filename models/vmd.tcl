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
