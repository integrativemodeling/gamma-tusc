#!/bash/bash

export SDIR=~/projects/gamma-tusc/src/
export PREFIX=closed_ins5/
export FN=tusc.pdb

for (( c1=$1; c1<=$2; c1++ ))
do
    cd $PREFIX/$c1/
    if [ ! -f tusc_autopsf.pdb ]
    then
        echo "processing" $c1
        python $SDIR/fix_het.py tusc.pdb tusc_tmp.pdb
        mv tusc_tmp.pdb tusc.pdb
        vmd -dispdev text -eofexit < ../vmd.tcl > vmd.log
        $IMPFAST python $SDIR/add_symmetry_entries.py tusc_autopsf.pdb sym.pdb 6
    fi
    cd ../../
done
