## Overview
This module creates models of the gamma tubulin small complex using MODELLER and IMP. It has config files to do flexible fitting in MDFF. And it has a set of analysis and plotting tools.

## Creating alignments for GCP2 and GCP3
#### Manual steps to prepare alignments
- get TUBGCP family from uniprot. align using PROMALS3D
- create "options" files with to mark what should be forced to SSE (and kept in model). In the paper this was done as an iterative process, changing the SSE options and repeating all the modeling steps below.
- get representative proteins of tubulin family (HHPRED at ss15, plus gamma tubulin human and yeast) and align using PROMALS3D

#### Extract alignment to human GCP4, but limit GCP4 to what's in the PDB file
```
> python src/align_to_pdb.py alignments/TUBGCP_promals.fasta data/gcp4_pos2.pdb GCP4_HUMAN alignments/gcp4_pdblim
> python src/align_to_pdb.py alignments/gtub_seed_promals.fasta data/gtub_human_3CB2.pdb GTUB_HUMAN alignments/gtub_pdblim
> python src/align_to_pdb.py alignments/gtub_seed_promals.fasta data/gtub_human_3CB2_edit.pdb GTUB_HUMAN alignments/gtub_edit_pdblim
```

#### Create a few types of alignments
"orig" alignments - no SSE additions
```
> python src/remove_insertions.py alignments/TUBGCP_promals.fasta GCP4_HUMAN GCP2_YEAST 10000 alignments/gcp2/orig
> python src/remove_insertions.py alignments/TUBGCP_promals.fasta GCP4_HUMAN GCP3_YEAST 10000 alignments/gcp3/orig
```

GCP2 - with and without SSEs
```
> python src/remove_insertions.py alignments/TUBGCP_promals.fasta GCP4_HUMAN GCP2_YEAST 5 alignments/gcp2/ins5_orig -p alignments/gcp4_pdblim.pir
> python src/remove_insertions.py alignments/TUBGCP_promals.fasta GCP4_HUMAN GCP2_YEAST 5 alignments/gcp2/ins5_sse -p alignments/gcp4_pdblim.pir -s alignments/gcp2/GCP2.options -o alignments/gcp2/ins5_to_orig

```
GCP3 - with and without SSEs
```
> python src/remove_insertions.py alignments/TUBGCP_promals.fasta GCP4_HUMAN GCP3_YEAST 5 alignments/gcp3/ins5_orig -p alignments/gcp4_pdblim.pir
> python src/remove_insertions.py alignments/TUBGCP_promals.fasta GCP4_HUMAN GCP3_YEAST 5 alignments/gcp3/ins5_sse -p alignments/gcp4_pdblim.pir -s alignments/gcp3/GCP3.options -o alignments/gcp3/ins5_to_orig

```

Gamma tubulin w/pdb and insertion removal but no SSEs
```
> python src/remove_insertions.py alignments/gtub_rp15_promals.fasta GTUB_HUMAN GTUB_YEAST 6 alignments/gtub/gtub_rp15_promals_ins6 -p alignments/gtub_pdblim.pir -o alignments/gtub/ins6_to_orig
> python src/remove_insertions.py alignments/gtub_rp15_promals.fasta GTUB_HUMAN GTUB_YEAST 10 alignments/gtub_edit/gtub_rp15_promals_ins10 -p alignments/gtub_edit_pdblim.pir -o alignments/gtub/ins10_to_orig
```

#### Combine the various alignments to setup modeling
From the alignment directory, call this "ad hoc" script to combine alignments and options files, adjusting numbering for everything.
```
> python create_combo.py
```

## Modeling
#### Homology modeling with MODELLER
Back to the root directory, call model.py to start building models.
```
> python src/model.py alignments/combo/ins5.pir 1 -d data/ -s alignments/combo/ins5.options -o models/closed_ins5/
```
Clean up the model by restoring correct chain names and residue indexes
```
> ~/imp_fast_clean/setup_environment.sh python src/fix_model.py -i models/closed_ins5/TARGET.B99990001.pdb -o models/closed_ins5/fixed.pdb -f alignments/gcp2/ins5_to_orig.pir alignments/gcp3/ins5_to_orig.pir alignments/gtub/ins6_to_orig.pir data/spc110_fragment.fasta -n 123344123344
```

#### EM flexible fitting with MDFF
- fixing chain names and res indexes. adding ytubulin + spc110.
- generate psf, add SSE and symmetry restraints
- run in MDFF
- run cleanup files

## Analysis
#### finding best models
Use VMD scoring script (in src). You can split this into multiple runs/files if you want
```
vmd -dispdev text -eofexit < src/score_cc.tcl -args data/ringmasked_crop.situs 6.9 models/closed_ins5 1 300 models/ccs_closed.dat
vmd -dispdev text -eofexit < src/score_cc.tcl -args data/open_map.situs 8.0 models/open_ins5 1 300 models/ccs_open.dat
```
Seperately run k-means clustering (MPI aware)
```
mpirun -n p8 $IMPFAST python src/cluster.py models/closed_ins5/ 300 tusc_flex.pdb 10 models/closed_clust10/
```
If you've already computed a distance matrix, just add the prefix to that file:
```
mpirun -n p8 $IMPFAST python src/cluster.py models/closed_ins5/ 300 tusc_flex.pdb 5 models/closed_clust5/ models/closed_clust10/cdist
```

#### calculating RMSF
Here I'm analyzing the RMSF of a particular cluster. This outputs both RMSF as a text file and also as the bfactors in a pdb file.
```
$IMPFAST python src/measure_rmsf.py "models/closed_clust4/cluster.1/*.pdb" models/closed_clust4/cluster.1/284.pdb models/closed_clust4/cluster.1/rmsf
$IMPFAST python src/measure_rmsf.py "models/open_ins5/*/tusc_flex.pdb" models/open_ins5/133/tusc_flex.pdb models/results/rmsf_open
```

#### calculating local similarity
calc_local.py

```
~/imp_fast_clean/setup_environment.sh python src/calc_local_similarity.py key_models/closed.pdb key_models/open.pdb BDGI local_scores/closed_d6_n2.txt -d 6

```
## Open state fitting
- setting up best models


## Some useful classes defined in sequence_tools.py
- SequenceShifter: input an alignment of a "final" modeling sequence to the raw original one. Then use this class to shift sequence numbering (from original to final).
- ModelOptions: store, read, and write options for modeling, including secondary structure restraints, insertions, and symmetry.
- InsertionRemover: input an alignment and this class will help you remove insertions above a specified length. Can convert columns to insertions when the PDB file is missing residues, and convert insertions back to non-insertions if they're on the SSE modeling "whitelist"