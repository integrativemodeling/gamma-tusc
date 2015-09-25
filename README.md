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
> python src/align_to_pdb.py alignments/gtub_seed_promals.fasta data/human_ytub_3CB2.pdb GTUB_HUMAN alignments/gtub_pdblim
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
> python src/model.py alignments/combo/ins5.pir 1 models/closed_ins5/ -d data/ -s alignments/combo/ins5.options -l 1
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
#### calculating RMSF
#### calculating local similarity
```
~/imp_fast_clean/setup_environment.sh python src/calc_local_similarity.py key_models/closed.pdb key_models/open.pdb BDGI local_scores/closed_d6_n2.txt -d 6

```
## Open state fitting
- setting up best models


## Some useful classes defined in sequence_tools.py
- SequenceShifter: input an alignment of a "final" modeling sequence to the raw original one. Then use this class to shift sequence numbering (from original to final).
- ModelOptions: store, read, and write options for modeling, including secondary structure restraints, insertions, and symmetry.
- InsertionRemover: input an alignment and this class will help you remove insertions above a specified length. Can convert columns to insertions when the PDB file is missing residues, and convert insertions back to non-insertions if they're on the SSE modeling "whitelist"