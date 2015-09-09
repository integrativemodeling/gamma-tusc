# Overview
This module creates models of the gamma tubulin small complex using MODELLER and IMP. It has config files to do flexible fitting in MDFF. And it has a set of analysis and plotting tools.

# Creating alignments for GCP2 and GCP3
## Some manual steps to prepare alignments
- get TUBGCP family from uniprot. align using PROMALS3D
- create "options" files with to mark what should be forced to SSE (and kept in model)
- get representative proteins of tubulin family (HHPRED at ss15, plus gamma tubulin human and yeast) and align using PROMALS3D

## Extract alignment to human GCP4, but limit GCP4 to what's in the PDB file
`> python src/align_to_pdb.py alignments/TUBGCP_promals.fasta data/gcp4_pos2.pdb GCP4_HUMAN alignments/gcp4_pdblim`
`> python src/align_to_pdb.py alignments/gtub_seed_promals.fasta data/human_ytub_3CB2.pdb GTUB_HUMAN alignments/gtub_pdblim`

## Create a few types of alignments
- "orig" alignments - no SSE additions
`> python src/remove_insertions.py alignments/TUBGCP_promals.fasta GCP4_HUMAN GCP2_YEAST 10000 alignments/gcp2/orig -p alignments/gcp4_pdblim.pir`
`> python src/remove_insertions.py alignments/TUBGCP_promals.fasta GCP4_HUMAN GCP3_YEAST 10000 alignments/gcp3/orig -p alignments/gcp4_pdblim.pir`

- "ins5" alignment: preserve requested SSEs, remove remaining insertions.
`> python src/remove_insertions.py alignments/TUBGCP_promals.fasta GCP4_HUMAN GCP2_YEAST 5 alignments/gcp2/ins5 -p alignments/gcp4_pdblim.pir -s alignments/gcp2/GCP2.options -o alignments/gcp2/ins5_to_orig`
`> python src/remove_insertions.py alignments/TUBGCP_promals.fasta GCP4_HUMAN GCP3_YEAST 5 alignments/gcp3/ins5 -p alignments/gcp4_pdblim.pir -s alignments/gcp3/GCP3.options -o alignments/gcp3/ins5_to_orig`

- the gtub alignment w/pdb and insertion removal but no SSEs
`> python src/remove_insertions.py alignments/gtub_rp15_promals.fasta GTUB_HUMAN GTUB_YEAST 6 alignments/gtub/gtub_rp15_promals_ins6 -p alignments/gtub_pdblim.pir`

### Combine the various alignments to setup modeling
`> create_combo.py`

# Homology modeling
- create a ModelOptions file for modeller containing:
  - insertion points and lengths, for creating distance restraints
  - symmetry breakpoints
  - domain restraints (if g-tubulin?)
  - all SSEs
  -(use SequenceShifter to help this out!)
- build model

# EM flexible fitting
- fixing chain names and res indexes. adding ytubulin + spc110.
- generate psf, add SSE and symmetry restraints
- run in MDFF
- run cleanup files

# Analysis
- finding best models
- calculating RMSF

# Open state fitting
- setting up best models


# Some useful classes defined in the tools file
SequenceShifter:
         __init__: read in an alignment-to-orig file.
         get_adjusted(raw sequence spec): return some sequence as an adjusted sequence
ModelOptions:
        just store SSE info, insertions, symmetry, domain restraints
        access functions like get_sses, get_symmetry_breaks, ...
        read/write command (possibly use YAML?)
