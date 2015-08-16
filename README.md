# gtusc
##### sequence tools file for some abstraction
SequenceShifter:
         __init__: read in an alignment-to-orig file.
         get_adjusted(raw sequence spec): return some sequence as an adjusted sequence
ModelOptions:
        just store SSE info, insertions, symmetry, domain restraints
        access functions like get_sses, get_symmetry_breaks, ...
        read/write command (possibly use YAML?)

##### creating alignments for GCP2 and GCP3
>>> (get TUBGCP family from uniprot. align using PROMALS3D)
>>> (create SSE files: mark what should be forced to SSE or just kept in model)
>>> (extract alignment to human GCP4, but limit GCP4 to what's in the PDB file)
>>> python ../src/align_to_pdb.py TUBGCP_promals.fasta ../data/gcp4_pos2.pdb GCP4_HUMAN gcp4_pdblim

for each of GCP2/GCP3:
    "orig" alignment: remove insertions
    "sse"  alignment: preserve requested SSEs, remove remaining insertions.
    insertion removal: create both the new alignment AND an alignment to the original model sequence.
                       do this with a class so that it's transparent!

##### homology modeling (GCP2+3 with restraints inc. symmetry)
combine final alignments of GCP2/GCP3, two copies each. (plus g-tubulin?)
create a ModelOptions file for modeller containing:
       insertion points and lengths, for creating distance restraints
       symmetry breakpoints
       domain restraints (if g-tubulin?)
       all SSEs
       (use SequenceShifter to help this out!)
build model

##### EM prep 1: fixing chain names and res indexes. adding ytubulin + spc110.

##### EM prep 2: generate psf, add solvent, symmetry restraints

##### flexible fitting

##### EM cleanup: stripping hydrogen and solvent
