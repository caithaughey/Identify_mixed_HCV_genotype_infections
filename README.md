# Detect mixed genotype HCV infections in PacBio sequence samples
HCV samples can contain so called mixed infections of more than one HCV genotype. This Python script can detect mixed infections in HCV samples sequenced with Pacific Biosciences (PacBio). It utilizes NCBI’s local command-line version of BLAST, multiple sequence alignment using MUSCLE, and Python module AlignInfo. For each sample the reads are matched against HCV genotype references from NCBI and grouped by which reference the read has closest similarity to. The reads are then filtered, aligned and a consensus sequence is created for the sample.

## Requirements

* Python >3
* BLAST+ command line tools
* MUSCLE (3.8.31)
