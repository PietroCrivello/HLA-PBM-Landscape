# HLA-PBM-Landscape
R scripts for the analysis of HLA peptide-binding motifs (PBM) to characterize the functional landscape of HLA antigen presentation across class I and class II. This repository is associated with the manuscript "Peptide-binding motifs define the functional architecture of human leukocyte antigens”, submitted to Genome Medicine

## Class I: Peptide_to_motifs

This workflow processes HLA class I peptide datasets retrieved from the MHC Motif Atlas (), and generates peptide-binding motifs for individual HLA alleles using MixMHCp ().

### Workflow

**Step 1. Prepare allele-specific peptide lists**

`01_Prepare_Allele_Specific_Peptide_Lists.R`

This script takes a tab-delimited master file containing peptide sequences and their corresponding HLA allele assignments. It removes the peptide-length column and duplicate rows, separates peptides by allele, and exports individual CSV files containing the peptide sequences for each allele.

The resulting files are named using the format `<allele>_unnested.csv` and are used as input for the next step.

**Input:** Master TSV file containing at least the columns `Peptide`, `Allele`, `Locus`, and `length`.

**Output:** One CSV file per HLA allele, containing peptide sequences without column headers.

**Step 2. Generate peptide-binding motifs using MixMHCp**

`02_Run_MixMHCp.R`

This script processes the allele-specific peptide-list files generated in Step 1. Each file is processed independently using MixMHCp to generate peptide-binding motif predictions.

**Input:** CSV files generated in Step 1.

**Output:** MixMHCp output files for each input peptide list.

### Requirements

- R and the R packages `readr` and `tidyverse`.
- MixMHCp installed and accessible to the script. Follow the installation instructions provided by the [MixMHCp developers](https://github.com/GfellerLab/MoDec) only if that is the actual software distribution you used; otherwise, replace this link with the appropriate MixMHCp documentation.

### Notes

The peptide lists must already contain allele assignments. The preprocessing script reorganizes the input data but does not infer HLA restriction or generate peptide sequences.


