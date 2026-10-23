# HLA-PBM-Landscape
R scripts for the analysis of HLA peptide-binding motifs (PBM) to characterize the functional landscape of HLA antigen presentation across class I and class II. This repository is associated with the manuscript "Peptide-binding motifs define the functional architecture of human leukocyte antigens”, submitted to Genome Medicine

## Class I: Peptide_to_motifs

This workflow processes HLA class I peptide datasets retrieved from the MHC Motif Atlas (http://mhcmotifatlas.org/class1, DOI: 10.1093/nar/gkac965), and generates peptide-binding motifs for individual HLA alleles using [MixMHCp (version 2.1)](https://github.com/GfellerLab/MixMHCp).

### Workflow

**Step 1. Prepare allele-specific peptide lists**

`1_Prepare_Allele_Specific_Peptide_Lists.R`

This script takes a tab-delimited master file containing peptide sequences and their corresponding HLA allele assignments. It removes duplicate rows and exports peptide sequences by allele to individual CSV files.

The resulting output are files named using the format `<allele>_unnested.csv`, which are used as input for the next step.

**Input:** Tab-delimited master file containing the columns `Peptide`, `Allele`, `Locus`.

**Output:** One CSV file per HLA allele, containing peptide sequences without column headers.

**Step 2. Generate peptide-binding motifs using MixMHCp**

`2_Run_MixMHCp.R`

This script processes the allele-specific peptide lists generated in Step 1. Each input file is processed independently using MixMHCp to generate peptide-binding motif outputs.

**Input:** CSV files generated in Step 1.

**Output:** MixMHCp output files for each input peptide list.

### Requirements

- MixMHCp must be installed separately. Follow the installation instructions provided by the [MixMHCp developers](https://github.com/GfellerLab/MixMHCp).
Before running the script, configure the paths to the MixMHCp executable, input folder and output folder. The script also requires the bias_phospho.txt file included in the MixMHCp installation, with its path configured accordingly.
 
### Notes
The peptide lists must already contain allele assignments. This workflow reorganizes the input data but does not infer HLA restriction from unassigned master peptide list.

## Class I: PBM_clustering

This workflow combines position weight matrices (PWMs) describing the PBM of each class I allotype, performs hierarchical clustering based on motif similarity, and visualizes the relationships between allotypes and their PBM.

### Workflow

**Step 3. Merging PWM**

`1_Merging_of_PWMs.R`

This script imports individual motif files and combines their position-specific values into a single matrix, with one column per selected allotype.

**Input:**
- `1_Input.txt`. A tab-delimited file containing two columns `Folder` and `Name` reporting the path to the directory containing the motif files generated in step 2, and the corresponding allotype name (matching the file name, without .txt extension).
- Individual motif files named Name.txt within the corresponding folder.

**Output:**
- Merged PWMs.txt. A tab-delimited matrix containing the amino-acid position and symbol in the first two columns, followed by one column for each motif reporting the amino acid frequencies.

**Step 4. Merging PWM**
