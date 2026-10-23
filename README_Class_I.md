# HLA-PBM-Landscape
R scripts for the analysis of HLA peptide-binding motifs (PBM) to characterize the functional landscape of HLA antigen presentation across class I and class II. This repository is associated with the manuscript "Peptide-binding motifs define the functional architecture of human leukocyte antigens”, submitted to Genome Medicine

## Class I: Peptide_to_motifs

This workflow processes HLA class I peptide datasets retrieved from the MHC Motif Atlas (http://mhcmotifatlas.org/class1, DOI: 10.1093/nar/gkac965), and generates peptide-binding motifs for individual HLA alleles using [MixMHCp (version 2.1)](https://github.com/GfellerLab/MixMHCp).
The list of peptides used in this analysis and their HLA assignment is provided under the folder `0-Peptide_datatsets`.

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

`3_Merging_of_PWMs.R`

This script imports individual motif files and combines their position-specific values into a single matrix, with one column per selected allotype.

**Input:**
- `3_Input.txt`. A tab-delimited file containing two columns `Folder` and `Name` reporting the path to the directory containing the motif files generated in step 2, and the corresponding allotype name (matching the file name, without .txt extension).
- Individual motif files named Name.txt within the corresponding folder.

**Output:**
- `Merged PWMs.txt`. A tab-delimited matrix containing the amino-acid position and symbol in the first two columns, followed by one column for each motif reporting the amino acid frequencies.

**Step 4. Perform Hierarchical clustering**

`4_PBM_Hierarchical_Clustering.R`

This script calculates correlation-based distances between HLA motifs and performs hierarchical clustering using Ward.D2 linkage. It generates dendrograms and exports the distance matrix and computational cluster assignments.

Defining and visualizing PBM groups:
The workflow involves two runs of the hierarchical-clustering script.
Initial clustering and group inspection. Run the script to generate the dendrograms and computational cluster assignments. These clusters provide an initial framework for inspecting similarities between motifs. PBM groups are then reviewed and finalized by manual inspection of the dendrograms and direct comparison of the corresponding motifs. The computational clusters are not necessarily identical to the final PBM groups.
Visualization of finalized PBM groups. Update the Group column in the annotation file with the manually finalized group assignments and rerun the script. The dendrogram labels and circular group ring are then colored according to these assignments.

**Input:**
- `Merged PWMs.txt`, generated in Step 3.
-  `4-Input.txt`, a tab-delimited annotation file containing: 1) `PWM`: motif identifier matching the corresponding column name in the merged matrix. `Group`: PBM group assignment; use 1 for allotypes not yet assigned to a group. `Name`: allele label used in the visualizations.

**Output:**
- Dendograms: 1) Dendrogram.jpeg, 2)Top_down_Dendrogram.jpeg, 3) combined_circular_dendrogram_ring.png, 4) group_ring_only.png
- Distance.txt: pairwise correlation-based distances.
- Alleles_with_clusters.txt: allele annotations with computational cluster assignments.
- Ordered_Alleles.txt: allele annotations ordered according to the dendrogram.

### Requirements and configuration
Set path_output to the desired output directory. When set to NULL, outputs are saved to the current working directory.

## Class I: tSNE Visualization

**Step 5. tSNE analysis**

`5_tSNE.R`
This script generates a two-dimensional t-distributed stochastic neighbor embedding (t-SNE) visualization of HLA peptide-binding motifs. Samples can be colored according to either their assigned peptide-binding motif (PBM) group or their computationally determined clusters.
Before running, configure the plot_type variable (PBM or Cluster) for the color assignment.

**Input:**
- `Merged PWMs.txt`, generated in Step 3.
- `5-Input.txt`, a tab-delimited annotation file containing: 1)`PWM`, identifier of the PWM included in the `Merged PWMs.txt` file, 2) `PBM`, manually annotated PBM group assignement, 3) `Cluster`, computationally assigned clusters from the analysis in step 4.

**Output:**
The script generates a JPEG figure named according to the selected group type: `tsne_plot_by_sample_PBM.jpeg` or `tsne_plot_by_sample_Cluster.jpeg`

## Class I: Immunopeptidome overlap

This workflow uses the peptide analysis performed with mixMHCp (Step 2) to retrieve peptide lists and compare immunopeptidome across allotypes. The first script merges the peptide lists, while the second script performs pairwise overlap comparisons.

**Step 6a. Immunopeptidome overlap analysis**

`6a_Merging_peptide_lists.R`
This script extracts peptide sequences from the PWM-specific output files generated by mixMHCp in step 2, applies predefined filtering criteria, and merges the extracted sequences into a single tab-delimited table for downstream immunopeptidome-overlap analysis.

**Input:**
- `6_Input.txt`, Tab-delimited summary table containing the PWM identifiers (column `PWM`) used to retrieve the corresponding MixMHCp output files (Responsability output).
- `Responsability output`, files generated by MixMHCp in step 2 and containing peptide sequences and their associated PWM assignment.

**Output:**
-`peptides_lists_merged.txt`, Tab-delimited table containing the peptide lists grouped by HLA allotype.

**Step 6b. Immunopeptidome overlap analysis**

`6b_Immunopeptidome_overlap_analysis.R`

This script calculates pairwise peptide overlap between HLA allotypes based on the peptide lists generated in the previous step (6a). It also compares overlap scores within and between predefined groups.

**Input:**
- `6_Input.txt`, same file used in step 6a, including `Name` of the allotypes and their corresponding `Group` assignment.
- `peptides_lists_merged.txt`, generated in step 6a

**Output:**
- `Analyzed peptides.txt`: Peptide sequences included in the analysis, optionally restricted by peptide length.
- `Combinations.txt`: Binary matrix indicating the presence or absence of each unique peptide in each allotype.
- `Overlap.txt`: Pairwise peptide-overlap matrix and the total number of unique peptides for each allotype.
- `Overlap comparison.txt`: Pairwise overlap scores separated into within-group and across-group comparisons.
