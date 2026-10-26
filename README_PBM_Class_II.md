# HLA-PBM-Landscape
R scripts for the analysis of HLA peptide-binding motifs (PBM) to characterize the functional landscape of HLA antigen presentation across class I and class II. This repository is associated with the manuscript "Peptide-binding motifs define the functional architecture of human leukocyte antigens”, submitted to Genome Medicine

## Class II: Peptide_to_motifs

This workflow analyzes HLA class II peptide datasets derived from the supplementary materials of Stražar et al. (DOI: 10.1016/j.immuni.2023.05.009) and Racle et al. (DOI: 10.1016/j.immuni.2023.03.009), as well as the [NetMHCIIpan-4.3 webserver](https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.3/) (DOI: 10.1126/sciadv.adj6367). The scripts generate peptide-binding motifs for individual HLA alleles using [MoDec v1.2](https://github.com/GfellerLab/MoDec).

The peptide lists used in the analyses, together with their HLA assignments, are provided in the folder `0-Class_II_datasets`. These files correspond to the input datasets used to generate the results presented in the manuscript.

### Workflow

**Step 1. Prepare allele-specific peptide lists**

`1_Peptides_to_allele_peptide_lists.R`

This script takes a tab-delimited master file containing peptide sequences and their corresponding HLA allele assignments (datasets in `0-Class_II_datasets`). It removes duplicate rows and exports peptide sequences by allele to individual CSV files.

The resulting output are files named using the format `<allele>_unnested.csv`, which are used as input for the next step.

**Input:** Tab-delimited master file containing the columns `Peptide`, `Allele`.

**Output:** One CSV file per HLA allele, containing peptide sequences without column headers.

**Step 2. Generate peptide-binding motifs using MoDec**

`2_Run_MODEC.R`

This script runs MoDec v1.2 on the allele-specific peptide lists generated in Step 1. Each input file is processed independently, and the output is saved in separate output folders.

**Input:**
- CSV files named <allele>_unnested.csv, generated in Step 1.

**Analysis parameters:**
--Kmax 5: sets the maximum number of clusters to 5.
--MHC2: enables HLA class II analysis.
--nruns 10: specifies 10 runs.
--makeReport: requests report generation.
--outAdd: enables the corresponding MoDec output option.

**Output:**
- MODEC output files saved in the specified output directories. The output filename is based on the input filename, with _output.txt appended to its basename.

**Parallel processing:** 
The script uses the R packages future and furrr to process multiple input files in parallel, with up to six workers configured by plan(multisession, workers = 6).

**Requirements:** 
- MoDec v1.2 executable (MoDec.exe), to be installed separately. Follow the installation instructions provided by the [MoDec developers](https://github.com/GfellerLab/MoDec).

Before running the scripts, update the input and output directory paths in 2_Run_MODEC.R to match your local environment. Ensure that the MoDec executable is accessible to the script.

## Class II: PBM_clustering

This workflow combines position weight matrices (PWMs) describing the PBM of each class II allotype, performs hierarchical clustering based on motif similarity, and visualizes the relationships between allotypes and their PBM.

### Workflow

**Step 3a. Collecting selected PWM**

`3a_Collect_selected_PWMs.R`

This script collects selected PWMs from MoDec output folders for HLA class II heterodimers. The input table is generated following manual review of the MoDec outputs and specifies the motifs selected for downstream analysis.

For each row, the script uses the column `MoDec_output`, and the `N of motifs` (K), and the selected motif index (i) to identify the corresponding PWM file. It copies the file to a dedicated output directory and renames it using the corresponding ID from the input table.

**Input:**

- `3a_PWM_selection.txt`: A tab-delimited table generated after manual inspection of the MoDec outputs, and containing the following columns: 1) `ID`: identifier used to name the extracted PWM file, 2) `MoDec_output`: identifies the corresponding MoDec output folder, 3) `N_of_motifs`: corresponding to the parameter K included in the MoDec analysis, 4) `Selected_for_analysis`: reporting the selected index (i) of the PWM.

**Output:**

- A directory containing the selected PWM files, named according to the ID column in the input file. The output directory is specified by the output_folder variable in the script.

**Step 3b. Merging PWM**

`3b_Merging_of_PWMs.R`

This script imports individual motif files and combines their position-specific values into a single matrix, with one column per selected allotype.

**Input:**
- `3b_Input.txt`. A tab-delimited file containing the following columns: 1) `Folder` with the path to the folder containing the selected PWM files, 2) `Name` with the name of the selected PWM file (without .txt extension). This correspond to the `ID` column from the previous input in Step 3a.
- Individual motif files named Name.txt within the corresponding folder.

**Output:**
- `Merged PWMs.txt`. A tab-delimited matrix containing the amino-acid position and symbol in the first two columns, followed by one column for each selected PWM containing the corresponding position-specific values.

**Step 4. Perform Hierarchical clustering**

`4_PBM_Hierarchical_Clustering.R`

This script calculates correlation-based distances between HLA motifs and performs hierarchical clustering using Ward.D2 linkage. It generates dendrograms and exports the distance matrix and computational cluster assignments.

Defining and visualizing PBM groups:
The workflow involves two runs of the hierarchical-clustering script.
Initial clustering and group inspection. Run the script to generate the dendrograms and computational cluster assignments. These clusters provide an initial framework for inspecting similarities between motifs. PBM groups are then reviewed and finalized by manual inspection of the dendrograms and direct comparison of the corresponding motifs. The computational clusters are not necessarily identical to the final PBM groups.
Visualization of finalized PBM groups. Update the Group column in the annotation file with the manually finalized group assignments and rerun the script. The dendrogram labels and circular group ring are then colored according to these assignments.

**Input:**
- `Merged PWMs.txt`, generated in Step 3.
-  `4-Input.txt`, a tab-delimited annotation file containing: 1) `PWM`: motif identifier matching the corresponding column name in the merged matrix. 2) `Group`: PBM group assignment; use "-" for allotypes not yet assigned to a group. 3) `Name`: allele label used in the visualizations.

**Output:**
- Dendograms: 1) Dendrogram.jpeg, 2)Top_down_Dendrogram.jpeg, 3) combined_circular_dendrogram_ring.png, 4) group_ring_only.png
- Distance.txt: pairwise correlation-based distances.
- Alleles_with_clusters.txt: allele annotations with computational cluster assignments.
- Ordered_Alleles.txt: allele annotations ordered according to the dendrogram.

### Requirements and configuration
Set path_output to the desired output directory. When set to NULL, outputs are saved to the current working directory.

## Class II: tSNE Visualization

**Step 5. tSNE analysis**

`5_tSNE.R`
This script generates a two-dimensional t-distributed stochastic neighbor embedding (t-SNE) visualization of HLA peptide-binding motifs. Samples can be colored according to either their assigned peptide-binding motif (PBM) group or their computationally determined clusters.
Before running, configure the plot_type variable (PBM or Cluster) for the color assignment.

**Input:**
- `Merged PWMs.txt`, generated in Step 3.
- `5-Input.txt`, a tab-delimited annotation file containing: 1)`PWM`, identifier of the PWM included in the `Merged PWMs.txt` file, 2) `PBM`, manually annotated PBM group assignement, 3) `Cluster`, computationally assigned clusters from the analysis in step 4.

Important: The rows in 5-Input.txt must be in exactly the same order as the PWM columns in Merged PWMs.txt, after the first two columns have been removed. The script currently uses the row order to associate annotations with motifs; it does not match entries using the PWM column.

**Output:**
The script generates a JPEG and SVG figures named according to the selected group type.

## Class II: Immunopeptidome overlap

This workflow uses the peptide analysis performed with MoDec (Step 2) to retrieve peptide lists and compare immunopeptidome across allotypes. The first scripts (6a and 6b) collect and merges the peptide lists, while the third script (6c) performs pairwise overlap comparisons.

**Step 6a. Retrieve peptide-level MoDec output files**

`6a_Retrieve_peptide_files.R`

This script retrieves the peptide-level output files corresponding to the selected MoDec motifs. Using the same input table as in Step 3a,  it locates the corresponding bestPepResp_K*.txt file in each MoDec output folder. The file is copied to a designated output directory and renamed using the corresponding ID.

**Input:**
- `3a_PWM_selection.txt`: A tab-delimited table generated after manual inspection of the MoDec outputs, and containing the following columns: 1) `ID`: identifier used to name the extracted PWM file, 2) `MoDec_output`: identifies the corresponding MoDec output folder, 3) `N_of_motifs`: corresponding to the parameter K included in the MoDec analysis, 4) `Selected_for_analysis`: reporting the selected index (i) of the PWM.
- The corresponding MoDec output folders resulting from Step 2.

**Output:**
- One retrieved and renamed .txt file per successfully located output file, containing the peptide-level results and motif assignments from the MoDec analysis.

Note: The script uses N_of_motifs to determine which output file to retrieve. The Selected_for_analysis column is read but is not directly used by the script.

**Step 6b. Merge peptide lists from selected MoDec motifs**

`6b_Merging_peptide_lists.R`

This script extracts and merges peptide sequences from the MoDec output files generated in step 2, and corresponding to the manually selected motifs recorded in the Step 3a input table. For each allotype, sequences are retained according to their motif assignment in the Motif_1st column. The extracted sequences are grouped according to a string derived from their ID (string length set to 10, to be adjusted for longer ID) and combined into a single tab-delimited table.

**Input:**
- `3a_PWM_selection.txt`
-  The retrieved MoDec output files from step 6a. Each file must contain the Motif_1st column and the sequence column to be extracted.

**Output:**
-`peptides_lists_merged.txt`, Tab-delimited table containing the peptide lists grouped by HLA allotype.

**Step 6c. Immunopeptidome overlap analysis**

`6c_Immunopeptidome_overlap_analysis.R`

This script calculates pairwise peptide overlap between HLA allotypes based on the merged peptide lists generated in the previous step (6a). It also compares overlap scores within and between assigned groups.

Pairwise overlap is calculated using the Sørensen–Dice coefficient. Peptide-length restriction can optionally be enabled to include only peptides within a specified length range.

**Input:**
- `peptide_lists_merged.txt` generated in step 6b.
- `6c_Input.txt`, overview table containing the following columns: 1) `Name`: allotype identifiers matching the column names in the merged peptide table. 2) `Group`: group assignments used to classify pairwise comparisons as within-group or across-group.

Configuration:
- path_output: output directory. Set to NULL to save the output files in the current working directory.
- Enable_restriction: set to TRUE to restrict the analysis by peptide length, or FALSE to include peptides of all lengths.
- Min_pep_length and Max_pep_length: minimum and maximum peptide lengths to include when length restriction is enabled.

**Output:**
- `Analyzed peptides.txt`: Peptide sequences included in the analysis, optionally restricted by peptide length.
- `Combinations.txt`: Binary matrix indicating the presence or absence of each unique peptide in each allotype.
- `Overlap.txt`: Pairwise peptide-overlap matrix and the total number of unique peptides for each allotype.
- `Overlap comparison.txt`: Pairwise overlap scores separated into within-group and across-group comparisons.
