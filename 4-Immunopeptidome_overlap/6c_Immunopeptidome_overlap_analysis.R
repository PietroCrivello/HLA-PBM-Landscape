######################################################################################################################################################
########################################################     !!!!!     Settings    !!!!!     #########################################################
######################################################################################################################################################

# Set desired pathway for output files (Folder Name!, e.g. "C:/Users/Anonym/Desktop/Test"). 
# Set to "NULL" for automatic saving in folder of this script.
path_output <- NULL

# Specify whether only peptides with a restricted length should be considered for overlap calculation
# Values for "Enable_restriction": TRUE (if peptides should be restricted) or FALSE
# Values for "Min_pep_length": any integer
# Values for "Max_pep_length": any integer > or = "Min_pep_length"
Enable_restriction <- FALSE
Min_pep_length <- 8
Max_pep_length <- 11


######################################################################################################################################################
###################################################          Import of data          #################################################################
######################################################################################################################################################

#Import data file containing  peptides for all alleles, for which the overlap should be computed
data <- read.table(file.choose(), sep = "\t", header = T, check.names = F)
Overview <- read.table(file.choose(), sep = "\t", header = T, check.names = F)
data <- data[,match(Overview$Name, names(data))]


######################################################################################################################################################
##################################################          Calculate overlap          ###############################################################
######################################################################################################################################################

# If enabled, restrict peptides by length
if(Enable_restriction == TRUE){
  
  for(i in 1:ncol(data)){
    
    TempPeptide <- data[,i][(nchar(data[,i]) >= Min_pep_length) & (nchar(data[,i]) <= Max_pep_length)]
    data[,i] <- c(TempPeptide, rep("", nrow(data) - length(TempPeptide)))
    
  }
  
}

# Export Analyzed peptide sequences
if(is.null(path_output)){
  
  path_output <- getwd()
  
}
write.table(data, paste0(path_output, "/Analyzed peptides.txt"), quote = F, sep = "\t", row.names = F, col.names = T)

# Determine unique peptides across all alleles
seq <-  c()
for(i in 1:ncol(data)){
  
  seq <- c(seq, data[,i])

}
seq <- unique(seq)
seq <- seq[!seq==""]

# Create combination dataframe
dataCombination <- matrix(c(seq, rep(0, ncol(data)*length(seq))), ncol = ncol(data)+1)
dataCombination <- as.data.frame(dataCombination)
names(dataCombination) <- c("Peptides", names(data))
for(i in 1:ncol(data)){
  
  dataCombination[,i+1] <- seq %in% data[,i]
  
}
dataCombination[dataCombination == FALSE] <- 0
dataCombination[dataCombination == TRUE] <- 1

# Export combination table
write.table(dataCombination, paste0(path_output, "/Combinations.txt"), quote = F, sep = "\t", row.names = F, col.names = T)


# Create overlap matrix
dataOverlap <- matrix(rep(0, ncol(data)*(ncol(data)+1)), ncol = ncol(data) + 1)
dataOverlap <- as.data.frame(dataOverlap)
names(dataOverlap) <- c(names(data), "total number of peptides")
row.names(dataOverlap) <- names(data)

for(i in 1:ncol(data)){
  
  for(j in 1:ncol(data)){
    
    dataOverlap[j,i] <- 2*sum(dataCombination[,i+1] == 1 & dataCombination[,j+1] == 1)/(sum(as.numeric(dataCombination[,i+1])) + sum(as.numeric(dataCombination[,j+1])))
    dataOverlap[i,j] <- 2*sum(dataCombination[,i+1] == 1 & dataCombination[,j+1] == 1)/(sum(as.numeric(dataCombination[,i+1])) + sum(as.numeric(dataCombination[,j+1])))
  
  }
  
  dataOverlap[i,ncol(data)+1] <- sum(as.numeric(dataCombination[,i+1]))
}

# Export overlap matrix
write.table(dataOverlap, paste0(path_output, "/Overlap.txt"), quote = F, sep = "\t", row.names = T, col.names = NA)



######################################################################################################################################################
########################################          Comparison of overlaps between groups          #####################################################
######################################################################################################################################################

# Remove total peptide count column for cluster comparison
dataOverlap[, names(dataOverlap) == "total number of peptides"] <- NULL

# Create cluster label matrix (upper triangle only, excluding self-overlaps)
Cluster_df <- matrix(NA, nrow = nrow(dataOverlap), ncol = ncol(dataOverlap))
rownames(Cluster_df) <- rownames(dataOverlap)
colnames(Cluster_df) <- colnames(dataOverlap)

for (i in 1:(nrow(Cluster_df) - 1)) {
  for (j in (i + 1):ncol(Cluster_df)) {
    g1 <- Overview$Group[which(Overview$Name == rownames(Cluster_df)[i])]
    g2 <- Overview$Group[which(Overview$Name == colnames(Cluster_df)[j])]
    Cluster_df[i, j] <- paste0(g1, "-", g2)
  }
}

# --------------------------------------------
# Extract within-group overlaps (no self-overlap)
# --------------------------------------------
temp_within <- c()
unique_groups <- sort(unique(Overview$Group))
for (grp in unique_groups) {
  # Find all pairwise combinations *within* the same group (excluding self-pairs)
  members <- Overview$Name[Overview$Group == grp]
  if (length(members) < 2) next
  combos <- combn(members, 2)
  for (k in 1:ncol(combos)) {
    i <- combos[1, k]
    j <- combos[2, k]
    if (!is.na(dataOverlap[i, j])) {
      temp_within <- c(temp_within, dataOverlap[i, j])
    }
  }
}

# --------------------------------------------
# Extract across-group overlaps
# --------------------------------------------
temp_across <- c()
group_combinations <- combn(unique_groups, 2)
for (k in 1:ncol(group_combinations)) {
  g1 <- group_combinations[1, k]
  g2 <- group_combinations[2, k]
  
  members1 <- Overview$Name[Overview$Group == g1]
  members2 <- Overview$Name[Overview$Group == g2]
  
  for (m1 in members1) {
    for (m2 in members2) {
      # Ensure we're reading only from upper triangle
      if (!is.na(dataOverlap[m1, m2])) {
        temp_across <- c(temp_across, dataOverlap[m1, m2])
      } else if (!is.na(dataOverlap[m2, m1])) {
        temp_across <- c(temp_across, dataOverlap[m2, m1])
      }
    }
  }
}

# --------------------------------------------
# Combine results into output table
# --------------------------------------------
Output_Cluster <- as.data.frame(matrix("", ncol = 2, nrow = max(length(temp_within), length(temp_across))))
Output_Cluster[, 1] <- c(temp_within, rep("", max(0, length(temp_across) - length(temp_within))))
Output_Cluster[, 2] <- c(temp_across, rep("", max(0, length(temp_within) - length(temp_across))))
names(Output_Cluster) <- c("within groups", "across groups")

# Export overlap comparison table
write.table(Output_Cluster, file = paste0(path_output, "/Overlap comparison.txt"),
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
