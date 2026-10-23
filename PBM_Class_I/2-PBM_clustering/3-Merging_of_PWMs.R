######################################################################################################################################################
########################################################     !!!!!     Settings    !!!!!     #########################################################
######################################################################################################################################################

# Set desired pathway for output files (Folder Name!, e.g. "C:/Users/Anonym/Desktop/Test"). 
# Set to "NULL" to save in the current working directory
path_output <- NULL



######################################################################################################################################################
###################################################          Import of data          #################################################################
######################################################################################################################################################

# Import of tab-separated file containing the path to the desired motif and the name of the allele
destClusters <- read.table(file.choose(), sep = "\t", header = T, check.names = F)



######################################################################################################################################################
########################################          Merging of position-weight matrices          #######################################################
######################################################################################################################################################

# Generation of output file containing in the first two columns the AA position and the AA of the corresponding values
Output <- data.frame("AA Position" = rep(1:9, times = 23),
                     "AA" = rep(c("A", "C",	"D",	"E",	"F",	"G",	"H",	"I",	"K",	"L",	"M",	"N",	"P",	"Q",	"R",	"S",	"T",	"V",	"W",	"Y", "s", "t", "y"), each = 9),
                     check.names = F)


# Filling of the output dataframe
for(i in 1:nrow(destClusters)){
  
  # Import of the position-weight matrix file of the desired peptide binding motifs
  TempFile <- readLines(paste0(destClusters$Folder[i], "/", 
                               destClusters$Name[i], ".txt"))
  
  # Transformation of position-weight matrix into a vector
  Tempdata <- c()
  for(j in 2:24){
    Tempdata <- c(Tempdata,  strsplit(TempFile[j], split="\t")[[1]][2:10])
  }
  
  # add the motif vector as a new column to the output data frame
  Output <- cbind(Output, Tempdata)
  names(Output)[i+2] <- destClusters$Name[i]
}

# Export of the output dataframe
if(is.null(path_output)){
  
  path_output <- getwd()
  
}
write.table(Output, paste0(path_output, "/Merged PWMs.txt"), sep = "\t", row.names = F, col.names = T, quote = F)




