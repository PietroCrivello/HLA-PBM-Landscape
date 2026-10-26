library(amap)
library(data.table)
library(circlize)
library(dendextend)
library(stats)

######################################################################################################################################################
########################################################     !!!!!     Settings    !!!!!     #########################################################
######################################################################################################################################################

path_output <- NULL
Post_order <- NULL
label_size <- 1.3

# Fixed list of colors for groups (up to 11)
group_colors <- c(
  "#999999",  # Grey
  "#2A5B88",  # Dark Blue
  "#317132",  # Dark Green
  "#B35900",  # Dark Orange
  "#6B3368",  # Dark Purple
  "#B2B229",  # Dark Yellow
  "#6BAED6",  # instead of "#377EB8",  # Blue
  "#E41A1C",  # Red
  "#666666",  # Dark Grey
  "#4D8372",  # Dark Teal
  "#B35745"   # Dark Salmon
)

#group_colors <- c(
# "#6B4620",  # Dark Brown
# "#A25487",  # Dark Pink
# "#4DAF4A",  # Green
# "#984EA3",  # Purple
# "#FF7F00",  # Orange
# "#FFFF33",  # Yellow
# "#A65628",  # Brown
# "#F781BF",  # Pink
# "#A31414",  # Dark Red
# "#66C2A5",  # Teal
# "#FC8D62"   # Salmon
#)

######################################################################################################################################################
###################################################          Import of data          #################################################################
######################################################################################################################################################

PWM <- read.table(file.choose(), sep = "\t", check.names = F, header = T)
PWM <- PWM[,3:ncol(PWM)]
PWM <- PWM[,order(names(PWM))]
t_PWM <- transpose(PWM)
colnames(t_PWM) <- rownames(PWM)
rownames(t_PWM) <- colnames(PWM)

Name_alleles <- read.table(file.choose(), sep = "\t", header = T, check.names = F)
Name_alleles <- Name_alleles[order(Name_alleles$Name),]

######################################################################################################################################################
##############################################          Hierarchical Clustering          #############################################################
######################################################################################################################################################

PWM_dist <- Dist(t_PWM, method = "correlation")
PWM_clust <- hclust(PWM_dist, method = "ward.D2")
PWM_clust$labels <- Name_alleles$Name
PWM_dend <- as.dendrogram(PWM_clust)

Name_alleles$Group[Name_alleles$Group == "-"] <- 0
Name_alleles$Group <- as.numeric(Name_alleles$Group)
num_groups <- max(Name_alleles$Group)

# Check if enough colors are provided
if(num_groups + 1 > length(group_colors)) {
  stop("Number of groups exceeds the number of provided colors. Please provide more colors.")
}

# Assign colors from the fixed list to groups
Name_alleles$Color <- group_colors[Name_alleles$Group + 1]

col_to_use <- Name_alleles$Color[match(labels(PWM_dend), Name_alleles$Name)]
PWM_dend <- color_labels(PWM_dend, col = col_to_use)
PWM_dend <- set(PWM_dend, "labels_cex", label_size)
PWM_dend_circular <- PWM_dend

if(is.null(path_output)){
  path_output <- getwd()
}

# Save horizontal dendrogram
jpeg(paste0(path_output, "/Dendrogram.jpeg"), width = 15, height = 30, unit = "cm", res = 500)
par(mar = c(5,5,5,20))
if(is.null(Post_order)){
  plot(PWM_dend, main = "", sub = "", axes = F, xlab = "", ylab = "", horiz = T)
}else{
  PWM_dend %>% rotate(Post_order) %>% plot(main = "", sub = "", axes = F, xlab = "", ylab = "", horiz = T)
}
dev.off()

# Save top-down dendrogram
PWM_dend_inverted <- rev(PWM_dend)
jpeg(paste0(path_output, "/Top_down_Dendrogram.jpeg"), width = 30, height = 15, unit = "cm", res = 500)
par(mar = c(20,5,5,5))
if(is.null(Post_order)){
  plot(PWM_dend_inverted, main = "", sub = "", axes = F, xlab = "", ylab = "", horiz = F)
}else{
  PWM_dend_inverted %>% rotate(Post_order) %>% plot(main = "", sub = "", axes = F, xlab = "", ylab = "", horiz = F)
}
dev.off()

######################################################################################################################################################
##############################################          Combined Circular Dendrogram + Group Ring         ###########################################
######################################################################################################################################################

png(paste0(path_output, "/combined_circular_dendrogram_ring.png"),
    width = 2000, height = 2000, res = 300, bg = "transparent")  # transparent bg

par(mar = c(2, 2, 2, 2))
par(family = "Arial")
par(cex = 0.7)
circos.clear()
circos.par(track.margin = c(0.01, 0.01), start.degree = 90, gap.degree = 0.5)

ordered_labels <- labels(PWM_dend_circular)
group_info <- Name_alleles[match(ordered_labels, Name_alleles$Name), ]
factors <- ordered_labels

# Initialize circular layout
circos.initialize(factors = factors, xlim = cbind(rep(0, length(factors)), rep(1, length(factors))))

# Track 1: Group ring (outer ring)
circos.track(ylim = c(0, 1), track.height = 0.05, bg.border = NA,
             panel.fun = function(x, y) {
               sector_index <- CELL_META$sector.index
               color <- group_info$Color[group_info$Name == sector_index]
               circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1,
                           col = color, border = NA)
             })

# Track 2: Dendrogram
circlize_dendrogram(PWM_dend_circular,
                    dend_track_height = 0.3,
                    labels_track_height = 0.4,
                    dend_line_width = 2)

dev.off()
circos.clear()

png("group_ring_only.png", width = 2000, height = 2000, res = 300, bg = "transparent")

circos.clear()
circos.par(start.degree = 90, gap.degree = 0.5)

# Initialize sectors (one per factor)
circos.initialize(factors = factors, xlim = cbind(rep(0, length(factors)), rep(1, length(factors))))

# Draw the colored group ring as one track
circos.track(ylim = c(0, 1), track.height = 0.05, bg.border = NA,
             panel.fun = function(x, y) {
               sector_index <- CELL_META$sector.index
               col_sector <- group_info$Color[group_info$Name == sector_index]
               circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1,
                           col = col_sector, border = NA)
             })

dev.off()
circos.clear()

######################################################################################################################################################
##############################################          Export distance matrix          ##############################################################
######################################################################################################################################################

dist_output <- matrix(NA, ncol = ncol(PWM), nrow = ncol(PWM))
dist_temp <- PWM_dist
for(i in 1:(ncol(PWM)-1)){
  dist_output[(i+1):ncol(PWM),i] <- dist_temp[1:(ncol(PWM)-i)]
  dist_output[i,(i+1):ncol(PWM)] <- dist_temp[1:(ncol(PWM)-i)]
  dist_temp <- dist_temp[(ncol(PWM)-i+1):length(dist_temp)]
}
for(i in 1:ncol(PWM)){
  dist_output[i,i] <- 0
}
dist_output <- as.data.frame(dist_output)
names(dist_output) <- Name_alleles$Name
row.names(dist_output) <- Name_alleles$Name

write.table(dist_output, file = paste0(path_output, "/Distance.txt"), sep = "\t", quote = F)

order <- order.dendrogram(PWM_dend)
ordered_items <- Name_alleles[order, ]

num_clusters <- 7
clusters <- cutree(PWM_dend, k = num_clusters)
Names_alleles_with_clusters <- cbind(Name_alleles, cluster = clusters)

write.table(Names_alleles_with_clusters, file = paste0(path_output, "/Alleles_with_clusters.txt"), sep = "\t", quote = F)
write.table(ordered_items, file = paste0(path_output, "/Ordered_Alleles.txt"), sep = "\t", quote = F)
