# making heatmap of immune profiling differences
# written by GK Fragiadakis
# February 20th, 2020
# all ppts included except 8000 and 8012 (see methods)
# set working directory to repo

library(devtools)
library(tidyverse)

feature_set <- read_csv(file = "data/combined_immune/combined_immune_feature_differences_set.csv")
feature_key <- read_csv(file = "data/combined_immune/immune_feature_key.csv")

# filter down to fiber only

fiber_set <- dplyr::filter(feature_set, Group == "Fiber")

# filter to only T1_T6 and remove that label

fiber_T6 <- fiber_set %>% 
  select(-contains('T1_T4'), -contains('T1_T5'))

colnames(fiber_T6) <- gsub("_T1_T6", "", colnames(fiber_T6))

# sort columns based on feature type and remove superfluous phosphoflow features

sorted_feature_key <- feature_key %>% 
  dplyr::filter(Time == "T6") %>% 
  arrange(Type)

  
sorted_feature_set <- select(fiber_T6, Participant, Group, !!(sorted_feature_key$Master_feature))

# function rescale onto -1 to 1 while maintaining directionality:

one2one <- function(x){
  
  newX <- x/max(abs(x), na.rm = TRUE)
  
  return(newX)
  
}

rescaled_values <- sorted_feature_set %>% dplyr::select(-Participant, -Group) %>% apply(., 2, one2one)

rescaled_feature_set <- sorted_feature_set %>% select(Participant) %>% bind_cols(., data.frame(rescaled_values))

rownames(rescaled_feature_set) <- rescaled_feature_set$Participant

rescaled_feature_set <- select(rescaled_feature_set, -Participant)

# color bars: 

Freq <- dplyr::filter(sorted_feature_key,Type == "Cell frequencies")
Cyt <- dplyr::filter(sorted_feature_key,Type == "Cytokine")
EndSig <- dplyr::filter(sorted_feature_key,Type == "Endogenous signaling")
SigCap <- dplyr::filter(sorted_feature_key,Type == "Signaling capacity")

immune_feature_colors = c(rep("red",dim(Freq)[1]),rep("blue",dim(Cyt)[1]),rep("cyan",dim(EndSig)[1]),rep("pink",dim(SigCap)[1]))
ppt_feature_colors <- c(rep("dark green", dim(rescaled_feature_set)[1]))

library(RColorBrewer)

main_title="Immune Feature Difference from Baseline to Maintenance: Fiber Only"
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
myColors_RedBlue <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(20))
fiberheatmap <-heatmap.3(
  x=as.matrix(rescaled_feature_set), 
  na.rm = TRUE, scale="none",  dendrogram = "row",
  Colv=colnames(rescaled_feature_set),
  Rowv = TRUE,
  ColSideColors=as.matrix(immune_feature_colors),
  RowSideColors = t(as.matrix(ppt_feature_colors)),
  symbreaks=FALSE, key=TRUE, symkey=FALSE,
  density.info="none", trace="none", 
  main=NULL, labCol=FALSE,
  cexRow=1, col=myColors_RedBlue,
  ColSideColorsSize=1, RowSideColorsSize=1, 
  KeyValueName="1 to 1 scale")


# new groups based on clustering (annotated from heatmap)
# need to check, ones not commented out are correct
#group1 <- c("8002","8001","8013","8009","8007","8006") # just missing 8000
#group2 <- c("8003","8022","8041", "8017","8029","8023") # just missing 8012
#group3 <- c("8038","8035","8037","8036","8018","8039") # now has 8036 (don't know where it was before)
group3 <- rownames(rescaled_feature_set)[fiberheatmap$rowInd][1:6]
group2 <- rownames(rescaled_feature_set)[fiberheatmap$rowInd][7:12]
group1 <- rownames(rescaled_feature_set)[fiberheatmap$rowInd][13:18]

group_df <- data.frame(group1 = group1, group2 = group2, group3 = group3) %>% gather(., key = Immune_group, value = Participant, group1:group3)
write_csv(group_df, "data/combined_immune/fiber_immune_groups.csv")



# experimenting with distances: 

main_title="Immune Feature Difference from Baseline to Maintenance for Fiber participants"
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
myColors_RedBlue <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(20))
fiberheatmapTEST <-heatmap.3(
  x=as.matrix(rescaled_feature_set), 
  distfun = function(x) dist(x,method = 'euclidean'),
  hclustfun = function(x) hclust(x, method = "ward.D"),
  na.rm = TRUE, scale="none",  dendrogram = "row",
  Colv=colnames(rescaled_feature_set),
  Rowv = TRUE,
  ColSideColors=as.matrix(immune_feature_colors),
  RowSideColors = t(as.matrix(ppt_feature_colors)),
  symbreaks=FALSE, key=TRUE, symkey=FALSE,
  density.info="none", trace="none", 
  main=NULL, labCol=FALSE,
  cexRow=1, col=myColors_RedBlue,
  ColSideColorsSize=1, RowSideColorsSize=1, 
  KeyValueName="1 to 1 scale")

group3 <- rownames(rescaled_feature_set)[fiberheatmapTEST$rowInd][1:6]
group2 <- rownames(rescaled_feature_set)[fiberheatmapTEST$rowInd][7:13]
group1 <- rownames(rescaled_feature_set)[fiberheatmapTEST$rowInd][14:18]

group_df <- bind_rows(data.frame(group1 = group1) %>% gather(., key = Immune_group, value = Participant), 
          data.frame(group2 = group2) %>% gather(., key = Immune_group, value = Participant), 
          data.frame(group3 = group3) %>% gather(., key = Immune_group, value = Participant))

write_csv(group_df, "data/combined_immune/fiber_immune_groups.csv")
