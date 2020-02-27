## Features for Combined Immune Profiling Analysis
# this script generates a feature key to be used in multi-modal immune profiling analysis
# to rerun script set working directory to repo directory

library(tidyverse)

olink_differences <- read_csv("data/Olink/cleaned/Olink_data_differences_with_imputed.csv")

cytof_signaling_differences <- read_csv("data/CyTOF/cleaned/cytof_signaling_differences_with_imputed.csv")

cytof_frequency_differences <- read_csv("data/CyTOF/cleaned/cytof_frequency_differences_with_imputed.csv")

phospho_differences <- read_csv("data/phosphoflow/cleaned/phospho_data_differences_restricted.csv")

make_feature_key <- function(data, feature_type){
  
  gathered_data <- gather(data, key = "Feature", value = "value", -c(Participant, Group)) %>% select(Feature)
  
  feature_data <- tibble(Feature = gathered_data$Feature, Type = feature_type) %>%
                    separate(col = Feature, into = c("Master_feature", "Time"), sep = "_T1_", remove = FALSE) %>%
                    select(Feature, Master_feature, Time, Type)
                  
  
}

feature_sets <- list("Cytokine" = olink_differences, 
                     "Cell frequencies" = cytof_frequency_differences, 
                     "Endogenous signaling" = cytof_signaling_differences,
                     "Signaling capacity" = phospho_differences)

type_list <- list()

for (feature_type in names(feature_sets)){
  
  feature_info <- make_feature_key(data = feature_sets[[feature_type]], feature_type = feature_type)
  
  type_list[[feature_type]] <- feature_info
  
}

feature_key <- reduce(type_list, bind_rows)

feature_key <- unique(feature_key)

combined_immune_feature_differences_set <- reduce(feature_sets, left_join)

write_csv(feature_key, path = "data/combined_immune/immune_feature_key.csv")

write_csv(combined_immune_feature_differences_set, path = "data/combined_immune/combined_immune_feature_differences_set.csv")

