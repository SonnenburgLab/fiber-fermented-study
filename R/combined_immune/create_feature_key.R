## Features for Combined Immune Profiling Analysis
# this script generates a feature key to be used in multi-modal immune profiling analysis
# to rerun script set working directory to repo directory

library(tidyverse)

olink_differences <- read_csv("data/Olink/cleaned/Olink_data_differences_with_imputed.csv")


make_feature_key <- function(data, feature_type){
  
  gathered_data <- gather(data, key = "Feature", value = "value", -c(Participant, Group)) %>% select(Feature)
  
  feature_data <- tibble(Feature = gathered_data$Feature, Type = feature_type) %>%
                    separate(col = Feature, into = c("Master_feature", "Time"), sep = "_T1_", remove = FALSE) 
                  
  
}