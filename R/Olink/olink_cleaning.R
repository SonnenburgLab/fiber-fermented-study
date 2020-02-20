### Olink data processing
# this script cleans up the output received from Olink for use in our analysis
# all participants are used in the analysis including the the two fiber participants who were not randomized (8037, 8038, see methods), with
# the exception of 2 ppts that were excluded (8000, 8012, see methods)
# note: to run script, set current directory to repo directory

library(tidyverse)

data_directory <- "data/Olink/raw/"

metadata_directory <- "data/metadata/"
  
raw_data <- read_csv(paste(data_directory, "olink_data_raw.csv", sep = ""))

assay_key <- read_csv(paste(data_directory, "olink_assay_key.csv", sep = ""))

# we will not include proteins that were at the limit of detection for > 25% of samples

extracted_percents <- as.character(assay_key$`Missing Data freq.`) %>% 
                          gsub(pattern = "%", replacement = "", x = . ) %>%
                          as.numeric(.)

proteins_passed <- assay_key %>%
                      mutate(Missing_data_percent = extracted_percents) %>%
                      dplyr::filter(Missing_data_percent < 25) %>% 
                      select(Assay)

# filtering for detected proteins, removing control samples, plasma samples, and parsing for metadata in tube labels:

diet_key <- read.csv(paste(metadata_directory, "diet_key.csv", sep = "")) %>% mutate(Participant = as.character(Participant))

filtered_data <- raw_data %>% 
                    gather(key = "Protein", value = "NPX", -TubeID, -`Plate ID`, -`QC Warning`) %>%
                    mutate(NPX = as.numeric(NPX)) %>%
                    dplyr::filter(Protein %in% proteins_passed$Assay) %>%
                    dplyr::filter(!str_detect(TubeID, 'CONTROL')) %>%
                    separate(col = TubeID, into = c("Study", "Participant", "Timepoint", "Aliquot"), 
                             sep = "-", remove = FALSE) %>%
                    separate(col = Aliquot, into = c("Sample_type", "Aliquot_no"), 
                             sep = "0", remove = TRUE) %>% 
                    left_join(., diet_key) %>%
                    dplyr::filter(Sample_type == "S") %>% 
                    select(Participant, Timepoint, Group, Protein, NPX) %>%
                    spread(key = Protein, value = NPX)

# remove samples that we are excluding (n = 2), where we are missing intervention timepoint (n = 1), and samples from 
# time points where we only have one participant (i.e. timepoints 3 and 7)

filtered_data_remove_missing <- filtered_data %>% 
                    dplyr::filter(!Participant %in% c("8035","8000","8012")) %>%
                    dplyr::filter(Timepoint != "03") %>%
                    dplyr::filter(Timepoint != "07")

write.csv(filtered_data_remove_missing, file = "data/Olink/cleaned/olink_data_cleaned.csv", row.names = FALSE)

## Fold change data

# we will still exclude the 2 excluded ppts, but we can keep the participant with the missing intervention timepoint and impute

all_participants <- levels(factor(filtered_data$Participant))

participant_list <- all_participants[!all_participants %in% c("8000", "8012")]

# to do: adapt code below to olink

create_imputed_row <- function(feature_means, participant, timepoint, diet){
  
  imputed_row <- tibble(Participant = participant,
                            Timepoint = timepoint,
                            Group = diet) %>%
                  bind_cols(feature_means)
  
  return(imputed_row)
  
}

add_imputed_data <- function(data, feature_means, participant_list, selected_timepoints){
  
  added_rows <- c()
  
  for (participant in participant_list){
    
    diet <- diet_key %>% dplyr::filter(Participant == participant) %>% select(Group) %>% mutate_if(is.factor, as.character) %>% .[1,1]
    
    for (time in selected_timepoints){
      
      entry_check <- data %>% 
        dplyr::filter(Participant == participant) %>% 
        dplyr::filter(Timepoint == time)  %>% 
        mutate_if(is.factor, as.character)
      
      if (dim(entry_check)[1] == 0) {
        
        imputed_row <- create_imputed_row(feature_means = feature_means, participant = participant, timepoint = time, diet = diet)
        
        added_rows <- bind_rows(added_rows, imputed_row)
        
      }
      
    }
  }
  
  data_with_imputed <- bind_rows(data, added_rows)
  
  return(data_with_imputed)
}

feature_means <- filtered_data %>% select(-c(Participant:Group)) %>% summarise_all(mean)

filtered_data_imputed <- add_imputed_data(data = filtered_data, 
                                                feature_means, 
                                                participant_list,
                                                selected_timepoints = c("01", "04", "05", "06"))


find_differences <- function(data, participant_list, reference_time, end_time_set){
  
  reference_data <- data %>% 
    dplyr::filter(Timepoint == reference_time) %>% 
    dplyr::filter(Participant %in% participant_list) %>% 
    arrange(Participant) %>% 
    select(-c(Participant:Group))
  
  difference_list <- list()
  
  for (time in end_time_set){
    
    
    end_data <- data %>% 
      dplyr::filter(Timepoint == time) %>% 
      dplyr::filter(Participant %in% participant_list) %>% 
      arrange(Participant)
    end_features <- end_data %>% select(-c(Participant:Group))
    
    difference_features <- end_features - reference_data
    
    colnames(difference_features) <- paste(colnames(difference_features), "_T", 
                                           substr(reference_time, 2, 2), "_T", 
                                           substr(time, 2, 2), sep = "")
    difference_data <- bind_cols(select(end_data, Participant, Group), difference_features)
    
    difference_list[[time]] <- difference_data           
  }
  
  difference_data <- difference_list %>% reduce(left_join)
  
  return(difference_data)
  
}

olink_data_differences <- find_differences(data = filtered_data_imputed, 
                                                  participant_list, 
                                                  reference_time = "01", 
                                                  end_time_set = c("04", "05", "06"))

write.csv(olink_data_differences, file = "data/Olink/cleaned/Olink_data_differences_with_imputed.csv", row.names = FALSE)
