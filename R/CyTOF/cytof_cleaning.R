# This script cleans up the output exported from immuneatlas.org including:
# * event counts (for frequencies)
# * medians (for steady-state signaling)
# the version of the exported data used for the publication's analysis can be found in the stated data directory

# to do's:
# finalized what participants are being included (donor removal, hospital person?) Right now its everyone

data_directory <- "data/CyTOF/exported/"
metadata_directory <- "data/metadata/"

### counts data

library(tidyverse)

counts_df <- read.csv(paste(data_directory, "cytof-exported-event-counts.csv", sep = "")) %>%
                    separate(., col = SAMPLE_NAME, into = c("Participant", "Timepoint_temp"), sep = "-", remove = FALSE) %>%
                    separate(., col = Timepoint_temp, into = c("Timepoint", "temp"), sep = "\\.", remove = TRUE) %>%
                    dplyr::select(-temp)

# ID samples with > 10K CD45+ cells

CD45pos_cells <- counts_df %>% 
                      dplyr::filter(population == "CD45+CD66-") %>% 
                      arrange(desc(eventCount)) %>% 
                      mutate(Greater_10K = (eventCount > 10000))

ggplot(CD45pos_cells, aes(x = eventCount)) + geom_histogram(bins = 50, aes(fill = (Greater_10K))) + theme_classic()

counts_df <- left_join(counts_df, select(CD45pos_cells, -population, -eventCount))

# derive cell frequencies

frequency_df <- counts_df %>% 
  spread(., key = population, value = eventCount) %>%
  mutate(B_cell_freq = `B cells`/`CD45+CD66-`) %>% 
  mutate(CD4_T_cell_freq = `CD4 T cells`/`CD45+CD66-`) %>%
  mutate(CM_CD4_T_cell_freq = `CM_CD4`/`CD45+CD66-`) %>%
  mutate(CM_CD8_T_cell_freq = `CM_CD8`/`CD45+CD66-`) %>%
  mutate(Effector_CD4_T_cell_freq = `Effector_CD4`/`CD45+CD66-`) %>%
  mutate(Effector_CD48_T_cell_freq = `Effector_CD8`/`CD45+CD66-`) %>%
  mutate(EM_CD4_T_cell_freq = `EM_CD4`/`CD45+CD66-`) %>%
  mutate(EM_CD8_T_cell_freq = `EM_CD8`/`CD45+CD66-`) %>%
  mutate(Naive_CD4_T_cell_freq = `Naive_CD4`/`CD45+CD66-`) %>%
  mutate(Naive_CD8_T_cell_freq = `Naive_CD8`/`CD45+CD66-`) %>%
  mutate(CD8_T_cell_freq = `CD8 T cells`/`CD45+CD66-`) %>%
  mutate(Classical_monocytes_freq = `Classical monocytes`/`CD45+CD66-`) %>%
  mutate(Nonclassical_monocytes_freq = `Nonclassical monocytes`/`CD45+CD66-`) %>%
  mutate(mDCs_freq = `mDCs`/`CD45+CD66-`) %>%
  mutate(NK_cell_freq = `NK Cells`/`CD45+CD66-`) %>%
  mutate(pDCs_freq = `pDCs`/`CD45+CD66-`) %>%
  mutate(Neutrophils_freq = `Neutrophils`/`Cells`) %>%
  dplyr::select(-(`B cells`:`pDCs`))

# adding diet annotation, restricting to baseline and maintenance phases, restricting to samples with > 10K

diet_key <- read.csv(paste(metadata_directory, "diet_key.csv", sep = ""))

diet_key$Participant <- as.character(diet_key$Participant)

cleaned_frequency_data <- frequency_df %>%
          left_join(., diet_key) %>%
          dplyr::filter(Timepoint != "03") %>% 
          dplyr::filter(Timepoint != "07") %>%
          mutate(Phase = car::recode(Timepoint, "c('01','02') = 'Baseline'; c('04','05','06') = 'Intervention'")) %>%
          dplyr::filter(Greater_10K == TRUE) %>%
          select(filename:Timepoint, Group, Phase, Greater_10K, everything())

write.csv(cleaned_frequency_data, file = "data/CyTOF/cleaned/cytof_frequencies_cleaned.csv")

### differences from baseline to intervention frequency data

# if we are missing baseline 1 (missing or low cell count), replace with baseline 2

# participants we are using for this: 

cleaned_frequency_data$Participant <- as.factor(cleaned_frequency_data$Participant)
participant_list <- levels(cleaned_frequency_data$Participant)

replace_baselines <- function(data, participant_list){
  
  updated_baselines <- data %>% mutate_if(is.factor, as.character)
  
  for (participant in participant_list){
    
    baseline_1 <- data %>% 
      dplyr::filter(Participant == participant) %>% 
      dplyr::filter(Timepoint == "01")  %>% 
      mutate_if(is.factor, as.character)
    
    baseline_2 <- data %>% 
      dplyr::filter(Participant == participant) %>% 
      dplyr::filter(Timepoint == "02") %>% 
      mutate_if(is.factor, as.character)
    
    diet <- diet_key %>% dplyr::filter(Participant == participant) %>% select(Group) %>% mutate_if(is.factor, as.character)
    
    if (dim(baseline_1)[1] == 0 & dim(baseline_2)[1] > 0){
      
      replaced_baseline <- baseline_2
      replaced_baseline[, "Timepoint"] <- "01"
      
      updated_baselines <- bind_rows(updated_baselines, replaced_baseline)
      
    }
  
  }
  return(updated_baselines)
}

cleaned_frequencies_updated_baselines <- replace_baselines(data = cleaned_frequency_data, participant_list)


# adding imputed values for low or missing data:

create_imputed_row <- function(feature_means, participant, timepoint, diet){
  
  imputed_row <- data.frame(filename = NA,
                            BATCH = "missing_sample",
                            SAMPLE_NAME = NA,
                            Participant = participant,
                            Timepoint = timepoint,
                            Group = diet,
                            Phase = "Baseline",
                            Greater_10K = FALSE,
                            feature_means,
                            stringsAsFactors = FALSE)
  
  return(imputed_row)
  
}

add_imputed_data <- function(data, feature_means, participant_list, selected_timepoints){
  
  added_rows <- c()
  
  for (participant in participant_list){
    
    diet <- diet_key %>% dplyr::filter(Participant == participant) %>% select(Group) %>% mutate_if(is.factor, as.character)
    
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

feature_means <- cleaned_frequency_data %>% select(-c(filename:Greater_10K)) %>% summarise_all(mean)

cleaned_frequencies_imputed <- add_imputed_data(data = cleaned_frequencies_updated_baselines, 
                                                       feature_means, 
                                                       participant_list,
                                                       selected_timepoints = c("01", "04", "05", "06"))


find_differences <- function(data, participant_list, reference_time, end_time_set){
  
  reference_data <- data %>% 
                      dplyr::filter(Timepoint == reference_time) %>% 
                      dplyr::filter(Participant %in% participant_list) %>% 
                      arrange(Participant) %>% 
                      select(-c(filename:Greater_10K))
  
  for (time in end_time_set){
    
    difference_list <- list()
    
    end_data <- data %>% 
                    dplyr::filter(Timepoint == time) %>% 
                    dplyr::filter(Participant %in% participant_list) %>% 
                    arrange(Participant)
    end_features <- end_data %>% select(-c(filename:Greater_10K))
    
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

cleaned_frequency_differences <- find_differences(data = cleaned_frequencies_imputed, 
                                                  participant_list, 
                                                  reference_time = "01", 
                                                  end_time_set = c("04", "05", "06"))

write.csv(cleaned_frequency_differences, file = "data/CyTOF/cleaned/cytof_frequency_differences_with_imputed.csv")

## signaling data

signaling_data <- read.csv(paste(data_directory, "cytof-exported-medians.csv", sep = "")) %>% 
                        separate(., col = reagent, into = c("Channel", "Protein"), sep = "_") %>% 
                        separate(., col = SAMPLE_NAME, into = c("Participant", "temp"), sep = "-", remove = FALSE) %>% 
                        separate(., col = temp, into = c("Timepoint", "temp"), sep = "\\.") %>% 
                        select(-temp) %>% 
                        mutate(median_transformed = asinh(median/5)) %>% 
                        dplyr::left_join(., diet_key) %>%
                        dplyr::left_join(select(CD45pos_cells, filename, Participant, Timepoint, Greater_10K)) %>% 
                        dplyr::filter(Timepoint != "03") %>% 
                        dplyr::filter(Timepoint != "07") %>%
                        mutate(Phase = car::recode(Timepoint, "c('01','02') = 'Baseline'; c('04','05','06') = 'Intervention'")) %>%
                        unite(col = "Feature", population, Protein, remove = FALSE) %>%
                        select(filename:Timepoint, Group, Phase, Greater_10K, everything())


# the set we use in our significance analysis is restricted to four major cell types and j
# excludes samples with low cell count

signaling_data_restricted <- signaling_data %>% 
                      dplyr::filter(population == "B cells" | 
                                    population == "CD4 T cells" |
                                    population == "CD8 T cells" | 
                                    population == "Classical monocytes") %>% 
                      dplyr::filter(Greater_10K == TRUE) %>% 
                      select(filename, BATCH, SAMPLE_NAME, Participant, Timepoint, Group, Phase, Greater_10K, Feature, median_transformed) %>%
                      spread(., key = Feature, value = median_transformed)

# up next: 
# decide what want to save
# run imputation/differences code on signaling and save
# then siggenes (new script)
