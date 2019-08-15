### Olink data processing
# this script cleans up the output received from Olink for use in our analysis
# all participants are used in the analysis with the exception of the two participants who were not randomized (see methods)
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

# remove samples that are not randomized (n = 2), where we are missing intervention timepoint (n = 1), and samples from 
# time points where we only have one participant (i.e. timepoints 3 and 7)

filtered_data_remove_missing <- filtered_data %>% 
                    dplyr::filter(!Participant %in% c("8035","8037","8038")) %>%
                    dplyr::filter(Timepoint != "03") %>%
                    dplyr::filter(Timepoint != "07")

write.csv(filtered_data_remove_missing, file = "data/Olink/cleaned/olink_data_cleaned.csv")



