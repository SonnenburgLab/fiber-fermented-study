# This script cleans up the output exported from immuneatlas.org including:
# * event counts (for frequencies)
# * medians (for steady-state signaling)
# the version of the exported data used for the publication's analysis can be found in the stated data directory

data_directory <- "data/CyTOF/exported/"
metadata_directory <- "data/metadata/"

### counts data

library(dplyr)
library(tidyr)
library(ggplot2)

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

frequency_df$Participant <- as.numeric(frequency_df$Participant)


# adding diet annotation, restricting to baseline and maintenance phases, restricting to samples with > 10K

diet_key <- read.csv(paste(metadata_directory, "diet_key.csv", sep = ""))

cleaned_frequency_data <- frequency_df %>%
          left_join(., diet_key) %>%
          dplyr::filter(Timepoint != "03") %>% 
          dplyr::filter(Timepoint != "07") %>%
          mutate(Phase = car::recode(Timepoint, "c('01','02') = 'Baseline'; c('04','05','06') = 'Intervention'")) %>%
          dplyr::filter(Greater_10K == TRUE) %>%
          select(filename:Timepoint, Group, Phase, Greater_10K, everything())

write.csv(cleaned_frequency_data, file = "data/CyTOF/cleaned/cytof_frequencies_cleaned.csv")

# differences from baseline to intervention data

# if we are missing baseline 1 (missing or low cell count), replace with baseline 2

cleaned_frequency_data$Participant <- as.factor(cleaned_frequency_data$Participant)

# if missing bl1
# if have bl2
# get that row and bind
# else impute

feature_means <- cleaned_frequency_data %>% select(-c(filename:Greater_10K)) %>% summarise_all(mean)

cleaned_frequencies_updated_baselines <- cleaned_frequency_data %>% mutate_if(is.factor, as.character)

for (participant in levels(cleaned_frequency_data$Participant)){
  
  baseline_1 <- cleaned_frequency_data %>% 
    dplyr::filter(Participant == participant) %>% 
    dplyr::filter(Timepoint == "01")  %>% 
    mutate_if(is.factor, as.character)
  
  baseline_2 <- cleaned_frequency_data %>% 
    dplyr::filter(Participant == participant) %>% 
    dplyr::filter(Timepoint == "02") %>% 
    mutate_if(is.factor, as.character)
  
  diet <- diet_key %>% dplyr::filter(Participant == participant) %>% select(Group) %>% mutate_if(is.factor, as.character)
  
  if (dim(baseline_1)[1] == 0 & dim(baseline_2)[1] > 0){
    
    replaced_baseline <- baseline_2
    replaced_baseline[, "Timepoint"] <- "01"
    
    cleaned_frequencies_updated_baselines <- bind_rows(cleaned_frequencies_updated_baselines, replaced_baseline)
    
    
  } else if (dim(baseline_1)[1] == 0 & dim(baseline_2)[1] == 0){
    
    replaced_baseline <- data.frame(filename = NA,
                                    BATCH = "missing_sample",
                                    SAMPLE_NAME = NA,
                                    Participant = participant,
                                    Timepoint = "01",
                                    Group = diet,
                                    Phase = "Baseline",
                                    Greater_10K = FALSE,
                                    feature_means,
                                    stringsAsFactors = FALSE)
    
    cleaned_frequencies_updated_baselines <- bind_rows(cleaned_frequencies_updated_baselines, replaced_baseline)
    
  }
}



# to do: 
# for most of them we can replace with baseline 2, for the ones we don't, we need to impute
# also need to impute for time 04, 05, 06
# convert these to functions
# then can grab from the other to merge/take differences etc