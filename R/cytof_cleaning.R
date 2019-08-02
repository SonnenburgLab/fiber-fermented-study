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
