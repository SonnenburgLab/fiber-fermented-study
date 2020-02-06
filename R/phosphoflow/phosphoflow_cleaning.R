### Phosphoflow data processing
# this script cleans up the output phosphoflow output received from the HIMC for use in our analysis
# all participants are used in the analysis with the exception of the two participants who were not randomized (see methods)
# this happened in 2 runs
# note: to run script, set current directory to repo directory

library(readxl)
library(dplyr)
library(tidyr)

metadata_directory <- "data/metadata/"

# format is different enough in the 2 runs that we will just parse them separately with a lot of redundant code then combine them. 

###### RUN 1 #######

file_path <- "data/phosphoflow/raw/run1/"
sheet_names <- excel_sheets(paste(file_path, "gardner study_8000_SDM .xlsx", sep = ""))
run = 1

file_list <- list.files(path = file_path, pattern = ".xlsx")


file_dfs <- list()
k = 1
for (selected_file in file_list){
  print(selected_file)
  tabs_list <- list()
  
  for (sheet in sheet_names){
    print(sheet)
    df <- read_excel(paste(file_path, selected_file, sep = ""), sheet = sheet)
    tabs_list[[sheet]] <- df
  }
  
  joined_file <- plyr::join_all(tabs_list, by = 'Sample', type = 'left')
  joined_file[["Plate"]] <- k
  file_dfs[[selected_file]] <- joined_file
  k <- k + 1
}

full_df <- bind_rows(file_dfs)

clean_df <- full_df %>% 
  dplyr::filter(Sample != "Mean") %>% 
  dplyr::filter(Sample != "SD") %>% 
  dplyr::filter(!is.na(Sample)) %>% 
  mutate(Control = grepl("0199|Control", Sample, ignore.case = TRUE))

# separate the controls from the samples and deal with two separate data frames:

sample_df <- dplyr::filter(clean_df, Control == FALSE)
controls_df <- dplyr::filter(clean_df, Control == TRUE)

# creating naming key with cleaned up names and sorting them into cell types and proteins

name_key <- data.frame(original_name = colnames(sample_df))
name_key <- name_key %>% dplyr::filter(original_name != "Plate") %>% 
  dplyr::filter(original_name != "Control") %>% 
  dplyr::filter(original_name != "Sample")

name_key[["Protein"]] <- "p"
phospho_list <- c("pSTAT3", "pSTAT1", "pSTAT5", "pERK", "pP38", "pPLCg", "Freq")

for (phospho in phospho_list){
  ind <- grep(phospho, name_key$original_name)
  name_key[ind, "Protein"] <- phospho
}

name_key[["Celltype"]] <- "c"
ct_list_regex <- c("B cells", "CD4\\+", "CD4\\+CD45RA\\+", "CD4\\+CD45RA\\-", "CD8", "CD8\\+CD45RA\\+", "CD8\\+CD45RA\\-", "Non BT", "Monocytes")
ct_list <- c("Bcells", "CD4Tcells", "CD4TcellsCD45RApos", "CD4TcellsCD45RAneg", "CD8Tcells", "CD8TcellsCD45RApos", "CD8TcellsCD45RAneg", "NonBT", "Monocytes")

for (i in 1:length(ct_list)){
  ind <- grep(ct_list_regex[i], name_key$original_name)
  name_key[ind, "Celltype"] <- ct_list[i]
}
name_key <- tidyr::unite(name_key, col = new_name, Celltype:Protein, sep = "_", remove = FALSE)
name_key$original_name <- droplevels(name_key$original_name)

# sorting into conditons, individuals, timepoints, etc

sample_df <- sample_df %>% 
  tidyr::separate(col = Sample, into = c("temp1", "Condition", "temp2"), sep = " ", remove = FALSE) %>%
  tidyr::separate(col = temp1, into = c("Study", "Participant", "Timepoint", "temp3"), sep = "-", remove = FALSE) %>%
  dplyr::select(-temp1, -temp2, -temp3, -Study)

sample_df$Condition <- gsub(pattern = "erk", replacement = "ERK", sample_df$Condition)

# lets clean up the column names now: 

sample_df_long <- gather(sample_df, key = "original_name", value = "value", 'singlet/Lymphocytes/B cells,Freq. of Parent':'singlet/Monocytes,90th %ile,<PE-A>,pPLCg')
sample_df_long$original_name <- factor(sample_df_long$original_name)
name_switch <- dplyr::select(name_key, original_name, new_name)
sample_df_long <- left_join(sample_df_long, name_switch)

updated_sample_df <- sample_df_long %>% dplyr::select(-original_name) %>% spread(key = new_name, value = value)

# separating STATs condition and the ERK/P38/PLCg conditions (LPS and US-ERK)

main_df <- updated_sample_df %>% gather(key = "new_name", value = "value", Bcells_Freq:NonBT_pSTAT5) %>%
  dplyr::filter(Condition != "US-ERK") %>% dplyr::filter(Condition != "LPS") %>%
  left_join(., name_key) %>% dplyr::filter(Protein != "pERK") %>% dplyr::filter(Protein != "pP38") %>% 
  dplyr::filter(Protein != "pPLCg") %>% 
  unite(col = "Participant_Timepoint", Participant:Timepoint, sep = "_", remove = FALSE) %>%
  dplyr::select(-original_name, -Protein, -Celltype) %>%
  spread(key = new_name, value = value) %>% mutate(DataSet = run)

LPS_df <- updated_sample_df %>% gather(key = "new_name", value = "value", Bcells_Freq:NonBT_pSTAT5) %>%
  dplyr::filter(Condition == "US-ERK" | Condition == "LPS") %>%
  left_join(., name_key) %>% dplyr::filter(Protein == "pERK" | Protein == "pP38" | Protein == "pPLCg") %>% 
  unite(col = "Participant_Timepoint", Participant:Timepoint, sep = "_", remove = FALSE) %>%
  dplyr::select(-original_name, -Protein, -Celltype) %>%
  spread(key = new_name, value = value) %>% mutate(DataSet = run)

# fold changes:
# STATs ("main")

basal_main_df <- main_df %>% dplyr::filter(Condition == "US") %>% arrange(Participant_Timepoint)
basal_main_numeric <- basal_main_df %>% select(-(Sample:Control))


conditions <- c("IFNa", "IFNg", "IL10", "IL2", "IL6")
fold_df_list <- list()

for (condition in conditions){
  cond_df <- main_df %>% dplyr::filter(Condition == condition) %>% 
    arrange(Participant_Timepoint) %>%
    select(-(Sample:Control))
  metadata <- main_df %>% dplyr::filter(Condition == condition) %>% 
    arrange(Participant_Timepoint) %>%
    select(Sample:Control)
  fold_df <- cond_df/basal_main_numeric
  fold_df <- bind_cols(metadata, fold_df)
  fold_df_list[[condition]] <- fold_df
}

fold_main_df <- bind_rows(fold_df_list) %>% mutate(DataSet = run)

# LPS 

basal_LPS_df <- LPS_df %>% dplyr::filter(Condition == "US-ERK") %>%
  arrange(Participant_Timepoint)
basal_LPS_numeric <- basal_LPS_df %>% select(-(Sample:Control))

cond_df <- LPS_df %>% dplyr::filter(Condition == "LPS") %>%
  arrange(Participant_Timepoint) %>%
  select(-(Sample:Control))
metadata <- LPS_df %>% dplyr::filter(Condition == "LPS") %>% 
  arrange(Participant_Timepoint) %>%
  select(Sample:Control)
fold_df <- cond_df/basal_LPS_numeric
fold_LPS_df <- bind_cols(metadata, fold_df) %>% mutate(DataSet = run)

# we want each feature per condition
# filter by condition, remove "Sample", add name to column name, be compiling that as a list, then merge

df_list <- list()

for (condition in conditions){
  df <- fold_main_df %>% dplyr::filter(Condition == condition) %>% dplyr::select(-Sample, -Condition)
  colnames(df)[colnames(df) %in% name_key$new_name] <- paste(condition, colnames(df)[colnames(df) %in% name_key$new_name], sep = "_")
  df_list[[condition]] <- df
}

fold_main_conditions_df <- plyr::join_all(df_list, type = 'left') %>% mutate(DataSet = run)

compiled_structures_run1 <- list(main_df, basal_main_df, main_name_key, lps_name_key, fold_main_df, fold_main_conditions_df, fold_LPS_df)

###### RUN 2 #######

file_path <- "data/phosphoflow/raw/run2/"
sheet_names <- excel_sheets(paste(file_path, "Plate1_Gardner study 2018_SDM .xlsx", sep = ""))
run = 2

file_list <- list.files(path = file_path, pattern = ".xlsx")

# read in all xlsx files and combine into a single data frame

file_dfs <- list()

k = 1
for (selected_file in file_list){
  print(selected_file)
  tabs_list <- list()
  
  for (sheet in sheet_names){
    print(sheet)
    df <- read_excel(paste(file_path, selected_file, sep = ""), sheet = sheet)
    tabs_list[[sheet]] <- df
    if (colnames(tabs_list[[sheet]])[1] == 'Sample:'){
      colnames(tabs_list[[sheet]])[1] <- 'Sample'
    }
  }
  
  joined_file <- plyr::join_all(tabs_list, by = 'Sample', type = 'left')
  joined_file[["Plate"]] <- k
  file_dfs[[selected_file]] <- joined_file
  k <- k + 1
}

full_df <- bind_rows(file_dfs)

# removing extraneous rows

clean_df <- full_df %>% 
  dplyr::filter(Sample != "Mean") %>% 
  dplyr::filter(Sample != "SD") %>% 
  dplyr::filter(!is.na(Sample)) %>% 
  mutate(Control = grepl("088|control|ctrl|0199|Control", Sample, ignore.case = TRUE))

#check for NAs in the data COME BACK TO THIS
na_loc <- apply(clean_df, 1, function(x) sum(is.na(x)) > 0)
clean_df$Sample[na_loc]

# let's separate the controls from the samples and work with them as separate data frames: 

sample_df <- dplyr::filter(clean_df, Control == FALSE)
controls_df <- dplyr::filter(clean_df, Control == TRUE)

# we need to separate the LPS samples from the cytokine samples because the channels correspond to different proteins:

main_sample_df <- mutate(sample_df,LPS=grepl("ERK|LPS",Sample)) %>% dplyr::filter(LPS == FALSE)
lps_sample_df <- mutate(sample_df,LPS=grepl("ERK|LPS",Sample)) %>% dplyr::filter(LPS == TRUE)

# note: one duplicated sample entry needed to be removed: 

main_sample_df <- main_sample_df[!duplicated(main_sample_df$Sample),]

# Create a key that maps the original name from the raw data to a cleaned up name, with a Cell type and Protein field

main_name_key <- data.frame(original_name = colnames(sample_df))


main_name_key <- main_name_key %>% dplyr::filter(original_name != "Plate") %>% 
  dplyr::filter(original_name != "Control") %>% 
  dplyr::filter(original_name != "Sample")

main_name_key[["Protein"]] <- "p"
phospho_list_ind <- c("APC","FITC","PE","Freq")
phospho_list <- c("pSTAT3", "pSTAT1", "pSTAT5","Freq")

for (i in 1:length(phospho_list)){
  
  ind <- grep(phospho_list_ind[i], main_name_key$original_name)
  main_name_key[ind, "Protein"] <- phospho_list[i]
  
}

main_name_key[["Celltype"]] <- "c"
ct_list_regex <- c("B cells", "CD4\\+", "CD4\\+CD45RA\\+", "CD4\\+CD45RA\\-", "CD8", "CD8\\+CD45RA\\+", "CD8\\+CD45RA\\-", "Non BT", "Monocytes")
ct_list <- c("Bcells", "CD4Tcells", "CD4TcellsCD45RApos", "CD4TcellsCD45RAneg", "CD8Tcells", "CD8TcellsCD45RApos", "CD8TcellsCD45RAneg", "NonBT", "Monocytes")

for (i in 1:length(ct_list)){
  
  ind <- grep(ct_list_regex[i], main_name_key$original_name)
  main_name_key[ind, "Celltype"] <- ct_list[i]
  
}

main_name_key <- tidyr::unite(main_name_key, col = new_name, Celltype:Protein, sep = "_", remove = FALSE)
main_name_key$original_name <- droplevels(main_name_key$original_name)

## creating naming key for LPS condition

lps_name_key <- data.frame(original_name = colnames(lps_sample_df))

lps_name_key <- lps_name_key %>% dplyr::filter(original_name != "Plate") %>% 
  dplyr::filter(original_name != "Control") %>% 
  dplyr::filter(original_name != "Sample") %>% 
  dplyr::filter(original_name != "LPS")

lps_name_key[["Protein"]] <- "p"
phospho_list_ind <- c("APC","FITC","PE","Freq")
phospho_list <- c("pERK", "pP38", "pPLCg","Freq")

for (i in 1:length(phospho_list)){
  
  ind <- grep(phospho_list_ind[i], lps_name_key$original_name)
  lps_name_key[ind, "Protein"] <- phospho_list[i]
  
}

lps_name_key[["Celltype"]] <- "c"
ct_list_regex <- c("B cells", "CD4\\+", "CD4\\+CD45RA\\+", "CD4\\+CD45RA\\-", "CD8", "CD8\\+CD45RA\\+", "CD8\\+CD45RA\\-", "Non BT", "Monocytes")
ct_list <- c("Bcells", "CD4Tcells", "CD4TcellsCD45RApos", "CD4TcellsCD45RAneg", "CD8Tcells", "CD8TcellsCD45RApos", "CD8TcellsCD45RAneg", "NonBT", "Monocytes")

for (i in 1:length(ct_list)){
  
  ind <- grep(ct_list_regex[i], lps_name_key$original_name)
  lps_name_key[ind, "Celltype"] <- ct_list[i]
  
}
lps_name_key <- tidyr::unite(lps_name_key, col = new_name, Celltype:Protein, sep = "_", remove = FALSE)
lps_name_key$original_name <- droplevels(lps_name_key$original_name)


# separating sample names into conditons, individuals, and timepoints fields

main_sample_df <- main_sample_df %>% 
                      tidyr::separate(col = Sample, into = c("temp1", "Condition", "temp2"), sep = " ", remove = FALSE) %>%
                      tidyr::separate(col = temp1, into = c("Participant", "Timepoint"), sep = "-", remove = FALSE) %>%
                      dplyr::select(-temp1, -temp2)

lps_sample_df <- lps_sample_df %>% 
                      tidyr::separate(col = Sample, into = c("temp1", "Condition", "temp2"), sep = " ", remove = FALSE) %>%
                      tidyr::separate(col = temp1, into = c("Participant", "Timepoint"), sep = "-", remove = FALSE) %>%
                      dplyr::select(-temp1, -temp2)

lps_sample_df$Condition <- gsub(pattern = "erk", replacement = "ERK", lps_sample_df$Condition)

#### Clean up column names (main): 

main_sample_df_long <- gather(main_sample_df, key = "original_name", value = "value", "Singlets/Lymphocytes/B cells | Freq. of Parent":"singlet/Monocytes,90th %ile,<PE-A>,pSTAT5")
main_sample_df_long$original_name <- factor(main_sample_df_long$original_name)
main_name_switch <- dplyr::select(main_name_key, original_name, new_name)
main_sample_df_long <- left_join(main_sample_df_long, main_name_switch)

main_updated_sample_df <- main_sample_df_long %>% dplyr::select(-original_name) %>% spread(key = new_name, value = value)

# LPS

lps_sample_df_long <- gather(lps_sample_df, key = "original_name", value = "value", "Singlets/Lymphocytes/B cells | Freq. of Parent":"singlet/Monocytes,90th %ile,<PE-A>,pSTAT5")
lps_sample_df_long$original_name <- factor(lps_sample_df_long$original_name)
lps_name_switch <- dplyr::select(lps_name_key, original_name, new_name)
lps_sample_df_long <- left_join(lps_sample_df_long, lps_name_switch)

lps_updated_sample_df <- lps_sample_df_long %>% dplyr::select(-original_name) %>% spread(key = new_name, value = value)

# separating STATs condition and the ERK/P38/PLCg conditions (LPS and US-ERK)

main_sample_df <- main_updated_sample_df %>% 
                            gather(key = "new_name", value = "value", Bcells_Freq:NonBT_pSTAT5) %>%
                            left_join(., main_name_key) %>%
                            unite(col = "Participant_Timepoint", Participant:Timepoint, sep = "_", remove = FALSE) %>%
                            dplyr::select(-original_name, -Protein, -Celltype) %>%
                            spread(key = new_name, value = value) %>% select(-LPS) %>% mutate(DataSet = run)

lps_sample_df <- lps_updated_sample_df %>% 
                            gather(key = "new_name", value = "value", Bcells_Freq:NonBT_pPLCg) %>%
                            left_join(., lps_name_key) %>% dplyr::filter(Protein == "pERK" | Protein == "pP38" | Protein == "pPLCg") %>%
                            unite(col = "Participant_Timepoint", Participant:Timepoint, sep = "_", remove = FALSE) %>%
                            dplyr::select(-original_name, -Protein, -Celltype) %>%
                            spread(key = new_name, value = value) %>% select(-LPS) %>% mutate(DataSet = run)

# fold changes:
# STATs ("main")

## unstimulated data

basal_main_df <- main_sample_df %>% dplyr::filter(Condition == "US") %>% arrange(Participant_Timepoint)
basal_main_numeric <- basal_main_df %>% select(-(Sample:Control))

conditions <- c("IFNa", "IFNg", "IL10", "IL2", "IL6")
fold_df_list <- list()

for (condition in conditions){
  cond_df <- main_sample_df %>% dplyr::filter(Condition == condition) %>% 
    arrange(Participant_Timepoint) %>%
    select(-(Sample:Control))
  metadata <- main_sample_df %>% dplyr::filter(Condition == condition) %>% 
    arrange(Participant_Timepoint) %>%
    select(Sample:Control)
  fold_df <- cond_df/basal_main_numeric
  fold_df <- bind_cols(metadata, fold_df)
  fold_df_list[[condition]] <- fold_df
}

fold_main_df <- bind_rows(fold_df_list) %>% mutate(DataSet = run)

# LPS fold change:

basal_LPS_df <- lps_sample_df %>% dplyr::filter(Condition == "US-ERK") %>%
  arrange(Participant_Timepoint)

basal_LPS_numeric <- basal_LPS_df %>% select(-(Sample:Control))

cond_df <- lps_sample_df %>% dplyr::filter(Condition == "LPS") %>%
  arrange(Participant_Timepoint) %>%
  select(-(Sample:Control))

metadata <- lps_sample_df %>% dplyr::filter(Condition == "LPS") %>% 
  arrange(Participant_Timepoint) %>%
  select(Sample:Control)

fold_df <- cond_df/basal_LPS_numeric

fold_LPS_df <- bind_cols(metadata, fold_df) %>% mutate(DataSet = run) %>% filter()

# turns out there's another thing we want: have each as a feature per condition
# filter by condition, remove "Sample", add name to column name, be compiling that as a list, then merge
# only done for main data frame and not LPS data frame because only stim with LPS for LPS plate 

df_list <- list()

for (condition in conditions){
  df <- fold_main_df %>% dplyr::filter(Condition == condition) %>% dplyr::select(-Sample, -Condition)
  colnames(df)[colnames(df) %in% main_name_key$new_name] <- paste(condition, colnames(df)[colnames(df) %in% main_name_key$new_name], sep = "_")
  df_list[[condition]] <- df
}

fold_main_conditions_df <- plyr::join_all(df_list, type = 'left') %>% mutate(DataSet = run)

compiled_structures_run2 <- list(main_sample_df, basal_main_df, main_name_key, lps_name_key, fold_main_df, fold_main_conditions_df, fold_LPS_df)


saveRDS(compiled_structures_run1, file = "data/phosphoflow/cleaned/compiled_phosphoflow_structures_run1.rds")
saveRDS(compiled_structures_run2, file = "data/phosphoflow/cleaned/compiled_phosphoflow_structures_run2.rds")

# bind rows of the first and second runs together and save as new objects

main_df_all <- rbind(compiled_structures_run1[[1]],compiled_structures_run2[[1]])
basal_main_df_all <- rbind(compiled_structures_run1[[2]],compiled_structures_run2[[2]])
fold_main_df_all <- rbind(compiled_structures_run1[[5]],compiled_structures_run2[[5]])
fold_main_conditions_df_all <- rbind(compiled_structures_run1[[6]],compiled_structures_run2[[6]])
fold_LPS_df_all <- rbind(compiled_structures_run1[[7]],compiled_structures_run2[[7]])

compiled_structures_all <- list(main_df_all, basal_main_df_all, fold_main_df_all, fold_main_conditions_df_all, fold_LPS_df_all)

# saving as R object for reference but will ignore in repo
# saving as .csvs for reference on repo

saveRDS(compiled_structures_all, file = "data/phosphoflow/cleaned/compiled_phosphoflow_structures_combined.rds")

# we are using a single set for the analysis that is the fold change of all conditions including LPS
# we just need to annotate the LPS condition to include LPS in the feature name: 

colnames(fold_LPS_df_all)[8:34] <- paste("LPS_", colnames(fold_LPS_df_all)[8:34], sep = "")

phosphoflow_feature_set <- left_join(fold_main_conditions_df_all, fold_LPS_df_all)

# for our final feature set for testing and for combination, differences, etc we must exclude our non-randomized ppts:

library(tidyverse)

diet_key <- read_csv(paste(metadata_directory, "diet_key.csv", sep = "")) %>% mutate(Participant = as.character(Participant))

phosphoflow_feature_set <- phosphoflow_feature_set %>% 
                                dplyr::filter(!(Participant %in% c("8037", "8038"))) %>% 
                                dplyr::select(-Participant_Timepoint, -Plate, -Control, -Sample, -Condition) %>%
                                left_join(., diet_key) %>%
                                select(Participant, Timepoint, Group, everything())

write_csv(phosphoflow_feature_set, path = "data/phosphoflow/cleaned/phosphoflow_feature_set.csv")

saveRDS(phosphoflow_feature_set, file = "data/phosphoflow/cleaned/phosphoflow_feature_set.rds")

# difference analysis for combined immune: 

feature_means <- phosphoflow_feature_set %>% select(-c(Participant:Group)) %>% summarise_all(mean)

participant_list <- levels(factor(phosphoflow_feature_set$Participant))

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


phospho_data_imputed <- add_imputed_data(data = phosphoflow_feature_set, 
                                          feature_means, 
                                          participant_list,
                                          selected_timepoints = c("01", "05", "06"))
# lools like we didn't actually need to impute anyting (all data there)

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

phospho_data_differences <- find_differences(data = phospho_data_imputed, 
                                           participant_list, 
                                           reference_time = "01", 
                                           end_time_set = c("05", "06"))

write_csv(phospho_data_differences, path = "data/phosphoflow/cleaned/phospho_data_differences_with_imputed.csv")


