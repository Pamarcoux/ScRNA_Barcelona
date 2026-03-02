#test
if (!require("pacman")) install.packages("pacman")

pacman::p_load(
  tidyverse, here, patchwork, ggpubr, Seurat, reshape2, paletteer, 
  RColorBrewer, gridExtra, cowplot,
  magick, grid, rio, stringr, pheatmap, ggnewscale, rstatix
)


CRCT_Share <- TRUE
Origin_list <- c("PBMC", "LN")

data_dir <- here("..", "SC_Flex_TALYIES_July2024")

source(here("../GlofiResistance/Code_figures/Replacement_table.R"))
Good_responder_irl <- c("FL4_LN1", "tFL5_PB1", "DLBCL6_LN1", "DLBCL16_LN1", "DLBCL17_LN1", "tFL8_LN1", "tFL10_LN1")
Bad_responder_irl <- c("FL10_LN2", "DLBCL9_LN1")

Pas_de_B_D6 <- c("tFL9_LN1", "DLBCL8_LN1", "DLBCL12_LN1", "tMZL3_LN1", "DLBCL10_LN1")
Removed <- c("DLBCL15_LN1") # DLBCL15_LN1 = Richter

data_talyies_full <- read_csv(file.path(data_dir, "Metadata/Data_TALYIES_grouped.csv"))
data_talyies_full <- data_talyies_full %>%
  relocate(c(B_cell_depletion_total), .before = Day) %>%
  mutate(Irl_Response = case_when(
    Sample_code %in% Good_responder_irl ~ "CR",
    Sample_code %in% Bad_responder_irl ~ "PD"
  )) %>%
  mutate(Irl_Response = factor(Irl_Response, levels = c("PD", "CR"))) |> 
  mutate(Review = FALSE)


Sample_order <- c(
  "FL1_PB1", "FL1_LN1_1", "FL1_LN1_2", "FL2_LN1", "FL1_LN2", "FL2_PB1", "FL3_PB1", "FL3_LN1", "FL4_LN1", "FL4_BM1", "FL5_PB1_1", "FL5_PB1_2", "FL5_PB1_3", "FL6_PB1_1", "FL6_PB1_2", "FL7_LN1", "FL8_LN1", "FL9_PB1", "FL10_PB1", "FL10_PB2", "FL10_LN2", "FL10_BM1", "FL11_PB1", "FL11_LN1", "FL12_PB1",
  "FL13_PB1", "FL14_LN1", "FL15_PB1", "FL16_PB1_1", "FL16_PB1_2", "FL17_LN1", "FL18_LN1", "FL19_LN1", "tFL1_LN1", "tFL2_LN1", "tFL3_PB1", "tFL3_LN1", "tFL4_PB1", "tFL5_PB1", "tFL6_PB1", "tFL7_PB1", "tFL8_LN1", "tFL9_LN1", "tFL10_LN1", "DLBCL1_LN1", "DLBCL2_LN1", "DLBCL3_PB1", "DLBCL4_PB1", "DLBCL5_PB1", "DLBCL6_LN1",
  "DLBCL7_LN1", "DLBCL8_LN1", "DLBCL9_LN1", "DLBCL10_LN1", "DLBCL11_LN1", "DLBCL12_LN1", "DLBCL13_LN1", "DLBCL14_LN1", "DLBCL15_LN1", "DLBCL16_LN1", "DLBCL17_LN1", "DLBCL18_LN1", "DLBCL19_LN1",
  "MZL1_S1", "MZL2_LN1", "MZL3_PB1", "MZL3_PB2", "MZL4_S1", "MZL5_PB1", "tMZL2_LN1", "tMZL3_LN1"
)

data_table_SC_population <- read.csv(file.path(data_dir, "Metadata/data_table_SC_population.csv"))


data_talyies_full <- data_talyies_full %>%
  mutate(
    Screening = as.logical(Screening),
    ScRNA_flex = as.logical(ScRNA_flex),
    ScRNA_cite = as.logical(ScRNA_cite),
    D_R = as.factor(D_R),
    Origin = factor(Origin, levels = c("PBMC", "LN", "BM", "Spleen")),
    Treatment = as.factor(Treatment),
    Response_type_total = factor(Response_type_total, levels = c("Low", "Medium", "High")),
    Pop = as.factor(Pop),
    Day = as.factor(Day)
  ) %>%
  mutate(Response_Origin = paste(Response_type_total, Origin, sep = "_")) %>%
  mutate(Disease = ifelse(Disease == "FL_DLBCL", "tFL", Disease)) %>%
  mutate(Disease = ifelse(Disease == "FLt", "tFL", Disease)) %>%
  mutate(Disease = as.factor(Disease)) %>%
  mutate(Sample_code_paper = Sample_code) %>%
  mutate(Sample_code_paper = recode(Sample_code_paper, !!!Replacement_table)) %>%
  mutate(Sample_code_paper = factor(Sample_code_paper, levels = Replacement_table)) %>%
  mutate(Sample_code = factor(Sample_code, levels = Sample_order)) %>%
  filter(!Sample_code %in% Removed)

# data_talyies_CD20_neg#
data_talyies_CD20_neg <- data_talyies_full %>%
  dplyr::filter(Day == "D0" & Screening == T & Pop == "CD3_CD4") %>%
  mutate(ratio_CD20 = Per_CD20 / Per_CD22_total) %>%
  relocate(ratio_CD20, .after = Per_CD20) %>%
  dplyr::filter(ratio_CD20 < 0.09 | Patient_code %in% c("201326988_20221028", "19T005587-1", "199601089", "202200610", "202100777", "201807661"))

## 19T005587-1 DLBCL1-LN 10% CD22 UT D6
## 201326988_20221028 FL10_PB1 Pas de B D6
## 202100777 / FL11_PB1 CD20 neg
## 202200610 / FL13_PB1 CD20neg
## 201807661 / tFL7 CD20 neg ?
## 025508129 / tMZL1_PB tMZL
## 19T058974 / tMZL2_LN1
# 23T019286 / tFL9_LN1 Pas de B D6

CD20_neg_data <- setdiff(c(data_talyies_CD20_neg$Patient_code, "19T005587-1", "199601089", "201807661", "025508129", "19T058974", "23T019286"), "201908255")
### 201908255 FL1_PB très peu de B CD20 mais reponse quand même

data_talyies_filtered <- data_talyies_full %>%
  dplyr::filter(Screening == T & ScRNA_flex == T & !Patient_code %in% CD20_neg_data & Disease %in% c("FL", "tFL", "DLBCL") & Day == "D6") %>%
  filter(!Sample_code %in% Pas_de_B_D6)

Sample_list <- data_talyies_filtered %>%
  distinct(Sample_code) %>%
  pull(Sample_code)

Sample_list_extended <- data_talyies_full %>%
  dplyr::filter(ScRNA_flex == T & !Patient_code %in% CD20_neg_data & Disease %in% c("FL", "tFL", "DLBCL")) %>%
  distinct(Sample_code) %>%
  pull(Sample_code)

Sample_list_CD20_pos <- data_talyies_full %>%
  dplyr::filter(Screening == T & Disease %in% c("FL", "tFL", "DLBCL") & !Patient_code %in% CD20_neg_data) %>%
  mutate(Sample_code = factor(Sample_code, levels = Sample_order)) %>%
  distinct(Sample_code) %>%
  pull(Sample_code)

#### Annotation Response_type ####
data_talyies_depletion <- data_talyies_full %>%
  filter(Treatment == "αCD20-TCB 0,1 nM") %>%
  select(Patient_code, B_cell_depletion_total, Sample_code, Response_type_total) %>%
  rename(B_cell_depletion_total_glofi = B_cell_depletion_total) %>%
  mutate(Response_type_total = case_when(
    Patient_code %in% CD20_neg_data ~ NA,
    Patient_code %in% c("19T005587-1") ~ "High",
    B_cell_depletion_total_glofi < 45 ~ "Low",
    B_cell_depletion_total_glofi > 45 ~ "High"
  )) %>%
  distinct() %>%
  mutate(Response_type_total = factor(Response_type_total, levels = c("Low", "High")))


data_talyies_full <- data_talyies_full %>%
  select(-Response_type_total) %>%
  full_join(data_talyies_depletion) %>%
  relocate(B_cell_depletion_total_glofi, .before = B_cell_depletion_total)

data_sample_info_complete <- data_talyies_full %>%
  select(
    Patient_code, Date_prelev, Cevi_code, Disease,
    Origin, D_R, Date_D0, Screening,
    ScRNA_flex, ScRNA_cite, Name_cite, Sample_code,
    Response_type_total, B_cell_depletion_total_glofi,
    Sample_code_paper
  ) %>%
  distinct(Patient_code, .keep_all = TRUE)

#### Data Papier Glofi #####
colors_response_type <- c("lightskyblue2", "lightseagreen")
colors_response_type <- setNames(colors_response_type, c("Low", "High"))
sample_colors_LN_PBMC <- setNames(c("#a00000", "#1a80bb"), c("LN", "PBMC"))
sample_colors_disease <- setNames(c("#D2691E", "palegreen4"), c("FL", "tFL\nDLBCL"))
colors_response_type_irl <- setNames(c("lightskyblue2", "lightseagreen"), c("PD", "CR"))

data_talyies <- data_talyies_full %>%
  dplyr::filter(Screening == T & Disease %in% c("FL", "tFL", "DLBCL") & !Patient_code %in% CD20_neg_data) %>%
  filter(!Sample_code %in% Pas_de_B_D6) %>%
  mutate(Sample_code = factor(Sample_code, levels = Sample_order))

data_talyies_glofi <- data_talyies_full %>%
  dplyr::filter(Screening == T & Disease %in% c("FL", "tFL", "DLBCL")) %>%
  dplyr::filter(!is.na(B_cell_depletion_total_glofi)) %>%
  filter(!Sample_code %in% Pas_de_B_D6) %>%
  mutate(Sample_code = factor(Sample_code, levels = Sample_order))

data_sample_info_filter <- data_talyies %>%
  mutate(Sample_code = factor(Sample_code, levels = Sample_order)) %>%
  select(Sample_code) %>%
  distinct() %>%
  left_join(data_sample_info_complete)

Sample_list_paper <- data_talyies_glofi %>%
  distinct(Sample_code_paper) %>%
  pull(Sample_code_paper)

Sample_list_flex <- data_talyies_filtered %>%
  distinct(Sample_code_paper) %>%
  pull(Sample_code_paper)

Sample_left <- setdiff(Sample_list_paper, Sample_list_flex)
colors_sample_left_list <- rev(paletteer_d("ggsci::default_igv", length(Sample_left)))
sample_color_left <- setNames(colors_sample_left_list, Sample_left)

colors_sample_flex <- rev(paletteer_d("rcartocolor::Bold"))[1:length(Sample_list_flex)]
sample_colors_flex <- setNames(colors_sample_flex, Sample_list_flex[order(factor(Sample_list_flex))])


sample_colors_paper <- c(sample_colors_flex, sample_color_left)


#### Data for Paper Plaftform #####
data_talyies_LN_PBMC <- data_talyies_full %>%
  filter(Review == FALSE) %>%
  dplyr::filter(Disease %in% c("FL", "DLBCL", "tFL") & Origin %in% c("LN", "PBMC"))

##### FL####
data_talyies_FL <- data_talyies_full %>%
  filter(Review == FALSE) %>%
  dplyr::filter(Disease %in% c("FL") & Origin %in% c("LN", "PBMC"))

Sample_list_FL <- data_talyies_FL %>%
  distinct(Sample_code) %>%
  arrange(Sample_code) %>%
  pull(Sample_code)

# Générer une palette de couleurs
# colors_FL <- paletteer_c("grDevices::Temps", length(Sample_list_FL))
colors_FL <- paletteer_d("ggsci::default_igv", length(Sample_list_FL))


# Associer chaque échantillon à une couleur
sample_colors_FL <- setNames(colors_FL, Sample_list_FL)

##### DLBCL_tFL####
data_talyies_tFL_DLBCL <- data_talyies_full %>%
  filter(Review == FALSE) %>%
  dplyr::filter(Disease %in% c("tFL", "DLBCL") & Origin %in% c("PBMC", "LN"))

Sample_list_tFL_DLBCL <- data_talyies_tFL_DLBCL %>%
  distinct(Sample_code) %>%
  pull(Sample_code)

# Générer une palette de couleurs
colors_tFL_DLBCL <- paletteer_c("grDevices::Spectral", length(Sample_list_tFL_DLBCL))
# Associer chaque échantillon à une couleur
sample_colors_tFL_DLBCL <- setNames(colors_tFL_DLBCL, Sample_list_tFL_DLBCL)

## Combine les 2 palettes###
sample_colors_all <- c(sample_colors_FL, sample_colors_tFL_DLBCL)

Tcells_Sample_metadata <- read.csv(file.path(data_dir, "Metadata/Tcells_Sample_metadata.csv")) %>% mutate(
  Screening = as.logical(Screening),
  ScRNA_flex = as.logical(ScRNA_flex),
  ScRNA_cite = as.logical(ScRNA_cite),
  D_R = as.factor(D_R),
  Origin = factor(Origin, levels = c("PBMC", "LN", "BM", "Spleen")),
  Disease = as.factor(Disease),
  T_clusters = as.factor(T_clusters),
  Response_type_total = factor(Response_type_total, levels = c("Low", "Medium", "High"))
)


#### Treatment Table #####
df_Cluster_response_cat <- read_csv(file.path(data_dir, "Metadata/Sample_Response_cluster.csv"))


sample_list <- data_talyies_full %>%
  dplyr::filter(Day == "D6", Disease %in% c("FL", "DLBCL", "tFL")) %>%
  distinct(Sample_code)

table_treatment <- read_csv(file.path(data_dir, "Metadata/liste_combo.csv")) %>%
  mutate(Treatment = apply(.[c("A", "B", "C")], 1, function(row) {
    paste(na.omit(row), collapse = " + ")
  })) %>%
  mutate(Treatment_reorder = apply(.[c("E", "F", "G")], 1, function(row) {
    paste(na.omit(row), collapse = " + ")
  })) %>%
  select(Treatment_reorder, Treatment, Treatment_type, Treatment_Cat) %>%
  # dplyr::filter(!grepl("TCB 10 nM", Treatment_reorder)) %>% #dplyr::filter les Treatments
  dplyr::filter(!grepl("GA101", Treatment_reorder) & !grepl("IL2v", Treatment_reorder) & !grepl("TCB 10 nM", Treatment_reorder) &
                  !grepl("ZB2", Treatment_reorder)) %>% # dplyr::filter les Treatments
  mutate(Treatment_reorder = factor(Treatment_reorder, levels = sort(unique(Treatment_reorder))))

all_combinations <- crossing(sample_list$Sample_code, table_treatment$Treatment_reorder) %>%
  rename(Sample_code = "sample_list$Sample_code", Treatment_reorder = "table_treatment$Treatment_reorder") %>%
  left_join(table_treatment, relationship = "many-to-many") %>%
  left_join(data_sample_info_complete) %>%
  dplyr::filter(Disease %in% c("FL", "DLBCL", "tFL") & Screening == TRUE) %>%
  select(Treatment, Sample_code, Treatment_reorder, Treatment_type, Treatment_Cat) %>%
  mutate(
    Treatment_Cat = factor(Treatment_Cat, levels = (c("UT", "αCD20-TCB", "Inhibiteur_CP", "Co_Activator", "ADC"))),
    Treatment_reorder = gsub(",", ".", Treatment_reorder),
    Treatment_reorder = gsub("αCD19-41BBL 0.1795 µg/mL", "αCD19-4-1BBL", Treatment_reorder),
    Treatment_reorder = gsub("αCD19-CD28 0.1464 µg/mL", "αCD19-CD28", Treatment_reorder),
    Treatment_reorder = gsub("αPDL1 10 µg/mL", "αPD-L1", Treatment_reorder),
    Treatment_reorder = gsub("αCD79-ct 1 µg/mL", "αCD79-ct", Treatment_reorder),
    Treatment_reorder = gsub("αCD79-MMAE 1 µg/mL", "αCD79-MMAE", Treatment_reorder),
    Treatment_reorder = gsub("αPD1-TIM3 1 µg/mL", "αPD-1-TIM-3", Treatment_reorder),
    Treatment_reorder = gsub("LAG3", "LAG-3", Treatment_reorder),
    Treatment_reorder = gsub("TIM3", "TIM-3", Treatment_reorder),
    Treatment_reorder = gsub("PD1", "PD-1", Treatment_reorder),
    Treatment_reorder = gsub("PDL", "PD-L", Treatment_reorder),
    Treatment_reorder = gsub("41BB ", "4-1BB", Treatment_reorder),
    Treatment_reorder = gsub("41BBL ", "4-1BBL", Treatment_reorder)
  )


#### Themes ####
theme_blood <- function() {
  theme_classic() + # Fond classique et épuré
    theme(
      axis.text.x = element_text(size = 6, face = "bold", color = "black"), # Labels des axes X
      axis.text.y = element_text(size = 6, face = "bold", color = "black"),
      axis.title = element_text(size = 8, face = "bold"),
      axis.title.y = element_text(size = 8, face = "bold", color = "black"), # Labels des axes X
      legend.title = element_text(size = 7, face = "bold"), # Titre de la légende
      legend.text = element_text(size = 6, face = "bold"),
      legend.box.spacing = unit(0.1, "cm"),
      legend.key.size = unit(0.4, "cm"),
      plot.margin = margin(
        t = 0.6, # Top margin
        r = 0.5, # Right margin
        b = 0, # Bottom margin
        l = 1
      ), # Left margin
      plot.title = element_text(hjust = 0.5, size = 9, face = "bold"), # Titre du graphique centré
      strip.placement = "outside", # Strips à l'extérieur
      strip.background = element_blank(), # Suppression du cadre autour des strips
      strip.text = element_text(size = 6, face = "bold", angle = 0)
    )
}

