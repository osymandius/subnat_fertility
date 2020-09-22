library(pointr)
library(tidyverse)
library(here)

source(here("R/inputs.R"))

iso3_current <- "TGO"

naomi_data_path <- "~/Imperial College London/HIV Inference Group - WP - Documents/Analytical datasets/naomi-data/"

# int <- read_sf("~/Imperial College London/HIV Inference Group - WP - Documents/Analytical datasets/naomi-data/TGO/data/tgo_areas.geojson") %>%
#   mutate(naomi_level = ifelse(area_level == 2, TRUE, FALSE))
# 
# st_write(int, "~/Imperial College London/HIV Inference Group - WP - Documents/Analytical datasets/naomi-data/TGO/data/tgo_areas.geojson", delete_dsn = TRUE)

list2env(make_areas_population(iso3_current, naomi_data_path, full = FALSE, population=FALSE, wide = FALSE, boundaries = FALSE), globalenv())

mics_rds_path <- "~/Imperial College London/HIV Inference Group - WP - Documents/Data/household surveys/MICS/datasets/archive/mics_rds" %>%
  list.files(full.names = TRUE) %>%
  sort() %>%
  head(1)

mics_rds_path <- paste0(mics_rds_path, "/mics_datasets_rds")

paths <- file.path(mics_rds_path, grep(tolower(iso3_current), list.files(mics_rds_path, full.names = FALSE), value = TRUE))

survey_ids <- toupper(grep(tolower(iso3_current), list.files(mics_rds_path, full.names = FALSE), value = TRUE)) %>%
  str_sub(end=-5)

df <- lapply(paths, readRDS)
names(df) <- survey_ids


### 2000

colnames(df[["TGO2000MICS"]]$wmto) <- tolower(colnames(df[["TGO2000MICS"]]$wmto))
colnames(df[["TGO2000MICS"]]$chto) <- tolower(colnames(df[["TGO2000MICS"]]$chto))
colnames(df[["TGO2000MICS"]]$hhto) <- tolower(colnames(df[["TGO2000MICS"]]$hhto))

df[["TGO2000MICS"]]$wm <- df[["TGO2000MICS"]]$wm %>%
  select("wiclno", "wihhno", "wilnno", "cmcdoi", "wdob", "wmweight") %>%
  rename(cluster = wiclno, hh_number = wihhno, line_number = wilnno, doi = cmcdoi) %>%
  filter(!is.na(wdob), !is.na(cluster), !is.na(hh_number), !is.na(line_number), !is.na(doi)) %>%
  arrange(cluster, hh_number, line_number) %>%
  mutate(unique_id = group_indices(., cluster, hh_number, line_number))

df[["TGO2000MICS"]]$ch <- df[["TGO2000MICS"]]$ch %>%
  select(chclno, chhhno, chctno, cdob) %>%
  rename(cluster = chclno, hh_number = chhhno, line_number = chctno)

df[["TGO2000MICS"]]$hh <- df[["TGO2000MICS"]]$hh %>%
  select(hi1, hi2, hi7) %>%
  rename(cluster = hi1, hh_number = hi2, mics_area_name = hi7)

# ????

### 2006

colnames(df[["TGO2006MICS"]]$wm) <- tolower(colnames(df[["TGO2006MICS"]]$wm))
colnames(df[["TGO2006MICS"]]$ch) <- tolower(colnames(df[["TGO2006MICS"]]$ch))
colnames(df[["TGO2006MICS"]]$hh) <- tolower(colnames(df[["TGO2006MICS"]]$hh))

df[["TGO2006MICS"]]$wm <- df[["TGO2006MICS"]]$wm %>%
  select("hh1", "hh2", "wm4", "cmcdoiw", "wdob", "wmweight") %>%
  rename(cluster = hh1, hh_number = hh2, line_number = wm4, doi = cmcdoiw) %>%
  filter(!is.na(wdob), !is.na(cluster), !is.na(hh_number), !is.na(line_number), !is.na(doi)) %>%
  arrange(cluster, hh_number, line_number) %>%
  mutate(unique_id = group_indices(., cluster, hh_number, line_number))

df[["TGO2006MICS"]]$ch <- df[["TGO2006MICS"]]$ch %>%
  select(hh2, hh1, uf6, cdob) %>%
  rename(cluster = hh1, hh_number = hh2, line_number = uf6)

df[["TGO2006MICS"]]$hh <- df[["TGO2006MICS"]]$hh %>%
  select(hh1, hh2, hh7) %>%
  rename(cluster = hh1, hh_number = hh2, mics_area_name = hh7)

lab_df <- data.frame("mics_area_name" = attr(df[["TGO2006MICS"]]$hh$mics_area_name, "labels"),
                     "mics_area_name_label" = str_to_title(names(attr(df[["TGO2006MICS"]]$hh$mics_area_name, "labels")))
)

df[["TGO2006MICS"]]$hh %>%
  left_join(lab_df) %>%
  left_join(areas_long %>% filter(area_level == 1),by=c("mics_area_name_label" = "area_name")) %>%
  filter(is.na(area_id)) %>%
  .$mics_area_name_label %>%
  unique

df[["TGO2006MICS"]]$hh <- df[["TGO2006MICS"]]$hh %>%
  left_join(lab_df) %>%
  mutate(
    mics_area_name_label = case_when(
      mics_area_name_label == "Lomé" ~ "Lome",
      mics_area_name_label == "Maritime (Sans Lomé)" ~ "Maritime",
      TRUE ~ mics_area_name_label
    )
  ) %>%
  left_join(areas_long %>% filter(area_level == 1), by=c("mics_area_name_label" = "area_name"))

df[["TGO2006MICS"]]$hh %>%
  filter(is.na(area_id)) %>%
  .$mics_area_name_label %>%
  unique

df[["TGO2006MICS"]]$wm <- df[["TGO2006MICS"]]$hh %>%
  select(cluster, hh_number, area_id) %>%
  left_join(df[["TGO2006MICS"]]$wm) %>%
  filter(!is.na(unique_id))

bh_df <- df[["TGO2006MICS"]]$wm %>%
  left_join(df[["TGO2006MICS"]]$ch) %>%
  filter(!is.na(cdob)) %>%
  select(unique_id, cdob)

### 2016

colnames(df[["TGO2010MICS"]]$wm) <- tolower(colnames(df[["TGO2010MICS"]]$wm))
colnames(df[["TGO2010MICS"]]$ch) <- tolower(colnames(df[["TGO2010MICS"]]$ch))
colnames(df[["TGO2010MICS"]]$hh) <- tolower(colnames(df[["TGO2010MICS"]]$hh))

df[["TGO2010MICS"]]$wm <- df[["TGO2010MICS"]]$wm %>%
  select("hh1", "hh2", "ln", "wdoi", "wdob", "wmweight") %>%
  rename(cluster = hh1, hh_number = hh2, line_number = ln, doi = wdoi) %>%
  filter(!is.na(wdob), !is.na(cluster), !is.na(hh_number), !is.na(line_number), !is.na(doi)) %>%
  arrange(cluster, hh_number, line_number) %>%
  mutate(unique_id = group_indices(., cluster, hh_number, line_number))

df[["TGO2010MICS"]]$ch <- df[["TGO2010MICS"]]$ch %>%
  select(hh2, hh1, uf6, cdob) %>%
  rename(cluster = hh1, hh_number = hh2, line_number = uf6)

df[["TGO2010MICS"]]$hh <- df[["TGO2010MICS"]]$hh %>%
  select(hh1, hh2, hh7) %>%
  rename(cluster = hh1, hh_number = hh2, mics_area_name = hh7)

lab_df <- data.frame("mics_area_name" = attr(df[["TGO2010MICS"]]$hh$mics_area_name, "labels"),
                     "mics_area_name_label" = str_to_title(names(attr(df[["TGO2010MICS"]]$hh$mics_area_name, "labels")))
)

df[["TGO2010MICS"]]$hh %>%
  left_join(lab_df) %>%
  left_join(areas_long %>% filter(area_level == 1),by=c("mics_area_name_label" = "area_name")) %>%
  filter(is.na(area_id)) %>%
  .$mics_area_name_label %>%
  unique

df[["TGO2010MICS"]]$hh <- df[["TGO2010MICS"]]$hh %>%
  left_join(lab_df) %>%
  mutate(
    mics_area_name_label = case_when(
      mics_area_name_label == "Lomé" ~ "Lome",
      TRUE ~ mics_area_name_label
    )
  ) %>%
  left_join(areas_long %>% filter(area_level == 1), by=c("mics_area_name_label" = "area_name"))

df[["TGO2010MICS"]]$hh %>%
  filter(is.na(area_id)) %>%
  .$mics_area_name_label %>%
  unique

df[["TGO2010MICS"]]$wm <- df[["TGO2010MICS"]]$hh %>%
  select(cluster, hh_number, area_id) %>%
  left_join(df[["TGO2010MICS"]]$wm) %>%
  filter(!is.na(unique_id))

bh_df <- df[["TGO2010MICS"]]$wm %>%
  left_join(df[["TGO2010MICS"]]$ch) %>%
  filter(!is.na(cdob)) %>%
  select(unique_id, cdob)
     