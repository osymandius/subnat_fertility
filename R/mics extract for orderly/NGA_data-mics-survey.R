library(pointr)
library(tidyverse)
library(here)

source(here("R/inputs.R"))

iso3_current <- "NGA"

naomi_data_path <- "~/Imperial College London/HIV Inference Group - WP - Documents/Analytical datasets/naomi-data/"

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


### 2007

colnames(df[["NGA2007MICS"]]$wm) <- tolower(colnames(df[["NGA2007MICS"]]$wm))
colnames(df[["NGA2007MICS"]]$ch) <- tolower(colnames(df[["NGA2007MICS"]]$ch))
colnames(df[["NGA2007MICS"]]$hh) <- tolower(colnames(df[["NGA2007MICS"]]$hh))

df[["NGA2007MICS"]]$wm <- df[["NGA2007MICS"]]$wm %>%
  select("hh1", "hh2", "wm4", "cmcdoiw", "wdob", "wmweight") %>%
  rename(cluster = hh1, hh_number = hh2, line_number = wm4, doi = cmcdoiw) %>%
  filter(!is.na(wdob), !is.na(cluster), !is.na(hh_number), !is.na(line_number), !is.na(doi)) %>%
  arrange(cluster, hh_number, line_number) %>%
  mutate(unique_id = group_indices(., cluster, hh_number, line_number))

df[["NGA2007MICS"]]$ch <- df[["NGA2007MICS"]]$ch %>%
  select(hh2, hh1, uf6, cdob) %>%
  rename(cluster = hh1, hh_number = hh2, line_number = uf6)

df[["NGA2007MICS"]]$hh <- df[["NGA2007MICS"]]$hh %>%
  select(hh1, hh2, hh7) %>%
  rename(cluster = hh1, hh_number = hh2, mics_area_name = hh7)

lab_df <- data.frame("mics_area_name" = attr(df[["NGA2007MICS"]]$hh$mics_area_name, "labels"),
                     "mics_area_name_label" = str_to_title(names(attr(df[["NGA2007MICS"]]$hh$mics_area_name, "labels")))
)

df[["NGA2007MICS"]]$hh %>%
  left_join(lab_df) %>%
  left_join(areas_long %>% filter(area_level == 2),by=c("mics_area_name_label" = "area_name")) %>%
  filter(is.na(area_id)) %>%
  .$mics_area_name_label %>%
  unique

df[["NGA2007MICS"]]$hh <- df[["NGA2007MICS"]]$hh %>%
  left_join(lab_df) %>%
  mutate(
    mics_area_name_label = case_when(
      mics_area_name_label == "Plataeu" ~ "Plateau",
      mics_area_name_label == "Abuja Fct" ~ "FCT",
      mics_area_name_label == "Cross-Rivers" ~ "Cross River",
      mics_area_name_label == "Akwa-Ibom" ~ "Akwa Ibom",
      TRUE ~ mics_area_name_label
    )
  ) %>%
  left_join(areas_long %>% filter(area_level == 2), by=c("mics_area_name_label" = "area_name"))

df[["NGA2007MICS"]]$hh %>%
  filter(is.na(area_id)) %>%
  .$mics_area_name_label %>%
  unique

df[["NGA2007MICS"]]$wm <- df[["NGA2007MICS"]]$hh %>%
  select(cluster, hh_number, area_id) %>%
  left_join(df[["NGA2007MICS"]]$wm) %>%
  filter(!is.na(unique_id))

bh_df <- df[["NGA2007MICS"]]$wm %>%
  left_join(df[["NGA2007MICS"]]$ch) %>%
  filter(!is.na(cdob)) %>%
  select(unique_id, cdob)

### 2011

colnames(df[["NGA2011MICS"]]$wm) <- tolower(colnames(df[["NGA2011MICS"]]$wm))
colnames(df[["NGA2011MICS"]]$ch) <- tolower(colnames(df[["NGA2011MICS"]]$ch))
colnames(df[["NGA2011MICS"]]$hh) <- tolower(colnames(df[["NGA2011MICS"]]$hh))

df[["NGA2011MICS"]]$wm <- df[["NGA2011MICS"]]$wm %>%
  select("hh1", "hh2", "wm4", "wdoi", "wdob", "wmweight") %>%
  rename(cluster = hh1, hh_number = hh2, line_number = wm4, doi = wdoi) %>%
  filter(!is.na(wdob), !is.na(cluster), !is.na(hh_number), !is.na(line_number), !is.na(doi)) %>%
  arrange(cluster, hh_number, line_number) %>%
  mutate(unique_id = group_indices(., cluster, hh_number, line_number))

df[["NGA2011MICS"]]$ch <- df[["NGA2011MICS"]]$ch %>%
  select(hh2, hh1, uf6, cdob) %>%
  rename(cluster = hh1, hh_number = hh2, line_number = uf6)

df[["NGA2011MICS"]]$hh <- df[["NGA2011MICS"]]$hh %>%
  select(hh1, hh2, hh7) %>%
  rename(cluster = hh1, hh_number = hh2, mics_area_name = hh7)

lab_df <- data.frame("mics_area_name" = attr(df[["NGA2011MICS"]]$hh$mics_area_name, "labels"),
                     "mics_area_name_label" = str_to_title(names(attr(df[["NGA2011MICS"]]$hh$mics_area_name, "labels")))
)

df[["NGA2011MICS"]]$hh %>%
  left_join(lab_df) %>%
  left_join(areas_long %>% filter(area_level == 2),by=c("mics_area_name_label" = "area_name")) %>%
  filter(is.na(area_id)) %>%
  .$mics_area_name_label %>%
  unique

df[["NGA2011MICS"]]$hh <- df[["NGA2011MICS"]]$hh %>%
  left_join(lab_df) %>%
  mutate(
    mics_area_name_label = case_when(
      mics_area_name_label == "Fct (Abuja)" ~ "FCT",
      TRUE ~ mics_area_name_label
    )
  ) %>%
  left_join(areas_long %>% filter(area_level == 2), by=c("mics_area_name_label" = "area_name"))

df[["NGA2011MICS"]]$hh %>%
  filter(is.na(area_id)) %>%
  .$mics_area_name_label %>%
  unique

df[["NGA2011MICS"]]$wm <- df[["NGA2011MICS"]]$hh %>%
  select(cluster, hh_number, area_id) %>%
  left_join(df[["NGA2011MICS"]]$wm) %>%
  filter(!is.na(unique_id))

bh_df <- df[["NGA2011MICS"]]$wm %>%
  left_join(df[["NGA2011MICS"]]$ch) %>%
  filter(!is.na(cdob)) %>%
  select(unique_id, cdob)

### 2016

colnames(df[["NGA2016MICS"]]$wm) <- tolower(colnames(df[["NGA2016MICS"]]$wm))
colnames(df[["NGA2016MICS"]]$ch) <- tolower(colnames(df[["NGA2016MICS"]]$ch))
colnames(df[["NGA2016MICS"]]$hh) <- tolower(colnames(df[["NGA2016MICS"]]$hh))

df[["NGA2016MICS"]]$wm <- df[["NGA2016MICS"]]$wm %>%
  select("hh1", "hh2", "ln", "wdoi", "wdob", "wmweight") %>%
  rename(cluster = hh1, hh_number = hh2, line_number = ln, doi = wdoi) %>%
  filter(!is.na(wdob), !is.na(cluster), !is.na(hh_number), !is.na(line_number), !is.na(doi)) %>%
  arrange(cluster, hh_number, line_number) %>%
  mutate(unique_id = group_indices(., cluster, hh_number, line_number))

df[["NGA2016MICS"]]$ch <- df[["NGA2016MICS"]]$ch %>%
  select(hh2, hh1, uf6, cdob) %>%
  rename(cluster = hh1, hh_number = hh2, line_number = uf6)

df[["NGA2016MICS"]]$hh <- df[["NGA2016MICS"]]$hh %>%
  select(hh1, hh2, hh7) %>%
  rename(cluster = hh1, hh_number = hh2, mics_area_name = hh7)

lab_df <- data.frame("mics_area_name" = attr(df[["NGA2016MICS"]]$hh$mics_area_name, "labels"),
                     "mics_area_name_label" = str_to_title(names(attr(df[["NGA2016MICS"]]$hh$mics_area_name, "labels")))
)

df[["NGA2016MICS"]]$hh %>%
  left_join(lab_df) %>%
  left_join(areas_long %>% filter(area_level == 2),by=c("mics_area_name_label" = "area_name")) %>%
  filter(is.na(area_id)) %>%
  .$mics_area_name_label %>%
  unique

df[["NGA2016MICS"]]$hh <- df[["NGA2016MICS"]]$hh %>%
  left_join(lab_df) %>%
  mutate(
    mics_area_name_label = case_when(
      mics_area_name_label == "Fct Abuja" ~ "FCT",
      TRUE ~ mics_area_name_label
    )
  ) %>%
  left_join(areas_long %>% filter(area_level == 2), by=c("mics_area_name_label" = "area_name"))

df[["NGA2016MICS"]]$hh %>%
  filter(is.na(area_id)) %>%
  .$mics_area_name_label %>%
  unique

df[["NGA2016MICS"]]$wm <- df[["NGA2016MICS"]]$hh %>%
  select(cluster, hh_number, area_id) %>%
  left_join(df[["NGA2016MICS"]]$wm) %>%
  filter(!is.na(unique_id))

bh_df <- df[["NGA2016MICS"]]$wm %>%
  left_join(df[["NGA2016MICS"]]$ch) %>%
  filter(!is.na(cdob)) %>%
  select(unique_id, cdob)
     