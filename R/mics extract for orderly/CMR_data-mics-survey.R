library(pointr)
library(tidyverse)
library(here)

source(here("R/inputs.R"))

mics_indicators <- read_csv(here::here("input_data/MICS_indicators.csv")) %>%
  pivot_longer(-c(label, id, filetype))

iso3_current <- "CMR"

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

### 2000

colnames(df[["CMR2000MICS"]]$wmca) <- tolower(colnames(df[["CMR2000MICS"]]$wmca))
colnames(df[["CMR2000MICS"]]$chca) <- tolower(colnames(df[["CMR2000MICS"]]$chca))
colnames(df[["CMR2000MICS"]]$hhca) <- tolower(colnames(df[["CMR2000MICS"]]$hhca))

df[["CMR2000MICS"]]$wmca <- df[["CMR2000MICS"]]$wmca %>%
  select("wihhno", "wiclno", "wilnno", "cmcdoi", "wdob", "wmweight") %>%
  rename(cluster = wiclno, hh_number = wihhno, line_number = wilnno, doi = cmcdoi) %>%
  filter(!is.na(wdob), !is.na(cluster), !is.na(hh_number), !is.na(line_number), !is.na(doi)) %>%
  arrange(cluster, hh_number, line_number) %>%
  mutate(unique_id = group_indices(., cluster, hh_number, line_number))

df[["CMR2000MICS"]]$chca <- df[["CMR2000MICS"]]$chca %>%
  select(hid, chclno, ln, cdob) %>%
  rename(cluster = chclno, hh_number = hid, line_number = ln)

df[["CMR2000MICS"]]$hhca <- df[["CMR2000MICS"]]$hhca %>%
  select(hi1, hi2, hi7) %>%
  rename(cluster = hi1, hh_number = hi2, mics_area_name = hi7)

### 2006

colnames(df[["CMR2006MICS"]]$wm) <- tolower(colnames(df[["CMR2006MICS"]]$wm))
colnames(df[["CMR2006MICS"]]$ch) <- tolower(colnames(df[["CMR2006MICS"]]$ch))
colnames(df[["CMR2006MICS"]]$hh) <- tolower(colnames(df[["CMR2006MICS"]]$hh))

df[["CMR2006MICS"]]$wm <- df[["CMR2006MICS"]]$wm %>%
  select("hh1", "hh2", "ln", "cmcdoiw", "wdob", "wmweight") %>%
  rename(cluster = hh1, hh_number = hh2, line_number = ln, doi = cmcdoiw) %>%
  filter(!is.na(wdob), !is.na(cluster), !is.na(hh_number), !is.na(line_number), !is.na(doi)) %>%
  arrange(cluster, hh_number, line_number) %>%
  mutate(unique_id = group_indices(., cluster, hh_number, line_number))

df[["CMR2006MICS"]]$ch <- df[["CMR2006MICS"]]$ch %>%
  select(hh2, hh1, uf6, cdob) %>%
  rename(cluster = hh1, hh_number = hh2, line_number = uf6)

df[["CMR2006MICS"]]$hh <- df[["CMR2006MICS"]]$hh %>%
  select(hh1, hh2, hh7) %>%
  rename(cluster = hh1, hh_number = hh2, mics_area_name = hh7)

lab_df <- data.frame("mics_area_name" = attr(df[["CMR2006MICS"]]$hh$mics_area_name, "labels"),
                     "mics_area_name_label" = str_to_title(names(attr(df[["CMR2006MICS"]]$hh$mics_area_name, "labels")))
)

df[["CMR2006MICS"]]$hh %>%
  left_join(lab_df) %>%
  left_join(areas_long,by=c("mics_area_name_label" = "area_name")) %>%
  filter(is.na(area_id)) %>%
  .$mics_area_name_label %>%
  unique

df[["CMR2006MICS"]]$hh <- df[["CMR2006MICS"]]$hh %>%
  left_join(lab_df) %>%
  mutate(
    mics_area_name_label = case_when(
      mics_area_name_label == "Centre" ~ "Centre (sans Yaoundé)",
      mics_area_name_label == "Littoral" ~ "Littoral (sans Douala)",
      mics_area_name_label == "Extreme Nord" ~ "Extreme-Nord",
      mics_area_name_label == "Yaounde" ~ "Yaoundé",
      mics_area_name_label == "Nord Ouest" ~ "Nord-Ouest",
      mics_area_name_label == "Sud Ouest" ~ "Sud-Ouest",
      TRUE ~ mics_area_name_label
    )
  ) %>%
  left_join(areas_long, by=c("mics_area_name_label" = "area_name"))

df[["CMR2006MICS"]]$hh %>%
  filter(is.na(area_id)) %>%
  .$mics_area_name_label %>%
  unique

df[["CMR2006MICS"]]$wm <- df[["CMR2006MICS"]]$hh %>%
  select(cluster, hh_number, area_id) %>%
  left_join(df[["CMR2006MICS"]]$wm) %>%
  filter(!is.na(unique_id))

bh_df <- df[["CMR2006MICS"]]$wm %>%
  left_join(df[["CMR2006MICS"]]$ch) %>%
  filter(!is.na(cdob)) %>%
  select(unique_id, cdob)

### 2014

df["CMR2014MICS"] <- df["CMR2014MICS"] %>%
  lapply("[", c("wm", "bh", "hh"))

colnames(df[["CMR2014MICS"]]$wm) <- tolower(colnames(df[["CMR2014MICS"]]$wm))
colnames(df[["CMR2014MICS"]]$bh) <- tolower(colnames(df[["CMR2014MICS"]]$bh))
colnames(df[["CMR2014MICS"]]$hh) <- tolower(colnames(df[["CMR2014MICS"]]$hh))

df[["CMR2014MICS"]]$wm <- df[["CMR2014MICS"]]$wm %>%
  select("wm1", "wm2", "wm4", "wdoi", "wdob", "wmweight") %>%
  rename(cluster = wm1, hh_number = wm2, line_number = wm4, doi = wdoi) %>%
  filter(!is.na(wdob), !is.na(cluster), !is.na(hh_number), !is.na(line_number), !is.na(doi)) %>%
  arrange(cluster, hh_number, line_number) %>%
  mutate(unique_id = group_indices(., cluster, hh_number, line_number))

df[["CMR2014MICS"]]$bh <- df[["CMR2014MICS"]]$bh %>%
  select(hh1, hh2, ln, bh4c) %>%
  rename(cluster = hh1, hh_number = hh2, line_number = ln, cdob = bh4c)

df[["CMR2014MICS"]]$hh <- df[["CMR2014MICS"]]$hh %>%
  select(hh1, hh2, hh7) %>%
  rename(cluster = hh1, hh_number = hh2, mics_area_name = hh7)
  
lab_df <- data.frame("mics_area_name" = attr(df[["CMR2014MICS"]]$hh$mics_area_name, "labels"),
           "mics_area_name_label" = str_to_title(names(attr(df[["CMR2014MICS"]]$hh$mics_area_name, "labels")))
)

df[["CMR2014MICS"]]$hh %>%
  left_join(lab_df) %>%
  left_join(areas_long,by=c("mics_area_name_label" = "area_name")) %>%
  filter(is.na(area_id)) %>%
  .$mics_area_name_label %>%
  unique

df[["CMR2014MICS"]]$hh <- df[["CMR2014MICS"]]$hh %>%
  left_join(lab_df) %>%
  mutate(
    mics_area_name_label = case_when(
      mics_area_name_label == "Centre (Sans Yaoundã©)" ~ "Centre (sans Yaoundé)",
      mics_area_name_label == "Littoral (Sans Douala)" ~ "Littoral (sans Douala)",
      mics_area_name_label == "Extrãªme-Nord" ~ "Extreme-Nord",
      mics_area_name_label == "Yaoundã©" ~ "Yaoundé",
      TRUE ~ mics_area_name_label
    )
  ) %>%
  left_join(areas_long, by=c("mics_area_name_label" = "area_name"))

df[["CMR2014MICS"]]$hh %>%
  filter(is.na(area_id)) %>%
  .$mics_area_name_label %>%
  unique

df[["CMR2014MICS"]]$wm <- df[["CMR2014MICS"]]$hh %>%
  select(cluster, hh_number, area_id) %>%
  left_join(df[["CMR2014MICS"]]$wm) %>%
  filter(!is.na(unique_id))

bh_df <- df[["CMR2014MICS"]]$wm %>%
  left_join(df[["CMR2014MICS"]]$bh) %>%
  filter(!is.na(cdob)) %>%
  select(unique_id, cdob)


     