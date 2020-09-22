library(pointr)
library(tidyverse)
library(here)

source(here("R/inputs.R"))

iso3_current <- "KEN"

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

colnames(df[["KEN2000MICS"]]$wmKE) <- tolower(colnames(df[["KEN2000MICS"]]$wmKE))
colnames(df[["KEN2000MICS"]]$chKE) <- tolower(colnames(df[["KEN2000MICS"]]$chKE))
colnames(df[["KEN2000MICS"]]$hhKE) <- tolower(colnames(df[["KEN2000MICS"]]$hhKE))

# No women weight??
df[["KEN2000MICS"]]$wmKE <- df[["KEN2000MICS"]]$wmKE %>%
  select("hi1", "hi2", "hl", "cmcdoi", "wdob", "wmweight") %>%
  rename(cluster = hi1, hh_number = hi2, line_number = hl, doi = cmcdoi) %>%
  filter(!is.na(wdob), !is.na(cluster), !is.na(hh_number), !is.na(line_number), !is.na(doi)) %>%
  arrange(cluster, hh_number, line_number) %>%
  mutate(unique_id = group_indices(., cluster, hh_number, line_number))

df[["KEN2000MICS"]]$chKE <- df[["KEN2000MICS"]]$chKE %>%
  select(hid, chclno, ln, cdob) %>%
  rename(cluster = chclno, hh_number = hid, line_number = ln)

df[["KEN2000MICS"]]$hhKE <- df[["KEN2000MICS"]]$hhKE %>%
  select(hi1, hi2, hi7) %>%
  rename(cluster = hi1, hh_number = hi2, mics_area_name = hi7)

### 2006

colnames(df[["KEN2005MICS"]]$wm) <- tolower(colnames(df[["KEN2005MICS"]]$wm))
colnames(df[["KEN2005MICS"]]$ch) <- tolower(colnames(df[["KEN2005MICS"]]$ch))
colnames(df[["KEN2005MICS"]]$hh) <- tolower(colnames(df[["KEN2005MICS"]]$hh))

df[["KEN2005MICS"]]$wm <- df[["KEN2005MICS"]]$wm %>%
  select("hh1", "hh2", "wm4", "cmcdoiw", "wdob", "wmweight") %>%
  rename(cluster = hh1, hh_number = hh2, line_number = wm4, doi = cmcdoiw) %>%
  filter(!is.na(wdob), !is.na(cluster), !is.na(hh_number), !is.na(line_number), !is.na(doi)) %>%
  arrange(cluster, hh_number, line_number) %>%
  mutate(unique_id = group_indices(., cluster, hh_number, line_number))

df[["KEN2005MICS"]]$ch <- df[["KEN2005MICS"]]$ch %>%
  select(hh2, hh1, uf6, cdob) %>%
  rename(cluster = hh1, hh_number = hh2, line_number = uf6)

df[["KEN2005MICS"]]$hh <- df[["KEN2005MICS"]]$hh %>%
  select(hh1, hh2, hh1a) %>%
  rename(cluster = hh1, hh_number = hh2, mics_area_name = hh1a)

lab_df <- data.frame("mics_area_name" = attr(df[["KEN2005MICS"]]$hh$mics_area_name, "labels"),
                     "mics_area_name_label" = str_to_title(names(attr(df[["KEN2005MICS"]]$hh$mics_area_name, "labels")))
)

df[["KEN2005MICS"]]$hh %>%
  left_join(lab_df) %>%
  left_join(areas_long %>% filter(area_level == 1),by=c("mics_area_name_label" = "area_name")) %>%
  filter(is.na(area_id)) %>%
  .$mics_area_name_label %>%
  unique

df[["KEN2005MICS"]]$hh <- df[["KEN2005MICS"]]$hh %>%
  left_join(lab_df) %>%
  mutate(
    mics_area_name_label = case_when(
      mics_area_name_label == "Karusi" ~ "Karuzi",
      mics_area_name_label == "Bururi Rural" ~ "Bururi",
      TRUE ~ mics_area_name_label
    )
  ) %>%
  left_join(areas_long, by=c("mics_area_name_label" = "area_name"))

df[["KEN2005MICS"]]$hh %>%
  filter(is.na(area_id)) %>%
  .$mics_area_name_label %>%
  unique

df[["KEN2005MICS"]]$wm <- df[["KEN2005MICS"]]$hh %>%
  select(cluster, hh_number, area_id) %>%
  left_join(df[["KEN2005MICS"]]$wm) %>%
  filter(!is.na(unique_id))

bh_df <- df[["KEN2005MICS"]]$wm %>%
  left_join(df[["KEN2005MICS"]]$ch) %>%
  filter(!is.na(cdob)) %>%
  select(unique_id, cdob)


     