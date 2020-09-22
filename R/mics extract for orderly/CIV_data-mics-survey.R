library(pointr)
library(tidyverse)
library(here)

source(here("R/inputs.R"))

iso3_current <- "CIV"

naomi_data_path <- "~/Imperial College London/HIV Inference Group - WP - Documents/Analytical datasets/naomi-data/"

int <- read_sf("~/Imperial College London/HIV Inference Group - WP - Documents/Analytical datasets/naomi-data/CIV/data/civ_areas.geojson") %>%
  mutate(naomi_level = ifelse(area_level %in% c(0,1), FALSE, TRUE))

st_write(int, "~/Imperial College London/HIV Inference Group - WP - Documents/Analytical datasets/naomi-data/CIV/data/civ_areas.geojson", delete_dsn = TRUE)

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

colnames(df[["CIV2000MICS"]]$CIwm) <- tolower(colnames(df[["CIV2000MICS"]]$CIwm))
colnames(df[["CIV2000MICS"]]$CIch) <- tolower(colnames(df[["CIV2000MICS"]]$CIch))
colnames(df[["CIV2000MICS"]]$CIhh) <- tolower(colnames(df[["CIV2000MICS"]]$CIhh))

df[["CIV2000MICS"]]$CIwm <- df[["CIV2000MICS"]]$CIwm %>%
  select("wihhno", "wiclno", "wilnno", "cmcdoi", "wdob", "wmweight") %>%
  rename(cluster = wiclno, hh_number = wihhno, line_number = wilnno, doi = cmcdoi) %>%
  filter(!is.na(wdob), !is.na(cluster), !is.na(hh_number), !is.na(line_number), !is.na(doi)) %>%
  arrange(cluster, hh_number, line_number) %>%
  mutate(unique_id = group_indices(., cluster, hh_number, line_number))

df[["CIV2000MICS"]]$CIch <- df[["CIV2000MICS"]]$CIch %>%
  select(hid, chclno, ln, cdob) %>%
  rename(cluster = chclno, hh_number = hid, line_number = ln)

df[["CIV2000MICS"]]$CIhh <- df[["CIV2000MICS"]]$CIhh %>%
  select(hh1, hh2, hh7) %>%
  rename(cluster = hh1, hh_number = hh2, mics_area_name = hh7)

lab_df <- data.frame("mics_area_name" = attr(df[["CIV2000MICS"]]$CIhh$mics_area_name, "labels"),
                     "mics_area_name_label" = str_to_title(names(attr(df[["CIV2000MICS"]]$CIhh$mics_area_name, "labels")))
)

df[["CIV2000MICS"]]$CIhh %>%
  left_join(lab_df) %>%
  left_join(areas_long %>% filter(area_level == 1),by=c("mics_area_name_label" = "area_name")) %>%
  filter(is.na(area_id)) %>%
  .$mics_area_name_label %>%
  unique

df[["CIV2000MICS"]]$CIhh <- df[["CIV2000MICS"]]$CIhh %>%
  left_join(lab_df) %>%
  mutate(
    mics_area_name_label = case_when(
      mics_area_name_label == "Hauts-Bassins" ~ "Hauts Bassins",
      mics_area_name_label == "Cascade" ~ "Cascades",
      mics_area_name_label == "Plateau-Central" ~ "Plateau Central",
      TRUE ~ mics_area_name_label
    )
  ) %>%
  left_join(areas_long, by=c("mics_area_name_label" = "area_name"))

df[["CIV2000MICS"]]$CIhh %>%
  filter(is.na(area_id)) %>%
  .$mics_area_name_label %>%
  unique

df[["CIV2000MICS"]]$CIwm <- df[["CIV2000MICS"]]$CIhh %>%
  select(cluster, hh_number, area_id) %>%
  left_join(df[["CIV2000MICS"]]$CIwm) %>%
  filter(!is.na(unique_id))

bh_df <- df[["CIV2000MICS"]]$CIwm %>%
  left_join(df[["CIV2000MICS"]]$CIch) %>%
  filter(!is.na(cdob)) %>%
  select(unique_id, cdob)

### 2006

colnames(df[["CIV2006MICS"]]$wm) <- tolower(colnames(df[["CIV2006MICS"]]$wm))
colnames(df[["CIV2006MICS"]]$ch) <- tolower(colnames(df[["CIV2006MICS"]]$ch))
colnames(df[["CIV2006MICS"]]$hh) <- tolower(colnames(df[["CIV2006MICS"]]$hh))

df[["CIV2006MICS"]]$wm <- df[["CIV2006MICS"]]$wm %>%
  select("hh1", "hh2", "wm4", "cmcdoiw", "wdob", "wmweight") %>%
  rename(cluster = hh1, hh_number = hh2, line_number = wm4, doi = cmcdoiw) %>%
  filter(!is.na(wdob), !is.na(cluster), !is.na(hh_number), !is.na(line_number), !is.na(doi)) %>%
  arrange(cluster, hh_number, line_number) %>%
  mutate(unique_id = group_indices(., cluster, hh_number, line_number))

df[["CIV2006MICS"]]$ch <- df[["CIV2006MICS"]]$ch %>%
  select(hh2, hh1, uf6, cdob) %>%
  rename(cluster = hh1, hh_number = hh2, line_number = uf6)

df[["CIV2006MICS"]]$hh <- df[["CIV2006MICS"]]$hh %>%
  select(hh1, hh2, hh7) %>%
  rename(cluster = hh1, hh_number = hh2, mics_area_name = hh7)

lab_df <- data.frame("mics_area_name" = attr(df[["CIV2006MICS"]]$hh$mics_area_name, "labels"),
                     "mics_area_name_label" = str_to_title(names(attr(df[["CIV2006MICS"]]$hh$mics_area_name, "labels")))
)

df[["CIV2006MICS"]]$hh %>%
  left_join(lab_df) %>%
  left_join(areas_long %>% filter(area_level == 1),by=c("mics_area_name_label" = "area_name")) %>%
  filter(is.na(area_id)) %>%
  .$mics_area_name_label %>%
  unique

df[["CIV2006MICS"]]$hh <- df[["CIV2006MICS"]]$hh %>%
  left_join(lab_df) %>%
  mutate(
    mics_area_name_label = case_when(
      mics_area_name_label == "Hauts-Bassins" ~ "Hauts Bassins",
      mics_area_name_label == "Cascade" ~ "Cascades",
      mics_area_name_label == "Plateau-Central" ~ "Plateau Central",
      TRUE ~ mics_area_name_label
    )
  ) %>%
  left_join(areas_long, by=c("mics_area_name_label" = "area_name"))

df[["CIV2006MICS"]]$hh %>%
  filter(is.na(area_id)) %>%
  .$mics_area_name_label %>%
  unique

df[["CIV2006MICS"]]$wm <- df[["CIV2006MICS"]]$hh %>%
  select(cluster, hh_number, area_id) %>%
  left_join(df[["CIV2006MICS"]]$wm) %>%
  filter(!is.na(unique_id))

bh_df <- df[["CIV2006MICS"]]$wm %>%
  left_join(df[["CIV2006MICS"]]$ch) %>%
  filter(!is.na(cdob)) %>%
  select(unique_id, cdob)


     