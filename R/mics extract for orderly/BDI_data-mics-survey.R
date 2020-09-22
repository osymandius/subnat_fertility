library(pointr)
library(tidyverse)
library(here)

source(here("R/inputs.R"))

iso3_current <- "BDI"

naomi_data_path <- "~/Imperial College London/HIV Inference Group - WP - Documents/Analytical datasets/naomi-data/"

int <- read_sf("~/Imperial College London/HIV Inference Group - WP - Documents/Analytical datasets/naomi-data/BDI/data/bdi_areas.geojson") %>%
  mutate(naomi_level = ifelse(area_level %in% c(0,1), FALSE, TRUE))

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

colnames(df[["BDI2000MICS"]]$BUwm) <- tolower(colnames(df[["BDI2000MICS"]]$BUwm))
colnames(df[["BDI2000MICS"]]$BUch) <- tolower(colnames(df[["BDI2000MICS"]]$BUch))
colnames(df[["BDI2000MICS"]]$BUHH) <- tolower(colnames(df[["BDI2000MICS"]]$BUHH))

# No women weight??
df[["BDI2000MICS"]]$BUwm <- df[["BDI2000MICS"]]$BUwm %>%
  select("wihhno", "wiclno", "wilnno", "cmcdoi", "wdob", "wmweight") %>%
  rename(cluster = wiclno, hh_number = wihhno, line_number = wilnno, doi = cmcdoi) %>%
  filter(!is.na(wdob), !is.na(cluster), !is.na(hh_number), !is.na(line_number), !is.na(doi)) %>%
  arrange(cluster, hh_number, line_number) %>%
  mutate(unique_id = group_indices(., cluster, hh_number, line_number))

df[["BDI2000MICS"]]$BUch <- df[["BDI2000MICS"]]$BUch %>%
  select(hid, chclno, ln, cdob) %>%
  rename(cluster = chclno, hh_number = hid, line_number = ln)

df[["BDI2000MICS"]]$BUHH <- df[["BDI2000MICS"]]$BUHH %>%
  select(hi1, hi2, hi7) %>%
  rename(cluster = hi1, hh_number = hi2, mics_area_name = hi7)

### 2006

colnames(df[["BDI2005MICS"]]$wm) <- tolower(colnames(df[["BDI2005MICS"]]$wm))
colnames(df[["BDI2005MICS"]]$ch) <- tolower(colnames(df[["BDI2005MICS"]]$ch))
colnames(df[["BDI2005MICS"]]$hh) <- tolower(colnames(df[["BDI2005MICS"]]$hh))

df[["BDI2005MICS"]]$wm <- df[["BDI2005MICS"]]$wm %>%
  select("hh1", "hh2", "wm4", "cmcdoiw", "wdob", "wmweight") %>%
  rename(cluster = hh1, hh_number = hh2, line_number = wm4, doi = cmcdoiw) %>%
  filter(!is.na(wdob), !is.na(cluster), !is.na(hh_number), !is.na(line_number), !is.na(doi)) %>%
  arrange(cluster, hh_number, line_number) %>%
  mutate(unique_id = group_indices(., cluster, hh_number, line_number))

df[["BDI2005MICS"]]$ch <- df[["BDI2005MICS"]]$ch %>%
  select(hh2, hh1, uf6, cdob) %>%
  rename(cluster = hh1, hh_number = hh2, line_number = uf6)

df[["BDI2005MICS"]]$hh <- df[["BDI2005MICS"]]$hh %>%
  select(hh1, hh2, hh1a) %>%
  rename(cluster = hh1, hh_number = hh2, mics_area_name = hh1a)

lab_df <- data.frame("mics_area_name" = attr(df[["BDI2005MICS"]]$hh$mics_area_name, "labels"),
                     "mics_area_name_label" = str_to_title(names(attr(df[["BDI2005MICS"]]$hh$mics_area_name, "labels")))
)

df[["BDI2005MICS"]]$hh %>%
  left_join(lab_df) %>%
  left_join(areas_long %>% filter(area_level == 1),by=c("mics_area_name_label" = "area_name")) %>%
  filter(is.na(area_id)) %>%
  .$mics_area_name_label %>%
  unique

df[["BDI2005MICS"]]$hh <- df[["BDI2005MICS"]]$hh %>%
  left_join(lab_df) %>%
  mutate(
    mics_area_name_label = case_when(
      mics_area_name_label == "Karusi" ~ "Karuzi",
      mics_area_name_label == "Bururi Rural" ~ "Bururi",
      TRUE ~ mics_area_name_label
    )
  ) %>%
  left_join(areas_long, by=c("mics_area_name_label" = "area_name"))

df[["BDI2005MICS"]]$hh %>%
  filter(is.na(area_id)) %>%
  .$mics_area_name_label %>%
  unique

df[["BDI2005MICS"]]$wm <- df[["BDI2005MICS"]]$hh %>%
  select(cluster, hh_number, area_id) %>%
  left_join(df[["BDI2005MICS"]]$wm) %>%
  filter(!is.na(unique_id))

bh_df <- df[["BDI2005MICS"]]$wm %>%
  left_join(df[["BDI2005MICS"]]$ch) %>%
  filter(!is.na(cdob)) %>%
  select(unique_id, cdob)


     