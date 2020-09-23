library(pointr)
library(tidyverse)
library(here)

source(here("R/inputs.R"))

iso3_current <- "COD"

naomi_data_path <- "~/Imperial College London/HIV Inference Group - WP - Documents/Analytical datasets/naomi-data/"

# int <- read_sf("~/Imperial College London/HIV Inference Group - WP - Documents/Analytical datasets/naomi-data/COD/data/cod_areas.geojson") %>%
#   mutate(naomi_level = ifelse(area_level == 3, TRUE, FALSE))
# 
# st_write(int, "~/Imperial College London/HIV Inference Group - WP - Documents/Analytical datasets/naomi-data/COD/data/cod_areas.geojson", delete_dsn = TRUE)

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


### 2001

colnames(df[["COD2001MICS"]]$wmto) <- tolower(colnames(df[["COD2001MICS"]]$wmto))
colnames(df[["COD2001MICS"]]$chto) <- tolower(colnames(df[["COD2001MICS"]]$chto))
colnames(df[["COD2001MICS"]]$hhto) <- tolower(colnames(df[["COD2001MICS"]]$hhto))


# ????

### 2006

colnames(df[["COD2010MICS"]]$wm) <- tolower(colnames(df[["COD2010MICS"]]$wm))
colnames(df[["COD2010MICS"]]$ch) <- tolower(colnames(df[["COD2010MICS"]]$ch))
colnames(df[["COD2010MICS"]]$hh) <- tolower(colnames(df[["COD2010MICS"]]$hh))

df[["COD2010MICS"]]$wm <- df[["COD2010MICS"]]$wm %>%
  select("hh1", "hh2", "wm4", "wdoi", "wdob", "wmweight") %>%
  rename(cluster = hh1, hh_number = hh2, line_number = wm4, doi = wdoi) %>%
  filter(!is.na(wdob), !is.na(cluster), !is.na(hh_number), !is.na(line_number), !is.na(doi)) %>%
  arrange(cluster, hh_number, line_number) %>%
  mutate(unique_id = group_indices(., cluster, hh_number, line_number))

df[["COD2010MICS"]]$ch <- df[["COD2010MICS"]]$ch %>%
  select(hh2, hh1, uf6, cdob) %>%
  rename(cluster = hh1, hh_number = hh2, line_number = uf6)

df[["COD2010MICS"]]$hh <- df[["COD2010MICS"]]$hh %>%
  select(hh1, hh2, hh7) %>%
  rename(cluster = hh1, hh_number = hh2, mics_area_name = hh7)

lab_df <- data.frame("mics_area_name" = attr(df[["COD2010MICS"]]$hh$mics_area_name, "labels"),
                     "mics_area_name_label" = str_to_title(names(attr(df[["COD2010MICS"]]$hh$mics_area_name, "labels")))
)

df[["COD2010MICS"]]$hh %>%
  left_join(lab_df) %>%
  left_join(areas_long %>% filter(area_level == 1),by=c("mics_area_name_label" = "area_name")) %>%
  filter(is.na(area_id)) %>%
  .$mics_area_name_label %>%
  unique

df[["COD2010MICS"]]$hh <- df[["COD2010MICS"]]$hh %>%
  left_join(lab_df) %>%
  mutate(
    mics_area_name_label = case_when(
      mics_area_name_label == "Bas Congo" ~ "Kongo-Central",
      mics_area_name_label == "Katanga" ~ "Tanganyika",
      mics_area_name_label == "Katanga" ~ "Haut-Lomami",
      mics_area_name_label == "Katanga" ~ "Lualaba",
      mics_area_name_label == "Katanga" ~ "Haut-Katanga",
      mics_area_name_label == "Sud Kivu" ~ "Sud-Kivu",
      mics_area_name_label == "Nord Kivu" ~ "Nord-Kivu",
      mics_area_name_label == "Kasai Oriental" ~ "Kasaï-Oriental",
      mics_area_name_label == "Kasai Oriental" ~ "Sankuru",
      mics_area_name_label == "Kasai Oriental" ~ "Kasaï-Oriental",
      
      mics_area_name_label == "Bandundu" ~ "Maï-Ndombe",
      mics_area_name_label == "Bandundu" ~ "Kwilu",
      mics_area_name_label == "Bandundu" ~ "Kwango",
      
      TRUE ~ mics_area_name_label
    )
  ) %>%
  left_join(areas_long %>% filter(area_level == 1), by=c("mics_area_name_label" = "area_name"))

df[["COD2010MICS"]]$hh %>%
  filter(is.na(area_id)) %>%
  .$mics_area_name_label %>%
  unique

df[["COD2010MICS"]]$wm <- df[["COD2010MICS"]]$hh %>%
  select(cluster, hh_number, area_id) %>%
  left_join(df[["COD2010MICS"]]$wm) %>%
  filter(!is.na(unique_id))

bh_df <- df[["COD2010MICS"]]$wm %>%
  left_join(df[["COD2010MICS"]]$ch) %>%
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
     