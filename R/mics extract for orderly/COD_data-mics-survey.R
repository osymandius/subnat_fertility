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

mics_dat <- list()
mics_dat$wm <- list()
mics_dat$bh <- list()


### 2001

colnames(df[["COD2001MICS"]]$wmto) <- tolower(colnames(df[["COD2001MICS"]]$wmto))
colnames(df[["COD2001MICS"]]$chto) <- tolower(colnames(df[["COD2001MICS"]]$chto))
colnames(df[["COD2001MICS"]]$hhto) <- tolower(colnames(df[["COD2001MICS"]]$hhto))


# ????

### 2010

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
  mutate(area_name = mics_area_name_label)
  # mutate(
  #   mics_area_name_label = case_when(
  #     mics_area_name_label == "Bas Congo" ~ "Kongo-Central",
  #     mics_area_name_label == "Katanga" ~ "Tanganyika",
  #     mics_area_name_label == "Katanga" ~ "Haut-Lomami",
  #     mics_area_name_label == "Katanga" ~ "Lualaba",
  #     mics_area_name_label == "Katanga" ~ "Haut-Katanga",
  #     mics_area_name_label == "Sud Kivu" ~ "Sud-Kivu",
  #     mics_area_name_label == "Nord Kivu" ~ "Nord-Kivu",
  #     mics_area_name_label == "Kasai Oriental" ~ "Kasaï-Oriental",
  #     mics_area_name_label == "Kasai Oriental" ~ "Sankuru",
  #     mics_area_name_label == "Kasai Oriental" ~ "Kasaï-Oriental",
  #     
  #     mics_area_name_label == "Bandundu" ~ "Maï-Ndombe",
  #     mics_area_name_label == "Bandundu" ~ "Kwilu",
  #     mics_area_name_label == "Bandundu" ~ "Kwango",
  #     
  #     TRUE ~ mics_area_name_label
  #   )
  # ) %>%
  # left_join(areas_long %>% filter(area_level == 1), by=c("mics_area_name_label" = "area_name"))

df[["COD2010MICS"]]$hh %>%
  filter(is.na(area_name)) %>%
  .$mics_area_name_label %>%
  unique

df[["COD2010MICS"]]$wm <- df[["COD2010MICS"]]$hh %>%
  select(cluster, hh_number, area_name) %>%
  left_join(df[["COD2010MICS"]]$wm) %>%
  mutate(survey_id = "COD2010MICS") %>%
  filter(!is.na(unique_id))

# Error in survey data - cluster x hh_number x line_number is not unique for unique id #2448. Remove one of the women
df[["COD2010MICS"]]$wm <- df[["COD2010MICS"]]$wm %>%
  add_count(unique_id) %>%
  filter(n == 1)

bh_df <- df[["COD2010MICS"]]$wm %>%
  left_join(df[["COD2010MICS"]]$ch) %>%
  filter(!is.na(cdob)) %>%
  select(unique_id, cdob)

mics_dat$wm$COD2010MICS <- df[["COD2010MICS"]]$wm
mics_dat$bh$COD2010MICS <- bh_df

### 2017

colnames(df[["COD2017MICS"]]$wm) <- tolower(colnames(df[["COD2017MICS"]]$wm))
colnames(df[["COD2017MICS"]]$bh) <- tolower(colnames(df[["COD2017MICS"]]$bh))
colnames(df[["COD2017MICS"]]$hh) <- tolower(colnames(df[["COD2017MICS"]]$hh))

df[["COD2017MICS"]]$wm <- df[["COD2017MICS"]]$wm %>%
  select("hh1", "hh2", "ln", "wdoi", "wdob", "wmweight") %>%
  rename(cluster = hh1, hh_number = hh2, line_number = ln, doi = wdoi) %>%
  filter(!is.na(wdob), !is.na(cluster), !is.na(hh_number), !is.na(line_number), !is.na(doi)) %>%
  arrange(cluster, hh_number, line_number) %>%
  mutate(unique_id = group_indices(., cluster, hh_number, line_number))

df[["COD2017MICS"]]$bh <- df[["COD2017MICS"]]$bh %>%
  select(hh2, hh1, wm3, bh4c) %>%
  rename(cluster = hh1, hh_number = hh2, line_number = wm3, cdob = bh4c)

df[["COD2017MICS"]]$hh <- df[["COD2017MICS"]]$hh %>%
  select(hh1, hh2, hh7) %>%
  rename(cluster = hh1, hh_number = hh2, mics_area_name = hh7)

lab_df <- data.frame("mics_area_name" = attr(df[["COD2017MICS"]]$hh$mics_area_name, "labels"),
                     "mics_area_name_label" = str_to_title(names(attr(df[["COD2017MICS"]]$hh$mics_area_name, "labels")))
)

df[["COD2017MICS"]]$hh %>%
  left_join(lab_df) %>%
  left_join(areas_long %>% filter(area_level == 1),by=c("mics_area_name_label" = "area_name")) %>%
  filter(is.na(area_id)) %>%
  .$mics_area_name_label %>%
  unique

df[["COD2017MICS"]]$hh <- df[["COD2017MICS"]]$hh %>%
  left_join(lab_df) %>%
  mutate(
    mics_area_name_label = case_when(
      mics_area_name_label == "Kongo Central" ~ "Kongo-Central",
      mics_area_name_label == "Sud Ubangi" ~ "Sud-Ubangi",
      mics_area_name_label == "Nord Ubangi" ~ "Nord-Ubangi",
      mics_area_name_label == "Bas Uele" ~ "Bas-Uele",
      mics_area_name_label == "Haut Uele" ~ "Haut-Uele",
      mics_area_name_label == "Nord Kivu" ~ "Nord-Kivu",
      mics_area_name_label == "Sud Kivu" ~ "Sud-Kivu",
      mics_area_name_label == "Haut Lomami" ~ "Haut-Lomami",
      mics_area_name_label == "Haut Katanga" ~ "Haut-Katanga",
      mics_area_name_label == "Kasai Oriental" ~ "Kasaï-Oriental",
      mics_area_name_label == "Kasai Central" ~ "Kasaï-Central",
      mics_area_name_label == "Maindombe" ~ "Maï-Ndombe",
      mics_area_name_label == "Kasai" ~ "Kasaï",
      TRUE ~ mics_area_name_label
    )
  ) %>%
  left_join(areas_long %>% filter(area_level == 1), by=c("mics_area_name_label" = "area_name"))

df[["COD2017MICS"]]$hh %>%
  filter(is.na(area_id)) %>%
  .$mics_area_name_label %>%
  unique

df[["COD2017MICS"]]$wm <- df[["COD2017MICS"]]$hh %>%
  select(cluster, hh_number, area_id) %>%
  left_join(df[["COD2017MICS"]]$wm) %>%
  mutate(survey_id = "COD2017MICS") %>%
  filter(!is.na(unique_id))

bh_df <- df[["COD2017MICS"]]$wm %>%
  left_join(df[["COD2017MICS"]]$bh) %>%
  filter(!is.na(cdob)) %>%
  select(unique_id, cdob)

mics_dat$wm$COD2017MICS <- df[["COD2017MICS"]]$wm
mics_dat$bh$COD2017MICS <- bh_df

saveRDS(mics_dat, "~/Documents/GitHub/subnat_fertility/countries/COD/data/COD_mics_dat.rds")
     