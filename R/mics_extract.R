library(pointr)
library(tidyverse)
library(sf)
devtools::load_all("~/Documents/GitHub/naomi")
source(here("R/inputs.R"))

iso3_current <- c("NGA")
naomi_data_path <- "~/Imperial College London/HIV Inference Group - Documents/Analytical datasets/naomi-data/"
list2env(make_areas_population(iso3_current, naomi_data_path, full = FALSE, population=FALSE, wide = FALSE, boundaries = FALSE), globalenv())

mics_indicators <- read_csv(here::here("input_data/MICS_indicators.csv")) %>%
  pivot_longer(-c(label, id, filetype))

mics_rds_path <- "~/Imperial College London/HIV Inference Group - Documents/Data/household surveys/MICS/datasets/archive/mics_rds" %>%
  list.files(full.names = TRUE) %>%
  sort() %>%
  head(1)

catalogue <- read_csv(file.path(mics_rds_path, "mics_survey_catalogue.csv"))
#' Identifying national surveys only based on location_prefix and manual filtering. 
#' Would be better to update the table with a national vs. subnational flag...
surveys <- catalogue %>%
  filter(region %in% c("Eastern and Southern Africa", "West and Central Africa"),
         datasets == "Available",
         round != "MICS2",
         location_prefix %in% iso3_current,
         !country %in% c("Kenya (Mombasa Informal Settlements)",
                         "Kenya (Nyanza Province)"),
         !survey_id %in% c("CIV2016MICS", "MOZ2008MICS")) %>%
  arrange(survey_id)

dataset_paths <- file.path(mics_rds_path, "mics_datasets_rds",
                           paste0(tolower(surveys$survey_id), ".rds")) %>%
  setNames(surveys$survey_id)

raw <- lapply(dataset_paths, readRDS) %>%
  lapply("[", c("wm", "bh", "hh"))

svy_df <- Map(function(names, raw) {
  
  indicators <- mics_indicators %>%
    filter(name == names)
    
    wm <- raw$wm
    colnames(wm) <- tolower(colnames(wm))
    wm <- wm %>%
      select(filter(mics_indicators, name == names, filetype == "wm")$value)
    colnames(wm) <- filter(mics_indicators, name == names, filetype == "wm")$id
    
    bh <- raw$bh
    colnames(bh) <- tolower(colnames(bh))
    bh <- bh %>%
      select(filter(mics_indicators, name == names, filetype == "bh")$value)
    colnames(bh) <- filter(mics_indicators, name == names, filetype == "bh")$id

  
  hh <- raw$hh
  colnames(hh) <- tolower(colnames(hh))
  hh <- hh %>%
    select(filter(mics_indicators, name == names, filetype == "hh")$value)
  colnames(hh) <- filter(mics_indicators, name == names, filetype == "hh")$id
  
  df <- list()
  df$wm <- wm %>%
    mutate(survey_id = names) %>%
    filter(!is.na(wdob), !is.na(cluster), !is.na(hh_number), !is.na(line_number), !is.na(doi)) %>%
    arrange(cluster, hh_number, line_number) %>%
    mutate(unique_id = group_indices(., cluster, hh_number, line_number))
  
  df$bh <- bh %>%
    mutate(survey_id = names)
  
  df$hh <- hh %>%
    mutate(survey_id = names)

  return(df)
    
}, names = names(raw), raw)

svy_wm <- svy_df %>%
  lapply("[[", "wm")

svy_hh <- svy_df %>%
  lapply("[[", "hh") %>%
  lapply(function(x) {
    x %>%
      left_join(data.frame(mics_area_name = attr(x$mics_area_name, "labels"), mics_area_name_label = str_to_title(names(attr(x$mics_area_name, "labels")))))
    
  })

svy_bh <- svy_df %>%
  lapply("[[", "bh")

## Check for mismatched names. Change names where area_id is NA to match names with area_ids
lapply(svy_hh, function(x) {
  
  lvl <- as.numeric(filter(mics_indicators, name == unique(x$survey_id), id == "mics_area_level")$value)
  iso <- substr(unique(x$survey_id), 0, 3)
  areas_long <- areas_long %>%
    filter(iso3 == iso,
           area_level == lvl)
  x <- x %>%
    full_join(areas_long, by=c("mics_area_name_label" = "area_name"))
  
  return(x)
  
}) %>%
  bind_rows %>%
  select(area_id, survey_id, mics_area_name_label) %>%
  unique() %>%
  filter(is.na(area_id) | is.na(survey_id))

svy_hh <- svy_hh %>%
  bind_rows() %>%
  mutate(
    mics_area_name_label = as.character(mics_area_name_label),
    mics_area_name_label = case_when(
      ## Lesotho
      mics_area_name_label == "Botha-Bothe" ~ "Butha-Buthe",
      mics_area_name_label == "Mohales Hoek" ~ "Mohale's Hoek",
      mics_area_name_label == "Qachas Nek" ~ "Qacha's Nek",
      ## Malawi
      mics_area_name_label == "Nkhata Bay" ~ "Nkhatabay",
      ## Mozambique
      mics_area_name_label == "Maputo Cidade" ~ "Cidade de Maputo",
      ## Nigeria
      mics_area_name_label == "Fct Abuja" ~ "FCT",
      TRUE ~ mics_area_name_label
    )
  ) %>%
  group_split(survey_id)

## Check again!
#' MWI: Likoma not sampled in either MICS survey
#' MWI: Mzimba present as merged Mzimba North & South. Area level 3 code - fix later on.
#' MWI: Neno not sampled - note that 'District Mwanza in 2003 was split into two districts, Neno and Mwanza'. Need spatial data from MWI to check definition of Mwanza in MWI2006MICS

lapply(svy_hh, function(x) {
  
  lvl <- as.numeric(filter(mics_indicators, name == unique(x$survey_id), id == "mics_area_level")$value)
  iso <- substr(unique(x$survey_id), 0, 3)
  areas_long <- areas_long %>%
    filter(iso3 == iso,
           area_level == lvl)
  x <- x %>%
    full_join(areas_long, by=c("mics_area_name_label" = "area_name"))
  return(x)
  
}) %>%
  bind_rows %>%
  select(area_id, survey_id, mics_area_name_label) %>%
  unique() %>%
  filter(is.na(area_id) | is.na(survey_id))

## Now left join to areas_long

svy_hh <- lapply(svy_hh, function(x) {
  
  lvl <- as.numeric(filter(mics_indicators, name == unique(x$survey_id), id == "mics_area_level")$value)
  iso <- substr(unique(x$survey_id), 0, 3)
  areas_long <- areas_long %>%
    filter(iso3 == iso,
           area_level == lvl)
  x <- x %>%
    left_join(areas_long, by=c("mics_area_name_label" = "area_name"))
  return(x)
  
}) %>%
  bind_rows() %>%
  mutate(area_id = case_when(
    mics_area_name_label == "Mzimba" ~ "MWI_3_05",
    TRUE ~ area_id
  )) %>%
  group_split(survey_id)

svy_wm <- Map(function(svy_wm, svy_hh) {
  
  svy_wm %>% 
    left_join(svy_hh %>% select(cluster, hh_number, area_id))
  
}, svy_wm, svy_hh)

bh_df <- Map(function(svy_bh, svy_wm, y) {
  
  if(y == "MOZ2008MICS") {
    
    
  }
  svy_wm %>%
    select(cluster, hh_number, line_number, unique_id) %>%
    left_join(svy_bh %>% select(cluster, hh_number, line_number, cdob)) %>%
    select(unique_id, cdob) %>%
    filter(!is.na(cdob))
}, svy_bh, svy_wm, y = names(svy_bh))

mics_dat <- list()

mics_dat$bh_df <- bh_df %>%
  setNames("NGA2016MICS")

mics_dat$wm <- svy_wm %>%
  setNames("NGA2016MICS")


saveRDS(dat, "input_data/mics_extract.rds")
