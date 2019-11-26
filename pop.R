library(tidyverse)

setwd("~/Documents/GitHub/subnat_fertility")

source("fertility_funs.R")

saveRDS(national_pred, file="pred_df_national.rds")

asfr_pred_country_subnat <- readRDS("asfr_pred_country_subnat.rds")

mod_list <- list(readRDS("MWI_poisson_mod.rds"), readRDS("LSO_poisson_mod.rds"), readRDS("RWA_poisson_mod.rds"), readRDS("ZWE_poisson_mod.rds"), readRDS("UGA_poisson_mod.rds"))

### WORLD POP. Filter areas on cluster area. Convert MWI from level 5 to level 4.

iso3_code <- c("MWI", "RWA", "ZWE", "LSO", "UGA")

areas_long <- readRDS("~/Documents/GitHub/naomi-data/data/area_hierarchy.rds")
areas_wide <- readRDS("~/Documents/GitHub/naomi-data/data/areas_wide.rds")
clusters <- readRDS("oli_clusters.rds")

get_age_groups()

quarter_ids <- convert_quarter_id(2, c(2000:2020))

population_age_female <- load_population_agesex("~/Documents/GitHub/naomi-data/data/population_agesex_wide.rds", "~/Documents/GitHub/naomi-data/data/area_hierarchy.rds") %>%
  filter(source == "worldpop_wpp19", iso3 %in% iso3_code) %>%
  group_split(iso3) %>%
  lapply(interpolate_population_agesex, quarter_ids) %>%
  lapply(function(x) {
    x %>%
      filter(sex=="female", age_group_id %in% 4:10) %>%
      mutate(labs = quarter_year_labels(quarter_id)) %>%
      separate(labs, into=c(NA, "period"), sep=-4) %>%
      select(-c(quarter_id, sex)) %>%
      type.convert %>%
      left_join(get_age_groups() %>%
                  select(age_group_id, age_group_label))
  }) %>%
  bind_rows %>%
  rename(agegr = age_group_label)

# areas_wide %>%
#   filter(iso3 != "MWI") %>%
#   bind_rows(
#     areas_wide %>%
#       filter(iso3 =="MWI") %>%
#       mutate(area_id = id4) %>%
#       select(-c(name5, id5)) %>%
#       distinct
#   )

pop_areas <- population_age_female %>%
  left_join(areas_wide %>%
              filter(iso3 != "MWI") %>%
              bind_rows(
                areas_wide %>%
                  filter(iso3 =="MWI") %>%
                  mutate(area_id = id4) %>%
                  select(-c(name5, id5)) %>%
                  distinct
              )) %>%
  group_by(period, agegr, id0) %>%
  mutate(id0_agepop = sum(population))  %>%
  group_by(period, agegr, id1) %>%
  mutate(id1_agepop = sum(population)) %>%
  group_by(period, agegr, id2) %>%
  mutate(id2_agepop = sum(population)) %>%
  group_by(period, agegr, id3) %>%
  mutate(id3_agepop = sum(population)) %>%
  group_by(period, agegr, id4) %>%
  mutate(id4_agepop = sum(population)) %>%
  mutate(id2_agepop = ifelse(is.na(id2), NA, id2_agepop),
         id3_agepop = ifelse(is.na(id3), NA, id3_agepop),
         id4_agepop = ifelse(is.na(id4), NA, id4_agepop),
  ) %>%
  select(-c(name5, id5)) %>%
  filter(!is.na(name0))

# pop_areas_mwi <- pop_areas %>%
#   filter(iso3 == "MWI") %>%
#   mutate(area_id = id4,
#          area_name = name4,
#          area_level = 4
#   ) %>%
#   group_by(period, agegr, id4) %>%
#   mutate(population = sum(population)) %>%
#   distinct()
# 
# pop_areas <- pop_areas %>%
#   filter(iso3 != "MWI") %>%
#   bind_rows(pop_areas_mwi)
# 
# 
# worldpop_rep <- worldpop %>%
#   left_join(areas_wide, by=c("iso3", "id" = "area_id")) %>%
#   group_split(year) %>%
#   lapply(function(x) {
#     n <- nrow(x)
#     year_start <- unique(x$year)
#     x <- x %>%
#       bind_rows(., ., ., ., .) %>%
#       mutate(year = rep(year_start:(year_start+4), each=n))
#     return(x)
#   })  %>%
#   bind_rows
# 
# pop_areas_wide <- worldpop_rep %>%
#   filter(id %in% unique(clusters$geoloc_area_id)) %>%
#   group_by(year, age, id0) %>%
#   mutate(id0_agepop = sum(population))  %>%
#   group_by(year, age, id1) %>%
#   mutate(id1_agepop = sum(population)) %>%
#   group_by(year, age, id2) %>%
#   mutate(id2_agepop = sum(population)) %>%
#   group_by(year, age, id3) %>%
#   mutate(id3_agepop = sum(population)) %>%
#   group_by(year, age, id4) %>%
#   mutate(id4_agepop = sum(population)) %>%
#   mutate(id2_agepop = ifelse(is.na(id2), NA, id2_agepop),
#          id3_agepop = ifelse(is.na(id3), NA, id3_agepop),
#          id4_agepop = ifelse(is.na(id4), NA, id4_agepop),
#   ) %>%
#   select(-c(name5, id5))
# 
# pop_areas_wide_mwi <- pop_areas_wide %>%
#   filter(iso3 == "MWI") %>%
#   mutate(id = id4,
#          name = name4,
#          level = 4
#   ) %>%
#   group_by(year, age, id4) %>%
#   mutate(population = sum(population)) %>%
#   distinct()
# 
# pop_areas <- pop_areas_wide %>%
#   filter(iso3 != "MWI") %>%
#   bind_rows(pop_areas_wide_mwi)

############ Get ASFR, TFR, Births, Births by age at all administrative levels from model.

mod_results <- Map(get_mod_results, mod_list, asfr_pred_country_subnat, list(pop_areas), list(areas_long), list(population_age_female))

########

quarter_ids <- convert_quarter_id(2, c(2000:2020))

worldpop_U1 <- read_csv("WorldPop_agesex.csv", col_types = cols(X1 = col_skip(), X = col_skip(), sex = col_character())) %>%
  filter(startage==0, iso3 %in% iso3_code) %>%
  group_by(iso3, id, level, name, source, year, age, startage, agespan) %>%
  summarise(population = sum(population)) %>%
  ungroup %>%
  mutate(variable = "U1_pop",
         quarter_id = convert_quarter_id(2, year),
         age_group_id = 999,
         sex = "both") %>%
  rename(area_id = id, area_level = level, area_name = name) %>%
  select(iso3, area_id, area_level, area_name, sex, source, quarter_id, age_group_id, population)%>%
  group_split(iso3) %>%
  lapply(interpolate_population_agesex, quarter_ids) %>%
  lapply(function(x) {
    x %>%
      mutate(labs = quarter_year_labels(quarter_id)) %>%
      separate(labs, into=c(NA, "period"), sep=-4) %>%
      select(-c(quarter_id, sex)) %>%
      type.convert
  }) %>%
  bind_rows %>%
  select(-age_group_id) %>%
  mutate(variable = "U1_pop") %>%
  rename(val = population)
 
######### National comparison of ASFR and TFR vs WPP 2019

wpp_asfr <- read_excel("WPP2019_FERT_F07_AGE_SPECIFIC_FERTILITY.xlsx")

wpp_asfr <- wpp_asfr %>%
  separate(period, into=c("period", NA), sep="-") %>%
  type.convert %>%
  mutate(iso3 = countrycode(country, "country.name", "iso3c"),
         area_name = countrycode(iso3, "iso3c", "country.name"),
         area_id = iso3
  ) %>%
  select(-c(countrycode, country)) %>%
  melt(id=c("iso3", "area_name", "area_id", "period"), variable.name="agegr", value.name = "val") %>%
  filter(iso3 %in% iso3_code) %>%
  group_split(period) %>%
  lapply(function(x) {
    n <- nrow(x)
    year_start <- unique(x$period)
    x <- x %>%
      bind_rows(., ., ., ., .) %>%
      mutate(period = rep(year_start:(year_start+4), each=n))
    return(x)
  })  %>%
  bind_rows %>%
  mutate(val = as.numeric(val)/1000,
         source = "WPP2019",
         variable = "asfr",
         area_level = 0)

wpp_tfr <- read_excel("WPP2019_FERT_F04_TOTAL_FERTILITY(1).xlsx")

wpp_tfr <- wpp_tfr %>%
  melt(id="area_name", variable.name="period", value.name = "val") %>%
  separate(period, into=c("period", NA), sep="-") %>%
  type.convert() %>%
  mutate(iso3 = countrycode(area_name, "country.name", "iso3c")) %>%
  filter(iso3 %in% iso3_code) %>%
  group_split(period) %>%
  lapply(function(x) {
    n <- nrow(x)
    year_start <- unique(x$period)
    x <- x %>%
      bind_rows(., ., ., ., .) %>%
      mutate(period = rep(year_start:(year_start+4), each=n))
    return(x)
  })  %>%
  bind_rows %>%
  mutate(source = "WPP2019",
         area_level = 0,
         area_id = iso3,
         variable = "tfr")


######### GBD 2017 national ASFR and TFR

gbd_asfr <- read_csv("IHME_GBD_2017_FERT_ESTIMATES_1950_2017_Y2018M11D08.CSV", col_types = cols(age_group_id = col_skip(), location_id = col_skip(), measure_id = col_skip(), measure_name = col_skip(), metric_name = col_skip(), sex_id = col_skip(), sex_name = col_skip()))

gbd_asfr <- gbd_asfr %>%
  mutate(iso3 = countrycode(location_name, "country.name", "iso3c")) %>%
  filter(iso3 %in% iso3_code) %>%
  separate(age_group_name, into=c("agegr", "up"), sep=" to ") %>%
  mutate(agegr = paste0(agegr, "-", up),
         area_id = iso3,
         area_level = 0,
         source = "GBD2017",
         variable = "asfr") %>%
  select(-up) %>%
  rename(area_name = location_name, period = year_id) 


gbd_tfr <- gbd_asfr %>%
  group_by(iso3, area_name, area_id, area_level,  period, source) %>%
  summarise(val = 5*sum(val),
            upper = 5*sum(upper),
            lower = 5*sum(lower),
  ) %>%
  mutate(variable = "tfr")


#### GBD Births

gbd_2000_pop <- read_csv("IHME_GBD_2017_POP_2000_2004_Y2018M11D08.CSV", col_types = cols(location_id = col_skip(), measure_id = col_skip(), measure_name = col_skip(), metric_name = col_skip()))

gbd_2005_pop <- read_csv("IHME_GBD_2017_POP_2005_2009_Y2018M11D08.CSV", col_types = cols(location_id = col_skip(), measure_id = col_skip(), measure_name = col_skip(), metric_name = col_skip()))

gbd_2010_pop <- read_csv("IHME_GBD_2017_POP_2010_2014_Y2018M11D08.CSV", col_types = cols(location_id = col_skip(), measure_id = col_skip(), measure_name = col_skip(), metric_name = col_skip()))

gbd_2015_pop <- read_csv("IHME_GBD_2017_POP_2015_2017_Y2018M11D08.CSV", col_types = cols(location_id = col_skip(), measure_id = col_skip(), measure_name = col_skip(), metric_name = col_skip()))

gbd_2000_pop <- gbd_2000_pop %>%
  filter(age_group_id == 164, sex_id ==3)

gbd_2005_pop <- gbd_2005_pop %>%
  filter(age_group_id == 164, sex_id ==3)

gbd_2010_pop <- gbd_2010_pop %>%
  filter(age_group_id == 164, sex_id ==3)

gbd_2015_pop <- gbd_2015_pop %>%
  filter(age_group_id == 164, sex_id ==3)

gbd_births <- gbd_2000_pop %>%
  bind_rows(gbd_2005_pop, gbd_2010_pop, gbd_2015_pop) %>%
  select(-c(sex_id, sex_name, age_group_id, age_group_name)) %>%
  rename(area_name = location_name, period = year_id) %>%
  mutate(iso3 = countrycode(area_name, "country.name", "iso3c"),
         area_level = 0,
         area_id = iso3,
         source = "GBD2017",
         variable = "births") %>%
  filter(iso3 %in% iso3_code)

## Spectrum TFR and ASFR


iso_cod <- list("MWI", "UGA", "LSO", "RWA")

zim_spec <- list("~/Documents/Data/Spectrum files/2018 final/SSA/zim_Bulawayo_inc.PJNZ", "~/Documents/Data/Spectrum files/2018 final/SSA/zim_Harare_inc.PJNZ", "~/Documents/Data/Spectrum files/2018 final/SSA/zim_Manicaland_inc.PJNZ", "~/Documents/Data/Spectrum files/2018 final/SSA/zim_Mashonaland Central_inc.PJNZ", "~/Documents/Data/Spectrum files/2018 final/SSA/zim_Mashonaland East_inc.PJNZ", "~/Documents/Data/Spectrum files/2018 final/SSA/zim_Mashonaland West_inc.PJNZ", "~/Documents/Data/Spectrum files/2018 final/SSA/zim_Masvingo_inc.PJNZ", "~/Documents/Data/Spectrum files/2018 final/SSA/zim_Matabeleland North_inc.PJNZ", "~/Documents/Data/Spectrum files/2018 final/SSA/zim_Matabeleland South_inc.PJNZ", "~/Documents/Data/Spectrum files/2018 final/SSA/zim_Midlands_inc.PJNZ")

zim_spec <- lapply(zim_spec, extract_pjnz_naomi)

names(zim_spec) <- c("Bulawayo", "Harare", "Manicaland", "Mashonaland Central", "Mashonaland East", "Mashonaland West", "Masvingo", "Matabeleland North", "Matabeleland South", "Midlands")

zim_province_tfr <- zim_spec %>%
  lapply(function(x) {x <- x %>%
    filter(sex=="female", age_group_id %in% c(4:10)) %>%
    group_by(year) %>%
    summarise(val = 5*sum(asfr))}) %>%
  map_df(~as.data.frame(.x), .id="area_name") %>%
  mutate(iso3 = "ZWE",
         area_level = 1,
         variable = "tfr",
         source = "Spectrum18") %>%
  rename("period" = "year") %>%
  left_join(areas_long %>% select(iso3, area_id, area_name))

moz_spec <- c("~/Documents/Data/Spectrum files/2018 final/SSA/1_MZ_Niassa_v5_63_updated census.pjnz", "~/Documents/Data/Spectrum files/2018 final/SSA/2_MZ_Cabo Delgado_v5_63_updated census.PJNZ", "~/Documents/Data/Spectrum files/2018 final/SSA/3_MZ_Nampula_v5_63_updated census_22_01_2018.PJNZ", "~/Documents/Data/Spectrum files/2018 final/SSA/4_MZ_Zambezia_v5_63_updated census22_1_2018.PJNZ", "~/Documents/Data/Spectrum files/2018 final/SSA/5_MZ_Tete_v5_63_updated census.PJNZ", "~/Documents/Data/Spectrum files/2018 final/SSA/6_MZ_Manica_v5_63_updated census.PJNZ", "~/Documents/Data/Spectrum files/2018 final/SSA/7_MZ_Sofala_v5_63_updated census_22_01_2018.PJNZ", "~/Documents/Data/Spectrum files/2018 final/SSA/8_MZ_Inhambane_v5_63_updated census22_01_2018.PJNZ", "~/Documents/Data/Spectrum files/2018 final/SSA/9_MZ_Gaza_v5_63_updated census.PJNZ", "~/Documents/Data/Spectrum files/2018 final/SSA/10_MZ_Maputo Provincia_v5_63_updated census.PJNZ", "~/Documents/Data/Spectrum files/2018 final/SSA/11_MZ_Maputo Cidade_v5_63_updated census_22_01_2018.PJNZ")

moz_spec <- lapply(moz_spec, extract_pjnz_naomi)

names(moz_spec) <- c("Niassa", "Cabo Delgado", "Nampula", "Zambezia", "Tete", "Manica", "Sofala", "Inhambane", "Gaza", "Maputo", "Cidade De Maputo")

moz_province_tfr <- moz_spec %>%
  lapply(function(x) {x <- x %>%
    filter(sex=="female", age_group_id %in% c(4:10)) %>%
    group_by(year) %>%
    summarise(val = 5*sum(asfr))}) %>%
  map_df(~as.data.frame(.x), .id="area_name") %>%
  mutate(iso3 = "MOZ",
         area_level = 1,
         variable = "tfr",
         source = "Spectrum18") %>%
  rename("period" = "year") %>%
  left_join(areas_long %>% select(iso3, area_id, area_name))


extract_pjnz_naomi("~/Documents/Data/Spectrum files/2018 final/SSA/zim_Bulawayo_inc.PJNZ")

foo <- spec %>%
  mutate(iso3 = "ZWE") %>% 
  bind_rows(
    
    extract_pjnz_naomi("~/Documents/Data/Spectrum files/2018 final/SSA/Malawi_2018_version_8.PJNZ") %>% mutate(iso3 = "MWI"),
    extract_pjnz_naomi("~/Documents/Data/Spectrum files/2018 final/SSA/Uganda  23May 2018.pjnz") %>% mutate(iso3 = "UGA"),
    extract_pjnz_naomi("~/Documents/Data/Spectrum files/2018 final/SSA/Lesotho 2018_v5_63_1 Feb.PJNZ") %>% mutate(iso3 = "LSO")
    
            )

plot_df$births_U1 <- plot_df$births_U1 %>%
  bind_rows(foo %>%
  filter(year %in% 2000:2020) %>%
  group_by(iso3, year) %>%
  summarise(births = sum(births)) %>%
  mutate(area_id = iso3,
         area_name = countrycode(iso3, "iso3c", "country.name"),
         area_level = 0,
         variable = "births",
         source = "Spectrum18") %>%
  rename(val = births, period = year)
  )

spec <- list()

spec <- extract_pjnz_naomi("~/Documents/Data/Spectrum files/2018 final/SSA/Malawi_2018_version_8.PJNZ")
spec$UGA <- extract_pjnz_naomi("~/Documents/Data/Spectrum files/2018 final/SSA/Uganda  23May 2018.pjnz")
spec$LSO <- extract_pjnz_naomi("~/Documents/Data/Spectrum files/2018 final/SSA/Lesotho 2018_v5_63_1 Feb.PJNZ")
spec$RWA <- extract_pjnz_naomi("~/Documents/Data/Spectrum files/2018 final/SSA/Rwanda _2018_final.PJNZ")

spec <- spec %>%
  select(year, age_group_id, age_group_label, asfr, sex) %>%
  filter(sex=="female", age_group_id>3, age_group_id<11) %>%
  select(-c(sex, age_group_id)) %>%
  rename(agegr = age_group_label, period=year) %>%
  mutate(iso3 = "ZWE",
         area_level = 0,
         area_name = "Zimbabwe",
         area_id = iso3,
         source = "Spectrum18")

spec_zwe_asfr <- spec %>%
  mutate(variable = "asfr") %>%
  rename(val=asfr)
    
spec_zwe_tfr <- spec %>%
  mutate(variable = "tfr") %>%
  group_by(iso3, area_id, area_name, area_level, period, source, variable) %>%
  summarise(tfr = 5*sum(asfr))

spec <- readRDS("spectrum_asfr_tfr.RDS")

spec$asfr <- spec$asfr %>%
  bind_rows(spec_zwe_asfr)

    spec_outputs <- list()
    
    spec$tfr <- spec$tfr %>%
      bind_rows(spec[[4]] %>%
      filter(age_group_id %in% 4:10, sex == "female") %>%
      mutate(iso3 = iso_cod[[4]]) %>%
      group_by(iso3, year) %>%
      summarise(val = 5*sum(asfr)) %>%
      ungroup %>%
      mutate(area_name = countrycode(iso3, "iso3c", "country.name"),
             area_id = iso3,
             area_level = 0,
             source = "Spectrum18",
             variable = "tfr") %>%
      rename(period = year)
      )
    
    spec$asfr <- spec$asfr %>%
      bind_rows(spec[[4]] %>%
      filter(age_group_id %in% 4:10, sex == "female") %>%
      select(year, age_group_label, asfr) %>%
    mutate(iso3 = iso_cod[[4]],
           area_name = countrycode(iso3, "iso3c", "country.name"),
           area_id = iso3,
           area_level = 0,
           source = "Spectrum18",
           variable = "asfr") %>%
      rename(agegr = age_group_label, period = year, val=asfr) %>%
      select(iso3, area_id, area_name, area_level, period, agegr, source, variable, val)
      )

    spec <- spec[c(5,6)]
    
    saveRDS(spec, file="spectrum_asfr_tfr.RDS")
    
plot_df <- list()
    
plot_df$tfr <- lapply(mod_results, "[[", "tfr") %>%
      bind_rows %>%
      rename(val=median) %>%
      bind_rows(wpp_tfr, gbd_tfr)
    
plot_df$asfr <- lapply(mod_results, "[[", "asfr") %>%
      bind_rows %>%
      rename(val=median) %>%
      bind_rows(wpp_asfr, gbd_asfr)
    
plot_df$births_U1 <- lapply(mod_results, "[[", "births") %>%
      bind_rows %>%
      bind_rows(gbd_births, worldpop_U1)
    
plot_df$births_by_age <- lapply(mod_results, "[[", "births_by_age") 
    
plot_df$tfr <- plot_df$tfr %>%
      bind_rows(zim_province_tfr, moz_province_tfr)
    
plot_df$asfr <- plot_df$asfr %>%
      bind_rows(spec$asfr)
    
plot_df$asfr <- plot_df$asfr %>%
      bind_rows(plot_df$asfr %>%
                  filter((source == "Spectrum18" | source == "Model"), area_level ==0 ) %>%
                  group_by(iso3, period, source) %>%
                  mutate(val = val/sum(val),
                  variable = "asfr_prop")
      )
    
# saveRDS(plot_df, file="plot_df.rds")
# plot_df <- readRDS("plot_df.rds")

