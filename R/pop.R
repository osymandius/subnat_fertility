library(tidyverse)

setwd("~/Documents/GitHub/subnat_fertility")

########

quarter_ids <- convert_quarter_id(2, c(2000:2020))

worldpop_U1 <- read_csv("WorldPop_agesex.csv", col_types = cols(X1 = col_skip(), X = col_skip(), sex = col_character())) %>%
  filter(startage==0, iso3 %in% iso3_codes) %>%
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
  filter(iso3 %in% iso3_codes) %>%
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

wpp_tfr <- read_excel("pop_data/WPP2019_FERT_F04_TOTAL_FERTILITY(1).xlsx")

wpp_tfr <- wpp_tfr %>%
  melt(id="area_name", variable.name="period", value.name = "val") %>%
  separate(period, into=c("period", NA), sep="-") %>%
  type.convert() %>%
  mutate(iso3 = countrycode(area_name, "country.name", "iso3c")) %>%
  filter(iso3 %in% c("UGA", "MWI")) %>%
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
  filter(iso3 %in% iso3_codes) %>%
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

zim_spec <- list("~/Documents/Data/Spectrum files/2018 final/SSA/zim_Bulawayo_inc.PJNZ", "~/Documents/Data/Spectrum files/2018 final/SSA/zim_Harare_inc.PJNZ", "~/Documents/Data/Spectrum files/2018 final/SSA/zim_Manicaland_inc.PJNZ", "~/Documents/Data/Spectrum files/2018 final/SSA/zim_Mashonaland Central_inc.PJNZ", "~/Documents/Data/Spectrum files/2018 final/SSA/zim_Mashonaland East_inc.PJNZ", "~/Documents/Data/Spectrum files/2018 final/SSA/zim_Mashonaland West_inc.PJNZ", "~/Documents/Data/Spectrum files/2018 final/SSA/zim_Masvingo_inc.PJNZ", "~/Documents/Data/Spectrum files/2018 final/SSA/zim_Matabeleland North_inc.PJNZ", "~/Documents/Data/Spectrum files/2018 final/SSA/zim_Matabeleland South_inc.PJNZ", "~/Documents/Data/Spectrum files/2018 final/SSA/zim_Midlands_inc.PJNZ")

zim_spec <- lapply(zim_spec, extract_pjnz_naomi)

get_age_groups()

names(zim_spec) <- c("Bulawayo", "Harare", "Manicaland", "Mashonaland Central", "Mashonaland East", "Mashonaland West", "Masvingo", "Matabeleland North", "Matabeleland South", "Midlands")

zim_province_tfr <- zim_spec %>%
  lapply(function(x) {x <- x %>%
    mutate(agegr = cut(age, c(0:16*5, Inf),c(paste0(0:15*5, "-", 0:15*5+4), "80+"), include.lowest = TRUE, right= FALSE)) %>%
    filter(sex=="female", age %in% c(15:49)) %>%
    group_by(year) %>%
    summarise(val = sum(asfr))}) %>%
  map_df(~as.data.frame(.x), .id="area_name") %>%
  mutate(iso3 = "ZWE",
         area_level = 1,
         variable = "tfr",
         source = "Spectrum18") %>%
  rename("period" = "year") %>%
  left_join(areas_long %>% filter(area_level == 1) %>% select(iso3, area_id, area_name))

moz_spec_paths <- c("~/Documents/Data/Spectrum files/2018 final/SSA/1_MZ_Niassa_v5_63_updated census.pjnz", "~/Documents/Data/Spectrum files/2018 final/SSA/2_MZ_Cabo Delgado_v5_63_updated census.PJNZ", "~/Documents/Data/Spectrum files/2018 final/SSA/3_MZ_Nampula_v5_63_updated census_22_01_2018.PJNZ", "~/Documents/Data/Spectrum files/2018 final/SSA/4_MZ_Zambezia_v5_63_updated census22_1_2018.PJNZ", "~/Documents/Data/Spectrum files/2018 final/SSA/5_MZ_Tete_v5_63_updated census.PJNZ", "~/Documents/Data/Spectrum files/2018 final/SSA/6_MZ_Manica_v5_63_updated census.PJNZ", "~/Documents/Data/Spectrum files/2018 final/SSA/7_MZ_Sofala_v5_63_updated census_22_01_2018.PJNZ", "~/Documents/Data/Spectrum files/2018 final/SSA/8_MZ_Inhambane_v5_63_updated census22_01_2018.PJNZ", "~/Documents/Data/Spectrum files/2018 final/SSA/9_MZ_Gaza_v5_63_updated census.PJNZ", "~/Documents/Data/Spectrum files/2018 final/SSA/10_MZ_Maputo Provincia_v5_63_updated census.PJNZ", "~/Documents/Data/Spectrum files/2018 final/SSA/11_MZ_Maputo Cidade_v5_63_updated census_22_01_2018.PJNZ")

moz_spec <- lapply(moz_spec_paths, extract_pjnz_naomi)

names(moz_spec) <- c("Niassa", "Cabo Delgado", "Nampula", "Zambézia", "Tete", "Manica", "Sofala", "Inhambane", "Gaza", "Maputo Província", "Cidade de Maputo")

moz_province_tfr <- moz_spec %>%
  lapply(function(x) {x <- x %>%
    mutate(agegr = cut(age, c(0:16*5, Inf),c(paste0(0:15*5, "-", 0:15*5+4), "80+"), include.lowest = TRUE, right= FALSE)) %>%
    filter(sex=="female", age %in% c(15:49)) %>%
    group_by(year) %>%
    summarise(val = sum(asfr))}) %>%
  map_df(~as.data.frame(.x), .id="area_name") %>%
  mutate(iso3 = "MOZ",
         area_level = 1,
         variable = "tfr",
         source = "Spectrum18") %>%
  rename("period" = "year") %>%
  left_join(areas_long %>% filter(area_level == 1) %>% select(iso3, area_id, area_name))

zmb_spec_paths <- c("~/Documents/Data/Spectrum files/2018 final/SSA/Zambia 2018_1_CENTRAL 15May2018.PJNZ", "~/Documents/Data/Spectrum files/2018 final/SSA/Zambia 2018_1_COPPERBELT 15May2018.PJNZ", "~/Documents/Data/Spectrum files/2018 final/SSA/Zambia 2018_1_EASTERN 15May2018.PJNZ", "~/Documents/Data/Spectrum files/2018 final/SSA/Zambia 2018_1_LUAPULA 15May2018.PJNZ", "~/Documents/Data/Spectrum files/2018 final/SSA/Zambia 2018_1_LUSAKA 15May2018.PJNZ", "~/Documents/Data/Spectrum files/2018 final/SSA/Zambia 2018_1_MUCHINGA may 2018.PJNZ", "~/Documents/Data/Spectrum files/2018 final/SSA/Zambia 2018_1_NORTHERN 15May2018.PJNZ", "~/Documents/Data/Spectrum files/2018 final/SSA/Zambia 2018_1_SOUTHERN 15 MAY 2018.PJNZ", "~/Documents/Data/Spectrum files/2018 final/SSA/Zambia 2018_1_WESTERN may 2018.PJNZ", "~/Documents/Data/Spectrum files/2018 final/SSA/Zambia 2018_1-NORTH-WESTERN 15May2018.PJNZ")

names(zmb_spec_paths) <- c("Central", "Copperbelt", "Eastern", "Luapula", "Lusaka", "Muchinga", "Northern", "Southern", "Western", "North-Western")

zmb_spec <- lapply(zmb_spec_paths, extract_pjnz_naomi)

zmb_province_tfr <- zmb_spec %>%
  lapply(function(x) {x <- x %>%
    mutate(agegr = cut(age, c(0:16*5, Inf),c(paste0(0:15*5, "-", 0:15*5+4), "80+"), include.lowest = TRUE, right= FALSE)) %>%
    filter(sex=="female", age %in% c(15:49)) %>%
    group_by(year) %>%
    summarise(val = sum(asfr))}) %>%
  map_df(~as.data.frame(.x), .id="area_name") %>%
  mutate(iso3 = "ZMB",
         area_level = 1,
         variable = "tfr",
         source = "Spectrum18") %>%
  rename("period" = "year") %>%
  left_join(areas_long %>% filter(area_level == 1) %>% select(iso3, area_id, area_name))

national_tfr <- bind_rows(
    extract_pjnz_naomi("~/Documents/Data/Spectrum files/2018 final/SSA/Malawi_2018_version_8.PJNZ") %>% mutate(iso3 = "MWI"),
    extract_pjnz_naomi("~/Documents/Data/Spectrum files/2018 final/SSA/Uganda  23May 2018.pjnz") %>% mutate(iso3 = "UGA"),
    extract_pjnz_naomi("~/Documents/Data/Spectrum files/2018 final/SSA/Lesotho 2018_v5_63_1 Feb.PJNZ") %>% mutate(iso3 = "LSO"),
    extract_pjnz_naomi("~/Documents/Data/Spectrum files/2018 final/SSA/TZ_Regions_May 2018.PJNZ") %>% mutate(iso3 = "TZA"),
    extract_pjnz_naomi("~/Documents/Data/Spectrum files/2018 final/SSA/Namibia_2018 national.PJNZ") %>% mutate(iso3 = "NAM"),
    extract_pjnz_naomi("~/Documents/Data/Spectrum files/2018 final/SSA/Swaziland File_2018_14May2018.PJNZ") %>% mutate(iso3 = "SWZ")
  ) %>%
  mutate(agegr = cut(age, c(0:16*5, Inf),c(paste0(0:15*5, "-", 0:15*5+4), "80+"), include.lowest = TRUE, right= FALSE)) %>%
  filter(sex=="female", age %in% c(15:49)) %>%
  group_by(year, iso3) %>%
  summarise(val = sum(asfr)) %>%
  mutate(area_level = 0,
         variable = "tfr",
         source = "Spectrum18",
         area_name = ifelse(iso3 ==  "SWZ", "Eswatini", countrycode(iso3, "iso3c", "country.name")),
         area_id = iso3) %>%
  rename("period" = "year") %>%
  filter(period < 2021)

spec_tfr <- national_tfr %>%
  bind_rows(zmb_province_tfr, zim_province_tfr, moz_province_tfr)

saveRDS(spec_tfr, file="spec_tfr_2019_12_3.rds")

plot_df <- list()
    
plot_df$tfr <- lapply(mod_results, "[[", "tfr") %>%
  bind_rows %>%
  rename(val=median) %>%
  select(-c(mean, sd, naomi_level)) %>%
  bind_rows(wpp_tfr, gbd_tfr, spec_tfr) %>%
  filter(period %in% 2000:2020) %>%
  mutate(area_name = ifelse(iso3 ==  "SWZ", "Eswatini", countrycode(iso3, "iso3c", "country.name"))) %>%
  left_join(areas_long %>% select(area_id, naomi_level))
    
plot_df$asfr <- lapply(mod_results, "[[", "asfr") %>%
  bind_rows %>%
  rename(val=median) %>%
  select(-c(mean, sd, naomi_level)) %>%
  bind_rows(wpp_asfr, gbd_asfr) %>%
  filter(period %in% 2000:2020) %>%
  mutate(area_name = ifelse(iso3 ==  "SWZ", "Eswatini", countrycode(iso3, "iso3c", "country.name"))) %>%
  left_join(areas_long %>% select(area_id, naomi_level))
    
# plot_df$births_U1 <- lapply(mod_results, "[[", "births") %>%
#       bind_rows %>%
#       bind_rows(gbd_births, worldpop_U1)
#     
# plot_df$births_by_age <- lapply(mod_results, "[[", "births_by_age") 

    
saveRDS(plot_df, file="plot_df_2019_12_4.rds")

plot_df$tfr %>%
  filter(area_level == 0) %>%
  ggplot(aes(x=period, y=val, color=source, group=source)) +
    geom_line() +
    facet_wrap(~area_id) +
    labs(title="TFR comparison | Admin0", y="TFR", x="Year")

plot_df$tfr %>%
  filter(area_level == 1, iso3 %in% c("MOZ")) %>%
  bind_rows(moz_res$tfr %>% filter(area_level == 1) %>% mutate(source = "Model no TIPS") %>% rename(val = median)) %>%
  ggplot(aes(x=period, y=val)) +
  geom_line(aes(color=source, group=source)) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=source), alpha=0.2) +
  facet_wrap(~area_id) +
  labs(title="TFR comparison | Admin1", y="TFR", x="Year")



