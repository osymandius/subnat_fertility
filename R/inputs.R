make_areas_population <- function(iso3_current, naomi_data_path, full=FALSE, wide=TRUE, boundaries=TRUE, population=TRUE, return_list = TRUE) {

paths <- file.path(naomi_data_path, iso3_current, "data")

files <- lapply(paths, function(paths) {
  
  files <- list.files(paths, full.names = TRUE)
  area <- files %>% str_subset(pattern = "areas.geojson") %>% str_subset(pattern = ".zip", negate=TRUE)
  
  if(population) {
    pop <- files %>% str_subset(pattern = "population")
    files <- c(area, pop)
    names(files) <- c("areas", "population")
  } else {
  files <- c(area)
  names(files) <- c("areas")
  }
  
  return(files)

})

names(files) <- iso3_current

areas_long <- lapply(files, "[[", "areas") %>%
  lapply(read_sf) %>% 
  lapply(function(x) {
    
    iso3_code <- x %>%
      filter(area_level == 0) %>%
      select(area_id) %>%
      unique %>%
      .$area_id
    
    x <- x %>%
      mutate(iso3 = iso3_code) %>%
      st_drop_geometry() %>%
      select(c("iso3", "area_id", "area_name", "area_level", "parent_area_id", "naomi_level"))
    
    return(x)
  }) %>% 
  bind_rows %>%
  arrange(iso3)

if(full)
  areas_full <- lapply(files, "[[", "areas") %>%
    lapply(st_read) %>% 
    bind_rows %>%
    arrange(iso3)

if(wide)
  areas_wide <- lapply(files, "[[", "areas") %>%
    lapply(read_sf) %>%
    lapply(function(x) {spread_areas(as.data.frame(x))}) %>%
    lapply(function(x) {x %>% mutate(iso3 = area_id0)}) %>%
    bind_rows %>%
    arrange(iso3)

if(boundaries)
  area_boundaries <- lapply(files, "[[", "areas") %>%
    lapply(read_sf) %>% 
    lapply(function(x) {
      
      iso3_code <- x %>%
        filter(area_level == 0) %>%
        select(area_id) %>%
        unique %>%
        .$area_id
      
      x <- x %>%
        mutate(iso3 = iso3_code) %>%
        select(-epp_level)
      
      return(x)
    }) %>% 
    bind_rows %>%
    arrange(iso3)

if(population)
  area_population <- lapply(files, "[[", "population") %>%
    lapply(read_csv) %>%
    lapply(left_join, areas_long) %>%
    bind_rows %>%
    mutate(period = year_labels(naomi:::calendar_quarter_to_quarter_id(calendar_quarter))) %>%
    select(iso3, "area_id" , "area_name", "source", "sex", "age_group", "population", "period") %>%
    arrange(iso3)



  df <- list()
  
  df$areas_long <- areas_long
  
  # if(length(iso3_current >1 & return_list == TRUE))
  #   df$areas_long <- df$areas_long %>% group_split(iso3)
  
  if(full) {
    df$areas_full <- areas_full
  
    # if(length(iso3_current >1 & return_list == TRUE))
    #   df$areas_full <- df$areas_full %>% group_split(iso3)
  }
  
  if(wide) {
    df$areas_wide <- areas_wide
  
    # if(length(iso3_current >1 & return_list == TRUE))
    #   df$areas_wide <- df$areas_wide %>% group_split(iso3)
  }
  
  if(boundaries) {
    df$boundaries <- area_boundaries
  
    # if(length(iso3_current >1 & return_list == TRUE))
    #   df$boundaries <- df$boundaries %>% group_split(iso3)
  }
  
  if(population) {
    df$population <- area_population
  
    # if(length(iso3_current >1 & return_list == TRUE))
      # df$population <- df$population %>% group_split(iso3)
  }
  
  return(df)


}

assign_cluster_area <- function(clusters, area_level) {

  iso3_current <- unique(clusters$iso3)

  if (area_level == "naomi") {

    int <- areas_long %>%
              filter(iso3 == iso3_current)

    area_level <- unique(int$area_level[int$naomi_level == TRUE])

  }
  
  areas <- clusters %>%
    rename(area_id = geoloc_area_id) %>%
    left_join(areas_wide %>% select(area_id, paste0("area_id", area_level))) %>%
    select(-area_id) %>%
    rename(area_id = paste0("area_id", area_level))
  
  area_list <- areas %>%
    group_by(survey_id) %>%
    group_split(keep=TRUE)
  
  names(area_list) <- area_list %>%
    lapply("[", "survey_id") %>% 
    lapply(unique) %>% 
    bind_rows %>% 
    .$survey_id
  
  return(area_list)
  
}

clusters_to_surveys <- function(surveys, cluster_areas, single_tips = TRUE) {

    level <- areas_long$area_level[areas_long$area_id == cluster_areas[[1]]$area_id[[1]]]


    ird <- dhs_datasets(fileType = "IR", fileFormat = "flat", surveyIds = surveys$SurveyId)

    ird <- ird %>%
      mutate(path = unlist(get_datasets(.))) %>%
      bind_rows()
      
    ir <- lapply(ird$path, readRDS) %>%
      lapply(function(x) {class(x) <- "data.frame"
      return(x)}) %>%
      Map(function(ir, surveys) {
        mutate(ir,
           surveyid = surveys$SurveyId,
           country = surveys$CountryName,
           survyear = surveys$SurveyYear,
           survtype = surveys$SurveyType)
      }, ., group_split(surveys, SurveyId))
    
    
    ## I think this mess was necessary because the order of cluster_areas was not always the same as the order for ir.. Double check. Otherwise the simple names(ir) <- names(cluster_areas) is sufficient
    # names(ir) <- ir %>% 
    #   lapply("[", "surveyid") %>% 
    #   lapply(unique) %>% 
    #   bind_rows %>% 
    #   left_join(clusters %>% 
    #               select(survey_id, DHS_survey_id) %>% 
    #               unique, 
    #             by=c("surveyid" = "DHS_survey_id")) %>% 
    #   select(-surveyid) %>% 
    #   .$survey_id
    
    names(ir) <- names(cluster_areas)

    if(unique(surveys$iso3) == "ETH") {
    
    cols_edit <- c("v008", "v011", "b3_01", "b3_02", "b3_03", "b3_04", "b3_05", "b3_06", "b3_07", "b3_08", "b3_09", "b3_10", "b3_11", "b3_12", "b3_13", "b3_14", "b3_15", "b3_16", "b3_17", "b3_18", "b3_19", "b3_20")
    
    ir <- ir %>%
      lapply(function(x) x %>% mutate_at(.vars = cols_edit, .funs = ~(.+92)))

    }
    
    
  if(level > 0) {
      
    ir <- Map(ir_by_area2, ir, cluster_areas[names(ir)], n=1:length(ir), total=length(ir)) %>%
      unlist(recursive = FALSE)
    
    survey_type <- ir %>%
      lapply("[", "survtype") %>%
      lapply(unique) %>%
      bind_rows
    
    if(!single_tips)
      tips_surv <- list("DHS" = c(0,15), "MIS" = c(0,5), "AIS" = c(0,5))[survey_type$survtype]
    else
      tips_surv <- list("DHS" = c(0:15), "MIS" = c(0:5), "AIS" = c(0:5))[survey_type$survtype]
    
  } else {

    ir <- lapply(ir, function(x) {
        mutate(x, area_id = unique(surveys$iso3))
      })
    
    if(!single_tips)
      tips_surv <- list("DHS" = c(0,15), "MIS" = c(0,5), "AIS" = c(0,5))[surveys$SurveyType]
    else  
      tips_surv <- list("DHS" = c(0:15), "MIS" = c(0:5), "AIS" = c(0:5))[surveys$SurveyType]
    
  }
    
  dat <- list()
  dat$ir <- ir
  dat$tips_surv <- tips_surv
  
  return(dat)

}

read_mics <- function(iso3_current, path_to_MICS = "~/Imperial College London/HIV Inference Group - Documents/Data/household surveys/MICS/datasets") {
  
  temp <- tempdir()
  path <- grep(countrycode(iso3_current, "iso3c", "country.name"), list.files(path_to_MICS, full.names=TRUE), value=TRUE)
  uz <- lapply(path, unzip, exdir = temp)
  
  check <- read.csv(here::here("input_data/MICS_list.csv")) %>%
    filter(status == "Completed") %>%
    mutate(iso3 = countrycode(country, "country.name", "iso3c")) %>%
    filter(iso3 == iso3_current)
  
  if(nrow(check) != length(path) ) warning(paste("Database has", nrow(check), "datasets for", toString(sort(type.convert(check$year))), "you have extracted", length(path)))
  
  mics_indicators <- read_csv(here::here("input_data/MICS_indicators.csv")) %>%
    pivot_longer(-c(label, id, filetype))
  
  df <- lapply(uz, function(x) {
    
    ##### HH
    
    hh_path <- grep("hh.sav", x, value = TRUE)
    hh_extract <- read_sav(hh_path)
    colnames(hh_extract) <- tolower(colnames(hh_extract))
    
    surv_year_id <- hh_extract %>%
      select(starts_with("hh")) %>%
      select(ends_with("y")) %>%
      colnames
    
    surv_year <- hh_extract %>%
      select(surv_year_id) %>%
      count(year = get(surv_year_id)) %>%
      slice(which.max(n)) %>%
      select(year) %>%
      as.integer()
    
    surv_round <- filter(check, year == surv_year)$round
    
    vars <- mics_indicators %>%
      filter(name == surv_round, filetype=="hh") %>%
      .$value
    
    names(vars) <- mics_indicators %>%
      filter(name == surv_round, filetype=="hh") %>%
      .$id
    
    hh_extract <- hh_extract %>%
      rename(province = vars[["province"]], cluster = vars[["cluster"]], hh_number = vars[["hh_number"]])
    
    prov_labels <- data.frame(province = attr(hh_extract$province, "labels"), area_name = str_to_title(names(attr(hh_extract$province, "labels")))) %>%
      left_join(areas_long %>% filter(area_level == 1))
    
    if(nrow(filter(prov_labels, is.na(area_id))) != 0)
      stop(paste(nrow(filter(prov_labels, is.na(area_id))), "province labels have not been matched with area ids"))
    
    hh_extract <- hh_extract %>%
      left_join(prov_labels %>% select(province, area_id))
    
    ##### WM
    
    wm_path <- grep("wm.sav", x, value = TRUE)
    wm_extract <- read_sav(wm_path)
    colnames(wm_extract) <- tolower(colnames(wm_extract))
    
    vars <- mics_indicators %>%
      filter(name == surv_round, filetype=="wm") %>%
      .$value
    
    names(vars) <- mics_indicators %>%
      filter(name == surv_round, filetype=="wm") %>%
      .$id
    
    wm_extract <- wm_extract %>%
      mutate(survyear = surv_year,
             surveyid = paste0(iso3_current, "MICS", survyear),
             survtype = "MICS") %>%
      rename(wdob = vars[["wdob"]], cluster = vars[["cluster"]], hh_number = vars[["hh_number"]], line_number = vars[["line_number"]], doi = vars[["doi"]]) %>%
      filter(!is.na(wdob), !is.na(cluster), !is.na(hh_number), !is.na(line_number), !is.na(doi)) %>%
      arrange(cluster, hh_number, line_number) %>%
      mutate(unique_id = group_indices(., cluster, hh_number, line_number))
    
    
    ###### BH
    
    bh_path <- grep("bh.sav", x, value=TRUE)
    bh_extract <- read_sav(bh_path)
    colnames(bh_extract) <- tolower(colnames(bh_extract))
    
    vars <- mics_indicators %>%
      filter(name == surv_round, filetype == "bh") %>%
      .$value
    
    names(vars) <- mics_indicators %>%
      filter(name == surv_round, filetype == "bh") %>%
      .$id
    
    bh_extract <- bh_extract %>%
      rename(cluster = vars[["cluster"]], hh_number = vars[["hh_number"]], line_number = vars[["line_number"]], cdob = vars[["cdob"]])
    
    df <- list()
    df$hh <- hh_extract
    df$wm <- wm_extract
    df$bh <- bh_extract
    
    
    return(df)
    
  })
  
  bh <- lapply(df, "[[", "bh")
  wm <- lapply(df, "[[", "wm")
  hh <- lapply(df, "[[", "hh")
  
  if(length(bh) != length(uz)) stop("Number of birth history files does not match number of surveys extracted")
  if(length(wm) != length(uz)) stop("Number of women files does not match number of surveys extracted")
  if(length(hh) != length(uz)) stop("Number of household files does not match number of surveys extracted")
  
  mics_dat <- list()
  
  wm <- Map(function(wm, hh) {
    
    wm %>% 
      left_join(hh %>% select(cluster, hh_number, area_id))
    
  }, wm, hh)
  
  bh_df <- Map(function(bh, wm) {
    wm %>%
      select(cluster, hh_number, line_number, unique_id) %>%
      left_join(bh %>% select(cluster, hh_number, line_number, cdob)) %>%
      select(unique_id, cdob) %>%
      filter(!is.na(cdob))
  }, bh, wm)
  
  mics_dat$bh_df <- bh_df
  mics_dat$wm <- wm
  
  return(mics_dat)
  
}

cmc_to_year <- function(cmc) {

  year <- 1900 + (cmc-1)/12

  return(year)

}


make_asfr_pred_df <- function(asfr, t2 = 2020) {
  
  area_df <- areas_long %>% 
    filter(iso3 == unique(asfr$iso3), area_level == areas_long$area_level[areas_long$area_id == asfr$area_id[1]])
  
  pred_df <- crossing(area_id = area_df$area_id, 
                      period = min(asfr$period):t2, 
                      age_group = asfr$age_group,
                      pys=1,
                      iso3 = unique(asfr$iso3))
  
  asfr_pred <- pred_df %>%
    bind_rows(asfr) %>%
    mutate(id.period = group_indices(., period),
           id.period2 = id.period,
           id.period3 = id.period,
           id.age_group = group_indices(., age_group),
           id.age_group2 = id.age_group,
           id.age_group3 = id.age_group,
           id.age_group.period = group_indices(., period, age_group),
           id.district = group_indices(., area_id),
           id.district2 = id.district,
           id.district3 = id.district,
           id.tips = (group_indices(., tips)),
           id.tips = ifelse(is.na(tips), NA, id.tips),
           id.tips = factor(id.tips),
           tips_dummy = ifelse(tips>5, 1, 0),
           id = 1:nrow(.)
    ) %>%
    mutate_if(is.factor, as.character)
  
  return(asfr_pred)
  
}

get_neighbourhood_structure <- function(asfr, areas_long, boundaries) {
  
  level <- areas_long$area_level[areas_long$area_id == asfr$area_id[1]]
  
  sh <- areas_long %>%
    filter(iso3 == unique(asfr$iso3), area_level == level) %>%
    mutate(area_idx = row_number())
  
  nb <- sh %>%
    left_join(boundaries) %>%
    st_as_sf %>%
    as("Spatial") %>%
    spdep::poly2nb() %>%
    `names<-`(sh$area_idx)
  
  nb2INLA(paste0("countries/", unique(asfr$iso3), "/adj/", unique(asfr$iso3), "_admin", level, ".adj"), nb)
  
}

get_asfr_pred_df <- function(iso3_current, area_level, areas_long, project) {

  if (area_level == "naomi") {

    int <- areas_long %>%
              filter(iso3 == iso3_current)

    area_level <- unique(int$area_level[int$naomi_level == TRUE])

  }


  if (project == FALSE) {

    dat <- readRDS(here("countries", paste0(iso3_current, "/data/", iso3_current, "_asfr_admin", area_level, ".rds")))

    year <- unique(max(filter(dat, !is.na(surveyid))$period))

    dat <- dat %>%
      filter(period <= year)

  } else {

    dat <- readRDS(here("countries", paste0(iso3_current, "/data/", iso3_current, "_asfr_admin", area_level, ".rds")))

  }
  
  return(dat)


}
