make_areas_population <- function(iso3_codes, path_to_naomi_data) {

paths <- paste0(path_to_naomi_data, iso3_codes, "/data")

files <- lapply(paths, function(paths) {
  
  files <- list.files(paths, full.names = TRUE)
  area <- files %>% str_subset(pattern = "areas.geojson")
  pop <- files %>% str_subset(pattern = "population")
  
  files <- c(area, pop)
  names(files) <- c("areas", "population")
  return(files)

})

names(files) <- iso3_codes

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
  bind_rows

areas_wide <- lapply(files, "[[", "areas") %>%
  lapply(read_sf) %>%
  lapply(function(x) {spread_areas(as.data.frame(x))}) %>%
  lapply(function(x) {x %>% mutate(iso3 = area_id0)}) %>%
  bind_rows

boundaries <- lapply(files, "[[", "areas") %>%
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
  bind_rows

population <- lapply(files, "[[", "population") %>%
  lapply(read_csv) %>%
  bind_rows %>%
  mutate(period = year_labels(calendar_quarter_to_quarter_id(calendar_quarter))) %>%
  select("area_id" , "area_name", "source", "sex", "age_group", "population", "period")



  df <- list()
  df$areas_long <- areas_long
  df$areas_wide <- areas_wide
  df$boundaries <- boundaries
  df$population <- population

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
      bind_rows

    ir <- lapply(ird$path, readRDS) %>%
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
    
    cols <- c("b3_01", "b3_02", "b3_03", "b3_04", "b3_05", "b3_06", "b3_07", "b3_08", "b3_09", "b3_10", "b3_11", "b3_12", "b3_13", "b3_14", "b3_15", "b3_16", "b3_17", "b3_18", "b3_19", "b3_20")
    
    ir[[1]]$v008 %<>% add(92)
    ir[[1]]$v011 %<>% add(92)
    ir[[1]][, cols] %<>% add(92)
    
    ir[[2]]$v008 %<>% add(92)
    ir[[2]]$v011 %<>% add(92)
    ir[[2]][, cols] %<>% add(92)
    
    ir[[3]]$v008 %<>% add(92)
    ir[[3]]$v011 %<>% add(92)
    ir[[3]][, cols] %<>% add(92)
    
    ir[[4]]$v008 %<>% add(92)
    ir[[4]]$v011 %<>% add(92)
    ir[[4]][, cols] %<>% add(92)
    
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

get_asfr_pred_df <- function(iso3_current, area_level, project) {

  if (area_level == "naomi") {

    int <- areas_long %>%
              filter(iso3 == iso3_current)

    area_level <- unique(int$area_level[int$naomi_level == TRUE])

  }


  if (project == FALSE) {

    dat <- readRDS(paste0("countries/", iso3_current, "/data/", iso3_current, "_asfr_admin", area_level, ".rds"))

    year <- unique(max(filter(dat, !is.na(surveyid))$period))

    dat <- dat %>%
      filter(period <= year)

  } else {

    dat <- readRDS(paste0("countries/", iso3_current, "/data/", iso3_current, "_asfr_admin", area_level, ".rds"))


  }
  
  return(dat)


}
