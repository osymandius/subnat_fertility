iso3_to_dhs <- function(iso3) {
  
  if (is.list(iso3)) {
    
    lapply(iso3, function(iso3) {
      
      dhs_countries(returnFields=c("CountryName", "DHS_CountryCode")) %>%
        filter(CountryName %in% sub("Swaziland", "Eswatini", countrycode(iso3, "iso3c", "country.name"))) %>%
        .$DHS_CountryCode
    })
    
  } else {
    
    dhs_countries(returnFields=c("CountryName", "DHS_CountryCode")) %>%
      filter(CountryName %in% sub("Swaziland", "Eswatini", countrycode(iso3, "iso3c", "country.name"))) %>% 
      .$DHS_CountryCode
    
  }
  
}

calc_asfr1 <- function(data,
                       y=NULL,
                       by = NULL,
                       agegr = NULL,
                       period = NULL,
                       cohort = NULL,
                       tips = NULL,
                       clusters=~v021,
                       strata=NULL,
                       id="caseid",
                       dob="v011",
                       intv = "v008",
                       weight= "v005",
                       varmethod = "none",
                       bvars = grep("^b3\\_[0-9]*", names(data), value=TRUE),
                       birth_displace = 1e-6,
                       origin=1900,
                       scale=12,
                       bhdata = NULL,
                       counts=FALSE,
                       clustcounts = FALSE){
  
  print(y)
  
  data$id <- data[[id]]
  data$dob <- data[[dob]]
  data$intv <- data[[intv]]
  data$weights <- data[[weight]] / mean(data[[weight]])
  
  if(is.null(by))
    by <- ~1
  
  vars <- unique(unlist(lapply(c(by, strata, clusters), all.vars)))
  f <- formula(paste("~", paste(vars, collapse = "+")))
  mf <- model.frame(formula = f, data = data, na.action = na.pass,
                    id = id, weights = weights, dob = dob, intv = intv)
  
  if(is.null(bhdata)) {
    births <- reshape(model.frame(paste("~", paste(bvars, collapse="+")),
                                  data, na.action=na.pass, id=id),
                      idvar="(id)", timevar="bidx",
                      varying=bvars, v.names="bcmc", direction="long")
  } else {
    if(length(bvars) > 1)
      stop("If `bhdata' is provided, bvars must provide a single variable name (length(bvars) = 1)")
    
    bhdata$id <- bhdata[[id]]
    bhdata$bcmc <- bhdata[[bvars]]
    births <- model.frame(~bcmc, data = bhdata, id = id)
    births$bidx <- ave(births$bcmc, births$`(id)`, FUN = seq_along)
  }
  births <- births[!is.na(births$bcmc), ]
  births$bcmc <- births$bcmc + births$bidx * birth_displace
  
  epis <- tmerge(mf, mf, id=`(id)`, tstart=`(dob)`, tstop=`(intv)`)
  epis <- tmerge(epis, births, id=`(id)`, birth = event(bcmc))
  
  aggr <- demog_pyears(f, epis, period=period, agegr=agegr, cohort=cohort, tips=tips,
                       event="birth", weights="(weights)", origin=origin, scale=scale)$data
  
  ## construct interaction of all factor levels that appear
  byvar <- intersect(c(all.vars(by), "agegr", "period", "cohort", "tips"),
                     names(aggr))
  aggr$byf <- interaction(aggr[byvar], drop=TRUE)
  
  ## prediction for all factor levels that appear
  pred <- data.frame(aggr[c(byvar, "byf")])[!duplicated(aggr$byf),]
  pred <- pred[order(pred$byf), ]
  
  if(counts || varmethod == "none"){
    mc <- model.matrix(~-1+byf, aggr)
    clong <- aggr[c("event", "pyears")]
    pred[c("births", "pys")] <- t(mc) %*% as.matrix(clong)
  }
  
  if(varmethod == "none") {
    
    pred$asfr <- pred$births / pred$pys
    pred$byf <- NULL
    if(!counts)
      pred[c("births", "pys")] <- NULL
    
  } else if(varmethod == "lin") {
    
    des <- survey::svydesign(ids=clusters, strata=strata, data=aggr, weights=~1)
    class(des) <- c("svypyears", class(des))
    
    ## fit model
    f <- if(length(levels(aggr$byf)) == 1)
      event ~ offset(log(pyears))
    else
      event ~ -1 + byf + offset(log(pyears))
    
    mod <- survey::svyglm(f, des, family=quasipoisson)
    
    ## prediction for all factor levels that appear
    pred$pyears <- 1
    
    asfr <- predict(mod, pred, type="response", vcov=TRUE)
    v <- vcov(asfr)
    dimnames(v) <- list(pred$byf, pred$byf)
    
    pred$asfr <- as.numeric(asfr)
    pred$se_asfr <- sqrt(diag(v))
    pred[c("byf", "pyears")] <- NULL
    attr(pred, "var") <- v
  } else if(varmethod %in% c("jkn", "jk1")) {
    
    ## Convert to array with events and PYs for each cluster
    ## reshape2::acast is MUCH faster than stats::reshape
    events_clust <- reshape2::acast(aggr, update(clusters, byf ~ .), value.var="event")
    pyears_clust <- reshape2::acast(aggr, update(clusters, byf ~ .), value.var="pyears")
    
    if(varmethod == "jkn"){
      aggr$strataid <- as.integer(interaction(aggr[all.vars(strata)], drop=TRUE))
      strataid <- drop(reshape2::acast(unique(aggr[c(all.vars(clusters), "strataid")]),
                                       update(clusters,  1 ~ .), value.var="strataid"))
    } else
      strataid <- NULL
    
    estdf <- jackknife(events_clust, pyears_clust, strataid)
    
    pred$asfr <- estdf$est
    pred$se_asfr <- estdf$se
    attr(pred, "var") <- vcov(estdf)
    pred$byf <- NULL
    if(clustcounts){
      attr(pred, "events_clust") <- events_clust
      attr(pred, "pyears_clust") <- pyears_clust
      attr(pred, "strataid") <- strataid
    }
  } else
    stop(paste0("varmethod = \"", varmethod, "\" is not recognized."))
  
  
  rownames(pred) <- NULL
  
  return(pred)
}

calc_asfr_mics <- function(data,
                       y=NULL,
                       by = NULL,
                       agegr = NULL,
                       period = NULL,
                       cohort = NULL,
                       tips = NULL,
                       clusters=~cluster,
                       strata=NULL,
                       id="unique_id",
                       dob="wdob",
                       intv = "doi",
                       weight= "weight",
                       varmethod = "none",
                       bvars = "cdob",
                       birth_displace = 1e-6,
                       origin=1900,
                       scale=12,
                       bhdata = bh_df,
                       counts=FALSE,
                       clustcounts = FALSE){
  
  print(y)
  
  data$id <- data[[id]]
  data$dob <- data[[dob]]
  data$intv <- data[[intv]]
  data$weights <- data[[weight]] / mean(data[[weight]])
  
  if(is.null(by))
    by <- ~1
  
  vars <- unique(unlist(lapply(c(by, strata, clusters), all.vars)))
  f <- formula(paste("~", paste(vars, collapse = "+")))
  mf <- model.frame(formula = f, data = data, na.action = na.pass,
                    id = id, weights = weights, dob = dob, intv = intv)
  
  if(is.null(bhdata)) {
    births <- reshape(model.frame(paste("~", paste(bvars, collapse="+")),
                                  data, na.action=na.pass, id=id),
                      idvar="(id)", timevar="bidx",
                      varying=bvars, v.names="bcmc", direction="long")
  } else {
    if(length(bvars) > 1)
      stop("If `bhdata' is provided, bvars must provide a single variable name (length(bvars) = 1)")
    
    bhdata$id <- bhdata[[id]]
    bhdata$bcmc <- bhdata[[bvars]]
    births <- model.frame(~bcmc, data = bhdata, id = id)
    births$bidx <- ave(births$bcmc, births$`(id)`, FUN = seq_along)
  }
  births <- births[!is.na(births$bcmc), ]
  births$bcmc <- births$bcmc + births$bidx * birth_displace
  
  epis <- tmerge(mf, mf, id=`(id)`, tstart=`(dob)`, tstop=`(intv)`)
  epis <- tmerge(epis, births, id=`(id)`, birth = event(bcmc))
  
  aggr <- demog_pyears(f, epis, period=period, agegr=agegr, cohort=cohort, tips=tips,
                       event="birth", weights="(weights)", origin=origin, scale=scale)$data
  
  ## construct interaction of all factor levels that appear
  byvar <- intersect(c(all.vars(by), "agegr", "period", "cohort", "tips"),
                     names(aggr))
  aggr$byf <- interaction(aggr[byvar], drop=TRUE)
  
  ## prediction for all factor levels that appear
  pred <- data.frame(aggr[c(byvar, "byf")])[!duplicated(aggr$byf),]
  pred <- pred[order(pred$byf), ]
  
  if(counts || varmethod == "none"){
    mc <- model.matrix(~-1+byf, aggr)
    clong <- aggr[c("event", "pyears")]
    pred[c("births", "pys")] <- t(mc) %*% as.matrix(clong)
  }
  
  if(varmethod == "none") {
    
    pred$asfr <- pred$births / pred$pys
    pred$byf <- NULL
    if(!counts)
      pred[c("births", "pys")] <- NULL
    
  } else if(varmethod == "lin") {
    
    des <- survey::svydesign(ids=clusters, strata=strata, data=aggr, weights=~1)
    class(des) <- c("svypyears", class(des))
    
    ## fit model
    f <- if(length(levels(aggr$byf)) == 1)
      event ~ offset(log(pyears))
    else
      event ~ -1 + byf + offset(log(pyears))
    
    mod <- survey::svyglm(f, des, family=quasipoisson)
    
    ## prediction for all factor levels that appear
    pred$pyears <- 1
    
    asfr <- predict(mod, pred, type="response", vcov=TRUE)
    v <- vcov(asfr)
    dimnames(v) <- list(pred$byf, pred$byf)
    
    pred$asfr <- as.numeric(asfr)
    pred$se_asfr <- sqrt(diag(v))
    pred[c("byf", "pyears")] <- NULL
    attr(pred, "var") <- v
  } else if(varmethod %in% c("jkn", "jk1")) {
    
    ## Convert to array with events and PYs for each cluster
    ## reshape2::acast is MUCH faster than stats::reshape
    events_clust <- reshape2::acast(aggr, update(clusters, byf ~ .), value.var="event")
    pyears_clust <- reshape2::acast(aggr, update(clusters, byf ~ .), value.var="pyears")
    
    if(varmethod == "jkn"){
      aggr$strataid <- as.integer(interaction(aggr[all.vars(strata)], drop=TRUE))
      strataid <- drop(reshape2::acast(unique(aggr[c(all.vars(clusters), "strataid")]),
                                       update(clusters,  1 ~ .), value.var="strataid"))
    } else
      strataid <- NULL
    
    estdf <- jackknife(events_clust, pyears_clust, strataid)
    
    pred$asfr <- estdf$est
    pred$se_asfr <- estdf$se
    attr(pred, "var") <- vcov(estdf)
    pred$byf <- NULL
    if(clustcounts){
      attr(pred, "events_clust") <- events_clust
      attr(pred, "pyears_clust") <- pyears_clust
      attr(pred, "strataid") <- strataid
    }
  } else
    stop(paste0("varmethod = \"", varmethod, "\" is not recognized."))
  
  
  rownames(pred) <- NULL
  
  return(pred)
}

ir_by_area2 <- function(ir, area_list, n, total) {
  
  print(paste(n, "of", total))
  
  ir_int <- ir %>%
    left_join(area_list, by=c("v001" = "cluster_id")) %>%
    filter(!is.na(area_id)) %>%
    group_split(area_id)
  
  return(ir_int)
  
}


get_pred <- function(mod_list, asfr_pred, asfr1, ...) {
  
  pred_size <- nrow(asfr_pred) - nrow(asfr1)
  
  if(subnational == FALSE) {
    
    if(multicountry==FALSE) {
    
    pred <- asfr_pred %>%
      filter(id<pred_size+1) %>%
      dplyr::select("country", "agegr", "period", "pys","id") %>%
      left_join(mod_list$summary.fitted.values[1:pred_size, ] %>%
                  mutate(id = 1:pred_size), by="id") %>%
      arrange(country, period, agegr) %>%
      mutate(agegr = factor(agegr))
    
    } else {
      
      pred <- asfr_pred %>%
        filter(id<pred_size+1) %>%
        dplyr::select("agegr", "period", "pys","id") %>%
        left_join(mod_list$summary.fitted.values[1:pred_size, ] %>%
                    mutate(id = 1:pred_size), by="id") %>%
        arrange(period, agegr) %>%
        mutate(agegr = factor(agegr))
      
    }
  } else {
    
    pred <- asfr_pred %>%
      filter(id<pred_size+1) %>%
      dplyr::select("area_name", "agegr", "period", "pys","id") %>%
      left_join(mod_list$summary.fitted.values[1:pred_size, ] %>%
                  mutate(id = 1:pred_size), by="id") %>%
      arrange(period, agegr) %>%
      mutate(agegr = factor(agegr))
    
  }
  
  return(pred)
  
}
#iso3_list, multicountry
run_mod_nat <- function(formula, asfr_pred) {
  mod <- inla(formula, family="poisson", data=asfr_pred, E=pys,
              control.family=list(link='log'),
              control.predictor=list(compute=TRUE, link=1),
              control.inla = list(strategy = "gaussian", int.strategy = "eb"),
              control.compute=list(config = TRUE, dic= TRUE, cpo=TRUE),
              verbose=TRUE)

  return(mod)
}


get_mod_results <- function(mod, asfr_pred_country_subnat, pop_areas, areas_long, population_age_female) {
  
  asfr1_country_subnat <- asfr_pred_country_subnat %>%
    filter(!is.na(surveyid))
  
  iso3 <- ifelse(unique(asfr_pred_country_subnat$country) == "Eswatini", "SWZ", countrycode(unique(asfr_pred_country_subnat$country), "country.name", "iso3c"))
  
  print("Sampling..")
  samples <- inla.posterior.sample(1000, mod)
  print("Done sampling")
  contents = mod$misc$configs$contents
  effect = "Predictor"
  id.effect = which(contents$tag==effect)
  ind.effect = contents$start[id.effect]-1 + (1:contents$length[id.effect])
  
  ind.effect <- 1:(nrow(asfr_pred_country_subnat) - nrow(asfr1_country_subnat))
  
  samples.effect = lapply(samples, function(x) x$latent[ind.effect] %>% exp)
  
  ident <- asfr_pred_country_subnat[ind.effect, ] %>%
    select(area_id, agegr, period) %>%
    mutate(iso3 = iso3) %>%
    left_join(pop_areas) %>%
    filter(period>1999) %>%
    mutate(ratio0 = population/id0_agepop,
           ratio1 = population/id1_agepop,
           ratio2 = population/id2_agepop,
           ratio3 = population/id3_agepop,
           ratio4 = population/id4_agepop,
           ratio5 = population/id5_agepop
    )
  
  samples_ident <- sapply(samples.effect, cbind) %>%
    data.frame %>%
    mutate(id = 1:nrow(.)) %>%
    left_join(asfr_pred_country_subnat[ind.effect , c(1:4, 25)], by="id") %>%
    select(-id) %>%
    left_join(ident) %>%
    filter(period >1999)
  
  #####
  
  admin0  <- samples_ident
  admin0[,1:1000] <- samples_ident[,1:1000] * samples_ident$ratio0
  
  admin0_asfr <- admin0 %>%
    select(1:1000, area_id0, period, agegr) %>%
    group_by(area_id0, period, agegr) %>%
    summarise_all(funs(sum))
  
  admin0_asfr$mean <- apply(admin0_asfr[, 4:1003], 1, mean)
  admin0_asfr$median <- apply(admin0_asfr[, 4:1003], 1, median)
  admin0_asfr$sd <- apply(admin0_asfr[, 4:1003], 1, sd)
  
  admin0_asfr <- admin0_asfr %>%
    mutate(lower = mean-(qnorm(0.95)*sd),
           upper = mean+(qnorm(0.95)*sd)
    ) %>%
    select(area_id0, period, agegr, mean, median, sd, lower, upper) %>%
    rename(area_id = area_id0)
  
  ###########
  
  admin1  <- samples_ident
  admin1[,1:1000] <- samples_ident[,1:1000] * samples_ident$ratio1
  
  admin1_asfr <- admin1 %>%
    select(1:1000, area_id1, period, agegr) %>%
    group_by(area_id1, period, agegr) %>%
    summarise_all(funs(sum))
  
  admin1_asfr$mean <- apply(admin1_asfr[, 4:1003], 1, mean)
  admin1_asfr$median <- apply(admin1_asfr[, 4:1003], 1, median)
  admin1_asfr$sd <- apply(admin1_asfr[, 4:1003], 1, sd)
  
  admin1_asfr <- admin1_asfr %>%
    mutate(lower = mean-(qnorm(0.95)*sd),
           upper = mean+(qnorm(0.95)*sd)
    ) %>%
    select(area_id1, period, agegr, mean, median, sd, lower, upper) %>%
    rename(area_id = area_id1)
  
  ##########
  admin2  <- samples_ident
  admin2[,1:1000] <- samples_ident[,1:1000] * samples_ident$ratio2
  
  admin2_asfr <- admin2 %>%
    select(1:1000, area_id2, period, agegr) %>%
    group_by(area_id2, period, agegr) %>%
    summarise_all(funs(sum))
  
  admin2_asfr$mean <- apply(admin2_asfr[, 4:1003], 1, mean)
  admin2_asfr$median <- apply(admin2_asfr[, 4:1003], 1, median)
  admin2_asfr$sd <- apply(admin2_asfr[, 4:1003], 1, sd)
  
  admin2_asfr  <- admin2_asfr %>%
    mutate(lower = mean-(qnorm(0.95)*sd),
           upper = mean+(qnorm(0.95)*sd)
    ) %>%
    select(area_id2, period, agegr, mean, median, sd, lower, upper) %>%
    rename(area_id = area_id2)
  
  
  ##############
  
  
  admin3  <- samples_ident
  admin3[,1:1000] <- samples_ident[,1:1000] * samples_ident$ratio3
  
  admin3_asfr <- admin3 %>%
    select(1:1000, area_id3, period, agegr) %>%
    group_by(area_id3, period, agegr) %>%
    summarise_all(funs(sum))
  
  admin3_asfr$mean <- apply(admin3_asfr[, 4:1003], 1, mean)
  admin3_asfr$median <- apply(admin3_asfr[, 4:1003], 1, median)
  admin3_asfr$sd <- apply(admin3_asfr[, 4:1003], 1, sd)
  
  admin3_asfr  <- admin3_asfr %>%
    mutate(lower = mean-(qnorm(0.95)*sd),
           upper = mean+(qnorm(0.95)*sd)
    ) %>%
    select(area_id3, period, agegr, mean, median, sd, lower, upper) %>%
    rename(area_id = area_id3)
  
  #############
  
  admin4  <- samples_ident
  admin4[,1:1000] <- samples_ident[,1:1000] * samples_ident$ratio4
  
  admin4_asfr <- admin4 %>%
    select(1:1000, area_id4, period, agegr) %>%
    group_by(area_id4, period, agegr) %>%
    summarise_all(funs(sum))
  
  admin4_asfr$mean <- apply(admin4_asfr[, 4:1003], 1, mean)
  admin4_asfr$median <- apply(admin4_asfr[, 4:1003], 1, median)
  admin4_asfr$sd <- apply(admin4_asfr[, 4:1003], 1, sd)
  
  admin4_asfr  <- admin4_asfr %>%
    mutate(lower = mean-(qnorm(0.95)*sd),
           upper = mean+(qnorm(0.95)*sd)
    ) %>%
    select(area_id4, period, agegr, mean, median, sd, lower, upper) %>%
    rename(area_id = area_id4)
  
  #############
  
  admin5  <- samples_ident
  admin5[,1:1000] <- samples_ident[,1:1000] * samples_ident$ratio5
  
  admin5_asfr <- admin5 %>%
    select(1:1000, area_id5, period, agegr) %>%
    group_by(area_id5, period, agegr) %>%
    summarise_all(funs(sum))
  
  admin5_asfr$mean <- apply(admin5_asfr[, 4:1003], 1, mean)
  admin5_asfr$median <- apply(admin5_asfr[, 4:1003], 1, median)
  admin5_asfr$sd <- apply(admin5_asfr[, 4:1003], 1, sd)
  
  admin5_asfr  <- admin5_asfr %>%
    mutate(lower = mean-(qnorm(0.95)*sd),
           upper = mean+(qnorm(0.95)*sd)
    ) %>%
    select(area_id5, period, agegr, mean, median, sd, lower, upper) %>%
    rename(area_id = area_id5)
  
  
  ####################################
  
  admin0_tfr <- admin0 %>%
    select(1:1000, area_id0, period) %>%
    group_by(area_id0, period) %>%
    summarise_all(function(x){5*sum(x)})
  
  admin0_tfr$mean <- apply(admin0_tfr[, 4:1002], 1, mean)
  admin0_tfr$median <- apply(admin0_tfr[, 4:1002], 1, median)
  admin0_tfr$sd <- apply(admin0_tfr[, 4:1002], 1, sd)
  
  admin0_tfr <- admin0_tfr %>%
    mutate(lower = mean-(qnorm(0.95)*sd),
           upper = mean+(qnorm(0.95)*sd)
    ) %>%
    select(area_id0, period, mean, median, sd, lower, upper) %>%
    rename(area_id = area_id0)
  
  ###########
  
  admin1_tfr <- admin1 %>%
    select(1:1000, area_id1, period) %>%
    group_by(area_id1, period) %>%
    summarise_all(function(x){5*sum(x)})
  
  
  admin1_tfr$mean <- apply(admin1_tfr[, 4:1002], 1, mean)
  admin1_tfr$median <- apply(admin1_tfr[, 4:1002], 1, median)
  admin1_tfr$sd <- apply(admin1_tfr[, 4:1002], 1, sd)
  
  admin1_tfr <- admin1_tfr %>%
    mutate(lower = mean-(qnorm(0.95)*sd),
           upper = mean+(qnorm(0.95)*sd)
    ) %>%
    select(area_id1, period, mean, median, sd, lower, upper) %>%
    rename(area_id = area_id1)
  
  #############
  
  admin2_tfr <- admin2 %>%
    select(1:1000, area_id2, period) %>%
    group_by(area_id2, period) %>%
    summarise_all(function(x){5*sum(x)})
  
  admin2_tfr$mean <- apply(admin2_tfr[, 4:1002], 1, mean)
  admin2_tfr$median <- apply(admin2_tfr[, 4:1002], 1, median)
  admin2_tfr$sd <- apply(admin2_tfr[, 4:1002], 1, sd)
  
  admin2_tfr  <- admin2_tfr %>%
    mutate(lower = mean-(qnorm(0.95)*sd),
           upper = mean+(qnorm(0.95)*sd)
    ) %>%
    select(area_id2, period, mean, median, sd, lower, upper) %>%
    rename(area_id = area_id2)
  
  ##############
  
  admin3_tfr <- admin3 %>%
    select(1:1000, area_id3, period) %>%
    group_by(area_id3, period) %>%
    summarise_all(function(x){5*sum(x)})
  
  admin3_tfr$mean <- apply(admin3_tfr[, 4:1002], 1, mean)
  admin3_tfr$median <- apply(admin3_tfr[, 4:1002], 1, median)
  admin3_tfr$sd <- apply(admin3_tfr[, 4:1002], 1, sd)
  
  admin3_tfr  <- admin3_tfr %>%
    mutate(lower = mean-(qnorm(0.95)*sd),
           upper = mean+(qnorm(0.95)*sd)
    ) %>%
    select(area_id3, period, mean, median, sd, lower, upper) %>%
    rename(area_id = area_id3)
  
  
  ##########
  
  admin4_tfr <- admin4 %>%
    select(1:1000, area_id4, period) %>%
    group_by(area_id4, period) %>%
    summarise_all(function(x){5*sum(x)})
  
  admin4_tfr$mean <- apply(admin4_tfr[, 4:1002], 1, mean)
  admin4_tfr$median <- apply(admin4_tfr[, 4:1002], 1, median)
  admin4_tfr$sd <- apply(admin4_tfr[, 4:1002], 1, sd)
  
  admin4_tfr  <- admin4_tfr %>%
    mutate(lower = mean-(qnorm(0.95)*sd),
           upper = mean+(qnorm(0.95)*sd)
    ) %>%
    select(area_id4, period, mean, median, sd, lower, upper) %>%
    rename(area_id = area_id4)
  
  
  #############
  
  
  admin5_tfr <- admin5 %>%
    select(1:1000, area_id5, period) %>%
    group_by(area_id5, period) %>%
    summarise_all(function(x){5*sum(x)})
  
  admin5_tfr$mean <- apply(admin5_tfr[, 4:1002], 1, mean)
  admin5_tfr$median <- apply(admin5_tfr[, 4:1002], 1, median)
  admin5_tfr$sd <- apply(admin5_tfr[, 4:1002], 1, sd)
  
  admin5_tfr  <- admin5_tfr %>%
    mutate(lower = mean-(qnorm(0.95)*sd),
           upper = mean+(qnorm(0.95)*sd)
    ) %>%
    select(area_id5, period, mean, median, sd, lower, upper) %>%
    rename(area_id = area_id5)
  
  
  #######
  
  mod_results <- list()
  
  mod_results$asfr <- admin0_asfr %>%
    bind_rows(admin1_asfr, admin2_asfr, admin3_asfr, admin4_asfr, admin5_asfr) %>%
    left_join(areas_long %>% select(-parent_area_id), by="area_id") %>%
    mutate(source = "Model",
           variable = "asfr") %>%
    #select(-c(mean, sd)) %>%
    filter(!is.na(area_id))
  
  mod_results$tfr <- admin0_tfr %>%
    bind_rows(admin1_tfr, admin2_tfr, admin3_tfr, admin4_tfr, admin5_tfr) %>%
    left_join(areas_long %>% select(-parent_area_id), by="area_id") %>%
    mutate(source = "Model",
           variable = "tfr") %>%
    #select(-c(mean, sd)) %>%
    filter(!is.na(area_id))
  
  
  mod_results$births_by_age <- mod_results$asfr %>%
    #select(-variable) %>%
    left_join(population_age_female %>% select(area_id, agegr, period, population)) %>%
    mutate(val = median*population,
           lower = lower*population,
           upper = upper*population,
           variable = "births",
           source = "Model"
    ) %>%
    select(-c(population, median))
  
  mod_results$births <- mod_results$births_by_age %>%
    group_by(iso3, area_id, area_name, area_level, period, source, variable) %>%
    summarise(val = sum(val),
              lower = sum(lower),
              upper = sum(upper))
  
  
  return(mod_results)
  
}


get_mod_results_test <- function(mod, asfr_pred_country_subnat) {

  
  asfr1_country_subnat <- asfr_pred_country_subnat %>%
    filter(!is.na(surveyid))
  
  print("Sampling..")
  samples <- inla.posterior.sample(1000, mod)
  print("Done sampling")
  contents = mod$misc$configs$contents
  effect = "Predictor"
  id.effect = which(contents$tag==effect)
  ind.effect = contents$start[id.effect]-1 + (1:contents$length[id.effect])
  
  ind.effect <- 1:(nrow(asfr_pred_country_subnat) - nrow(asfr1_country_subnat))
  
  samples.effect = lapply(samples, function(x) x$latent[ind.effect] %>% exp)
  
  ident <- asfr_pred_country_subnat[ind.effect, ] %>%
    select(area_id, age_group, period)
  
  qtls <- apply(sapply(samples.effect, cbind), 1, quantile, c(0.025, 0.5, 0.975))
  
  samples_ident <- ident %>%
    mutate(lower = qtls[1,],
           median = qtls[2,],
           upper = qtls[3,]
    )
  
  return(samples_ident)
  
}

get_boundaries <- function(iso3_code) {
  
  area_file <- paste0("~/Documents/GitHub/naomi-data/", iso3_code, "/data/", tolower(iso3_code), "_areas.geojson")
  
  boundary_sf <- read_sf(area_file) %>%
    mutate(iso3 = iso3_code
    ) %>%
    select(iso3, area_id, area_name, area_level, parent_area_id, naomi_level, geometry)
  
  return(boundary_sf)
}

get_wide_areas <- function(iso3_code) {
  
  area_file <- paste0("~/Documents/GitHub/naomi-data/", iso3_code, "/data/", tolower(iso3_code), "_areas.geojson")
  
  areas <- read_sf(area_file)
  
  spread_areas(as.data.frame(areas)) %>%
    left_join(select(areas, area_id))
  
}


load_population_agesex <- function(pop_agesex_wide_path, areas_path) {

    wide <- readRDS(pop_agesex_wide_path)

    long <- tidyr::gather(wide, area_id, population,
                          tidyselect::matches("^[A-Z]{3}[^a-zA-Z]*$",
                                              ignore.case = FALSE))

    if(!is.null(areas_path)) {
      areas <- readRDS(areas_path)
      long <- dplyr::left_join(long,
                               dplyr::select(areas, -parent_area_id),
                               by = "area_id")
      long <- dplyr::select(long, iso3, area_id, area_level, area_name,
                            dplyr::everything())
    }

    long
  }


interpolate_fertility_population <- function(population_agesex, calendar_quarters) {
  
  quarter_ids <- calendar_quarter_to_quarter_id(calendar_quarters)
  dfall <- dplyr::distinct(dplyr::select(population_agesex, -calendar_quarter, -population, -quarter_id))
  
  df <- dplyr::select(population_agesex, calendar_quarter, area_id, source, sex, age_group, population) %>%
    dplyr::mutate(quarter_id = calendar_quarter_to_quarter_id(calendar_quarter))
  
  tidyr::expand(df,
                tidyr::nesting(calendar_quarter = calendar_quarters, quarter_id = quarter_ids),
                tidyr::nesting(area_id, source, sex, age_group)) %>%
    dplyr::full_join(df, by = names(.)) %>%
    dplyr::group_by(area_id, source, sex, age_group) %>%
    dplyr::mutate(population = exp(zoo::na.approx(log(population), quarter_id, na.rm = FALSE)),
                  population = tidyr::replace_na(population, 0)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(quarter_id %in% quarter_ids) %>%
    mutate(period = year_labels(quarter_id)) %>%
    left_join(get_age_groups() %>% select(age_group, age_group_label)) %>%
    select(area_id, period, age_group_label, population) %>%
    rename(agegr = age_group_label)
  
}

sample_tmb_test <- function(fit, nsample = 1000, rng_seed = NULL,
                            random_only = FALSE, verbose = FALSE) {
  
  
  if(!random_only) {
    if(verbose) print("Calculating joint precision")
    hess <- sdreport_joint_precision(fit$obj, fit$par.fixed)
    
    if(verbose) print("Inverting precision for joint covariance")
    cov <- solve(hess)
    
    if(verbose) print("Drawing sample")
    ## TODO: write a version of rmvnorm that uses precision instead of covariance
    smp <- mvtnorm::rmvnorm(nsample, fit$par.full, as.matrix(cov))
    
  } else {
    r <- fit$obj$env$random
    par_f <- fit$par.full[-r]
    
    par_r <- fit$par.full[r]
    hess_r <- fit$obj$env$spHess(fit$par.full, random = TRUE)
    smp_r <- rmvnorm_sparseprec(nsample, par_r, hess_r)
    
    smp <- matrix(0, nsample, length(fit$par.full))
    smp[ , r] <- smp_r
    smp[ ,-r] <- matrix(par_f, nsample, length(par_f), byrow = TRUE)
    colnames(smp)[r] <- colnames(smp_r)
    colnames(smp)[-r] <- names(par_f)
  }
  
  if(verbose) print("Simulating outputs")
  sim <- apply(smp, 1, fit$obj$report)
  
  r <- fit$obj$report()
  
  if(verbose) print("Returning sample")
  fit$sample <- Map(vapply, list(sim), "[[", lapply(lengths(r), numeric), names(r))
  is_vector <- vapply(fit$sample, class, character(1)) == "numeric"
  fit$sample[is_vector] <- lapply(fit$sample[is_vector], as.matrix, nrow = 1)
  names(fit$sample) <- names(r)
  
  fit
}

make_adjacency_matrix <- function(iso3_current, areas_long, boundaries, exclude_districts = exc, level="naomi") {
  
  if (level == "naomi") {
    
    int <- areas_long %>%
      filter(iso3 == iso3_current)
    
    level <- unique(int$area_level[int$naomi_level == TRUE])
    
  }
  
  sh <- areas_long %>%
    filter(iso3 == iso3_current, area_level == level, !area_id %in% exclude_districts) %>%
    mutate(area_idx = row_number())
  
  #' Neighbor list
  nb <- sh %>%
    left_join(boundaries) %>%
    st_as_sf %>%
    as("Spatial") %>%
    spdep::poly2nb() %>%
    `names<-`(sh$area_idx)
  
  if(iso3_current == "MOZ") {
    #Make KaTembe adjacent to KaMpfumu and Nhlamankulu
    nb[[6]] <- c(nb[[6]], 1, 2)
    nb[[1]] <- c(nb[[1]], 6)
    nb[[2]] <- c(nb[[2]], 6)
    
    #Make KaNyaka adjacent to KaMpfumu	and KaMaxakeni
    nb[[7]] <- c(1, 3)
    nb[[1]] <- c(nb[[1]], 7)
    nb[[3]] <- c(nb[[3]], 7)
    
    nb <- lapply(nb, as.integer)
    class(nb) <- "nb"
  }
  
  adj <- nb2mat(nb, zero.policy=TRUE, style="B")
  R_spatial <- INLA::inla.scale.model(diag(rowSums(adj)) - 0.99*adj,
                                      constr = list(A = matrix(1, 1, nrow(adj)), e = 0))
  
  return(R_spatial)
  
}


make_rw_structure_matrix <- function(x, order, adjust_diagonal = TRUE) {
  
  
  
  D_mat <- diff(diag(x), differences = order)
  R_mat <- t(D_mat) %*% D_mat
  
  if(adjust_diagonal) {
    diag(R_mat) <- diag(R_mat) + 1E-6
  }
  
  R_mat <- as(R_mat, "dgCMatrix")
  
  return(R_mat)
  
}

area_populations <- function(population, areas_wide) {
  
  base_area_pop <- areas_wide %>%
    left_join(population)
  
  level_ids <- str_subset(colnames(areas_wide), "area_id[0-9]")
  
  group_area_pops <- function(level_ids, base_area_pop) {
    
    base_area_pop %>%
      group_by(.data[[level_ids]], sex, age_group, period) %>%
      summarise(population = sum(population)) %>%
      rename(area_id = .data[[level_ids]])
    
  }
  
  area_populations <- lapply(level_ids, group_area_pops, base_area_pop) %>%
    bind_rows
  
}

make_model_frames <- function(iso3_current, population, asfr, mics_asfr = NULL, exclude_districts = "", project = FALSE) {
  
  population <- area_populations(population, areas_wide) %>%
    filter(sex == "female") %>%
    ungroup %>%
    select(-sex)
  
  population <- crossing(area_id = unique(population$area_id),
           age_group = unique(population$age_group),
           period = 1995:2020
  ) %>%
    left_join(population) %>%
    group_by(area_id, age_group) %>%
    mutate(population = exp(zoo::na.approx(log(population), period, na.rm = FALSE))) %>%
    fill(population, .direction="up")
    
  
  # population <- crossing(area_id = unique(population$area_id),
  #          age_group = unique(population$age_group),
  #          period = 1995:(min(population$period)-1)
  # ) %>%
  #   left_join(population %>% filter(period == min(period)) %>% select(-period)) %>%
  #   bind_rows(population)
   
  if(!project) {
  
  if(!is.null(mics_asfr)) {
    df <- asfr %>%
      bind_rows(mics_asfr)
  } else {
    df <- asfr
    
  }
  
  max_year <- max(df$period)
    
  } else {
    
    max_year <- project
  }
  
  area_merged <- st_read(file.path(naomi_data_path, iso3_current, "data", paste0(tolower(iso3_current), "_areas.geojson")))
  areas <- create_areas(area_merged = area_merged)
  area_aggregation <- create_area_aggregation(area_merged$area_id[area_merged$naomi_level], areas) %>%
    filter(!model_area_id %in% exclude_districts)
  
  ## Make model frame.
  mf_model <- crossing(period = 1995:max_year,
                 age_group = c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49"),
                 # area_id = unique(area_aggregation$model_area_id)) %>%
                 # area_id = filter(areas_long, iso3 == iso3_current, area_level == 1)$area_id) %>%
                 area_id = iso3_current) %>%
    left_join(population %>%
                select(area_id, period, age_group, population)
    ) %>%
    # mutate(area_id = factor(area_id, levels = unique(area_aggregation$model_area_id)),
    # mutate(area_id = factor(area_id, levels = filter(areas_long, iso3 == iso3_current, area_level == 1)$area_id),
    mutate(area_id = factor(iso3_current),
           age_group = factor(age_group, levels = c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49")),
           period = factor(period)
           # restype = ifelse(area_id %in% c(
           #  filter(areas_long, parent_area_id == "ETH_1_10")$area_id,
           #  filter(areas_long, str_detect(area_name, "Town"))$area_id,
           #  filter(areas_long, area_name %in% c("Harari", "Fafen (Jijiga)"))$area_id),
           #  1, 0)
    ) %>%
    arrange(period, area_id, age_group) %>%
    mutate(idx = factor(row_number()),
           id.interaction_3d = factor(group_indices(., age_group, period, area_id)),
           # id.interaction_age_time = factor(group_indices(., age_group, period)),
           id.interaction1 = factor(group_indices(., age_group, period)),
           id.interaction2 = factor(group_indices(., period, area_id)),
           id.interaction3 = factor(group_indices(., age_group, area_id))
    ) %>%
   droplevels()
  
  obs <- asfr %>%
    mutate(period = factor(period, levels(mf_model$period))) %>%
    filter(!is.na(surveyid), !area_id %in% exclude_districts) %>%
    select(area_id, period, age_group, tips, births, pys) %>%
    left_join(mf_model) %>%
    mutate(tips_dummy = as.integer(tips > 5),
           tips_f = factor(tips),
           #####
           # urban_dummy = ifelse(area_id %in% c(
           #   filter(areas_long, parent_area_id == "ETH_1_10")$area_id,
           #   filter(areas_long, str_detect(area_name, "Town"))$area_id,
           #   filter(areas_long, area_name %in% c("Harari", "Fafen (Jijiga)"))$area_id),
           # 1, 0),
           #####
           age_group = factor(age_group, levels(mf_model$age_group)),
           area_id = factor(area_id, levels(mf_model$area_id)),
           period = factor(period, levels(mf_model$period)),
    )
  
  mf <- list()
  mf$mf_model <- mf_model
  mf$dist$obs <- obs
  mf$mics_toggle <- 0
  mf$out_toggle <- 0
  
  ## Outputs

  if(unique(asfr %>%left_join(areas_long) %>% .$naomi_level)) {  
    
    mf_out <- crossing(
      area_id = area_aggregation$area_id,
      age_group = unique(mf_model$age_group),
      period = unique(mf_model$period)
    ) %>%
      arrange(area_id, age_group, period) %>%
      mutate(out_idx = row_number()) %>%
      droplevels()
  
    join_out <- crossing(area_aggregation,
                         age_group = unique(mf_model$age_group),
                         period = unique(mf_model$period)) %>%
      full_join(mf_model %>%
                  select(area_id, age_group, period, idx), by = c("model_area_id" = "area_id",
                                                                  "age_group",
                                                                  "period")
      ) %>%
      full_join(mf_out) %>%
      # full_join(mf_out, by=c("area_id",
      #                        "period",
      #                        "age_group_out" = "age_group")
      #           ) %>%
      mutate(x=1) %>%
      filter(!is.na(model_area_id))
    
    A_out <- spMatrix(nrow(mf_out), nrow(mf_model), join_out$out_idx, as.integer(join_out$idx), join_out$x)
    
    # mf_out_restype <- crossing(
    #   age_group = unique(mf_model$age_group),
    #   period = unique(mf_model$period),
    #   restype = c(1, 0)
    #   ) %>%
    # arrange(age_group, period) %>%
    # mutate(out_idx = row_number()) %>%
    # droplevels()
    # 
    # join_out_restype <- crossing(
    #   age_group = unique(mf_model$age_group),
    #   period = unique(mf_model$period),
    #   restype = c(1, 0)
    # ) %>%
    # full_join(mf_model %>%
    #             select(age_group, period, restype, idx)) %>%
    # full_join(mf_out_restype) %>%
    # # full_join(mf_out, by=c("area_id",
    # #                        "period",
    # #                        "age_group_out" = "age_group")
    # #           ) %>%
    # mutate(x=1)
    # 
    # A_out_restype <- spMatrix(nrow(mf_out_restype), nrow(mf_model), join_out_restype$out_idx, as.integer(join_out_restype$idx), join_out_restype$x)
    # 
    # mf$out$mf_out_restype <- mf_out_restype
    # mf$out$A_out_restype <- A_out_restype
    
    mf$out$mf_out <- mf_out
    mf$out$A_out <- A_out
    mf$out_toggle <- 1
  
  }
  
  if(!is.null(mics_asfr)) {
    
    
    
    mf_mics <- crossing(area_id = unique(mics_asfr$area_id),
                       period = unique(mf_model$period),
                       age_group = unique(mf_model$age_group)
    ) %>%
      filter(!area_id %in% unique(filter(areas_wide, area_id %in% exclude_districts)$area_id1)) %>%
      mutate(idx = factor(row_number()))
    
    join_mics <- mf_mics %>%
      rename(idx_row = idx) %>%
      left_join(area_aggregation) %>%
      left_join(mf_model, by=c("age_group", "period", "model_area_id" = "area_id")) %>%
      mutate(idx_col = row_number(),
             x=1) %>%
      type.convert()
    
    A_mics <- sparseMatrix(i = join_mics$idx_row, j=join_mics$idx_col, x=join_mics$x, use.last.ij = TRUE)
    
    
    obs_mics <- mics_asfr %>%
      filter(!area_id %in% unique(filter(areas_wide, area_id %in% exclude_districts)$area_id1)) %>%
      mutate(period = factor(period, levels(mf_model$period))) %>%
      left_join(mf_mics) %>%
      select(area_id, period, age_group, tips, births, pys, idx) %>%
      mutate(tips_dummy = as.integer(tips > 2),
             tips_f = factor(tips, levels(obs$tips_f)),
             age_group = factor(age_group, levels(mf_model$age_group)),
             idx =factor(idx, levels(mf_mics$idx))
      )
    
    mf$mics$obs <- obs_mics
    mf$mics$A_mics <- A_mics
    mf$mics_toggle <- 1
    
  }
  
  return(mf)
}
