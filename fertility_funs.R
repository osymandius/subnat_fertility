calc_asfr1 <- function(data,
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


repeat_formulae <- function(formulae, formulae_number) {
  
  form1 <- list()
  
  form1[[1]] <- formulae[[1]]
  
  form2 <- list()
  
  form2[[1]] <- formulae[[2]]
  
  form3 <- list()
  
  form3[[1]] <- formulae[[3]]
  
  form4 <- list()
  
  form4[[1]] <- formulae[[4]]
    
  
  form1 <- rep(form1, formulae_number)
  form2 <- rep(form2, formulae_number)
  form3 <- rep(form3, formulae_number)
  form4 <- rep(form4, formulae_number)
  
  formulae <- c(form1, form2, form3, form4)
  
  return(formulae)
}

get_pred <- function(mod_list, asfr_pred, asfr1) {
  
  pred_size <- nrow(asfr_pred) - nrow(asfr1)
  
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
  
  return(pred)
  
}
#iso3_list, multicountry
run_mod <- function(formulae, asfr_pred) {
  mod <- inla(formulae, family="poisson", data=asfr_pred, E=pys,
              control.family=list(link='log'),
              control.predictor=list(compute=TRUE, link=1),
              control.inla = list(strategy = "gaussian", int.strategy = "eb"),
              control.compute=list(config = TRUE, dic= TRUE, cpo=TRUE),
              verbose=TRUE)
  
  # if(multicountry) {
  #   print(iso3_list)
  # }
    
  return(mod)
}