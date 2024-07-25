##Airborne spread

##All parameters are given for spread between each pair of production types
#curr_state = current disease state at each time: S/S_Infected/L/B/C/R
#max_distance = the maximum distance of airborne spread: n-by-n matrix (only linear dropoff)
#pop = main population data: dataframe(prod_type,unitsize,Lat,Lon,Init_state)
#dropoff = "linear"/"exp": (ntype*ntype) matrix
#deg_start = area at risk of exposure start(degree): n-by-n matrix
#deg_end = area at risk of exposure end(degree): n-by-n matrix
#p = probability of spread 1 km apart: (ntype*ntype) matrix
#delay_pdf = names of pdf for airborne transport delay (days) or "fixed": (ntype*ntype) matrix
#delay_par = pdf parameters for airborne transport delay (days) or days(in case delay_pdf=="fixed"): (ntype*ntype*d) array
#prev_chart_use = a vector of length ntype: "Yes"/"No"
#prev = ntype-to-(total number of days spent in disease states) data.frame: prevalence chart by production types
#status = N-by-sim_day data.frame: disease status from day 1 up to current day
airspread <- function(curr_state,max_distance,pop,dropoff,deg_start,deg_end,p,delay_pdf,delay_par,prev_chart_use,prev,status){
  disease_status <- curr_state
  source_index <- which(disease_status=="B" | disease_status=="C") #source_index
  if(is_empty(source_index)){
    infected_source <- NA
    exposure_target <- NA
    infected_target <- NA
    delays <- NA
    source_prod_type <- NA
    target_prod_type <- NA
    results <- data.frame(infected_source,exposure_target,infected_target,delays,source_prod_type,target_prod_type)
    outputs <- list(Infection = results, Status = curr_state)
    return(outputs)
  }else{ #at least one source unit
    infected_source <- c()
    exposure_target <- c()
    infected_target  <- c()
    delays <- c()
    source_type <- as.vector(as.matrix(pop[source_index,1])) #source production type
    for (i in 1:length(source_index)) {  #for each source unit i
      sus_ind  <- which(disease_status=="S")
      if(!is_empty(sus_ind)){
        target_type <- as.vector(as.matrix(pop[sus_ind,1])) #target production type
        #compute the distance between source i and all possible targets
        d <- rep(NA,length(sus_ind))
        for (k in 1:length(sus_ind)) {
          x <- (as.numeric(pop[source_index[i],4])-as.numeric(pop[sus_ind[k],4])) * cos(as.numeric(pop[sus_ind[k],3]))
          y <- as.numeric(pop[source_index[i],3]) - as.numeric(pop[sus_ind[k],3])
          d[k] <- abs(24901/360 * sqrt(x^2 + y^2))
        } #some calculated distance might be NA
        #calculate the direction of target from the source
        direction <- directionfromsource(as.numeric(pop[source_index[i],3]),as.numeric(pop[source_index[i],4]),as.vector(as.matrix(pop[sus_ind,3])),
                                         as.vector(as.matrix(pop[sus_ind,4])))
        for (j in 1:length(sus_ind)) { #for each pair of source-target
          if(dropoff[source_type[i],target_type[j]]=="linear"){
            #check if recipient unit with distance is within desired max distance
            if(!is.na(d[j]) & d[j] < max_distance[source_type[i],target_type[j]]){
              #check if target j shares the same location with the source and/or unknown direction from source and/or outside the range of area at risk
              if(!is.na(deg_start[source_type[i],target_type[j]])){
                if(!is.na(deg_end[source_type[i],target_type[j]])){ #start=known,end=known as indicated
                  if(!is.na(direction[j]) & direction[j]>deg_start[source_type[i],target_type[j]] & direction[j]<deg_end[source_type[i],target_type[j]]){
                    #calculate prob of infection and rejection sampling
                    dist_fac <- (max_distance[source_type[i],target_type[j]]-d[j])/(max_distance[source_type[i],target_type[j]]-1)
                    source_fac <- SizeFactor(as.vector(as.matrix(pop[,2])),source_index[i])
                    target_fac <- SizeFactor(as.vector(as.matrix(pop[,2])),sus_ind[j])
                    if(prev_chart_use[source_type[i]]=="Yes"){
                      st  <- status[source_index[i],] #disease status of source unit from day 1 up to current day
                      sum <- 0 #the number of days since source was infected
                      idx <- 0
                      while (st[length(st)-idx]=="L" | st[length(st)-idx]=="B" | st[length(st)-idx]=="C") {
                        sum <- sum + 1
                        idx <- idx + 1
                      }
                      P <-  prev[source_type[i],sum]*dist_fac*source_fac*target_fac*p[source_type[i],target_type[j]]
                    }else{P <-  dist_fac*source_fac*target_fac*p[source_type[i],target_type[j]]}
                    if(!is.na(P) & P!=Inf){
                      r <- runif(1,0,1)
                      if(r < P){ #infected
                        infected_source <- append(infected_source,source_index[i])
                        exposure_target <- append(exposure_target,sus_ind[j])
                        infected_target <- append(infected_target,sus_ind[j])
                        #updated disease status after the delays
                        if(delay_pdf[source_type[i],target_type[j]]=="fixed" | is.na(delay_pdf[source_type[i],target_type[j]])){
                          num_delay <- delay_par[source_type[i],target_type[j],1]
                          delays <- append(delays,num_delay)
                        }else{
                          num_delay <- delays(delay_pdf[source_type[i],target_type[j]],delay_par[source_type[i],target_type[j],])
                          delays <- append(delays,num_delay)
                        }
                        #update the disease status indicated that the unit is no longer considered
                        disease_status[sus_ind[j]] <- "S_Infected"
                      }else{ #not infected
                        infected_source <- append(infected_source,source_index[i])
                        exposure_target <- append(exposure_target,NA)
                        infected_target <- append(infected_target,NA)
                        delays <- append(delays,NA)
                      }
                    }else{
                      infected_source <- append(infected_source,source_index[i])
                      exposure_target <- append(exposure_target,NA)
                      infected_target <- append(infected_target,NA)
                      delays <- append(delays,NA)
                    }
                  }else{ #source-target is not inside the area at risk
                    infected_source <- append(infected_source,source_index[i])
                    exposure_target <- append(exposure_target,NA)
                    infected_target <- append(infected_target,NA)
                    delays <- append(delays,NA)
                  }
                }else{ #By default, start = known, end=360
                  if(!is.na(direction[j]) & direction[j]>deg_start[source_type[i],target_type[j]]){
                    #calculate prob of infection and rejection sampling
                    dist_fac <- (max_distance[source_type[i],target_type[j]]-d[j])/(max_distance[source_type[i],target_type[j]]-1)
                    source_fac <- SizeFactor(as.vector(as.matrix(pop[,2])),source_index[i])
                    target_fac <- SizeFactor(as.vector(as.matrix(pop[,2])),sus_ind[j])
                    if(prev_chart_use[source_type[i]]=="Yes"){
                      st  <- status[source_index[i],] #disease status of source unit from day 1 up to current day
                      sum <- 0 #the number of days since source was infected
                      idx <- 0
                      while (st[length(st)-idx]=="L" | st[length(st)-idx]=="B" | st[length(st)-idx]=="C") {
                        sum <- sum + 1
                        idx <- idx + 1
                      }
                      P <-  prev[source_type[i],sum]*dist_fac*source_fac*target_fac*p[source_type[i],target_type[j]]
                    }else{P <-  dist_fac*source_fac*target_fac*p[source_type[i],target_type[j]]}
                    if(!is.na(P) & P!=Inf){
                      r <- runif(1,0,1)
                      if(r < P){ #infected
                        infected_source <- append(infected_source,source_index[i])
                        exposure_target <- append(exposure_target,sus_ind[j])
                        infected_target <- append(infected_target,sus_ind[j])
                        #updated disease status after the delays
                        if(delay_pdf[source_type[i],target_type[j]]=="fixed" | is.na(delay_pdf[source_type[i],target_type[j]])){
                          num_delay <- delay_par[source_type[i],target_type[j],1]
                          delays <- append(delays,num_delay)
                        }else{
                          num_delay <- delays(delay_pdf[source_type[i],target_type[j]],delay_par[source_type[i],target_type[j],])
                          delays <- append(delays,num_delay)
                        }
                        #update the disease status indicated that the unit is no longer considered
                        disease_status[sus_ind[j]] <- "S_Infected"
                      }else{ #not infected
                        infected_source <- append(infected_source,source_index[i])
                        exposure_target <- append(exposure_target,NA)
                        infected_target <- append(infected_target,NA)
                        delays <- append(delays,NA)
                      }
                    }else{
                      infected_source <- append(infected_source,source_index[i])
                      exposure_target <- append(exposure_target,NA)
                      infected_target <- append(infected_target,NA)
                      delays <- append(delays,NA)
                    }
                  }else{ #source-target is not inside the area at risk
                    infected_source <- append(infected_source,source_index[i])
                    exposure_target <- append(exposure_target,NA)
                    infected_target <- append(infected_target,NA)
                    delays <- append(delays,NA)
                  }
                }
              }else{
                if(!is.na(deg_end[source_type[i],target_type[j]])){ #By default, start = 0, end=known
                  if(!is.na(direction[j]) & direction[j]<deg_end[source_type[i],target_type[j]]){
                    #calculate prob of infection and rejection sampling
                    dist_fac <- (max_distance[source_type[i],target_type[j]]-d[j])/(max_distance[source_type[i],target_type[j]]-1)
                    source_fac <- SizeFactor(as.vector(as.matrix(pop[,2])),source_index[i])
                    target_fac <- SizeFactor(as.vector(as.matrix(pop[,2])),sus_ind[j])
                    if(prev_chart_use[source_type[i]]=="Yes"){
                      st  <- status[source_index[i],] #disease status of source unit from day 1 up to current day
                      sum <- 0 #the number of days since source was infected
                      idx <- 0
                      while (st[length(st)-idx]=="L" | st[length(st)-idx]=="B" | st[length(st)-idx]=="C") {
                        sum <- sum + 1
                        idx <- idx + 1
                      }
                      P <-  prev[source_type[i],sum]*dist_fac*source_fac*target_fac*p[source_type[i],target_type[j]]
                    }else{P <-  dist_fac*source_fac*target_fac*p[source_type[i],target_type[j]]}
                    if(!is.na(P) & P!=Inf){
                      r <- runif(1,0,1)
                      if(r < P){ #infected
                        infected_source <- append(infected_source,source_index[i])
                        exposure_target <- append(exposure_target,sus_ind[j])
                        infected_target <- append(infected_target,sus_ind[j])
                        #updated disease status after the delays
                        if(delay_pdf[source_type[i],target_type[j]]=="fixed" | is.na(delay_pdf[source_type[i],target_type[j]])){
                          num_delay <- delay_par[source_type[i],target_type[j],1]
                          delays <- append(delays,num_delay)
                        }else{
                          num_delay <- delays(delay_pdf[source_type[i],target_type[j]],delay_par[source_type[i],target_type[j],])
                          delays <- append(delays,num_delay)
                        }
                        #update the disease status indicated that the unit is no longer considered
                        disease_status[sus_ind[j]] <- "S_Infected"
                      }else{ #not infected
                        infected_source <- append(infected_source,source_index[i])
                        exposure_target <- append(exposure_target,NA)
                        infected_target <- append(infected_target,NA)
                        delays <- append(delays,NA)
                      }
                    }else{
                      infected_source <- append(infected_source,source_index[i])
                      exposure_target <- append(exposure_target,NA)
                      infected_target <- append(infected_target,NA)
                      delays <- append(delays,NA)
                    }
                  }else{ #source-target is not inside the area at risk
                    infected_source <- append(infected_source,source_index[i])
                    exposure_target <- append(exposure_target,NA)
                    infected_target <- append(infected_target,NA)
                    delays <- append(delays,NA)
                  }
                }else{ #By default, start = 0, end=360 => no restriction on area at risk
                  if(!is.na(direction[j])){
                    #calculate prob of infection and rejection sampling
                    dist_fac <- (max_distance[source_type[i],target_type[j]]-d[j])/(max_distance[source_type[i],target_type[j]]-1)
                    source_fac <- SizeFactor(as.vector(as.matrix(pop[,2])),source_index[i])
                    target_fac <- SizeFactor(as.vector(as.matrix(pop[,2])),sus_ind[j])
                    if(prev_chart_use[source_type[i]]=="Yes"){
                      st  <- status[source_index[i],] #disease status of source unit from day 1 up to current day
                      sum <- 0 #the number of days since source was infected
                      idx <- 0
                      while (st[length(st)-idx]=="L" | st[length(st)-idx]=="B" | st[length(st)-idx]=="C") {
                        sum <- sum + 1
                        idx <- idx + 1
                      }
                      P <-  prev[source_type[i],sum]*dist_fac*source_fac*target_fac*p[source_type[i],target_type[j]]
                    }else{P <-  dist_fac*source_fac*target_fac*p[source_type[i],target_type[j]]}
                    if(!is.na(P) & P!=Inf){
                      r <- runif(1,0,1)
                      if(r < P){ #infected
                        infected_source <- append(infected_source,source_index[i])
                        exposure_target <- append(exposure_target,sus_ind[j])
                        infected_target <- append(infected_target,sus_ind[j])
                        #updated disease status after the delays
                        if(delay_pdf[source_type[i],target_type[j]]=="fixed" | is.na(delay_pdf[source_type[i],target_type[j]])){
                          num_delay <- delay_par[source_type[i],target_type[j],1]
                          delays <- append(delays,num_delay)
                        }else{
                          num_delay <- delays(delay_pdf[source_type[i],target_type[j]],delay_par[source_type[i],target_type[j],])
                          delays <- append(delays,num_delay)
                        }
                        #update the disease status indicated that the unit is no longer considered
                        disease_status[sus_ind[j]] <- "S_Infected"
                      }else{ #not infected
                        infected_source <- append(infected_source,source_index[i])
                        exposure_target <- append(exposure_target,NA)
                        infected_target <- append(infected_target,NA)
                        delays <- append(delays,NA)
                      }
                    }else{
                      infected_source <- append(infected_source,source_index[i])
                      exposure_target <- append(exposure_target,NA)
                      infected_target <- append(infected_target,NA)
                      delays <- append(delays,NA)
                    }
                  }else{ #source-target is not inside the area at risk
                    infected_source <- append(infected_source,source_index[i])
                    exposure_target <- append(exposure_target,NA)
                    infected_target <- append(infected_target,NA)
                    delays <- append(delays,NA)
                  }
                }
              }
            }else{ #outside max distance between units or no distance value
              infected_source <- append(infected_source,source_index[i])
              exposure_target <- append(exposure_target,NA)
              infected_target <- append(infected_target,NA)
              delays <- append(delays,NA)
            }
          }else{ #exponential dropoff
            #check if target j shares the same location with the source and/or unknown direction from source and/or outside the range of area at risk
            if(!is.na(deg_start[source_type[i],target_type[j]])){
              if(!is.na(deg_end[source_type[i],target_type[j]])){ #start=known,end=known as indicated
                if(!is.na(direction[j]) & direction[j]>deg_start[source_type[i],target_type[j]] & direction[j]<deg_end[source_type[i],target_type[j]]){
                  #calculate prob of infection and rejection sampling
                  source_fac <- SizeFactor(as.vector(as.matrix(pop[,2])),source_index[i])
                  target_fac <- SizeFactor(as.vector(as.matrix(pop[,2])),sus_ind[j])
                  if(prev_chart_use[source_type[i]]=="Yes"){
                    st  <- status[source_index[i],] #disease status of source unit from day 1 up to current day
                    sum <- 0 #the number of days since source was infected
                    idx <- 0
                    while (st[length(st)-idx]=="L" | st[length(st)-idx]=="B" | st[length(st)-idx]=="C") {
                      sum <- sum + 1
                      idx <- idx + 1
                    }
                    P <-  prev[source_type[i],sum]*source_fac*target_fac*(p[source_type[i],target_type[j]])^d[j]
                  }else{P <-  source_fac*target_fac*(p[source_type[i],target_type[j]])^d[j]}
                  if(!is.na(P) & P!=Inf){
                    r <- runif(1,0,1)
                    if(r < P){ #infected
                      infected_source <- append(infected_source,source_index[i])
                      exposure_target <- append(exposure_target,sus_ind[j])
                      infected_target <- append(infected_target,sus_ind[j])
                      #updated disease status after the delays
                      if(delay_pdf[source_type[i],target_type[j]]=="fixed" | is.na(delay_pdf[source_type[i],target_type[j]])){
                        num_delay <- delay_par[source_type[i],target_type[j],1]
                        delays <- append(delays,num_delay)
                      }else{
                        num_delay <- delays(delay_pdf[source_type[i],target_type[j]],delay_par[source_type[i],target_type[j],])
                        delays <- append(delays,num_delay)
                      }
                      #update the disease status indicated that the unit is no longer considered
                      disease_status[sus_ind[j]] <- "S_Infected"
                    }else{ #not infected
                      infected_source <- append(infected_source,source_index[i])
                      exposure_target <- append(exposure_target,NA)
                      infected_target <- append(infected_target,NA)
                      delays <- append(delays,NA)
                    }
                  }else{
                    infected_source <- append(infected_source,source_index[i])
                    exposure_target <- append(exposure_target,NA)
                    infected_target <- append(infected_target,NA)
                    delays <- append(delays,NA)
                  }
                }else{ #source-target is not inside the area at risk
                  infected_source <- append(infected_source,source_index[i])
                  exposure_target <- append(exposure_target,NA)
                  infected_target <- append(infected_target,NA)
                  delays <- append(delays,NA)
                }
              }else{ #By default, start = known, end=360
                if(!is.na(direction[j]) & direction[j]>deg_start[source_type[i],target_type[j]]){
                  #calculate prob of infection and rejection sampling
                  source_fac <- SizeFactor(as.vector(as.matrix(pop[,2])),source_index[i])
                  target_fac <- SizeFactor(as.vector(as.matrix(pop[,2])),sus_ind[j])
                  if(prev_chart_use[source_type[i]]=="Yes"){
                    st  <- status[source_index[i],] #disease status of source unit from day 1 up to current day
                    sum <- 0 #the number of days since source was infected
                    idx <- 0
                    while (st[length(st)-idx]=="L" | st[length(st)-idx]=="B" | st[length(st)-idx]=="C") {
                      sum <- sum + 1
                      idx <- idx + 1
                    }
                    P <-  prev[source_type[i],sum]*source_fac*target_fac*(p[source_type[i],target_type[j]])^d[j]
                  }else{P <-  source_fac*target_fac*(p[source_type[i],target_type[j]])^d[j]}
                  if(!is.na(P) & P!=Inf){
                    r <- runif(1,0,1)
                    if(r < P){ #infected
                      infected_source <- append(infected_source,source_index[i])
                      exposure_target <- append(exposure_target,sus_ind[j])
                      infected_target <- append(infected_target,sus_ind[j])
                      #updated disease status after the delays
                      if(delay_pdf[source_type[i],target_type[j]]=="fixed" | is.na(delay_pdf[source_type[i],target_type[j]])){
                        num_delay <- delay_par[source_type[i],target_type[j],1]
                        delays <- append(delays,num_delay)
                      }else{
                        num_delay <- delays(delay_pdf[source_type[i],target_type[j]],delay_par[source_type[i],target_type[j],])
                        delays <- append(delays,num_delay)
                      }
                      #update the disease status indicated that the unit is no longer considered
                      disease_status[sus_ind[j]] <- "S_Infected"
                    }else{ #not infected
                      infected_source <- append(infected_source,source_index[i])
                      exposure_target <- append(exposure_target,NA)
                      infected_target <- append(infected_target,NA)
                      delays <- append(delays,NA)
                    }
                  }else{
                    infected_source <- append(infected_source,source_index[i])
                    exposure_target <- append(exposure_target,NA)
                    infected_target <- append(infected_target,NA)
                    delays <- append(delays,NA)
                  }
                }else{ #source-target is not inside the area at risk
                  infected_source <- append(infected_source,source_index[i])
                  exposure_target <- append(exposure_target,NA)
                  infected_target <- append(infected_target,NA)
                  delays <- append(delays,NA)
                }
              }
            }else{
              if(!is.na(deg_end[source_type[i],target_type[j]])){ #By default, start = 0, end=known
                if(!is.na(direction[j]) & direction[j]<deg_end[source_type[i],target_type[j]]){
                  #calculate prob of infection and rejection sampling
                  source_fac <- SizeFactor(as.vector(as.matrix(pop[,2])),source_index[i])
                  target_fac <- SizeFactor(as.vector(as.matrix(pop[,2])),sus_ind[j])
                  if(prev_chart_use[source_type[i]]=="Yes"){
                    st  <- status[source_index[i],] #disease status of source unit from day 1 up to current day
                    sum <- 0 #the number of days since source was infected
                    idx <- 0
                    while (st[length(st)-idx]=="L" | st[length(st)-idx]=="B" | st[length(st)-idx]=="C") {
                      sum <- sum + 1
                      idx <- idx + 1
                    }
                    P <-  prev[source_type[i],sum]*source_fac*target_fac*(p[source_type[i],target_type[j]])^d[j]
                  }else{P <-  source_fac*target_fac*(p[source_type[i],target_type[j]])^d[j]}
                  if(!is.na(P) & P!=Inf){
                    r <- runif(1,0,1)
                    if(r < P){ #infected
                      infected_source <- append(infected_source,source_index[i])
                      exposure_target <- append(exposure_target,sus_ind[j])
                      infected_target <- append(infected_target,sus_ind[j])
                      #updated disease status after the delays
                      if(delay_pdf[source_type[i],target_type[j]]=="fixed" | is.na(delay_pdf[source_type[i],target_type[j]])){
                        num_delay <- delay_par[source_type[i],target_type[j],1]
                        delays <- append(delays,num_delay)
                      }else{
                        num_delay <- delays(delay_pdf[source_type[i],target_type[j]],delay_par[source_type[i],target_type[j],])
                        delays <- append(delays,num_delay)
                      }
                      #update the disease status indicated that the unit is no longer considered
                      disease_status[sus_ind[j]] <- "S_Infected"
                    }else{ #not infected
                      infected_source <- append(infected_source,source_index[i])
                      exposure_target <- append(exposure_target,NA)
                      infected_target <- append(infected_target,NA)
                      delays <- append(delays,NA)
                    }
                  }else{
                    infected_source <- append(infected_source,source_index[i])
                    exposure_target <- append(exposure_target,NA)
                    infected_target <- append(infected_target,NA)
                    delays <- append(delays,NA)
                  }
                }else{ #source-target is not inside the area at risk
                  infected_source <- append(infected_source,source_index[i])
                  exposure_target <- append(exposure_target,NA)
                  infected_target <- append(infected_target,NA)
                  delays <- append(delays,NA)
                }
              }else{ #By default, start = 0, end=360 => no restriction on area at risk
                if(!is.na(direction[j])){
                  #calculate prob of infection and rejection sampling
                  source_fac <- SizeFactor(as.vector(as.matrix(pop[,2])),source_index[i])
                  target_fac <- SizeFactor(as.vector(as.matrix(pop[,2])),sus_ind[j])
                  if(prev_chart_use[source_type[i]]=="Yes"){
                    st  <- status[source_index[i],] #disease status of source unit from day 1 up to current day
                    sum <- 0 #the number of days since source was infected
                    idx <- 0
                    while (st[length(st)-idx]=="L" | st[length(st)-idx]=="B" | st[length(st)-idx]=="C") {
                      sum <- sum + 1
                      idx <- idx + 1
                    }
                    P <-  prev[source_type[i],sum]*source_fac*target_fac*(p[source_type[i],target_type[j]])^d[j]
                  }else{P <-  source_fac*target_fac*(p[source_type[i],target_type[j]])^d[j]}
                  if(!is.na(P) & P!=Inf){
                    r <- runif(1,0,1)
                    if(r < P){ #infected
                      infected_source <- append(infected_source,source_index[i])
                      exposure_target <- append(exposure_target,sus_ind[j])
                      infected_target <- append(infected_target,sus_ind[j])
                      #updated disease status after the delays
                      if(delay_pdf[source_type[i],target_type[j]]=="fixed" | is.na(delay_pdf[source_type[i],target_type[j]])){
                        num_delay <- delay_par[source_type[i],target_type[j],1]
                        delays <- append(delays,num_delay)
                      }else{
                        num_delay <- delays(delay_pdf[source_type[i],target_type[j]],delay_par[source_type[i],target_type[j],])
                        delays <- append(delays,num_delay)
                      }
                      #update the disease status indicated that the unit is no longer considered
                      disease_status[sus_ind[j]] <- "S_Infected"
                    }else{ #not infected
                      infected_source <- append(infected_source,source_index[i])
                      exposure_target <- append(exposure_target,NA)
                      infected_target <- append(infected_target,NA)
                      delays <- append(delays,NA)
                    }
                  }else{
                    infected_source <- append(infected_source,source_index[i])
                    exposure_target <- append(exposure_target,NA)
                    infected_target <- append(infected_target,NA)
                    delays <- append(delays,NA)
                  }
                }else{ #source-target is not inside the area at risk
                  infected_source <- append(infected_source,source_index[i])
                  exposure_target <- append(exposure_target,NA)
                  infected_target <- append(infected_target,NA)
                  delays <- append(delays,NA)
                }
              }
            }
          }
        } #for each susceptible j
      }else{ #no pair of source unit i and susceptible target
        infected_source <- append(infected_source,source_index[i])
        exposure_target <- append(exposure_target,NA)
        infected_target <- append(infected_target,NA)
        delays <- append(delays,NA)
      }
    } #for each source unit i
    infected_results <- data.frame(infected_source,exposure_target,infected_target,delays)
    infected_results <- infected_results %>% filter(!is.na(exposure_target))
    if(nrow(infected_results)>0){
      source_prod_type <- as.vector(as.matrix(pop[infected_results[,1],1]))
      target_prod_type <- as.vector(as.matrix(pop[infected_results[,2],1]))
      results <- cbind(infected_results,source_prod_type,target_prod_type)
      outputs <- list(Infection = results, Status = disease_status)
    }else{
      source_prod_type <- NA
      target_prod_type <- NA
      infected_source  <- NA
      exposure_target  <- NA
      infected_target  <- NA
      delays  <- NA
      results <- data.frame(infected_source,exposure_target,infected_target,delays,source_prod_type,target_prod_type)
      outputs <- list(Infection = results, Status = curr_state)
    }
    return(outputs) # a list of a data.frame of infection information[[1]], disease status[[2]]
  }
}

