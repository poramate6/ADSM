#direct.indirect disease spread

#Arguments:
#movement_info = movement information: data.frame(path,source_ind,source_type,target_type,N)
#pop = population data file
#par = movement distance parameters: source-recipient (n-n-d) array
#pdf = movement distance distribution: source-recipient (n-by-n) matrix
#prob = probability of infection given exposure: source-recipient (n-by-n) matrix
#delay_pdf = shipping delay: "fixed" or names of distribution: source-recipient (n-by-n) matrix
#delay_par = shipping delay parameters: (source-recipient n-n-d array)
#curr_state = current disease state for all units #vector of character:S/S_Infected/L/B/C/R
#prev_chart_use = a vector of length ntype: "Yes"/"No"
#prev = ntype-to-(total number of days spent in disease states) data.frame: prevalence chart by production types
#status = N-by-sim_day data.frame: disease status from day 1 up to current day

#return: infection information and status

spread <- function(movement_info,pop,pdf,par,prob,curr_state,delay_pdf,delay_par,prev_chart_use,prev,status){
  path <- movement_info[1,1] #disease contact path
  movement_info <- movement_info %>% filter(!is.na(N) & N!=0) #remove zero or NA number of movements(N) from consideration
  movement_info <- movement_info[,2:5] #source_ind,source_type,target_type,N
  if(nrow(movement_info)>0){
    disease_status <- curr_state
    source_unit   <- c()
    exposure_unit <- c()
    infected_unit <- c()
    num_delays <- c()
    for (i in 1:nrow(movement_info)) { #for each source unit i
      for (j in 1:movement_info[i,4]) { #for each movement from source i
        all_index <- 1:nrow(pop)
        non_source_ind   <- all_index[-movement_info[i,1]] #remove source unit from a list of potential recipients
        non_source_ind   <- non_source_ind[!(non_source_ind %in% exposure_unit) & !(non_source_ind %in% infected_unit)]#remove exposed units from previous consideration
        non_source_type  <- as.vector(as.matrix(pop[non_source_ind,1])) #non-source unit production type
        desired_rec_type <- which(non_source_type==movement_info[i,3]) #desired target production type for source i
        if(!is_empty(desired_rec_type)){
          recipient_ind <- non_source_ind[desired_rec_type] #indexes of non-source with desired production type targets
          #sample a reference distance from a source unit i
          if(pdf[movement_info[i,2],movement_info[i,3]]=="fixed" | is.na(pdf[movement_info[i,2],movement_info[i,3]])){
            if(!is.na(par[movement_info[i,2],movement_info[i,3],1])){
              distance <- par[movement_info[i,2],movement_info[i,3],1]
              if(distance==0){ #no infection happens because no movement
                source_unit   <- append(source_unit,movement_info[i,1])
                exposure_unit <- append(exposure_unit,NA)
                infected_unit <- append(infected_unit,NA)
                num_delays <- append(num_delays,NA)
              }else{
                #compute the distance between source and potential recipient
                d <- rep(NA,length(recipient_ind))
                for (k in 1:length(recipient_ind)) {
                  x <- as.numeric((pop[movement_info[i,1],4]-pop[recipient_ind[k],4]) * cos(pop[recipient_ind[k],3]))
                  y <- as.numeric(pop[movement_info[i,1],3] - pop[recipient_ind[k],3])
                  d[k] <- abs(24901/360 * sqrt(x^2 + y^2) - distance)
                }
                rem_d <- which(is.na(d) & d==Inf) #distance can't be calculated or NA or Inf
                if(!is_empty(rem_d)){
                  d <- d[-rem_d]
                  recipient_ind <- recipient_ind[-rem_d] #update recipients with known distance from the source
                }
                smallest_ind <- which(d==min(d)) #index that represents a location of unit that has smallest distance from source i(not unit ID)
                target_ind <- recipient_ind[smallest_ind]
                if(length(target_ind) > 1){ #more than one recipient with the exact same distance
                  tar_ind <- sample(target_ind,1,prob = as.vector(as.matrix(pop[target_ind,2]))) #randomly select one unit w.r.t unit size
                  if(disease_status[tar_ind]=="S"){
                    source_unit   <- append(source_unit,movement_info[i,1])
                    exposure_unit <- append(exposure_unit,tar_ind) #record exposure
                    r <- runif(1,0,1)
                    if(path=="dir"){ #direct contact check if the prevalence chart is used
                      if(prev_chart_use[movement_info[i,2]]=="Yes"){
                        st  <- status[movement_info[i,1],] #disease status of source unit from day 1 up to current day
                        sum <- 0 #the number of days since source was infected
                        idx <- 0
                        while (st[length(st)-idx]=="L" | st[length(st)-idx]=="B" | st[length(st)-idx]=="C") {
                          sum <- sum + 1
                          idx <- idx + 1
                        }
                        P <- prev[movement_info[i,2],sum]
                      }else{P <- prob[movement_info[i,2],movement_info[i,3]]}
                    }else{
                      P <- prob[movement_info[i,2],movement_info[i,3]] #indirect contact
                    }
                    if(r < P & !is.na(P) & P!=Inf){ #potential recipient gets infected
                      infected_unit <- append(infected_unit,tar_ind)
                      if(delay_pdf[movement_info[i,2],movement_info[i,3]]=="fixed" | is.na(delay_pdf[movement_info[i,2],movement_info[i,3]])){
                        if(!is.na(delay_par[movement_info[i,2],movement_info[i,3],1])){
                          dy <- delay_par[movement_info[i,2],movement_info[i,3],1]
                          num_delays <- append(num_delays,dy)
                          disease_status[tar_ind] <- "S_Infected" #update disease status for infection
                        }else{ #by default, delay=0
                          num_delays <- append(num_delays,0)
                          disease_status[tar_ind] <- "S_Infected"
                        }
                      }else{ #known delay pdf
                        if(all(!is.na(delay_par[movement_info[i,2],movement_info[i,3],]))){ #known par
                          dy <- delays(delay_pdf[movement_info[i,2],movement_info[i,3]],delay_par[movement_info[i,2],movement_info[i,3],])
                          num_delays <- append(num_delays,dy)
                          disease_status[tar_ind] <- "S_Infected" #update disease status for infection
                        }else{ #unknown par: by default, delay=0
                          num_delays <- append(num_delays,0)
                          disease_status[tar_ind] <- "S_Infected"
                        }
                      }
                    }else{
                      infected_unit <- append(infected_unit,NA)
                      num_delays <- append(num_delays,NA)
                    }
                  }else{
                    source_unit   <- append(source_unit,movement_info[i,1])
                    exposure_unit <- append(exposure_unit,tar_ind) #record exposure but not infected
                    infected_unit <- append(infected_unit,NA)
                    num_delays <- append(num_delays,NA)
                  }
                }else if(length(target_ind) == 1){ #exactly one potential recipient
                  if(disease_status[target_ind]=="S"){
                    source_unit   <- append(source_unit,movement_info[i,1])
                    exposure_unit <- append(exposure_unit,target_ind) #record exposure
                    r <- runif(1,0,1)
                    if(path=="dir"){ #direct contact check if the prevalence chart is used
                      if(prev_chart_use[movement_info[i,2]]=="Yes"){
                        st  <- status[movement_info[i,1],] #disease status of source unit from day 1 up to current day
                        sum <- 0 #the number of days since source was infected
                        idx <- 0
                        while (st[length(st)-idx]=="L" | st[length(st)-idx]=="B" | st[length(st)-idx]=="C") {
                          sum <- sum + 1
                          idx <- idx + 1
                        }
                        P <- prev[movement_info[i,2],sum]
                      }else{P <- prob[movement_info[i,2],movement_info[i,3]]}
                    }else{
                      P <- prob[movement_info[i,2],movement_info[i,3]] #indirect contact
                    }
                    if(r < P & !is.na(P) & P!=Inf){ #potential recipient gets infected
                      infected_unit <- append(infected_unit,target_ind)
                      if(delay_pdf[movement_info[i,2],movement_info[i,3]]=="fixed" | is.na(delay_pdf[movement_info[i,2],movement_info[i,3]])){
                        if(!is.na(delay_par[movement_info[i,2],movement_info[i,3],1])){
                          dy <- delay_par[movement_info[i,2],movement_info[i,3],1]
                          num_delays <- append(num_delays,dy)
                          disease_status[target_ind] <- "S_Infected" #update disease status for infection
                        }else{ #by default, delay=0
                          num_delays <- append(num_delays,0)
                          disease_status[target_ind] <- "S_Infected"
                        }
                      }else{ #known delay pdf
                        if(all(!is.na(delay_par[movement_info[i,2],movement_info[i,3],]))){ #known par
                          dy <- delays(delay_pdf[movement_info[i,2],movement_info[i,3]],delay_par[movement_info[i,2],movement_info[i,3],])
                          num_delays <- append(num_delays,dy)
                          disease_status[target_ind] <- "S_Infected" #update disease status for infection
                        }else{ #unknown par: by default, delay=0
                          num_delays <- append(num_delays,0)
                          disease_status[target_ind] <- "S_Infected"
                        }
                      }
                    }else{
                      infected_unit <- append(infected_unit,NA)
                      num_delays <- append(num_delays,NA)
                    }
                  }else{
                    source_unit   <- append(source_unit,movement_info[i,1])
                    exposure_unit <- append(exposure_unit,target_ind) #record exposure but not infected
                    infected_unit <- append(infected_unit,NA)
                    num_delays <- append(num_delays,NA)
                  }
                }
              }
            }else{ #unspecified movement distance at risk from the source i
              source_unit   <- append(source_unit,movement_info[i,1])
              exposure_unit <- append(exposure_unit,NA)
              infected_unit <- append(infected_unit,NA)
              num_delays <- append(num_delays,NA)
            }

          }else{ #known movement distance pdf
            if(all(!is.na(par[movement_info[i,2],movement_info[i,3],]))){
              distance <- movement_distance(pdf[movement_info[i,2],movement_info[i,3]],par[movement_info[i,2],movement_info[i,3],])
              if(distance==0){ #no infection happens
                source_unit   <- append(source_unit,movement_info[i,1])
                exposure_unit <- append(exposure_unit,NA)
                infected_unit <- append(infected_unit,NA)
                num_delays <- append(num_delays,NA)
              }else{
                #compute the distance between source and potential recipient
                d <- rep(NA,length(recipient_ind))
                for (k in 1:length(recipient_ind)) {
                  x <- as.numeric((pop[movement_info[i,1],4]-pop[recipient_ind[k],4]) * cos(pop[recipient_ind[k],3]))
                  y <- as.numeric(pop[movement_info[i,1],3] - pop[recipient_ind[k],3])
                  d[k] <- abs(24901/360 * sqrt(x^2 + y^2) - distance)
                }
                rem_d <- which(is.na(d) & d==Inf) #distance can't be calculated or NA or Inf
                if(!is_empty(rem_d)){
                  d <- d[-rem_d]
                  recipient_ind <- recipient_ind[-rem_d] #update recipients with known distance from the source
                }
                smallest_ind <- which(d==min(d)) #index that represents a location of unit that has smallest distance from source i(not unit ID)
                target_ind <- recipient_ind[smallest_ind]
                if(length(target_ind) > 1){ #more than one recipient with the exact same distance
                  tar_ind <- sample(target_ind,1,prob = as.vector(as.matrix(pop[target_ind,2]))) #randomly select one unit w.r.t unit size
                  if(disease_status[tar_ind]=="S"){
                    source_unit   <- append(source_unit,movement_info[i,1])
                    exposure_unit <- append(exposure_unit,tar_ind) #record exposure
                    r <- runif(1,0,1)
                    if(path=="dir"){ #direct contact check if the prevalence chart is used
                      if(prev_chart_use[movement_info[i,2]]=="Yes"){
                        st  <- status[movement_info[i,1],] #disease status of source unit from day 1 up to current day
                        sum <- 0 #the number of days since source was infected
                        idx <- 0
                        while (st[length(st)-idx]=="L" | st[length(st)-idx]=="B" | st[length(st)-idx]=="C") {
                          sum <- sum + 1
                          idx <- idx + 1
                        }
                        P <- prev[movement_info[i,2],sum]
                      }else{P <- prob[movement_info[i,2],movement_info[i,3]]}
                    }else{
                      P <- prob[movement_info[i,2],movement_info[i,3]] #indirect contact
                    }
                    if(r < P & !is.na(P) & P!=Inf){ #potential recipient gets infected
                      infected_unit <- append(infected_unit,tar_ind)
                      if(delay_pdf[movement_info[i,2],movement_info[i,3]]=="fixed" | is.na(delay_pdf[movement_info[i,2],movement_info[i,3]])){
                        if(!is.na(delay_par[movement_info[i,2],movement_info[i,3],1])){
                          dy <- delay_par[movement_info[i,2],movement_info[i,3],1]
                          num_delays <- append(num_delays,dy)
                          disease_status[tar_ind] <- "S_Infected" #update disease status for infection
                        }else{ #by default, delay=0
                          num_delays <- append(num_delays,0)
                          disease_status[tar_ind] <- "S_Infected"
                        }
                      }else{ #known delay pdf
                        if(all(!is.na(delay_par[movement_info[i,2],movement_info[i,3],]))){ #known par
                          dy <- delays(delay_pdf[movement_info[i,2],movement_info[i,3]],delay_par[movement_info[i,2],movement_info[i,3],])
                          num_delays <- append(num_delays,dy)
                          disease_status[tar_ind] <- "S_Infected" #update disease status for infection
                        }else{ #unknown par: by default, delay=0
                          num_delays <- append(num_delays,0)
                          disease_status[tar_ind] <- "S_Infected"
                        }
                      }
                    }else{
                      infected_unit <- append(infected_unit,NA)
                      num_delays <- append(num_delays,NA)
                    }
                  }else{
                    source_unit   <- append(source_unit,movement_info[i,1])
                    exposure_unit <- append(exposure_unit,tar_ind) #record exposure but not infected
                    infected_unit <- append(infected_unit,NA)
                    num_delays <- append(num_delays,NA)
                  }
                }else if(length(target_ind) == 1){ #exactly one potential recipient
                  if(disease_status[target_ind]=="S"){
                    source_unit   <- append(source_unit,movement_info[i,1])
                    exposure_unit <- append(exposure_unit,target_ind) #record exposure
                    r <- runif(1,0,1)
                    if(path=="dir"){ #direct contact check if the prevalence chart is used
                      if(prev_chart_use[movement_info[i,2]]=="Yes"){
                        st  <- status[movement_info[i,1],] #disease status of source unit from day 1 up to current day
                        sum <- 0 #the number of days since source was infected
                        idx <- 0
                        while (st[length(st)-idx]=="L" | st[length(st)-idx]=="B" | st[length(st)-idx]=="C") {
                          sum <- sum + 1
                          idx <- idx + 1
                        }
                        P <- prev[movement_info[i,2],sum]
                      }else{P <- prob[movement_info[i,2],movement_info[i,3]]}
                    }else{
                      P <- prob[movement_info[i,2],movement_info[i,3]] #indirect contact
                    }
                    if(r < P & !is.na(P) & P!=Inf){ #potential recipient gets infected
                      infected_unit <- append(infected_unit,target_ind)
                      if(delay_pdf[movement_info[i,2],movement_info[i,3]]=="fixed" | is.na(delay_pdf[movement_info[i,2],movement_info[i,3]])){
                        if(!is.na(delay_par[movement_info[i,2],movement_info[i,3],1])){
                          dy <- delay_par[movement_info[i,2],movement_info[i,3],1]
                          num_delays <- append(num_delays,dy)
                          disease_status[target_ind] <- "S_Infected" #update disease status for infection
                        }else{ #by default, delay=0
                          num_delays <- append(num_delays,0)
                          disease_status[target_ind] <- "S_Infected"
                        }
                      }else{ #known delay pdf
                        if(all(!is.na(delay_par[movement_info[i,2],movement_info[i,3],]))){ #known par
                          dy <- delays(delay_pdf[movement_info[i,2],movement_info[i,3]],delay_par[movement_info[i,2],movement_info[i,3],])
                          num_delays <- append(num_delays,dy)
                          disease_status[target_ind] <- "S_Infected" #update disease status for infection
                        }else{ #unknown par: by default, delay=0
                          num_delays <- append(num_delays,0)
                          disease_status[target_ind] <- "S_Infected"
                        }
                      }
                    }else{
                      infected_unit <- append(infected_unit,NA)
                      num_delays <- append(num_delays,NA)
                    }
                  }else{
                    source_unit   <- append(source_unit,movement_info[i,1])
                    exposure_unit <- append(exposure_unit,target_ind) #record exposure but not infected
                    infected_unit <- append(infected_unit,NA)
                    num_delays <- append(num_delays,NA)
                  }
                }
              }
            }else{
              source_unit   <- append(source_unit,movement_info[i,1])
              exposure_unit <- append(exposure_unit,NA)
              infected_unit <- append(infected_unit,NA)
              num_delays <- append(num_delays,NA)
            }
          }
        }else{ #no desired production type targets
          source_unit   <- append(source_unit,movement_info[i,1])
          exposure_unit <- append(exposure_unit,NA)
          infected_unit <- append(infected_unit,NA)
          num_delays <- append(num_delays,NA)
        }
      } #iteration for potential recipients for source unit i
    } #iteration for source unit i
    infection <- data.frame(source_unit,exposure_unit,infected_unit,num_delays)
    if(nrow(infection)>0){
      source_prod_type <- as.vector(as.matrix(pop[infection[,1],1]))
      target_prod_type <- as.vector(as.matrix(pop[infection[,2],1]))
      results <- cbind(infection,source_prod_type,target_prod_type)
      outputs <- list(results,disease_status)
    }else{
      source_unit <- NA
      exposure_unit <- NA
      infected_unit <- NA
      num_delays <- NA
      source_prod_type <- NA
      target_prod_type <- NA
      results <- data.frame(source_unit,exposure_unit,infected_unit,num_delays,source_prod_type,target_prod_type)
      outputs <- list(Infection = results,Status = curr_state)
    }
    return(outputs) #return a data frame contains information about infection [[1]], disease status [[2]]
  }else{ #no shipping movements from the sources(N is 0 or NA)
    source_unit <- NA
    exposure_unit <- NA
    infected_unit <- NA
    num_delays <- NA
    source_prod_type <- NA
    target_prod_type <- NA
    results <- data.frame(source_unit,exposure_unit,infected_unit,num_delays,source_prod_type,target_prod_type)
    outputs <- list(Infection = results,Status = curr_state)
    return(outputs)
  }
}
