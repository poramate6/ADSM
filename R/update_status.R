#update disease status for the entire population for the remaining simulated days

#parameters specified for each production type

#pop = population data file: text comma-delimited (*.csv)
#global_status = updated pop status for entire time-frame: N-by-t
#spread = infection information: data.frame infection-by-6[source_ind,expose_ind,target_ind,delay,source_type,target_type])
#pdf = names of pdf for duration of state or "fixed" or NA or Inf(remaining for entire time-frame): #prod_type-by-4(dL,dB,dC,dR)
#par = pdf parameters or number(for "fixed" pdf) or NA or Inf(remaining for entire time-frame): #n-4-d array (dL,dB,dC,dR)
#curr_state = current status after infection at simulated day: N-by-1
#sim_day = simulated day
#time_frame = the entire time-frame

state_transition <- function(pop,global_status,spread,pdf,par,curr_state,sim_day,time_frame){
  if(sim_day==time_frame){ #we are not considering infection in the last day of simulation because it does not change any disease status
    update_state <- global_status
    update_state[,sim_day] <- curr_state #to get rid off more than one contacts occur on the same simulated day
    if(nrow(spread)>0){
      Delay<- spread[,4]
      LDur <- rep(0,nrow(spread))
      BDur <- rep(0,nrow(spread))
      CDur <- rep(0,nrow(spread))
      NDur <- rep(0,nrow(spread))
      Duration <- data.frame(Delay,LDur,BDur,CDur,NDur)
    }else{ #no infection
      Delay <- NA
      LDur  <- 0
      BDur  <- 0
      CDur  <- 0
      NDur  <- 0
      Duration <- data.frame(Delay,LDur,BDur,CDur,NDur)
    }
    #Duration$Delay[is.na(Duration$Delay)] <- 0
    Duration[is.na(Duration)] <- 0
    outputs <- list(Status = update_state, Duration = Duration)
    return(outputs) #return [[1]] N-by-t data.frame [[2]] data.frame duration
  }else{ #at least one simulated day
    update_state <- global_status
    spread_info <- spread
    #remove non-infected units from consideration
    spread_info <- spread_info[!is.na(spread_info[,3]),]
    if(nrow(spread_info)>0){
      Delay <- 0
      LDur  <- 0
      BDur  <- 0
      CDur  <- 0
      NDur  <- 0
      Duration <- data.frame(Delay,LDur,BDur,CDur,NDur)
      for (i in 1:nrow(spread_info)) { #each new infected unit i
        n_state <- rep(0,5) #duration of each state(S_Infected,L,B,C,N)
        if(!is.na(spread_info[i,4])){ #for S_Infected
          n_state[1] <- spread_info[i,4] #number of delays
        }else{ #by default delay=0 and Turn to "L" on the next day
          n_state[1] <- 0
        }
        for (j in 1:ncol(pdf)) { #for L,B,C,N #j+1 is added because [1]="S_Infected"
          if(pdf[spread_info[i,6],j]!=Inf){
            if(pdf[spread_info[i,6],j]=="fixed"){
              if(is.numeric(par[spread_info[i,6],j,1])){ #fixed number
                n_state[j+1] <- par[spread_info[i,6],j,1]
              }else if(par[spread_info[i,6],j,1]==Inf){
                n_state[j+1] <- Inf
              }else{n_state[j+1] <- 0}
            }else if(is.na(pdf[spread_info[i,6],j])){
              if(is.numeric(par[spread_info[i,6],j,1])){ #fixed number
                n_state[j+1] <- par[spread_info[i,6],j,1]
              }else if(par[spread_info[i,6],j,1]==Inf){
                n_state[j+1] <- Inf
              }else{n_state[j+1] <- 0}
            }else{ #specify pdf name
              if(par[spread_info[i,6],j,1]==Inf){
                n_state[j+1] <- Inf
              }else if(is.na(par[spread_info[i,6],j,1])){
                n_state[j+1] <- 0
              }else{ #simulate the duration based on pdf with par
                d <- delays(pdf[spread_info[i,6],j],par[spread_info[i,6],j,])
                d <- as.integer(d)
                n_state[j+1] <- d
              }
            }
          }else{n_state[j+1] <- Inf}
        } #we get all durations for infected unit i
        #update the disease status for unit  i
        day_remain <- time_frame - sim_day #start a day after infection occurs (at least one day remain)
        state <- c("S_Infected","L","B","C","N")
        unit_state <- c()
        check_inf <- which(n_state==Inf)
        if(is_empty(check_inf)){ #all non-Inf durations
          for (k in 1:length(n_state)) { #for each state: S_Infected,L,B,C,N
            unit_state <- append(unit_state,rep(state[k],n_state[k]))
          } #for each state
          day_left <- day_remain-length(unit_state)
          if(day_left > 0){unit_state <- append(unit_state,rep("S",day_left))
          }else if(day_left < 0){
            unit_state <- unit_state[1:day_remain] #record within the time-frame
          }else{unit_state <- unit_state}
        }else{ #at least one state with Inf
          state_inf <- min(check_inf) #last state can be reached and stays til the end of simulation
          state_inf <- as.integer(state_inf)
          if(state_inf==1){ #delays = Inf
            unit_state <- append(unit_state,rep(state[1],day_remain))
          }else{
            for (n in 1:(state_inf-1)) { #for each state up to last possible state
              unit_state <- append(unit_state,rep(state[n],n_state[n]))
            } #for each state
            day_left <- day_remain-length(unit_state)
            if(day_left > 0){unit_state <- append(unit_state,rep(state[state_inf],day_left))
            }else if(day_left < 0){
              unit_state <- unit_state[1:day_remain] #record within the time-frame
            }else{unit_state <- unit_state}
          }
        }
        #replace the updated states for infected unit i
        start_day <- sim_day+1 #convert in the next day after infection
        update_state[as.numeric(spread_info[i,3]),start_day:time_frame] <- unit_state
        Duration <- rbind(Duration,n_state) #record duration for infected unit i
      } #each infected unit i
      Duration <- Duration[-1,] #remove the first row
    }else{ #no infected units in spread_info
      update_state <- update_state
      Delay <- NA
      LDur  <- 0
      BDur  <- 0
      CDur  <- 0
      NDur  <- 0
      Duration <- data.frame(Delay,LDur,BDur,CDur,NDur)
    }
    update_state[,sim_day] <- curr_state #to get rid off more than one contacts occur on the same simulated day
    #Duration$Delay[is.na(Duration$Delay)] <- 0
    Duration[is.na(Duration)] <- 0
    outputs <- list(Status = update_state, Duration = Duration)
    return(outputs) #return [[1]] updated current pop status after infection at simulated time: N-by-t data.frame [[2]] data.frame duration
  } #at least one simulated day
}
