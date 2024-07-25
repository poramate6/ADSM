#update the initial status for source units before the simulation starts

initial_update <- function(pop,pdf,par,time_frame){
  init_state <- as.vector(as.matrix(pop[,5]))
  updated_status <- as.data.frame(create_status(pop,time_frame,1))
  updated_status[,1] <- init_state
  source_ind <- which(init_state=="L" | init_state=="B" | init_state=="C") #there must be at least one to simulate(by assumption)
  for (i in 1:length(source_ind)) { #for each source unit i
    source_type <- as.integer(pop[source_ind[i],1])
    length_day <- rep(0,4)
    if(init_state[source_ind[i]]=="L"){
      for (j in 1:ncol(pdf)) { #for L,B,C,N
        if(pdf[source_type,j]!=Inf){
          if(pdf[source_type,j]=="fixed"){
            if(is.numeric(par[source_type,j,1])){ #fixed number
              length_day[j] <- par[source_type,j,1]
            }else if(par[source_type,j,1]==Inf){
              length_day[j] <- Inf
            }else{length_day[j] <- 0}
          }else if(is.na(pdf[source_type,j])){
            if(is.numeric(par[source_type,j,1])){ #fixed number
              length_day[j] <- par[source_type,j,1]
            }else if(par[source_type,j,1]==Inf){
              length_day[j] <- Inf
            }else{length_day[j] <- 0}
          }else{ #specify pdf name
            if(par[source_type,j,1]==Inf){
              length_day[j] <- Inf
            }else if(is.na(par[source_type,j,1])){
              length_day[j] <- 0
            }else{ #simulate the duration based on pdf with par
              d <- delays(pdf[source_type,j],par[source_type,j,])
              d <- as.integer(d)
              length_day[j] <- d
            }
          }
        }else{length_day[j] <- Inf}
      }
    }else if(init_state[source_ind[i]]=="B"){
      for (j in 2:ncol(pdf)) { #for L,B,C,N
        if(pdf[source_type,j]!=Inf){
          if(pdf[source_type,j]=="fixed"){
            if(is.numeric(par[source_type,j,1])){ #fixed number
              length_day[j] <- par[source_type,j,1]
            }else if(par[source_type,j,1]==Inf){
              length_day[j] <- Inf
            }else{length_day[j] <- 0}
          }else if(is.na(pdf[source_type,j])){
            if(is.numeric(par[source_type,j,1])){ #fixed number
              length_day[j] <- par[source_type,j,1]
            }else if(par[source_type,j,1]==Inf){
              length_day[j] <- Inf
            }else{length_day[j] <- 0}
          }else{ #specify pdf name
            if(par[source_type,j,1]==Inf){
              length_day[j] <- Inf
            }else if(is.na(par[source_type,j,1])){
              length_day[j] <- 0
            }else{ #simulate the duration based on pdf with par
              d <- delays(pdf[source_type,j],par[source_type,j,])
              d <- as.integer(d)
              length_day[j] <- d
            }
          }
        }else{length_day[j] <- Inf}
      }
    }else{
      for (j in 3:ncol(pdf)) { #for L,B,C,N
        if(pdf[source_type,j]!=Inf){
          if(pdf[source_type,j]=="fixed"){
            if(is.numeric(par[source_type,j,1])){ #fixed number
              length_day[j] <- par[source_type,j,1]
            }else if(par[source_type,j,1]==Inf){
              length_day[j] <- Inf
            }else{length_day[j] <- 0}
          }else if(is.na(pdf[source_type,j])){
            if(is.numeric(par[source_type,j,1])){ #fixed number
              length_day[j] <- par[source_type,j,1]
            }else if(par[source_type,j,1]==Inf){
              length_day[j] <- Inf
            }else{length_day[j] <- 0}
          }else{ #specify pdf name
            if(par[source_type,j,1]==Inf){
              length_day[j] <- Inf
            }else if(is.na(par[source_type,j,1])){
              length_day[j] <- 0
            }else{ #simulate the duration based on pdf with par
              d <- delays(pdf[source_type,j],par[source_type,j,])
              d <- as.integer(d)
              length_day[j] <- d
            }
          }
        }else{length_day[j] <- Inf}
      }
    }
    #generate global status for source unit i
    #update the disease status for unit  i
    day_remain <- time_frame - 1
    state <- c("L","B","C","N")
    unit_state <- c()
    check_inf <- which(length_day==Inf)
    if(is_empty(check_inf)){ #all non-Inf durations
      for (k in 1:length(length_day)) { #for each state: L,B,C,N
        unit_state <- append(unit_state,rep(state[k],length_day[k]))
      } #for each state
      day_left <- day_remain-length(unit_state)
      if(day_left > 0){unit_state <- append(unit_state,rep("S",day_left))
      }else if(day_left < 0){
        unit_state <- unit_state[1:day_remain] #record within the time-frame
      }else{unit_state <- unit_state}
    }else{ #at least one state with Inf
      if(init_state[source_ind[i]]=="L"){num <- 1}else if(init_state[source_ind[i]]=="B"){num <- 2}else{num <- 3} #position of current status
      state_inf <- min(check_inf) #last state can be reached and stays til the end of simulation
      state_inf <- as.integer(state_inf)
      if(state_inf > num){
        for (n in 1:(state_inf-1)) { #for each state up to last possible state
          unit_state <- append(unit_state,rep(state[n],length_day[n]))
        }
        day_left <- day_remain-length(unit_state)
        if(day_left > 0){unit_state <- append(unit_state,rep(state[state_inf],day_left))
        }else if(day_left < 0){
          unit_state <- unit_state[1:day_remain] #record within the time-frame
        }else{unit_state <- unit_state}
      }else{
        unit_state <- append(unit_state,rep(state[num],day_remain))
      }
    }
    #replace the updated states for infected unit i
    updated_status[source_ind[i],2:time_frame] <- unit_state
  }
  return(updated_status)
}

