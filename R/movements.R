#All parameters are given for spread between each pair of production types

#Arguments
#path = disease transmission path #text: "dir","ind"
#prod_type = production type (1,...,n) #integer from the pop: data.frame
#curr_state = current disease state for all units #vector of character:S/S_Infected/L/B/C/R
#pdf = name of probability distribution for the number of movement or "fixed" #source-recipient n-by-n matrix: entry = pdf name(character)
#movement_rate = movement rate for a pair of source-recipient production type or some number(in case pdf=="fixed) #n-n-d array
#L_infect: n-n matrix,"Yes" / "No" Can latent unit spread disease?
#B_infect: n-n matrix,"Yes" / "No" Can sub-clinically infectious unit spread disease?

#Return = data.frame: source_index, source_type, target_type, N(number of movements)

nummovement <- function(path,prod_type,curr_state,pdf,movement_rate,L_infect,B_infect){
  if(path=="dir"){ #direct contact
    #search for source units
    source_ind <- which(curr_state=="L" | curr_state=="B" | curr_state=="C")
    if(!is_empty(source_ind)){
      source_type <- prod_type[source_ind] #source production type
      ind <- c()
      target_type <- c()
      for (i in 1:length(source_ind)) {
        if(curr_state[source_ind[i]]=="L"){
          potential_type <- which(!is.na(movement_rate[source_type[i],,1]) & L_infect[source_type[i],]=="Yes")
          if(length(potential_type)!=0){
            ind <- append(ind,rep(source_ind[i],length(potential_type)))
            target_type <- append(target_type,potential_type)
          }
        }else if(curr_state[source_ind[i]]=="B"){
          potential_type <- which(!is.na(movement_rate[source_type[i],,1]) & B_infect[source_type[i],]=="Yes")
          if(length(potential_type)!=0){
            ind <- append(ind,rep(source_ind[i],length(potential_type)))
            target_type <- append(target_type,potential_type)
          }
        }else{ #status = C
          potential_type <- which(!is.na(movement_rate[source_type[i],,1]))
          if(length(potential_type)!=0){
            ind <- append(ind,rep(source_ind[i],length(potential_type)))
            target_type <- append(target_type,potential_type)
          }
        }
      }
      type <- prod_type[ind] #source production types, include all rows
      if(!is_empty(target_type)){
        N <- rep(NA,length(ind))
        for (j in 1:length(ind)) {
          if(all(is.na(movement_rate[type[j],target_type[j],]))){ #No movement from source to target
            N[j] <- 0
          }else{
            N[j] <- sample_movements(pdf[type[j],target_type[j]],movement_rate[type[j],target_type[j],])
          }
        }
        path <- rep("dir",length(ind))
        N_movement <- data.frame(path,ind,type,target_type,N)
      }else{ #no pair of source-target production types
        path <- "dir"
        ind  <- NA
        type <- NA
        target_type <- NA
        N <- NA
        N_movement <- data.frame(path,ind,type,target_type,N)
      }
    }else{ #no source units on this day
      path <- "dir"
      ind  <- NA
      type <- NA
      target_type <- NA
      N <- NA
      N_movement <- data.frame(path,ind,type,target_type,N)
    }
  }else{ #indirect contact
    source_ind <- which(curr_state=="B" | curr_state=="C") #search for indirect source units
    if(!is_empty(source_ind)){
      source_type <- prod_type[source_ind] #source production type
      ind <- c()
      target_type <- c()
      for (i in 1:length(source_ind)) {
        if(curr_state[source_ind[i]]=="B"){
          potential_type <- which(!is.na(movement_rate[source_type[i],,1]) & B_infect[source_type[i],]=="Yes")
          if(length(potential_type)!=0){
            ind <- append(ind,rep(source_ind[i],length(potential_type)))
            target_type <- append(target_type,potential_type)
          }
        }else{ #status = C
          potential_type <- which(!is.na(movement_rate[source_type[i],,1]))
          if(length(potential_type)!=0){
            ind <- append(ind,rep(source_ind[i],length(potential_type)))
            target_type <- append(target_type,potential_type)
          }
        }
      }
      type <- prod_type[ind] #source production types, include all rows
      if(!is_empty(target_type)){
        N <- rep(NA,length(ind))
        for (j in 1:length(ind)) {
          if(all(is.na(movement_rate[type[j],target_type[j],]))){ #No movement from source to target
            N[j] <- 0
          }else{
            N[j] <- sample_movements(pdf[type[j],target_type[j]],movement_rate[type[j],target_type[j],])
          }
        }
        path <- rep("ind",length(ind))
        N_movement <- data.frame(path,ind,type,target_type,N)
      }else{ #no pair of source-target production types
        path <- "ind"
        ind  <- NA
        type <- NA
        target_type <- NA
        N <- NA
        N_movement <- data.frame(path,ind,type,target_type,N)
      }
    }else{ #no source units on this day
      path <- "ind"
      ind  <- NA
      type <- NA
      target_type <- NA
      N <- NA
      N_movement <- data.frame(path,ind,type,target_type,N)
    }
  }
  return(N_movement)
  #return dataframe: path, source_index, source_ProdType,target_ProdType,number of movements or shipments(can be NA or 0)
}
