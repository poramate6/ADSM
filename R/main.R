##main function

##All parameters are given for spread between each pair of production types
#N= population size
#n= number of population production size
#d= highest dimension of parameters
#datafile = population data file: text comma-delimited (*.csv): N-by-5
#ntype = the number of population production types

#Direct/Indirect
#L_infect_dir: n-n matrix,"Yes" / "No" Can latent unit spread disease by direct contact?
#B_infect_dir: n-n matrix,"Yes" / "No" Can sub-clinically infectious unit spread disease by direct contact?
#B_infect_ind: n-n matrix,"Yes" / "No" Can sub-clinically infectious unit spread disease by indirect contact?
#m= N-by-N matrix number of movements discrete pdf (name/"fixed"/NA): number of units receiving shipments form the source unit per day
#m_par = N-N-d array movement parameters (par/NA)
#d = N-by-N matrix movement distance pdf (name/"fixed"/NA) in km.
#d_par = N-N-d array distance parameters (par/NA)
#delay = N-by-N matrix shipping delays(days) pdf (name/"fixed"/NA) in days
#delay_par = N-N-d array shipping delay parameters (par/NA)
#p = N-by-N matrix probability of infection given exposure (0-to-1))

#Airborne
#max = N-by-N matrix maximum distance of spread in km (scalar values).
#dropoff = N-by-N matrix dropoffs ("linear"/"exp")
#deg_s = N-by-N matrix starting wind direction, area at risk of exposure (degree 0-to-360)
#deg_e = N-by-N matrix ending wind direction, area at risk of exposure (degree 0-to-360)
#p = N-by-N matrix probability of infection given exposure (0-to-1))
#delay_air = N-by-N matrix shipping delays(days) pdf (name/"fixed"/NA) in days
#delay_air_par = N-N-d array shipping delay parameters (par/NA)

#Duration
#dur = N-by-4 matrix of duration states discrete pdf: L/B/C/R in days
#dur_par = N-4-d array of duration parameters
#prev_chart_use = a vector of length ntype: "Yes"/"No"
#prev = ntype-to-(total number of days spent in disease states) data.frame: prevalence chart by production types

#t = finite time-frames for simulation (days)
#iteration = the number of iterations or reps in simulation (finite positive integer)

adsm <- function(datafile,ntype,
                 L_infect_dir,B_infect_dir,m_dir,m_dir_par,d_dir,d_dir_par,delay_dir,delay_dir_par,p_dir, #direct contact input parameters
                 B_infect_ind,m_ind,m_ind_par,d_ind,d_ind_par,delay_ind,delay_ind_par,p_ind, #indirect contact input parameters
                 max,dropoff,deg_s,deg_e,p_air,delay_air,delay_air_par,    #airborne contact input parameters
                 dur,dur_par,prev_chart_use,prev, #duration of states input parameters
                 t,iteration){ #simulation input parameters
  status <- create_status(datafile,t,iteration) #initial disease status before simulation starts # 3-dimensional array with dim=c(N,t,iteration)
  #initial disease status would be the same for each iteration
  Rep <- c()
  Day <- c()
  Path <- c()
  Source <- c()
  SourceType <- c()
  Infection <- c()
  Type <- c()
  Delay <- c()
  LDur <- c()
  BDur <- c()
  CDur <- c()
  RDur <- c()
  Expose <- data.frame(Rep,Day,Path,Source,SourceType,Infection,Type)
  Results <- data.frame(Rep,Day,Path,Source,SourceType,Infection,Type,Delay,LDur,BDur,CDur,RDur)
  for (i in 1:iteration) {
    itr_status <- initial_update(datafile,dur,dur_par,t) #global status for i iteration: N-by-t data.frame
    for (j in 1:t) { #simulate the spread on each day
      #spread independently for all transmission paths
      #the simulation is the sum of all processes:contact+state-transition: update the disease status before put them back to next time t
      route <- c("dir","ind","air")
      pick_path <- sample(route,3)
      for (k in 1:length(pick_path)) {
        print(paste("iteration ",i,"- day ",j," by ",pick_path[k]))
        if(pick_path[k]=="air"){ #airborne contact
          #determine infection
          sp_air <- airspread(itr_status[,j],max,datafile,dropoff,deg_s,deg_e,p_air,delay_air,delay_air_par,prev_chart_use,prev,itr_status[,1:j])
          #update the disease status after infection
          update_status <- state_transition(datafile,itr_status,sp_air[[1]],dur,dur_par,sp_air[[2]],j,t)
          itr_status <- update_status[[1]]
          #record exposure
          new_exposure <- sp_air[[1]] %>% filter(!is.na(exposure_target))
          if(nrow(new_exposure)>0){
            for (m in 1:nrow(new_exposure)) {
              add <- c(i,j,"Airborne",new_exposure[m,1],new_exposure[m,5],new_exposure[m,2],new_exposure[m,6])
              Expose <- rbind(Expose,add)
            }
          }
          #record new infection and duration
          new_infect <- sp_air[[1]] %>% filter(!is.na(infected_target))
          if(nrow(new_infect)>0){
            Dur <- update_status[[2]]
            for (n in 1:nrow(new_infect)) {
              add <- c(i,j,"Airborne",new_exposure[m,1],new_exposure[m,5],new_infect[n,3],new_infect[n,6],Dur[n,1],Dur[n,2],Dur[n,3],Dur[n,4],Dur[n,5])
              Results <- rbind(Results,add)
            }
          }
        }else if(pick_path[k]=="dir"){ #direct contact
          #sample the number of movements
          N_dir <- nummovement(pick_path[k],as.vector(as.matrix(datafile[,1])),itr_status[,j],
                               m_dir,m_dir_par,L_infect_dir,B_infect_dir)
          #determine infection
          sp_dir <- spread(N_dir,datafile,d_dir,d_dir_par,p_dir,itr_status[,j],delay_dir,delay_dir_par,prev_chart_use,prev,itr_status[,1:j])
          #update the disease status after infection
          update_status <- state_transition(datafile,itr_status,sp_dir[[1]],dur,dur_par,sp_dir[[2]],j,t)
          itr_status <- update_status[[1]]
          #record exposure
          new_exposure <- sp_dir[[1]] %>% filter(!is.na(exposure_unit))
          if(nrow(new_exposure)>0){
            for (m in 1:nrow(new_exposure)) {
              add <- c(i,j,"Direct",new_exposure[m,1],new_exposure[m,5],new_exposure[m,2],new_exposure[m,6])
              Expose <- rbind(Expose,add)
            }
          }
          #record new infection and duration
          new_infect <- sp_dir[[1]] %>% filter(!is.na(infected_unit))
          if(nrow(new_infect)>0){
            Dur <- update_status[[2]]
            for (n in 1:nrow(new_infect)) {
              add <- c(i,j,"Direct",new_exposure[m,1],new_exposure[m,5],new_infect[n,3],new_infect[n,6],Dur[n,1],Dur[n,2],Dur[n,3],Dur[n,4],Dur[n,5])
              Results <- rbind(Results,add)
            }
          }
        }else{ #indirect contact
          #sample the number of movements
          N_ind <- nummovement(pick_path[k],as.vector(as.matrix(datafile[,1])),itr_status[,j],
                               m_ind,m_ind_par,L_infect_dir,B_infect_ind)
          #determine infection
          sp_ind <- spread(N_ind,datafile,d_ind,d_ind_par,p_ind,itr_status[,j],delay_ind,delay_ind_par,prev_chart_use,prev,itr_status[,1:j])
          #update the disease status after infection
          update_status <- state_transition(datafile,itr_status,sp_ind[[1]],dur,dur_par,sp_ind[[2]],j,t)
          itr_status <- update_status[[1]]
          #record exposure
          new_exposure <- sp_ind[[1]] %>% filter(!is.na(exposure_unit))
          if(nrow(new_exposure)>0){
            for (m in 1:nrow(new_exposure)) {
              add <- c(i,j,"Indirect",new_exposure[m,1],new_exposure[m,5],new_exposure[m,2],new_exposure[m,6])
              Expose <- rbind(Expose,add)
            }
          }
          #record new infection and duration
          new_infect <- sp_ind[[1]] %>% filter(!is.na(infected_unit))
          if(nrow(new_infect)>0){
            Dur <- update_status[[2]]
            for (n in 1:nrow(new_infect)) {
              add <- c(i,j,"Indirect",new_exposure[m,1],new_exposure[m,5],new_infect[n,3],new_infect[n,6],Dur[n,1],Dur[n,2],Dur[n,3],Dur[n,4],Dur[n,5])
              Results <- rbind(Results,add)
            }
          }
        }
      } #for each contact route
    } #for each simulated day
    itr_status[itr_status=="S_Infected"] <- "S"
    status[,,i] <- as.matrix(itr_status) #update disease status at iteration i
  } #for each iteration

  Inf_info <- Results
  #Expose component
  if(nrow(Expose)>0){
    colnames(Expose) <- c("Rep","Day","Path","Source","SourceType","Infection","Type")
    Expose$Infection <- as.integer(Expose$Infection)
    Expose.AU <- datafile[Expose$Infection,2]  #S_size
    source.lat <- datafile[Expose$Source,3] #S_Lat
    source.lon <- datafile[Expose$Source,4] #S_Lon
    recipient.lat <- datafile[Expose$Infection,3] #R_Lat
    recipient.lon <- datafile[Expose$Infection,4] #R_Lon
    Expose <- cbind(Expose,Expose.AU,source.lat,source.lon,recipient.lat,recipient.lon)
    #Expose <- Expose %>% select(Rep,Day,Path,Source,SourceType,source.lat,source.lon,Infection,Type,recipient.lat,recipient.lon,Expose.AU)
    colnames(Expose) <- c("Rep","Day","Path","Source","SourceType","Infection","Type","Expose.AU","source.lat","source.lon","recipient.lat","recipient.lon")
    Expose <- Expose[,c("Rep","Day","Path","Source","SourceType","source.lat","source.lon","Infection","Type","recipient.lat","recipient.lon","Expose.AU")]
    colnames(Expose) <- c("Run","Day","Path","S_ID","S_Type","S_Lat","S_Lon","R_ID","R_Type","R_Lat","R_Lon","R_Size")
    Expose$Run  <- as.integer(Expose$Run)
    Expose$Day  <- as.integer(Expose$Day)
    Expose$S_ID <- as.integer(Expose$S_ID)
    Expose$R_ID <- as.integer(Expose$R_ID)
    Expose$R_Size <- as.integer(Expose$R_Size)
  }else{
    Run <- NA
    Day <- NA
    Path <- NA
    S_ID <- NA
    S_Type <- NA
    S_Lat <- NA
    S_Lon <- NA
    R_ID <- NA
    R_Type <- NA
    R_Lat <- NA
    R_Lon <- NA
    R_Size <- NA
    Expose <- data.frame(Run,Day,Path,S_ID,S_Type,S_Lat,S_Lon,R_ID,R_Type,R_Lat,R_Lon,R_Size)
  }

  #Infection component
  if(nrow(Results)>0){
    colnames(Results) <- c("Rep","Day","Path","Source","SourceType","Infection","Type","Delay","LDur","BDur","CDur","NDur")
    Results$Infection <- as.integer(Results$Infection)
    Infection.AU <- datafile[Results$Infection,2]
    s.lat <- datafile[Results$Source,3] #S_Lat
    s.lon <- datafile[Results$Source,4] #S_Lon
    r.lat <- datafile[Results$Infection,3] #R_Lat
    r.lon <- datafile[Results$Infection,4] #R_Lon
    Results <- cbind(Results,Infection.AU,s.lat,s.lon,r.lat,r.lon)
    #Results <- Results %>% select(Rep,Day,Path,Source,SourceType,s.lat,s.lon,Infection,Type,r.lat,r.lon,Infection.AU,Delay,LDur,BDur,CDur,NDur)
    colnames(Results) <- c("Rep","Day","Path","Source","SourceType","Infection","Type","Delay","LDur","BDur","CDur","NDur","Infection.AU","s.lat","s.lon",
                           "r.lat","r.lon")
    Results <- Results[,c("Rep","Day","Path","Source","SourceType","s.lat","s.lon","Infection","Type","r.lat","r.lon","Infection.AU","Delay","LDur","BDur",
                          "CDur","NDur")]
    colnames(Results) <- c("Run","Day","Path","S_ID","S_Type","S_Lat","S_Lon","R_ID","R_Type","R_Lat","R_Lon","R_Size","Delay","LDur","BDur","CDur","NDur")
    Results$Run <- as.integer(Results$Run)
    Results$Day <- as.integer(Results$Day)
    Results$S_ID <- as.integer(Results$S_ID)
    Results$R_ID <- as.integer(Results$R_ID)
    Results$R_Size <- as.integer(Results$R_Size)
    Results$Delay <- as.numeric(Results$Delay)
    Results$LDur <- as.numeric(Results$LDur)
    Results$BDur <- as.numeric(Results$BDur)
    Results$CDur <- as.numeric(Results$CDur)
    Results$NDur <- as.numeric(Results$NDur)
  }else{
    Run <- NA
    Day <- NA
    Path <- NA
    S_ID <- NA
    S_Type <- NA
    S_Lat <- NA
    S_Lon <- NA
    R_ID <- NA
    R_Type <- NA
    R_Lat <- NA
    R_Lon <- NA
    R_Size <- NA
    Delay <- NA
    LDur <- NA
    BDur <- NA
    CDur <- NA
    NDur <- NA
    Results <- data.frame(Run,Day,Path,S_ID,S_Type,S_Lat,S_Lon,R_ID,R_Type,R_Lat,R_Lon,R_Size,Delay,LDur,BDur,CDur,NDur)
  }

  #Progress
  Rep1 <- c()
  Rep2 <- c()
  Day1 <- c()
  Day2 <- c()
  ID  <-  c()
  State <- c()
  Lat <- c()
  Lon <- c()
  n_S <- c()
  n_L <- c()
  n_B <- c()
  n_C <- c()
  n_N <- c()
  for (i in 1:dim(status)[3]) {
    for (j in 1:dim(status)[2]) {
      ID.L <- which(status[,j,i]=="L")
      if(!is_empty(ID.L)){
        Rep1  <- append(Rep1,rep(i,length(ID.L)))
        Day1  <- append(Day1,rep(j,length(ID.L)))
        ID    <- append(ID,ID.L)
        State <- append(State,rep("L",length(ID.L)))
        Lat   <- append(Lat,as.numeric(unlist(datafile[ID.L,3])))
        Lon   <- append(Lon,as.numeric(unlist(datafile[ID.L,4])))
      }
      ID.B <- which(status[,j,i]=="B")
      if(!is_empty(ID.B)){
        Rep1  <- append(Rep1,rep(i,length(ID.B)))
        Day1  <- append(Day1,rep(j,length(ID.B)))
        ID    <- append(ID,ID.B)
        State <- append(State,rep("B",length(ID.B)))
        Lat   <- append(Lat,as.numeric(unlist(datafile[ID.B,3])))
        Lon   <- append(Lon,as.numeric(unlist(datafile[ID.B,4])))
      }
      ID.C <- which(status[,j,i]=="C")
      if(!is_empty(ID.C)){
        Rep1  <- append(Rep1,rep(i,length(ID.C)))
        Day1  <- append(Day1,rep(j,length(ID.C)))
        ID    <- append(ID,ID.C)
        State <- append(State,rep("C",length(ID.C)))
        Lat   <- append(Lat,as.numeric(unlist(datafile[ID.C,3])))
        Lon   <- append(Lon,as.numeric(unlist(datafile[ID.C,4])))
      }
      ID.N <- which(status[,j,i]=="N")
      if(!is_empty(ID.N)){
        Rep1  <- append(Rep1,rep(i,length(ID.N)))
        Day1  <- append(Day1,rep(j,length(ID.N)))
        ID    <- append(ID,ID.N)
        State <- append(State,rep("N",length(ID.N)))
        Lat   <- append(Lat,as.numeric(unlist(datafile[ID.N,3])))
        Lon   <- append(Lon,as.numeric(unlist(datafile[ID.N,4])))
      }
      Rep2 <- append(Rep2,i)
      Day2 <- append(Day2,j)
      n_S  <- append(n_S,length(status[,j,i][status[,j,i]=="S"]))
      n_L  <- append(n_L,length(status[,j,i][status[,j,i]=="L"]))
      n_B  <- append(n_B,length(status[,j,i][status[,j,i]=="B"]))
      n_C  <- append(n_C,length(status[,j,i][status[,j,i]=="C"]))
      n_N  <- append(n_N,length(status[,j,i][status[,j,i]=="N"]))
    }
  }
  Prog <- data.frame(Rep1,Day1,ID,State,Lat,Lon)
  colnames(Prog) <- c("Run","Day","ID","Status","Lat","Lon")

  #Compartments
  Comp <- data.frame(Rep2,Day2,n_S,n_L,n_B,n_C,n_N)
  colnames(Comp) <- c("Run","Day","n_Susceptible","n_Latent","n_Subclinical","n_Clinical","n_Recovery")

  #Variables
  Rep3 <- c()
  Day  <- c()
  Path <- c()
  R_Type <- c()
  ecU <- c()
  ecA <- c()
  enU <- c()
  enA <- c()
  icU <- c()
  icA <- c()
  inU <- c()
  inA <- c()
  for (i in 1:iteration) { #Run
    path <- c("Direct","Indirect","Airborne","All")
    for (j in 1:4) { #path + 1
      type <- c(1:ntype,"All")
      for (k in 1:(ntype+1)) { #R_Type+1
        n_ecU <- c(0)
        n_ecA <- c(0)
        n_enU <- c(0)
        n_enA <- c(0)
        n_icU <- c(0)
        n_icA <- c(0)
        n_inU <- c(0)
        n_inA <- c(0)
        for (m in 1:t) { #Day
          Rep3 <- append(Rep3,i)
          Path <- append(Path,path[j])
          R_Type <- append(R_Type,type[k])
          Day  <- append(Day,m)
          if(path[j]!="All"){
            if(type[k]!="All"){ #not cumulative path/type
              ep  <- Expose %>% filter(Run==i,Day==m,Path==path[j],R_Type==type[k])
              ife <- Results %>% filter(Run==i,Day==m,Path==path[j],R_Type==type[k])
            }else{ #not cumulative path + all types
              ep  <- Expose %>% filter(Run==i,Day==m,Path==path[j])
              ife <- Results %>% filter(Run==i,Day==m,Path==path[j])
            }
          }else{
            if(type[k]!="All"){ #all paths + not cumulative type
              ep  <- Expose %>% filter(Run==i,Day==m,R_Type==type[k])
              ife <- Results %>% filter(Run==i,Day==m,R_Type==type[k])
            }else{ #all paths/types
              ep  <- Expose %>% filter(Run==i,Day==m)
              ife <- Results %>% filter(Run==i,Day==m)
            }
          }
          c <- m+1
          n_ecU <- append(n_ecU,n_ecU[c-1] + nrow(ep))
          n_ecA <- append(n_ecA,n_ecA[c-1] + sum(ep$R_Size))
          n_enU <- append(n_enU,nrow(ep))
          n_enA <- append(n_enA,sum(ep$R_Size))
          n_icU <- append(n_icU,n_icU[c-1] + nrow(ife))
          n_icA <- append(n_icA,n_icA[c-1] + sum(ife$R_Size))
          n_inU <- append(n_inU,nrow(ife))
          n_inA <- append(n_inA,sum(ife$R_Size))
          ecU <- append(ecU,n_ecU[c])
          ecA <- append(ecA,n_ecA[c])
          enU <- append(enU,n_enU[c])
          enA <- append(enA,n_enA[c])
          icU <- append(icU,n_icU[c])
          icA <- append(icA,n_icA[c])
          inU <- append(inU,n_inU[c])
          inA <- append(inA,n_inA[c])
        }
      }
    }
  }
  V <- data.frame(Rep3,Path,R_Type,Day,ecU,ecA,enU,enA,icU,icA,inU,inA)
  colnames(V) <- c("Run","Path","R_Type","Day","expcU","expcA","expnU","expnA","infcU","infcA","infnU","infnA")
  V$Run <- as.integer(V$Run)
  V$Day <- as.integer(V$Day)


  #IterationSummary
  Rep <- c()
  expnUDir <- c()
  expnADir <- c()
  expnUInd <- c()
  expnAInd <- c()
  expnUAir <- c()
  expnAAir <- c()
  expnUAll <- c()
  expnAAll <- c()
  infnUDir <- c()
  infnADir <- c()
  infnUInd <- c()
  infnAInd <- c()
  infnUAir <- c()
  infnAAir <- c()
  infnUAll <- c()
  infnAAll <- c()
  summary  <- data.frame(Rep,expnUDir,expnADir,expnUInd,expnAInd,expnUAir,expnAAir,expnUAll,expnAAll,
                         infnUDir,infnADir,infnUInd,infnAInd,infnUAir,infnAAir,infnUAll,infnAAll)
  if(nrow(Expose)>0){
    for (i in 1:iteration) {
      exp <- Expose %>% filter(Run==i)
      if(nrow(exp)>0){
        expUdir <- nrow(exp %>% filter(Path=="Direct"))
        expAdir <- exp %>% filter(Path=="Direct")
        expUind <- nrow(exp %>% filter(Path=="Indirect"))
        expAind <- exp %>% filter(Path=="Indirect")
        expUair <- nrow(exp %>% filter(Path=="Airborne"))
        expAair <- exp %>% filter(Path=="Airborne")
        inf <- Results %>% filter(Run==i)
        if(nrow(inf)>0){
          infUdir <- nrow(inf %>% filter(Path=="Direct"))
          infAdir <- inf %>% filter(Path=="Direct")
          infUind <- nrow(inf %>% filter(Path=="Indirect"))
          infAind <- inf %>% filter(Path=="Indirect")
          infUair <- nrow(inf %>% filter(Path=="Airborne"))
          infAair <- inf %>% filter(Path=="Airborne")
          summary <- rbind(summary,c(i,expUdir,sum(expAdir$R_Size),expUind,sum(expAind$R_Size),expUair,sum(expAair$R_Size),nrow(exp),sum(exp$R_Size),
                                     infUdir,sum(infAdir$R_Size),infUind,sum(infAind$R_Size),infUair,sum(infAair$R_Size),nrow(inf),sum(inf$R_Size)))
        }else{ #exposed but no infection
          summary <- rbind(summary,c(i,expUdir,sum(expAdir$R_Size),expUind,sum(expAind$R_Size),expUair,sum(expAair$R_Size),nrow(exp),sum(exp$R_Size),rep(0,8)))
        }
      }else{ #no exposure and infection at Run=i
        summary <- rbind(summary,c(i,rep(0,16)))
      }
    }
  }else{ #no exposure and infection for all iterations
    for (i in 1:iteration) {
      summary <- rbind(summary,c(i,rep(0,16)))
    }
  }
  InfDur <- c()
  source_idx <- which(datafile[,5]=="L" | datafile[,5]=="B" | datafile[,5]=="C")
  s <- c("S","N")
  for (j in 1:iteration) {
    sum <- 0
    for (k in 1:t) {
      if(length(status[-source_idx,k,j])>0){
        if(all(status[-source_idx,k,j] %in% s)){ #true if no disease in population
          sum <- sum
        }else{sum <- sum + 1}
      }else{ #all units are source
        sum <- sum
      }
    }
    InfDur <- append(InfDur,sum)
  }
  summary <- cbind(summary,InfDur)
  colnames(summary) <- c("Run","expnUDir","expnADir","expnUInd","expnAInd","expnUAir","expnAAir","expnUAll","expnAAll",
                         "infnUDir","infnADir","infnUInd","infnAInd","infnUAir","infnAAir","infnUAll","infnAAll","InfDur")

  #All
  variable <- c("expnUDir","expnADir","expnUInd","expnAInd","expnUAir","expnAAir","expnUAll","expnAAll",
                "infnUDir","infnADir","infnUInd","infnAInd","infnUAir","infnAAir","infnUAll","infnAAll",
                "Delay","LDur","BDur","CDur","NDur","InfDur")
  Mean   <- c()
  SD     <- c()
  Min    <- c()
  LQ     <- c()
  Med    <- c()
  UQ     <- c()
  Max    <- c()
  for (i in 2:17) {
    Mean <- append(Mean,mean(summary[,i]))
    SD <- append(SD,sd(summary[,i]))
    Min <- append(Min,min(summary[,i]))
    LQ <- append(LQ,as.numeric(quantile(summary[,i],probs = c(0.25))))
    Med <- append(Med,median(summary[,i]))
    UQ <- append(UQ,as.numeric(quantile(summary[,i],probs = c(0.75))))
    Max <- append(Max,max(summary[,i]))
  }
  if(nrow(Inf_info)>0){
    Mean <- append(Mean,c(mean(Results$Delay),mean(Results$LDur),mean(Results$BDur),mean(Results$CDur),mean(Results$NDur),mean(summary[,18])))
    SD   <- append(SD,c(sd(Results$Delay),sd(Results$LDur),sd(Results$BDur),sd(Results$CDur),sd(Results$NDur),sd(summary[,18])))
    Min  <- append(Min,c(min(Results$Delay),min(Results$LDur),min(Results$BDur),min(Results$CDur),min(Results$NDur),min(summary[,18])))
    LQ   <- append(LQ,c(as.numeric(quantile(Results$Delay,probs = c(0.25))),as.numeric(quantile(Results$LDur,probs = c(0.25))),
                        as.numeric(quantile(Results$BDur,probs = c(0.25))),as.numeric(quantile(Results$CDur,probs = c(0.25))),
                        as.numeric(quantile(Results$NDur,probs = c(0.25))),as.numeric(quantile(summary[,18],probs = c(0.25)))))
    Med  <- append(Med,c(median(Results$Delay),median(Results$LDur),median(Results$BDur),median(Results$CDur),median(Results$NDur),median(summary[,18])))
    UQ   <- append(UQ,c(as.numeric(quantile(Results$Delay,probs = c(0.75))),as.numeric(quantile(Results$LDur,probs = c(0.75))),
                        as.numeric(quantile(Results$BDur,probs = c(0.75))),as.numeric(quantile(Results$CDur,probs = c(0.75))),
                        as.numeric(quantile(Results$NDur,probs = c(0.75))),as.numeric(quantile(summary[,18],probs = c(0.75)))))
    Max  <- append(Max,c(max(Results$Delay),max(Results$LDur),max(Results$BDur),max(Results$CDur),max(Results$NDur),max(summary[,18])))
  }else{
    Mean <- append(Mean,rep(0,6))
    SD <- append(SD,rep(0,6))
    Min <- append(Min,rep(0,6))
    LQ <- append(LQ,rep(0,6))
    Med <- append(Med,rep(0,6))
    UQ <- append(UQ,rep(0,6))
    Max <- append(Max,rep(0,6))
  }
  all  <- data.frame(variable,Mean,SD,Min,LQ,Med,UQ,Max)

  outputs <- list(Status = status, Expose = Expose, Infection = Results, Progress = Prog, Compartments = Comp,
                  Variables = V, IterationSummary = summary, All = all)
  return(outputs)
  # adsm$Status             N-t-iteration array: disease status for the entire population and time-frame
  # adsm$Expose             data.frame exposures: Run Day Path S_ID S_Type S_Lat S_Lon R_ID R_Type R_Lat R_Lon R_Size
  # adsm$Infection          data.frame: Run Day Path S_ID S_Type S_Lat S_Lon R_ID R_Type R_Lat R_Lon R_Size
  #Delay: Number of days for all infected units in Delay
  #LDur: Number of days for all infected units in latent state (exclude source units)
  #BDur: Number of days for all infected units in sub-clinically infectious state
  #CDur: Number of days for all infected units in clinically infectious state
  #NDur: Number of days for all infected units in recovery state before returning to susceptible
  # adsm$Progress           #Run Day ID Status Lat Lon
  # adsm$Compartments       #Run Day n_Susceptible n_Latent n_Subclinical n_Clinical n_Recovery
  # adsm$Variables          #Run Day Path R_Type expcU expcA expnU expnA infcU infcA infnU infnA
  # adsm$IterationSummary   itr-10 data.frame
  #iteration
  #expnUDir: number of new exposures by direct contact
  #expnADir: number of new animal exposures by direct contact
  #expnUInd: number of exposures by indirect contact
  #expnAInd: number of animal exposures by indirect contact
  #expnUAir: number of exposures by airborne spread
  #expnAAir: number of animal exposures by airborne spread
  #expnUAll: number of exposures, total for all routes of spread
  #expnAAll: number of animal exposures, total for all routes of spread
  #infnUDir: number of new infections caused by direct contact
  #infnADir: number of new animal infections caused by direct contact
  #infnUInd: number of new infections caused by indirect contact
  #infnAInd: number of new animal infections caused by indirect contact
  #infnUAir: number of new infections caused by airborne spread
  #infnAAir: number of new animal infections caused by airborne spread
  #infnUAll: number of new infections, total for all routes of spread
  #infnAAll: number of new animal infections, total for all routes of spread
  #InfDur: Number of days from first infection until no more disease is present in the population
  # adsm$All                #variable,Mean,SD,Min,LQ,Med,UQ,Max : summary statistics
}
