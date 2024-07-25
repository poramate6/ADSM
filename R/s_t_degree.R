#Find the degree between source and targets

#Arguments:
#source_lat = source latitude
#source_lon = source longitude
#Lat = list of target latitude
#Lon = list of target longitude

directionfromsource <- function(source_lat,source_lon,Lat,Lon){
  deg <- rep(NA,length(Lat))
  for (i in 1:length(deg)) {
    x <- Lon[i] - source_lon
    y <- Lat[i] - source_lat
    if(x==0){
      if(y>0){
        deg[i] <- 0
      }else if(y<0){
        deg[i] <- 180
      }else{deg[i] <- NA} #source
    }else{ #x!=0
      if(y==0){
        if(x>0){deg[i]<-90}else{deg[i]<-270}
      }else{ #x!=0,y!=0
        if(y/x==Inf){
          if(x>0 & y>0){
            deg[i] <- 0
          }else if(x>0 & y<0){
            deg[i] <- 180
          }else if(x<0 & y<0){
            deg[i] <- 180
          }else if(x<0 & y>0){deg[i]<-360}
        }else{if(x>0){deg[i] <- ((pi/2) - atan(y/x))*(180/pi)}else{deg[i] <- (3*(pi/2) - atan(y/x))*(180/pi)}}
      }
    }
  }
  return(deg) #return direction(degree) of targets from the source: it could be NA
}
