#create the array of population disease status: unit-day-iteration

create_status <- function(data,t,iteration){
  n <- nrow(data)*t*iteration
  st <- array(rep("S",n),dim = c(nrow(data),t,iteration))
  st[,1,] <- as.matrix(data[,5]) #all entries are character indicates the disease status
  return(st) #array
}


