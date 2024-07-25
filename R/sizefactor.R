#Calculate the SizeFactor for unit
#unitsize = pop unit size: total-by-1 vector
#index = unit of interest

SizeFactor <- function(unitsize,index){
  size <- append(unitsize,0)
  d <- density(size)
  xx <- d$x
  dx <- xx[2L] - xx[1L] #binwidth
  yy <- d$y
  C <- sum(yy)*dx #The area under the estimated density curve
  p <- (sum(yy[xx <= unitsize[index]])*dx) / C #scaled it by C for a proper probability estimation
  sf <- p*2
  return(sf)
}
