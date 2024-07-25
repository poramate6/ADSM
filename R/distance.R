#pdf for movement distance must have non-negative real-valued support: continuous pdf

movement_distance <- function(pdf,parameters){
  if(pdf=="Beta"){
    D <- rbeta(1,parameters[1],parameters[2])
  }else if(pdf=="Chi-squared"){
    D <- rchisq(1,parameters[1])
  }else if(pdf=="Exp"){
    D <- rexp(1,parameters[1])
  }else if(pdf=="F"){
    D <- rf(1,parameters[1],parameters[2])
  }else if(pdf=="Gamma"){
    D <- rgamma(1,parameters[1],parameters[2])
  }else if(pdf=="Uniform"){
    D <- runif(1,parameters[1],parameters[2])
  }else if(pdf=="Pert"){
    D <- rpert(1,parameters[1],parameters[2],parameters[3])
  }else if(pdf=="Triangular"){
    D <- rtriang(1,parameters[1],parameters[2],parameters[3])
  }else if(pdf=="Poisson"){
    D <- rpois(1,parameters[1])
  }else if(pdf=="Geometric"){
    D <- rgeom(1,parameters[1])
  }else if(pdf=="fixed"){
    D <- parameters[1]
  }else{print("Error: undefined pdf for the shipment distance from the source")}
  return(D)
}
