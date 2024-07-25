#pdf for number of movements must be discrete probability distribution with non-negative integer support

#Arguments:
#pdf = name of probability distribution
#parameters = vector of parameters for given pdf

sample_movements <- function(pdf,parameters){
  if(pdf=="Binomial"){
    N <- rbinom(1,parameters[1],parameters[2])
  }else if(pdf=="Geometric"){
    N <- rgeom(1,parameters[1])
  }else if(pdf=="Hypergeometric"){
    N <- rhyper(1,parameters[2],parameters[1]-parameters[2],parameters[3])
  }else if(pdf=="NegBin"){
    N <- rnbinom(1,parameters[1],parameters[2])
  }else if(pdf=="Poisson"){
    N <- rpois(1,parameters[1])
  }else if(pdf=="fixed"){
    N <- ceiling(parameters[1])
  }else if(pdf=="Beta"){
    N <- ceiling(rbeta(1,parameters[1],parameters[2]))
  }else if(pdf=="Chi-squared"){
    N <- ceiling(rchisq(1,parameters[1]))
  }else if(pdf=="Exp"){
    N <- ceiling(rexp(1,parameters[1]))
  }else if(pdf=="F"){
    N <- ceiling(rf(1,parameters[1],parameters[2]))
  }else if(pdf=="Gamma"){
    N <- ceiling(rgamma(1,parameters[1],parameters[2]))
  }else if(pdf=="Uniform"){
    N <- ceiling(runif(1,parameters[1],parameters[2]))
  }else if(pdf=="Pert"){
    N <- ceiling(rpert(1,parameters[1],parameters[2],parameters[3]))
  }else if(pdf=="Triangular"){
    N <- ceiling(rtriang(1,parameters[1],parameters[2],parameters[3]))
  }else{print("Error: undefined pdf for the number of animal shipment")}
  return(N)
}
