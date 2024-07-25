#pdf for the number of dalays must be discrete probability distribution with non-negative support

delays <- function(pdf,parameters){
  if(pdf=="Binomial"){
    D <- rbinom(1,parameters[1],parameters[2])
  }else if(pdf=="Geometric"){
    D <- rgeom(1,parameters[1])
  }else if(pdf=="Hypergeometric"){
    D <- rhyper(1,parameters[2],parameters[1]-parameters[2],parameters[3])
  }else if(pdf=="NegBin"){
    D <- rnbinom(1,parameters[1],parameters[2])
  }else if(pdf=="Poisson"){
    D <- rpois(1,parameters[1])
  }else if(pdf=="Triangular"){
    D <- rtriang(1,parameters[1],parameters[2],parameters[3])
  }else if(pdf=="Chi-squared"){
    D <- ceiling(rchisq(1,parameters[1]))
  }else if(pdf=="Exp"){
    D <- ceiling(rexp(1,parameters[1]))
  }else if(pdf=="F"){
    D <- ceiling(rf(1,parameters[1],parameters[2]))
  }else if(pdf=="Gamma"){
    D <- ceiling(rgamma(1,parameters[1],parameters[2]))
  }else if(pdf=="Uniform"){
    D <- ceiling(runif(1,parameters[1],parameters[2]))
  }else if(pdf=="Pert"){
    D <- ceiling(rpert(1,parameters[1],parameters[2],parameters[3]))
  }else if(pdf=="Beta"){
    D <- ceiling(rbeta(1,parameters[1],parameters[2]))
  }else if(pdf=="fixed"){
    D <- ceiling(parameters[1])
  }else{print("Error: undefined pdf for delays")}
  return(D)
}

