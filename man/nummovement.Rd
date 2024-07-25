\name{nummovement}
\alias{nummovement}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The number of movements: direct and indirect paths
}
\description{
This function generates the number of contacts or movements from a source unit to a target unit by either direct path or indirect path. It can be only determined when a movement rate for a pair of source-target production types is specified.

Note: Latent unit cannot be a source of infection by indirect contact.
}
\usage{
nummovement(path = "dir",
prod_type = NULL,
curr_state = NULL,
pdf = NULL,
movement_rate = NULL,
L_infect = NULL,
B_infect = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{path}{
      A route that disease spreads from a source unit to a target unit. There are       two options: "dir" and "ind"
  }
  \item{prod_type}{
     Production types for all units in the animal population, use integer value to indicate each production type
  }
  \item{curr_state}{
    A character value or text indicates current disease states from all units. If
    "S_Infected" is specified, the unit was infected but currently in the delay period.
  }
  \item{pdf}{
 A matrix, where each entry indicates the animal movement probability
 density function name (character value) from a source production type i
 to target production type j.
  }
  \item{movement_rate}{
  1 or higher-dimensional array containing a mean number of targets per source unit per day, where higher dimensions corresponding to the number of parameters for pdf. If specified and a pair of production types has a pdf with 1-dimensional parameter, put all NAs for higher dimension entries.
  }
  \item{L_infect}{
  A matrix, where each entry indicates a text either "Yes" or "No" for the latent unit can spread disease from a source production type i to target production type j.
  }
  \item{B_infect}{
  A matrix, where each entry indicates a text either "Yes" or "No" for the subclinical infectious unit can spread disease from a source production type i to target production type j.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{Movements}{
A data.frame containing path of disease spread, source index,
source production type, target production type, and the number of
movements. If no remaining source units in the population, it returns one
row of all NAs.
}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# helps you determine the number of animal-to-animal movements for swine population
nummovement(path = "dir",
prod_type = c(2,5,1,1,...,5),
curr_state = c("S","S","L","S",...,"S"),
...)

# It is also suitable for determining the movements of people among units
nummovement(path = "ind",
prod_type = c(1,2,3,1,...,2),
curr_state = c("C","S","L","N",...,"S"),
...)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
