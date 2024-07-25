\name{airspread}
\alias{airspread}

\title{
Airborne exposure and infection
}
\description{
This function determines if a susceptible unit for a desired target production type is exposed and infected by airborne contact. All pairs of source unit and susceptible unit will be considered based on dropoff type of spread, desired maximum distance of disease spread, and being within area at risk.
}
\usage{
airspread(curr_state = NULL,
max_distance = NULL,
pop = NULL,
dropoff = NULL,
deg_start = NULL,
deg_end = NULL,
p = NULL,
delay_pdf = NULL,
delay_par - NULL,
prev_chart_use = NULL,
prev = NULL,
status = NULL)
}

\arguments{
  \item{curr_state}{
A character value or text indicates current disease states from all units. If
"S_Infected" is specified, the unit was infected but currently in the delay period.
}
  \item{max_distance}{
A matrix, where each entry indicates the maximum distance of disease spread from a source production type i to target production type j (only linear dropoff type).
}
  \item{pop}{
Population data.file
}
  \item{dropoff}{
A matrix, either “exp” or “linear”, where each entry indicates the dropoff type of disease spread from a source production type i to target production type j
}
  \item{deg_start}{
A matrix, in degrees from 0 to 360, where each entry indicates the starting degree of area at risk of exposure between a source production type i and target production type j. If 0 specified, it is a north direction, and 90 is east direction. If put NA, by default, it starts from the north direction or deg = 0.
}
  \item{deg_end}{
A matrix, in degrees from 0 to 360, where each entry indicates the ending degree of area at risk of exposure between a source production type i and target production type j. If put NA, by default, it ends at north direction or deg = 360 or 0.
}
  \item{p}{
A matrix, 0 to 1, where each entry indicates the probability of infection at 1 km from the source given a source production type i and target production type j.
}
  \item{delay_pdf}{
A matrix, where each entry indicates the airborne transport delay probability density function name (character value) from a source production type i to target production type j.
}
  \item{delay_par}{
1 or higher-dimensional array containing the parameters for delay_pdf. If specified and a pair of production types has a pdf with 1-dimensional parameter, put all NAs for higher dimension entries.
}
  \item{prev_chart_use}{
A matrix, where each entry indicates if the prevalence vs. time is used to determine the infection from a source production type i to target production type j.
}
  \item{prev}{
Prevalence in the animal population, 0 to 1.
}
  \item{status}{
data.frame population disease status from first day to a desired day.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{Infection}{
A data.frame containing source index, exposed index, infected index, transport delay days, source production type, and target production type. If no exposure and infection, it returns one row of all NAs.
  }
  \item{Status}{
The updated current status for all units.
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
airspread(curr_state = c("S","S","C","L",...,"S"),
..., prev = 0.001,...)
}
