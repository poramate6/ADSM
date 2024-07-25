\name{spread}
\alias{spread}

\title{
Exposures and infections: direct and indirect paths
}
\description{
For each movement between source-to-target production types, spread() determines if susceptible units which is a desired target production type is exposed and infected.
}

\usage{
spread(movement_info = NULL,
pop = NULL,
pdf = NULL,
par = NULL,
prob = NULL,
curr_state = NULL,
delay_pdf = NULL,
delay_par = NULL,
prev_chart_use = NULL,
prev = NULL,
status = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{movement_info}{
A data.frame containing columns of disease spread, source index,
source production type, target production type, and the number of
movements. If no movement specified, put a row of all NAs.
}
  \item{pop}{
Population data.file
}
  \item{pdf}{
A matrix, where each entry indicates the animal movement distance probability
 density function name (character value) from a source production type i
 to target production type j.
}
  \item{par}{
3 or higher-dimensional array containing the parameters for pdf. If specified and a pair of production types has a pdf with 1-dimensional parameter, put all NAs for higher
dimension entries.
}
  \item{prob}{
A matrix, 0 to 1 value, where each entry indicates the probability of
infection given exposure for a source production type i and target
production type j.
}
  \item{curr_state}{
A character value or text indicates current disease states from all units. If
"S_Infected" is specified, the unit was infected but currently in the delay period.
}
  \item{delay_pdf}{
A matrix, where each entry indicates the shipping delay probability
 density function name (character value) from a source production type i
 to target production type j.
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
data.frame population disease status from first day to desired day.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{Infection}{
A data.frame containing source index, exposed index, infected index, shipping delay days, source production type, and target production type. If no exposure and infection, it returns one row of all NAs.
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

spread(...,
curr_state = c("S","S","C","L",...,"S"),...
prev = 0.05,...)
}
