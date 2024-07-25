\name{state_transition}
\alias{state_transition}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Disease transition
}
\description{
This function updates the disease status for the entire population starting from the desired day up to the end of specified timeframe. It is a sum of stochastic processes from the infections by direct contact, indirect contact, and airborne contact.
}
\usage{
state_transition(pop = NULL,
global_status = NULL,
spread = NULL,
pdf = NULL,
par = NULL,
curr_state = NULL,
sim_day = NULL,
time_frame = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pop}{
Population data.frame
}
  \item{global_status}{
A matrix of current population disease status, the columns is
the timeframe (day), and the row is unit.
}
  \item{spread}{
Data.frame containing source index, exposed index, infected index,
shipping/transport delay days, source production type, target production type.
}
  \item{pdf}{
A matrix, each entry indicates the state j period probability density
function name (character value) for a production type i.

If j = 1, it is latent.

If j = 2, it is sub-clinically infectious.

If j = 3, it is clinically infectious.

If j = 4, it is recovery.
}
  \item{par}{
1 or higher-dimensional containing the state j period probability density
function parameter for a production type i. The dimension must correspond to the input parameters. If a probability density function has only 1-dimensional parameter, put NAs for higher dimension entries.
}
  \item{curr_state}{
A character value or text indicates current disease states from all units. If
"S_Infected" is specified, the unit was infected but currently in the delay period.
}
  \item{sim_day}{
Positive integer value that indicates the starting day to update the disease
status.
}
  \item{time_frame}{
Positive integer value that indicates the number of simulated days.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{Status}{
A data.frame containing the updated disease status after infection at sim_day and progression through time.
  }
  \item{Duration}{
A data.frame containing the number of days in shipping/transport delay, latent, sub-clinically infectious, clinical infectious, and recovery.
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

#update disease status from day 3 up to day 30

state_transition(pop = NULL,...,
curr_state = c("S","S","C","L",...,"S"),
sim_day = 3,
time_frame = 30)
}
