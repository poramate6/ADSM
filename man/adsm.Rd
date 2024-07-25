\name{adsm}
\alias{adsm}

\title{
Simulate animal disease spread
}
\description{
It simulates animal disease spread for each scenario which is characterized by the input parameters and herd-level population datafile. The simulation is a stochastic, spatial time model. It is a sum of disease transmissions and transitions. At least one source of infected unit (L/B/C) in the population is required, otherwise the simulation will not run
}
\usage{
adsm(datafile = NULL,
ntype = NULL,
L_infect_dir = NULL,
B_infect_dir = NULL,
m_dir = NULL,
m_dir_par = NULL,
d_dir = NULL,
d_dir_par = NULL,
delay_dir = NULL,
delay_dir_par = NULL,
p_dir = NULL,
B_infect_ind = NULL,
m_ind = NULL,
m_ind_par = NULL,
d_ind = NULL,
d_ind_par = NULL,
delay_ind = NULL,
delay_ind_par = NULL,
p_ind = NULL,
max = NULL,
dropoff = NULL,
deg_s = NULL,
deg_e = NULL,
p_air = NULL,
delay_air = NULL,
delay_air_par = NULL,
dur = NULL,
dur_par = NULL,
prev_chart_use = NULL,
prev = NULL,
t = NULL,
iteration = NULL)
}

\arguments{
  \item{datafile}{
Population data.file containing production type, unit size, latitude, longitude, and initial disease status.
}
  \item{ntype}{
Number of population production types, numeric value
}
  \item{L_infect_dir}{
A matrix, “Yes” / “No” (character), where each entry indicates if the
latent units can spread disease by direct contact from a source
production type i to target production type j.
}
  \item{B_infect_dir}{
A matrix, “Yes” / “No” (character), where each entry indicates if the
sub-clinically infectious units can spread disease by direct contact
from a source production type i to target production type j.
}
  \item{m_dir}{

A matrix, where each entry indicates the animal movement by direct
contact probability density function name (character value) for a
source production type i and target production type j.
}
  \item{m_dir_par}{
1 or higher-dimensional array containing a mean number of targets per source unit per day by diect contact, where higher dimensions correspond to the number of parameters for pdf. If specified and a pair of production types has a pdf with 1-dimensional parameter, put all NAs for higher dimension entries.
}
  \item{d_dir}{
A matrix, where each entry indicates the animal movement distance by direct contact probability
 density function name (character value) from a source production type i
 to target production type j.
}
  \item{d_dir_par}{
1 or higher-dimensional array containing the animal movement distance by direct contact pdf parameters, where higher dimensions correspond to the number of parameters for pdf. If specified and a pair of production types has a pdf with 1-dimensional parameter, put all NAs for higher dimension entries.
}
  \item{delay_dir}{
A matrix, where each entry indicates the shipping delay by direct contact probability
 density function name (character value) from a source production type i
 to target production type j.
}
  \item{delay_dir_par}{
1 or higher-dimensional array containing the parameters for delay_dir. If specified and a pair of production types has a pdf with 1-dimensional parameter, put all NAs for higher dimension entries.
}
  \item{p_dir}{
A matrix, 0 to 1 value, where each entry indicates the probability of
infection given exposure by direct contact for a source production type i and target production type j.
}
  \item{B_infect_ind}{
A matrix, “Yes” / “No” (character), where each entry indicates if the
sub-clinically infectious units can spread disease by indirect contact
from a source production type i to target production type j.
}
  \item{m_ind}{
A matrix, where each entry indicates the animal movement by indirect
contact probability density function name (character value) for a
source production type i and target production type j.
}
  \item{m_ind_par}{
1 or higher-dimensional array containing a mean number of targets per source unit per day by indiect contact, where higher dimensions correspond to the number of parameters for pdf. If specified and a pair of production types has a pdf with 1-dimensional parameter, put all NAs for higher dimension entries.
}
  \item{d_ind}{
A matrix, where each entry indicates the animal movement distance by indirect contact probability
 density function name (character value) from a source production type i
 to target production type j.
}
  \item{d_ind_par}{
1 or higher-dimensional array containing the animal movement distance by indirect contact pdf parameters, where higher dimensions correspond to the number of parameters for pdf. If specified and a pair of production types has a pdf with 1-dimensional parameter, put all NAs for higher dimension entries.
}
  \item{delay_ind}{
A matrix, where each entry indicates the shipping delay by indirect contact probability
 density function name (character value) from a source production type i
 to target production type j.
}
  \item{delay_ind_par}{
1 or higher-dimensional array containing the parameters for delay_ind. If specified and a pair of production types has a pdf with 1-dimensional parameter, put all NAs for higher dimension entries.
}
  \item{p_ind}{
A matrix, 0 to 1 value, where each entry indicates the probability of
infection given exposure by indirect contact for a source production type i and target production type j.
}
  \item{max}{
A matrix, where each entry indicates the maximum distance of airborne spread from a source production type i to target production type j (only linear dropoff type).
}
  \item{dropoff}{
A matrix, either “exp” or “linear”, where each entry indicates the dropoff type of airborne spread from a source production type i to target production type j.
}
  \item{deg_s}{
A matrix, in degrees from 0 to 360, where each entry indicates the starting degree of area at risk of exposure between a source production type i and target production type j. If 0 specified, it is a north direction, and 90 is east direction. If put NA, by default, it starts from the north direction or deg = 0.
}
  \item{deg_e}{
A matrix, in degrees from 0 to 360, where each entry indicates the ending degree of area at risk of exposure between a source production type i and target production type j. If put NA, by default, it ends at north direction or deg = 360 or 0.
}
  \item{p_air}{
A matrix, 0 to 1, where each entry indicates the probability of airborne infection at 1 km from the source given a source production type i and target production type j.
}
  \item{delay_air}{
A matrix, where each entry indicates the airborne transport delay probability density function name (character value) from a source production type i to target production type j.
}
  \item{delay_air_par}{
1 or higher-dimensional array containing the parameters for delay_air. If specified and a pair of production types has a pdf with 1-dimensional parameter, put all NAs for higher dimension entries.
}
  \item{dur}{
A matrix, each entry indicates the state j period (duration) probability density
function name (character value) for a production type i.

If j = 1, it is latent.

If j = 2, it is sub-clinically infectious.

If j = 3, it is clinically infectious.

If j = 4, it is recovery.
}
  \item{dur_par}{
1 or higher-dimensional containing the state j period pdf parameter for a production type i. The dimension must correspond to the input parameters. If a probability density function has only 1-dimensional parameter, put NAs for higher dimension entries.
}
  \item{prev_chart_use}{
Whether or not a prevalence chart is used for each production type, "Yes"/"No".
}
  \item{prev}{
A relational function containing a prevalence chart by production types: prevalence vs. time. It is a data.frame with (ntype) rows and (number of days spent in disease states) columns.
}
  \item{t}{
Number of simulated days
}
  \item{iteration}{
Number of iterations
}
}
\details{
Each input parameter for disease contact is defined by a pair of source-target production types.

Available pdf names: Binomial, Geometric, Hypergeometric, NegBin, Poisson, Beta, Chi-squared, Exp , F, Gamma, Uniform, Pert, Triangular

*NegBin: Negative Binomial

*Exp: Exponential

*Pert: BetaPert
}
\value{
Let N be he total number of units.

\item{Status}{
N×t×iteration array of disease status for the entire population and timeframe
}
\item{Expose}{
A data.frame of exposures containing Iteration, Day, Path, Source ID, Source production type, Source Latitude,
Source Lontitude, Recipient ID, Recipient production type, Recipient Latitude, Recipient Lontitude, Recipient unit size.
}
\item{Infection}{
A data.frame of infections containing Iteration, Day, Path, Source ID, Source production type, Source Latitude,
Source Lontitude, Recipient ID, Recipient production type, Recipient Latitude, Recipient Lontitude, Recipient unit size,
durations in latent, sub-clinically infectious, clinically infectious, recovery.
}
\item{Progress}{
A data.frame of disease transition containing Iteration, Day, unit ID, Status, Latitude, Lontitude.
}
\item{Compartments}{
A data.frame of disease compartments containing Iterations, Day, Number of susceptible, latent, subclinical, clinical, recovery units.
}
\item{Variables}{
A data.frame containing Iteration, Day, Path, Recipient production type,

expcU is the number of cumulative exposure units,

expcA is the number of cumulative animal exposures,

expnU is the number of new exposure units,

expnA is the number of new animal exposures,

infcU is the number of cumulative infected units,

infcA is the number of cumulative infected animals,

infnU is the number of new infected units,

infnA is the number of new infected animals.
}
\item{IterationSummary}{
A data.frame containing Iteration,

expnUDir is the number of new exposures by direct contact,

expnADir is the number of new animal exposures by direct contact,

expnUInd is the number of new exposures by indirect contact,

expnAInd is the number of new animal exposures by indirect contact,

expnUAir is the number of new exposures by airborne spread,

expnAAir is the number of new animal exposures by airborne spread,

expnUAll is the number of new exposures from all routes of spread,

expnAAll is the number of new animal exposures from all routes of spread,

infnUDir is the number of new infections caused by direct contact,

infnADir is the number of new animal infections caused by direct contact,

infnUInd is the number of new infections caused by indirect contact,

infnAInd is the number of new animal infections caused by indirect contact,

infnUAir is the number of new infections caused by airborne spread,

infnAAir is the number of new animal infections caused by airborne spread,

infnUAll is the number of new infections from all routes of spread,

infnAAll is the number of new animal infections from all routes of spread,

InfDur is the number of days from first infection until no more disease is present in the population.
}
\item{All}{
A summary statistic for all variables.
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

##Simulate the disease spread in swine population which has four production types
##over 70 days
##repeat 200 iterations

adsm(datafile = swine,
ntype = 4,...,
prev_chart_use = rep("No",4),...,
t = 70,
iteration = 200)

}
