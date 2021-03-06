% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/magnifiedIV.R
\name{magnifiedIV}
\alias{magnifiedIV}
\title{Magnified Instrumental Variables}
\usage{
magnifiedIV(
  formula,
  data,
  ngroups = 4,
  grouping = "groupCF",
  groupformula,
  est = "group",
  p = 1/4,
  ivsubset = NULL,
  ivna.action = NULL,
  ivweights = NULL,
  ivcontrasts = NULL,
  silent = FALSE,
  ...
)
}
\arguments{
\item{formula}{A formula of the form \code{y ~ x1 + x2 | x1 + z2 + z3} where \code{x1} is an exogenous variable included in the first and second stage, \code{x2} is an endogenous variable being instrumented for, and \code{z2} and \code{z3} are instruments included only in the first stage.}

\item{data}{A data.frame.}

\item{ngroups}{Number of groups to split the data into.}

\item{grouping}{A string variable indicating either the name of a variable in \code{data} to be used as a grouping (or a vector of names if there's more than one excluded instrument), or \code{'groupCF'} or \code{'groupSearch'} to indicate that causal forest or GroupSearch should be used to create groups, respectively.}

\item{groupformula}{If using \code{grouping = 'groupCF'} or \code{grouping = 'groupSearch'}, the formula to pass to those functions as their \code{formula} argument. If there are multiple instruments, instead pass a list of functions, to be used in the same order the excluded instruments are listed in \code{formula}.}

\item{est}{A string variable equal to \code{'group'}, \code{'weight'}, or \code{'both'} to indicate whether the grouped, weighted, or group-and-weighted version of Magnified IV should be used. Note that weighted versions can only be used with one endogenous variable so there is a defined "first stage".}

\item{p}{If using \code{est = 'weight'} or \code{est = 'both'} with \code{grouping = 'groupCF'} or \code{grouping = 'groupSearch'}, weights will be constructed by taking the group effects, generating a first-stage F statistic for each observation as though every observation in the sample had that first-stage effect, and raising that F statistic to the power \code{p}. Set \code{p = 1/4} to get a similar weighting scheme as \code{est = 'group'} does, or \code{p = -1/4} to recover the average treatment effect among compliers (under monotonicity and very good estimates of first-stage treatment effect heterogeneity).}

\item{ivsubset, ivna.action, ivweights, ivcontrasts}{Additional options to be passed to \code{AER::ivreg}, corresponding to \code{subset}, \code{subset}, \code{na.action}, etc.. These are all preceded by 'iv' to avoid overlap with same-named options being passed in \code{...}. Note that \code{ivweights} will be ignored if \code{est = 'weight'} or \code{est = 'both'}. \code{ivweights} with \code{est = 'groups'} can be used to implement weighting schemes other than the power-of-an-F-statistic suggested in Huntington-Klein (2019).}

\item{silent}{Set to TRUE to suppress any messages related to group-creation slow times.}

\item{...}{Additional options to be passed to \code{groupSearch} or \code{groupCF}, as appropriate.}
}
\description{
This command runs the magnified instrumental variables estimator, of either the group-based or weighted variety.
}
\details{
This function will:

1. Use \code{groupSearch()}, \code{groupCF} or a supplied grouping variable, as requested, to find appropriate groupings over which the effects of the instruments vary.

2. If the weighted version of Magnified IV is requested, uses \code{factorPull()} to estimate effects within each group, and calculates sample weights for each observation, which are the first-stage F statistic that would be achieved if all observations in the sample had the same instrument effect as that observation, to the power \code{p}, then runs \code{AER::ivreg()} with those weights.

3. If the grouped version of Magnified IV is requested, runs \code{AER::ivreg()} with those groups included as controls and interacted with the instrument.

You can get more control over the estimation process (say, using a different command than \code{AER::ivreg()}) by making use of the constituent parts to get groups or weights and then setting the interaction or weights yourself.

One instance in which you may want to get more control over the estimation process is if you have multiple endogenous variables. With multiple endogenous variables and multiple instruments, there are lots of different ways to determine groupings. The \code{magnifiedIV} function will only allow one grouping per instrument, and if you want \code{groupSearch} or \code{groupCF} to determine that grouping for you, you must be able to specify a single formula like \code{x1 ~ z1 | w1 + w2} to create that grouping. Can't do \code{x1 + x2 ~ z1 + z2 | w1 + w2}, that doesn't work. So with multiple excluded instruments that are each "for" multiple endogenous variables, you will want to do something by hand.

It is recommended that you read Huntington-Klein (2019) before using this command to understand which versions appear to work well, and what the Super-Local Average Treatment Effect (SLATE) is.
}
\examples{
# This repeats the Stock-Watson textbook cigarette tax IV model
# from the AER::ivreg() help file.
# ... but it's magnified!

# Get data
data(CigarettesSW, package = 'AER')

CigarettesSW$rprice <- CigarettesSW$price/CigarettesSW$cpi
CigarettesSW$rincome <- CigarettesSW$income/CigarettesSW$population/CigarettesSW$cpi
CigarettesSW$tdiff <- (CigarettesSW$taxs - CigarettesSW$tax)/CigarettesSW$cpi

Cig95 <- CigarettesSW[CigarettesSW$year == "1995",]

# Run magnified IV, using causal forest to create effect quantiles
magIV <- magnifiedIV(log(packs) ~ log(rprice) + log(rincome) | log(rincome) + tdiff + I(tax/cpi),
                     data = Cig95,
                     ngroups = 5,
                     grouping = 'groupCF',
                     groupformula = list(log(rprice) ~ tdiff | I(tax/cpi) + population + taxs,
                                         log(rprice) ~ I(tax/cpi) | tdiff + population + taxs))
summary(magIV)
summary(magIV, vcov = sandwich, df = Inf, diagnostics = TRUE)

}
