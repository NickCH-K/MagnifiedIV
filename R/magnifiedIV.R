#' Magnified Instrumental Variables
#'
#' This command runs the magnified instrumental variables estimator, of either the group-based or weighted variety.
#'
#' This function will:
#'
#' 1. Use \code{groupSearch()}, \code{groupCF} or a supplied grouping variable, as requested, to find appropriate groupings over which the effects of the instruments vary.
#'
#' 2. If the weighted version of Magnified IV is requested, uses \code{factorPull()} to estimate effects within each group, and calculates sample weights for each observation, which are the first-stage F statistic that would be achieved if all observations in the sample had the same instrument effect as that observation, to the power \code{p}, then runs \code{AER::ivreg()} with those weights.
#'
#' 3. If the grouped version of Magnified IV is requested, runs \code{AER::ivreg()} with those groups included as controls and interacted with the instrument.
#'
#' You can get more control over the estimation process (say, using a different command than \code{AER::ivreg()}) by making use of the constituent parts to get groups or weights and then setting the interaction or weights yourself.
#'
#' One instance in which you may want to get more control over the estimation process is if you have multiple endogenous variables. With multiple endogenous variables and multiple instruments, there are lots of different ways to determine groupings. The \code{magnifiedIV} function will only allow one grouping per instrument, and if you want \code{groupSearch} or \code{groupCF} to determine that grouping for you, you must be able to specify a single formula like \code{x1 ~ z1 | w1 + w2} to create that grouping. Can't do \code{x1 + x2 ~ z1 + z2 | w1 + w2}, that doesn't work. So with multiple excluded instruments that are each "for" multiple endogenous variables, you will want to do something by hand.
#'
#' It is recommended that you read Huntington-Klein (2019) before using this command to understand which versions appear to work well, and what the Super-Local Average Treatment Effect (SLATE) is.
#'
#' @param formula A formula of the form \code{y ~ x1 + x2 | x1 + z2 + z3} where \code{x1} is an exogenous variable included in the first and second stage, \code{x2} is an endogenous variable being instrumented for, and \code{z2} and \code{z3} are instruments included only in the first stage.
#' @param data A data.frame.
#' @param ngroups Number of groups to split the data into.
#' @param grouping A string variable indicating either the name of a variable in \code{data} to be used as a grouping (or a vector of names if there's more than one excluded instrument), or \code{'groupCF'} or \code{'groupSearch'} to indicate that causal forest or GroupSearch should be used to create groups, respectively.
#' @param groupformula If using \code{grouping = 'groupCF'} or \code{grouping = 'groupSearch'}, the formula to pass to those functions as their \code{formula} argument. If there are multiple instruments, instead pass a list of functions, to be used in the same order the excluded instruments are listed in \code{formula}. Note that this does mean that each IV must have a single endogenous variable it is "for" when creating groupings. If you don't want that, you'll have to create groups by some other method.
#' @param est A string variable equal to \code{'group'}, \code{'weight'}, or \code{'both'} to indicate whether the grouped, weighted, or group-and-weighted version of Magnified IV should be used. Note that weighted versions can only be used with one endogenous variable so there is a defined "first stage".
#' @param p If using \code{est = 'weight'} or \code{est = 'both'} with \code{grouping = 'groupCF'} or \code{grouping = 'groupSearch'}, weights will be constructed by taking the group effects, generating a first-stage F statistic for each observation as though every observation in the sample had that first-stage effect, and raising that F statistic to the power \code{p}. Set \code{p = 1/4} to get a similar weighting scheme as \code{est = 'group'} does, or \code{p = -1/4} to recover the average treatment effect among compliers (under monotonicity and very good estimates of first-stage treatment effect heterogeneity).
#' @param ivsubset,ivna.action,ivweights,ivcontrasts Additional options to be passed to \code{AER::ivreg}, corresponding to \code{subset}, \code{subset}, \code{na.action}, etc.. These are all preceded by 'iv' to avoid overlap with same-named options being passed in \code{...}. Note that \code{ivweights} will be ignored if \code{est = 'weight'} or \code{est = 'both'}. \code{ivweights} with \code{est = 'groups'} can be used to implement weighting schemes other than the power-of-an-F-statistic suggested in Huntington-Klein (2019).
#' @param silent Set to TRUE to suppress any messages related to group-creation slow times.
#' @param ... Additional options to be passed to \code{groupSearch} or \code{groupCF}, as appropriate.
#' @examples
#' # This repeats the Stock-Watson textbook cigarette tax IV model
#' # from the AER::ivreg() help file.
#' # ... but it's magnified!
#'
#' # Get data
#' data(CigarettesSW, package = 'AER')
#'
#' CigarettesSW$rprice <- CigarettesSW$price/CigarettesSW$cpi
#' CigarettesSW$rincome <- CigarettesSW$income/CigarettesSW$population/CigarettesSW$cpi
#' CigarettesSW$tdiff <- (CigarettesSW$taxs - CigarettesSW$tax)/CigarettesSW$cpi
#'
#' Cig95 <- CigarettesSW[CigarettesSW$year == "1995",]
#'
#' # Run magnified IV, using causal forest to create effect quantiles
#' magIV <- magnifiedIV(log(packs) ~ log(rprice) + log(rincome) | log(rincome) + tdiff + I(tax/cpi),
#'                      data = Cig95,
#'                      ngroups = 5,
#'                      grouping = 'groupCF',
#'                      groupformula = list(log(rprice) ~ tdiff | I(tax/cpi) + population + taxs,
#'                                          log(rprice) ~ I(tax/cpi) | tdiff + population + taxs))
#' summary(magIV)
#' summary(magIV, vcov = sandwich, df = Inf, diagnostics = TRUE)
#'
#' @export

magnifiedIV <- function(formula, data, ngroups = 4, grouping = 'groupCF', groupformula = NULL, est = 'group', p = 1/4, ivsubset = NULL, ivna.action = getOption("na.action", default = "na.omit"), ivweights = NULL, ivcontrasts = NULL, silent = FALSE, ...) {
  if (nrow(data) <= ngroups) {
    stop("Not enough observations to run. You need at least as many observations as groups.")
  }
  if (nrow(data)/ngroups < 5) {
    warning("Each group will have fewer than five observations in it. While there isn't a specific value for too few observations per group, few observations per group may not reduce IV bias.")
  }
  if (!(as.integer(ngroups) == ngroups)) {
    stop("ngroups must be an integer.")
  }
  if (!is.character(grouping)) {
    stop('grouping must be a string.')
  }
  if (length(grouping) == 1) {
    if (!(grouping %in% c(names(data),'groupCF', 'groupSearch'))) {
      stop('grouping must be a vector of variable names in data, or "groupCF" or "groupSearch".')
    }
  } else {
    if (min(grouping %in% names(data)) == 0) {
      stop('grouping must be a vector of variable names in data, or "groupCF" or "groupSearch".')
    }
  }
  if (!is.character(est)) {
    stop('est must be a string.')
  }
  if (!(est %in% c('group','weight','both'))) {
    stop('est must be "group", "weight", or "both".')
  }
  if (!missing(ivweights) & est %in% c('weight', 'both')) {
    stop('User-provided ivweights cannot be mixed with est = \'weight\' or est = \'both\'.')
  }

  dots <- list(...)

  # Lifted from AER::ivreg
  cl <- match.call()
  if (missing(data))
    data <- environment(formula)
  mf <- match.call()
  m <- match(c("formula", "data"),
             names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  formula <- Formula::as.Formula(formula)
  stopifnot(length(formula)[1] == 1L, length(formula)[2] %in%
              1:2)
  has_dot <- function(formula) inherits(try(terms(formula),
                                            silent = TRUE), "try-error")
  if (has_dot(formula)) {
    f1 <- formula(formula, rhs = 1)
    f2 <- formula(formula, lhs = 0, rhs = 2)
    if (!has_dot(f1) & has_dot(f2))
      formula <- Formula::as.Formula(f1, update(formula(formula,
                                                        lhs = 0, rhs = 1), f2))
  }
  mf$formula <- formula
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  Y <- model.response(mf, "numeric")
  mt <- terms(formula, data = data)
  mtX <- terms(formula, data = data, rhs = 1)
  X <- model.matrix(mtX, mf)
  if (length(formula)[2] < 2L) {
    mtZ <- NULL
    Z <- NULL
  }
  else {
    mtZ <- delete.response(terms(formula, data = data, rhs = 2))
    Z <- model.matrix(mtZ, mf)
  }
  # End borrowed code

  # Figure out what our excluded instruments are
  fsnames <- colnames(Z)
  ssnames <- colnames(X)
  dvname <- as.character(formula)[2]
  exclnames <- fsnames[!(fsnames %in% ssnames) & !(fsnames == '(Intercept)')]
  endonames <- ssnames[!(ssnames %in% fsnames) & !(ssnames == '(Intercept)')]
  exonames <- ssnames[(ssnames %in% fsnames) & !(ssnames == '(Intercept)')]

  if (length(exclnames) < length(endonames)) {
    stop('There must be at least as many excluded instruments as endogenous variables.')
  }
  if (length(endonames) > 1 & est %in% c('weight','both')) {
    warning('est = \'weight\' and est = \'both\' can only be used with a single endogenous variable.')
  }
  if (length(endonames) > 1) {
    warning('See help file section on multiple endogenous variables, as behavior may not be expected.')
  }

  # Create groupings
  groupname <- c()
  # For single instruments
  if (exists('groupformula')) {
    if (!is.list(groupformula)) {
      groupformula <- list(groupformula)
    }
  }
  if (length(groupformula) != length(exclnames) & (identical(grouping,'groupCF') | identical(grouping,'groupsearch'))) {
    stop('There must be at exactly one formula in the groupformula list for each excluded instrument.')
  }
  # Go through each instrument and get groupings
  # Create separate dataset so if formula includes . it doesn't include the groups
  data_w_groups <- data

  for (ex in 1:length(exclnames)) {
    if (identical(grouping,'groupCF')) {

      if (!silent) {
        message(paste0('Starting causal forest at ',Sys.time(),'. This may take a moment.'))
      }

      data_w_groups[,ncol(data_w_groups)+1] <- groupCF(formula = groupformula[[ex]], data = data, ngroups = ngroups, ...)
      if (length(exclnames) == 1) {
        groupname[ex] <- 'groupCF.'
      } else {
        groupname[ex] <- paste0('groupCF',
                                stringr::str_remove_all(
                                  exclnames[ex],
                                  '[:punct:]'),
                                '.')
      }
      names(data_w_groups)[ncol(data_w_groups)] <- groupname[ex]
    } else if (identical(grouping,'groupSearch')) {
      data_w_groups[,ncol(data_w_groups)+1] <- groupSearch(formula = groupformula[[ex]], data = data, ngroups = ngroups, silent = silent, ...)
        if (length(exclnames) == 1) {
          groupname[ex] <- 'groupSearch.'
        } else {
          groupname[ex] <- paste0('groupSearch',
                                  stringr::str_remove_all(
                                    exclnames[ex],
                                    '[:punct:]'),
                                  '.')
        }
        names(data_w_groups)[ncol(data_w_groups)] <- groupname[ex]
    } else {
      groupname[ex] <- grouping[ex]
      if (!is.factor(data[[groupname[ex]]])) {
        data_w_groups[[groupname[ex]]] <- factor(data_w_groups[[groupname[ex]]])
      }
    }
  }

  exoform <- paste(exonames, collapse = '+')
  inxs <- paste(
    paste(exclnames,groupname,sep='*'),
    collapse = '+'
  )
  # IF WE NEED WEIGHTS GET EM
  if (est %in% c('weight','both')) {

    # Get SSE for first stage for personalized F statistic construction
    fsform <- formula(paste(endonames,
                      '~',
                      inxs,
                      '+',
                      exoform, sep = ''))

    fsmodel <- lm(fsform, data = data_w_groups)
    sse0 <- sum((fsmodel$residuals)^2)

    magIVcoef <- data.frame(row.names = 1:nrow(data_w_groups))
    for (ex in 1:length(exclnames)) {
      magIVcoef[[paste('V',ex,sep='')]] <- factorPull(fsmodel,
                                                      data = data_w_groups,
                                                      factor = groupname[ex],
                                                      interaction = exclnames[ex],
                                                      addterm = exclnames[ex])
    }

    # Create prediction for each observation with all exludeds set to 0
    mpred <- rowSums(t(coef(fsmodel)[c('(Intercept)',exonames)]*
      t(model.matrix(fsmodel)[,c('(Intercept)',exonames)])))
    # Get model matrix just for the excludeds alone
    excl <- model.matrix(fsmodel)[,exclnames]

    # For each obs, get the F stat we'd have if those were the coefs!
    n <- nrow(magIVcoef)
    k <- length(fsmodel$coefficients)

    myF <- sapply(1:nrow(magIVcoef),
                  function(x) ownPersonalSSEsus(magIVcoef[x,],
                                                excl,
                                                mpred,
                                                sse0, n, k))
    # Do this step in case p < 0
    myweight <- ifelse(myF == 0,
                       0,
                       myF^p)

  }

  ### Construct ivreg formula
  ivspec <- paste(inxs,exoform, sep = '+')
  ssspec <- paste(c(endonames,exoform,groupname), collapse = '+')
  fullspec <- formula(
    paste0(
      dvname,'~',
      ssspec,'|',
      ivspec))

  ### Call ivreg and return
  if (est == 'both') {
    iv <- AER::ivreg(fullspec,
                     data = data_w_groups,
                     subset = ivsubset,
                     na.action = ivna.action,
                     weights = myweight,
                     contrasts = ivcontrasts)
  } else if (est == 'group') {
    iv <- AER::ivreg(fullspec,
                     data = data_w_groups,
                     subset = ivsubset,
                     na.action = ivna.action,
                     weights = ivweights,
                     contrasts = ivcontrasts)
  } else if (est == 'weight') {
    # Remove group controls and interactions from weight-only version
    ivspec <- paste(c(exclnames,exoform), collapse = '+')
    ssspec <- paste(c(endonames,exoform), collapse = '+')
    fullspec <- formula(
      paste0(
        dvname,'~',
        ssspec,'|',
        ivspec))
    iv <- AER::ivreg(fullspec,
                     data = data_w_groups,
                     subset = ivsubset,
                     na.action = ivna.action,
                     weights = myweight,
                     contrasts = ivcontrasts)
  }

  return(iv)
}


#' Create Personalized F-Statistic
#'
#' This function creates an F-statistic of the personalized kind described in Huntington-Klein (2019).
#'
#' In particular, it calculates (N-K)*Var(x-hat | individual coefficients)/Var(x-hat).
#'
#' In Magnified IV this is raised to the power p to be used as a regression weight. This function is largely used internally for \code{magnifiedIV()} but it is exported here for use in case you'd like to use it to create your own weights. It is perhaps not as user-friendly as it could be, but it is largely an internal function.
#'
#' @param co This is a vector of coefficients unique to one individual, containing only the coefficients that vary across the sample. For example, if the grouped regression model contains \code{z*group}, then \code{co} would be a single value containing an individual's \code{z} effect.
#' @param excl This is a data frame containing just the variable(s) over which the effects vary. So if the regression model contains \code{z*group}, this would be a data frame containing only \code{z}.
#' @param mpred This is a vector containing the regression model prediction *if all variables in \code{excl} were set to 0*.
#' @param sse0,n,k These are the sum of squared errors, the number of observations, and the number of coefficients in the original regression model.
#' @examples
#'
#' df <- data.frame(w1 = rnorm(1000),
#'                  w2 = rnorm(1000),
#'                  e1 = rnorm(1000),
#'                  e2 = rnorm(1000),
#'                  z = rnorm(1000),
#'                  groups = factor(floor(0:999/100)))
#' df$x <- df$w1+df$w2+df$z+df$e1
#'
#' fsmodel <- lm(x ~ z*groups + w1 + w2, data = df)
#'
#' sse0 <- sum((fsmodel$residuals)^2)
#'
#' indivfx <- factorPull(fsmodel,
#'                       data = df,
#'                       factor = 'groups',
#'                       interaction = 'z',
#'                       addterm = 'z')
#'
#' # Create prediction for each observation with just the coefficients that are the
#' # same for everyone
#' mpred <- rowSums(t(coef(fsmodel)[c('(Intercept)','w1','w2')]*
#'                      t(model.matrix(fsmodel)[,c('(Intercept)','w1','w2')])))
#' # Get a data frame of just the variables with effects that vary
#' excl <- data.frame(z = df[['z']])
#'
#' # Get the relevant N and K
#' n <- nrow(df)
#' k <- length(fsmodel$coefficients)
#'
#' # And finally produce that individualized F statistic
#' indiv.F <- sapply(indivfx, function(x)
#'   ownPersonalSSEsus(x, excl, mpred, sse0, n, k))
#' @export

ownPersonalSSEsus <- function(co, excl, mpred, sse0, n, k) {

  # Construct prediction with the given coefficients
  pred <- mpred +
    rowSums(t(co*
              t(excl)))

  ss1 <- sum((pred)^2)

  Fstat <- (n-k)*ss1/sse0
  return(Fstat)
}
