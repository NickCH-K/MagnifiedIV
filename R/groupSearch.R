#' groupSearch
#'
#' This function performs the GroupSearch algorithm as in Huntington-Klein (2019) "Instruments with Heterogeneous Effects: Monotonicity, Bias, and Localness."
#'
#' The GroupSearch algorithm is naive. It simply tries a bunch of random groupings, and for each grouping attempts to predict \code{x} with \code{z}. It returns the grouping that produces the highest F-statistic as a factor vector. Be aware before using \code{groupSearch()} that in the original paper it only did a mediocre job at picking up effect heterogeneity. You may want to use \code{groupCF()} instead.
#'
#' This function is called by magnifiedIV. You can also run Magnified IV by yourself without the magnifiedIV function (with any estimator) by running groupSearch, then adding the resulting group variable as a control in both IV stages and also interacted with the instrument. Or use \code{factorPull()} to get the individual-level effects estimates and use those to construct a sample weight.
#'
#' @param formula A formula of the form \code{x ~ z | w1 + w2} where \code{x} is an endogenous variable in an instrumental variables model, being predicted here. \code{z} is one of the instruments, and \code{w1} and \code{w2}, etc., are covariates to be partialed out, if any.
#' @param data A data.frame.
#' @param weights Estimation weights.
#' @param ngroups Number of groups to split the data into.
#' @param ntries Number of groupings to attempt.
#' @param id A variable in \code{data} that indiates that observations with the same value of \code{id} should always be in the same group.
#' @param silent Suppress the progress report.
#' @param ... Additional arguments to be passed to \code{lm()}. Note that \code{na.action} will be ignored for partialling-out of the covariates, and if you prefer a different \code{na.action} for this purpose you should partial-out by hand before running \code{groupSearch}.
#'
#' @examples
#' # Get data
#' data(CPS1985, package = 'AER')
#'
#' # Split the data into 10 random groups 100 times, and each time see how the effect of
#' # education in predicting wages varies across the sample, after controlling
#' # for all the other variables in the data, plus a squared term on experience.
#' # Return the group with the largest resulting F statistic.
#' edeffect <- groupSearch(wage ~ education |
#'                           experience + I(experience^2) + age + ethnicity +
#'                           region + gender + occupation +
#'                           sector + union + married,
#'                         data = CPS1985, ngroups = 10)
#'
#' table(edeffect)
#' @export

groupSearch <- function(formula, data, weights, ngroups = 4L, ntries = 100L, id = NULL, silent = FALSE, ...) {
  if (nrow(data) <= ngroups) {
    stop("Not enough observations to run. You need at least as many observations as groups.")
  }
  if (nrow(data)/ngroups < 5) {
    warning("Each group will have fewer than five observations in it. While there isn't a specific value for too few observations per group, few observations per group may not reduce IV bias.")
  }
  if (!(as.integer(ngroups) == ngroups)) {
    stop("ngroups must be an integer.")
  }
  if (!(as.integer(ntries) == ntries)) {
    stop("ntries must be an integer.")
  }
  # Check for dots and create a version without na.action
  dots_no_na <- list(...)
  dots_no_na[['na.action']] <- NULL

  # This is mostly lifted from AER::ivreg()
  cl <- match.call()
  if (missing(data))
    data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "weights", "id"),
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
  X <- model.response(mf, "numeric")
  mt <- terms(formula, data = data)
  mtZ <- terms(formula, data = data, rhs = 1)
  Z <- model.matrix(mtZ, mf)
  if (length(formula)[2] < 2L) {
    mtW <- NULL
    W <- NULL
  }
  else {
    mtW <- delete.response(terms(formula, data = data, rhs = 2))
    W <- model.matrix(mtW, mf)
  }
  weights <- model.weights(mf)
  ID <- model.extract(mf, "id")
  # END HERE

  # Make sure we only have one x and one z, and if we have controls or weights
  if (is.matrix(X)) {
    stop('Only one x variable allowed.')
  }
  if (ncol(Z) > 2) {
    stop('Only one z variable allowed.')
  }
  Z <- Z[,2]
  haveW <- !is.null(W)
  haveWT <- !is.null(weights)

  # Partial out if we have controls
  if (haveW & haveWT) {
    X <- residuals(lm(X ~ W, weights = weights, dots_no_na, na.action = na.exclude))
    Z <- residuals(lm(Z ~ W, weights = weights, dots_no_na, na.action = na.exclude))
  } else if (haveW & !haveWT) {
    X <- residuals(lm(X ~ W, dots_no_na, na.action = na.exclude))
    Z <- residuals(lm(Z ~ W, dots_no_na, na.action = na.exclude))
  }

  # Now just the variables we need
  if (haveWT) {
    data <- data.frame(
      x = X,
      z = Z,
      wts = weights,
      groups = NA
    )
  } else {
    data <- data.frame(
      x = X,
      z = Z,
      groups = NA
    )
  }
  if (!is.null(ID)) {
    data$origorder <- 1:nrow(data)
    data$id <- ID
    data <- data[order(data$id),]

    # For expanding groups later
    idlens <- rle(data$id)$lengths
  }

  winnerF <- 0
  winnergroups <- 0

  # For progress report
  reports <- ceiling(ntries*(1:9)/10)

  for (b in 1:ntries) {

    # Progress report
    if (!silent) {
      if (b == 1) {
        message(paste('Starting groupSearch at',Sys.time()))
      } else if (b %in% reports) {
         pct = as.integer(which(b == reports)*10)
         message(paste0(pct,'% done at ',Sys.time()))
      }
    }

    # Build groups
    if (is.null(ID)) {
      data$groups <- sample(
        factor(1:ngroups),
        nrow(data),
        replace = T)
    } else {
      data$groups <- rep(
        sample(
          factor(1:ngroups),
          length(idlens),
          replace = T
        ),
        idlens)
    }

    if (haveWT) {
      f <- summary(lm(x~z*groups, weights = wts, data=data, ...))$fstatistic[1]
    } else{
      f <- summary(lm(x~z*groups, data=data, ...))$fstatistic[1]
    }


    if (f > winnerF) {
      winnerF <- f
      winnergroups <- data$groups
    }
  }

  if (is.null(ID)) {
    return(winnergroups)
  } else {
    data$groups <- winnergroups
    data <- data[order(data$origorder),]
    return(data$groups)
  }
}

