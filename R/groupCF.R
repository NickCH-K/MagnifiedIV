#' groupCF
#'
#' This is a wrapper for \code{grf::causal_forest()} which performs a causal forest estimation, pulls out the individual-level effect estimates (tau), and then returns a factor variable containing the quantiles of tau.
#'
#' This function is called by magnifiedIV. You can also run Magnified IV by yourself (with any estimator) by running groupSearch, then adding the resulting group variable as a control in both IV stages and also interacted with the instrument. Or use \code{grf::causal_forest()} directly to estimate tau on the individual level, and use that to construct a sample weight.
#'
#' @param formula A formula of the form \code{x ~ z | w1 + w2} where \code{x} is an endogenous variable in an instrumental variables model, being predicted here. \code{z} is one of the instruments, and \code{w1} and \code{w2}, etc., are covariates to be used to predict the effect of the instrument.
#' @param data A data.frame.
#' @param ngroups Number of quantiles to split tau into.
#' @param ... Additional arguments to be passed to \code{grf::causal_forest()}.
#'
#' @examples
#' # Get data
#' data(CPS1985, package = 'AER')
#'
#' # See how the effect of education on wage varies over all the variables in the data
#' # and then split the resulting individual coefficient estimates into 10 quantiles
#' # (note this example probably does not satisfy the theoretical unconfoundedness
#' # assumption of causal forest; this is just a code example)
#' edeffect <- groupCF(wage ~ education |
#'                       experience + age + ethnicity +
#'                       region + gender + occupation +
#'                       sector + union + married,
#'                     data = CPS1985, ngroups = 10)
#'
#' table(edeffect)
#' @export

groupCF <- function(formula, data, ngroups = 4L, ...) {
  if (nrow(data) < ngroups) {
    stop("Not enough observations to run. You need at least as many observations as groups.")
  }
  if (nrow(data)/ngroups < 5) {
    warning("Each group will have fewer than five observations in it. While there isn't a specific value for too few observations per group, few observations per group may not reduce IV bias.")
  }
  if (!(as.integer(ngroups) == ngroups)) {
    stop("ngroups must be an integer.")
  }

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
    stop('groupCF() requires covariates.')
  }
  else {
    mtW <- delete.response(terms(formula, data = data, rhs = 2))
    W <- model.matrix(mtW, mf)
  }
  # END HERE

  # Make sure we only have one x and one z, and if we have controls or weights
  if (is.matrix(X)) {
    stop('Only one x variable allowed.')
  }
  if (ncol(Z) > 2) {
    stop('Only one z variable allowed.')
  }
  Z <- Z[,2]

  cf <- grf::causal_forest(W, X, Z, ...)

  tau.hat <- predict(cf)$predictions

  groups <- cut(
      tau.hat,
      quantile(tau.hat,(0:ngroups)/ngroups),
      include.lowest = TRUE
    )

  return(groups)
}

