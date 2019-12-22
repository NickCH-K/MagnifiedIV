#' factorPull
#'
#' This function takes any estimated regression compatible with \code{broom::tidy()} that contains a factor variable, or a factor variable interacted with something else. It will pull out those coefficients and line them up with the original data so you can store them as a variable.
#'
#' If multiple coefficients related to the factor are dropped from the model, be sure to check the result.
#'
#' @param model A model object compatible with \code{broom::tidy()}.
#' @param data A data set containing the factor variable you want coefficients for.
#' @param factor A string with the name of the factor variable you'd like to get the coefficients of. Note if the variable is included in regression via the \code{factor()} function, include that. So if your model is \code{mpg~factor(cyl)+hp}, do \code{factor = 'factor(cyl)'}.
#' @param interaction A string with the name of a variable being interacted with \code{factor} for which you'd like to get the different interaction terms. Specify the name as it shows up in the regression table. So for example if \code{factor = 'X'} and your model has coefficients \code{Z1:X1, Z1:X2, Z2:X1, Z2:X2}, then \code{interaction = 'Z1'} will get the \code{Z1:X1, Z1:X2} coefficients.
#' @param specify For more complex sets of interactions, like three-way interactions, or if your \code{factor} variable is included via some function more complex than \code{factor()}, you can specify the exact naming structure of the coefficients, with \code{{value}} standing in for the factor value. So for example if you had coefficients \code{Z:X1:W, Z:X2:W} in your model, you could do \code{specify = 'Z:X{value}:W'}. This will override anything in \code{interaction}.
#' @param basevalue If \code{specify} is used, what value should the \code{factor} reference group (or any excluded terms) be given? By default \code{0}. Set to \code{NA} to leave omitted terms as \code{NA}.
#' @param addterm A string indicating a coefficient in the model to be added to all terms. Commonly this is the coefficient for the variable being interacted with, to give an overall effect of that variable, or \code{'(Intercept)'} to give a mean prediction.
#' @param value A string with the column name of the \code{broom::tidy()} result that you would like to extract. By default this is \code{value = 'estimate'} to get the coefficients.
#' @param na.predict Should observations in the data but not in the model, which do have nonmissing values of \code{factor}, be included in predictions? Setting this to \code{FALSE} requires that sending your \code{model} through \code{fitted()} produces a named vector.
#' @examples
#' df <- data.frame(w1 = rnorm(1000),
#'                  w2 = rnorm(1000),
#'                  e1 = rnorm(1000),
#'                  e2 = rnorm(1000),
#'                  z = rnorm(1000),
#'                  groups = factor(floor(0:999/100)))
#' df$x <- df$w1+df$w2+df$z+df$e1
#'
#' # Create a model with an interaction between z and groups
#' lm_inx <- lm(x ~ z*groups + w1 + w2, data = df)
#'
#' # Get the effect of z for each group
#' indivfx <- factorPull(lm_inx,
#'                       data = df,
#'                       factor = 'groups',
#'                       interaction = 'z',
#'                       addterm = 'z')
#'
#' # Create a model with a fixed effect for groups
#' lm_inx2 <- lm(x ~ groups + z + w1 + w2, data = df)
#'
#' # Get the fixed effect for each group, with 0 for the reference group
#' indivfe <- factorPull(lm_inx,
#'                       data = df,
#'                       factor = 'groups')
#' # Or add in the intercept to get the full fixed effect
#' indivfe2 <- factorPull(lm_inx,
#'                        data = df,
#'                        factor = 'groups',
#'                        addterm = '(Intercept)')
#' @export

factorPull <- function(model, data, factor, interaction = NULL, specify = NULL, basevalue = 0, addterm = NULL, value = 'estimate', na.predict = TRUE) {
  if (!is.character(value)) {
    stop('value must be a string.')
  }
  if (!is.logical(na.predict)) {
    stop('na.predict must be a string.')
  }
  if (!is.character(factor)) {
    stop('factor must be a string.')
  }
  if (!is.null(interaction)) {
    if (!is.character(interaction)) {
      stop('interaction must be a string.')
    }
  }
  if (!is.null(specify)) {
    if (!is.character(specify)) {
      stop('specify must be a string.')
    }
    if (!stringr::str_detect(specify,'\\{value\\}')) {
      stop('specify must contain {value} somewhere.')
    }
  }

  # Get tidied model
  mf <- broom::tidy(model)
  terms <- mf$term

  # Get the actual name of the variable
  if (stringr::str_detect(factor,'factor\\(')) {
    fname <- stringr::str_sub(stringr::str_extract(factor,'\\(.*\\)'),2,-2)
  } else {
    fname <- factor
  }

  # Create specify if not given
  if (is.null(specify)) {
    if (is.null(interaction)) {
      specify <- paste0(factor,'{value}')
    } else {
      # Determine whether it's factor:interaction or the other way around.
      escapestr <- escape_them(paste0(interaction,':',factor))
      inxfac <- sum(
        stringr::str_detect(
          terms,
          escapestr
          )) > 0
      if (inxfac) {
        specify <- paste0(interaction,':',factor,'{value}')
      } else {
        specify <- paste0(factor,'{value}:',interaction)
      }
    }
  }

  # Get the values of the factor variable
  vals <- unique(data[[fname]])
  vals <- vals[!is.na(vals)]

  # Get the relevant coefficient names
  cnames <- stringr::str_replace(specify,'\\{value\\}',as.character(vals))

  # Line up factor values, coefficient names, and values
  coefs <- data.frame(term = cnames, vals = vals)
  coefs <- merge(coefs, mf[,c('term',value)], by = 'term', all.x = TRUE)

  # For base category
  coefs[[value]][is.na(coefs[[value]])] <- basevalue

  # Adding-in
  if (!is.null(addterm)) {
    coefs[[value]] <- coefs[[value]] + mf[[value]][mf$term == addterm]
  }

  # Duplicate to match data
  row.names(coefs) <- coefs$vals

  # And bring in
  data[,ncol(data)+1] <- coefs[factor(data[[fname]]),value]
  varname <- names(data)[ncol(data)]

  # If we want those who weren't in the regression to get NA, then do that
  if (!na.predict) {
    obsin <- names(fitted(model))

    data[,ncol(data)+1] <- FALSE
    data[obsin,ncol(data)] <- TRUE

    data[[varname]] <- ifelse(data[[ncol(data)]],
                              data[[varname]],
                              NA)
  }

  return(data[[varname]])
}


escape_them <- function(s) {

  specs <- c('!','"','#','$','%','&','’','(',')',
             '*','+',',','-','.','/',':',';','<',
             '=','>','?','@','[',']','^','_','`',
             '{','|','}','~')

  for (c in specs) {
    s <- stringr::str_replace_all(s,
      paste0('\\',c),
      paste0('\\\\',c)
    )
  }

  return(s)
}
