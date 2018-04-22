#' Interaction Screening for Ultra-High Dimensional Data
#'
#' Extended variable selection approaches to jointly model main and interaction effects from high-dimensional data orignally proposed by Hao and Zhang (2014) and extended by Gosik and Wu (2016).
#' Based on a greedy forward approach, their model can identify all possible interaction effects through two algorithms, iFORT and iFORM, which have been proved to possess sure screening property in an ultrahigh-dimensional setting.
#'
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under ‘Details’.
#' @param data data.frame of your data with the response and all p predictors
#' @param strong  logical TRUE to use strong heredity or FALSE to use weak heredity (default TRUE)
#' @param higher_order  logical TRUE indicating to include order-3 interactions in the search (default FALSE)
#' @return a summary of the linear model returned after the selection procedure
#' @author Kirk Gosik
#' @details
#' Runs the iForm selection procedure on the dataset and returns a linear model
#' of the final selected model.
#' @seealso \code{lm}\
#' @seealso \code{model.matrix}
#' @export
#' @importFrom stats lm
#' @importFrom stats summary.lm
#' @importFrom stats model.matrix
#' @importFrom stats model.frame



iForm <- function(formula, data, strong = TRUE, higher_order = FALSE) {

  dat <- model.frame(formula, data)
  y <- dat[ , 1]
  x <- dat[ , -1]
  p <- ncol(x)
  n <- nrow(x)
  solution <- NULL
  model <- NULL
  step <- 1
  bic <- NULL

  candidate <- as.data.frame(x)
  fit <- iformselect(x, y, p, n, candidate, solution, model, bic, step, strong, higher_order)


  y <- fit$y
  solution <- fit$solution
  model <- fit$model
  bic <- fit$bic

  model <- data.frame(solution[ , 1 : which.min(bic), drop = FALSE])
  lm(y ~ . + 0 , data = model)

}



#' Inner workings for different selections under different higher-orders and strength of marginality.
#'
#' Extended variable selection approaches to jointly model main and interaction effects from high-dimensional data orignally proposed by Hao and Zhang (2014) and extended by Gosik and Wu (2016).
#' Based on a greedy forward approach, their model can identify all possible interaction effects through two algorithms, iFORT and iFORM, which have been proved to possess sure screening property in an ultrahigh-dimensional setting.
#'
#' @param x design matrix of predictors
#' @param y response variable
#' @param p the size of the predictor set
#' @param n the size of the number of observations
#' @param candidate the current candidate set of predictors to select from
#' @param solution the current set of predictors already selected
#' @param model the set of predictors to use in the final model
#' @param bic the cutoff value for determining the model set
#' @param step the step in the iteration currently on
#' @param strong indicator of the strength of marginality to be used
#' @param higher_order  logical TRUE indicating to include order-3 interactions in the search (default FALSE)
#' @return a summary of the linear model returned after the selection procedure
#' @author Kirk Gosik
#' @details
#' Runs the iForm selection procedure on the dataset and returns a linear model
#' of the final selected model.
#' @seealso \code{lm}\
#' @seealso \code{model.matrix}
#' @export
#' @importFrom stats lm
#' @importFrom stats summary.lm
#' @importFrom stats model.matrix
#' @importFrom stats formula
#' @importFrom utils combn


iformselect <- function( x, y, p, n, candidate, solution, model, bic, step, strong, higher_order ) {

  repeat{

    rss <- sapply(colnames(candidate), function(i){
      xx <- as.matrix(cbind(solution, candidate[ , i]))
      tryCatch({
        sum((y - xx %*% (solve(t(xx) %*% xx)) %*% (t(xx) %*% y)) ^ 2)
      }, error = function(e) Inf)
    })

    solution <- cbind(solution, candidate[ , which.min(rss)])
    colnames(solution)[step] <- colnames(candidate)[which.min(rss)]
    candidate <- candidate[ , -which.min(rss)]

    if( ncol(solution) > 1 ) {
      if( strong ){
        interaction_formula <- paste("~0+(", paste(colnames(solution)[colnames(solution) %in% colnames(x)], collapse= "+"),")^2")
      }else{
        interaction_formula <- paste("~0+(", paste(colnames(solution)[colnames(solution) %in% colnames(x)], collapse= "+"),")*(", paste(colnames(x),collapse = "+"),")")
      }

      if(higher_order){
        if(sum(!{colnames(solution) %in% colnames(x)}) > 2) {

          if(strong){
          tmp <- Filter(function(x) {length(x) == 2}, strsplit(colnames(solution), "[.]|[:]"))
          tmp <- combn(length(tmp), 3, function(x) Reduce(intersect, combn(x, 2, function(y) Reduce(union, tmp[y]), simplify=F)), simplify=F)
          tmp <- Map(function(x) {paste0(x, collapse=":")}, Filter(function(x){length(x)==3}, tmp))
            if(length(tmp) > 0) {interaction_formula <- paste(interaction_formula, paste0(unlist(tmp), collapse = "+"), sep = "+")}

          }else{
            tmp <- Filter(function(x) {length(x) == 2}, strsplit(colnames(solution), "[.]|[:]"))
            tmp <- Map(function(x) {paste0(x, collapse=":")}, Filter(function(x){length(x)==2}, tmp))
            if(length(tmp) > 0) {interaction_formula <- paste(interaction_formula, paste("(",paste(unlist(tmp), collapse="+"), ")*(", paste(colnames(x), collapse = "+"), ")"), sep = "+")}
          }

        }
      }

      interaction_formula <- as.formula(interaction_formula)
      interactions <- model.matrix(interaction_formula, x)
      interactions <- interactions[ , setdiff(colnames(interactions), colnames(solution)), drop = F]
      interactions <- interactions[ , setdiff(colnames(interactions), colnames(candidate)), drop = F]
      candidate <- cbind(candidate, interactions)
    }

    bic[step] <- log(as.numeric(min(rss))/n) + (dim(solution)[2]) * (log(n) + 2 * log(p))/n
    if((dim(solution)[2]) > (n/log(n, 2)))break
    step <- step + 1

  }

  list(y = y, solution = solution, model = model, bic = bic)
}


