#' Interaction Screening for Ultra-High Dimensional Data
#'
#' Extended variable selection approaches to jointly model main and interaction effects from high-dimensional data orignally proposed by Hao and Zhang (2014) and extended by Gosik and Wu (2016).
#' Based on a greedy forward approach, their model can identify all possible interaction effects through two algorithms, iFORT and iFORM, which have been proved to possess sure screening property in an ultrahigh-dimensional setting.
#'
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under ‘Details’.
#' @param data data.frame of your data with the response and all p predictors
#' @param strong  logical TRUE to use strong heredity or FALSE to use weak heredity (default TRUE)
#' @param higher_order  logical TRUE indicating to include order-3 interactions in the search (default FALSE)
#' @param use_ff  logical TRUE to use atomic data structur that are stored on disk via the ff package. (default FALSE)
#' @return a summary of the linear model returned after the selection procedure
#' @author Kirk Gosik
#' @details
#' Runs the iFORM selection procedure on the dataset and returns a linear model
#' of the final selected model.
#' @seealso \code{lm}\
#' @seealso \code{model.matrix}
#' @export
#' @importFrom stats lm
#' @importFrom stats summary.lm
#' @importFrom stats model.matrix
#' @importFrom ff as.ff


iForm <- function(formula, data, strong = TRUE, higher_order = FALSE, use_ff = FALSE) {

  dat <- model.frame(formula, data)
  y <- dat[ , 1]
  x <- dat[ , -1]
  p <- ncol(x)
  n <- nrow(x)
  solution <- NULL
  model <- NULL
  step <- 1
  bic <- NULL

  if( use_ff ) {

    candidate <- as.ff(as.matrix(x), colnames = colnames(x), overwrite = TRUE)
    fit <- iformselect_ff(x, y, p, n, candidate, solution, model, bic, step, strong, higher_order)


  }else{

    candidate <- as.data.frame(x)
    fit <- iformselect(x, y, p, n, candidate, solution, model, bic, step, strong, higher_order)

    }

  y <- fit$y
  solution <- fit$solution
  model <- fit$model
  bic <- fit$bic

  model <- data.frame(solution[ , 1 : which.min(bic), drop = FALSE])
  lm(y ~ . + 0 , data = model)

}





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






iformselect_ff <- function(x, y, p, n, candidate, solution, model, bic, step, strong, higher_order){

  repeat{

    rss <- sapply(colnames(candidate), function(i){
      xx <- as.matrix(cbind(solution, candidate[ , i]))
      tryCatch({
        sum((y - xx %*% (solve(t(xx) %*% xx)) %*% (t(xx) %*% y)) ^ 2)
      }, error = function(e) Inf)
    })

    rss <- as.numeric(rss)

    solution <- cbind(solution, candidate[ , which.min(rss)])
    colnames(solution)[step] <- colnames(candidate)[which.min(rss)]
    cnames <- setdiff(colnames(candidate), colnames(solution))

    candidate[, which.min(rss) : (ncol(candidate) - 1)] <- candidate[, (which.min(rss) + 1) : ncol(candidate)]
    dim(candidate) <- c(n, (ncol(candidate) - 1))
    colnames(candidate) <- cnames

    if( ncol(solution) > 1 ) {
      if( strong ){
        interaction_formula <- paste("~0+(", paste(colnames(solution)[colnames(solution) %in% colnames(x)], collapse= "+"),")^2")
      }else{
        interaction_formula <- paste("~0+(", paste(colnames(solution)[colnames(solution) %in% colnames(x)], collapse= "+"),")*(", paste(colnames(x),collapse = "+"),")")
      }

      if( higher_order ) {

        if( sum(!{colnames(solution) %in% colnames(x)}) > 2 ) {

          if( strong ) {

            tmp <- Filter(function(x) {length(x) == 2}, strsplit(colnames(solution), "[.]|[:]"))
            tmp <- combn(length(tmp), 3, function(x) Reduce(intersect, combn(x, 2, function(y) Reduce(union, tmp[y]), simplify=F)), simplify=F)
            tmp <- Map(function(x) {paste0(x, collapse=":")}, Filter(function(x){length(x)==3}, tmp))

            if( length(tmp) > 0 ) {

              interaction_formula <- paste(interaction_formula, paste0(unlist(tmp), collapse = "+"), sep = "+")

              }

          }else{

            tmp <- Filter(function(x) {length(x) == 2}, strsplit(colnames(solution), "[.]|[:]"))
            tmp <- Map(function(x) {paste0(x, collapse=":")}, Filter(function(x){length(x)==2}, tmp))
            if( length(tmp) > 0 ) {

              interaction_formula <- paste(interaction_formula, paste("(",paste(unlist(tmp), collapse="+"), ")*(", paste(colnames(x), collapse = "+"), ")"), sep = "+")

              }

          }

        }
      }

      interaction_formula <- as.formula(interaction_formula)
      interactions <- model.matrix(interaction_formula, x)
      interactions <- interactions[ , setdiff(colnames(interactions), colnames(solution)), drop = F]
      interactions <- interactions[ , setdiff(colnames(interactions), colnames(candidate)), drop = F]

      if( ncol(interactions) > 0 ){
        newdimstart <- ncol(candidate) + 1
        newdimend <- ncol(candidate) + ncol(interactions)
        dim(candidate) <- c(n, newdimend)
        candidate[, newdimstart : newdimend] <- interactions
        colnames(candidate) <- c(cnames, colnames(interactions))
      }

    }

    bic[step] <- log(as.numeric(min(rss))/n) + (dim(solution)[2]) * (log(n) + 2 * log(p))/n
    if((dim(solution)[2]) > (n / log(n, 2))) break
    step <- step + 1

  }

  list(y = y, solution = solution, model = model, bic = bic)
}

