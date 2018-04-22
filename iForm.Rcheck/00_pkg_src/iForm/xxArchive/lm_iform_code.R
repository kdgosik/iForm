## need to add , strong = TRUE, higher_order = FALSE, use_ff = FALSE
## change fit to iform select

iform_lm <- function (formula, data, subset, na.action, method = "qr",
          model = TRUE, x = FALSE, y = FALSE, qr = TRUE, singular.ok = TRUE,
          contrasts = NULL, offset, ...)
{
  ret.x <- x
  ret.y <- y
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action",
               "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  if (method == "model.frame")
    return(mf)
  else if (method != "qr")
    warning(gettextf("method = '%s' is not supported. Using 'qr'",
                     method), domain = NA)
  mt <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  w <- as.vector(model.weights(mf))
  if (!is.null(w) && !is.numeric(w))
    stop("'weights' must be a numeric vector")
  offset <- as.vector(model.offset(mf))
  if (!is.null(offset)) {
    if (length(offset) != NROW(y))
      stop(gettextf("number of offsets is %d, should equal %d (number of observations)",
                    length(offset), NROW(y)), domain = NA)
  }
  if (is.empty.model(mt)) {
    x <- NULL
    z <- list(coefficients = if( is.matrix(y) ) {
      matrix(, 0, 3)
      }else {
        numeric()
        },
      residuals = y,
      fitted.values = 0 * y,
      rank = 0L,
      df.residual = length(y))

    if (!is.null(offset)) {
      z$fitted.values <- offset
      z$residuals <- y - offset
    }
  }
  else {
    x <- model.matrix(mt, mf, contrasts)
    z <- if (is.null(w))
      lm.fit(x, y, offset = offset, singular.ok = singular.ok,
             ...)
    else lm.wfit(x, y, w, offset = offset, singular.ok = singular.ok,
                 ...)
  }
  class(z) <- c(if (is.matrix(y)) "mlm", "lm")
  z$na.action <- attr(mf, "na.action")
  z$offset <- offset
  z$contrasts <- attr(x, "contrasts")
  z$xlevels <- .getXlevels(mt, mf)
  z$call <- cl
  z$terms <- mt
  if (model)
    z$model <- mf
  if (ret.x)
    z$x <- x
  if (ret.y)
    z$y <- y
  if (!qr)
    z$qr <- NULL
  z
}




## iformselect ######################################

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




## iformselect_ff ######################################

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

