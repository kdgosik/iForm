data <- data.frame(matrix(rbinom(10000,1,.5),nrow=200) - .5)
response <- 3*data[,1] + 4*data[,2] + 2*data[,4] +
  3*data[,2]*data[,3] + 3*data[,2]*data[,4] + rnorm(200,0,1)

data <- data.frame(response,data)
formula <- response ~ .


iForm.weak <- function(data, response){

  y <- data[ , response]
  x <- data[ , -which(colnames(data) == response)]
  p <- dim(x)[2]
  n <- dim(x)[1]
  candidate <- data.frame(x)
  solution <- NULL
  model <- NULL
  step <- 1
  bic <- NULL

  repeat{

    rss <- sapply(colnames(candidate), function(i){
      xx <- as.matrix(cbind(solution, candidate[ , i]))
      try(sum((y - xx %*% (solve(t(xx) %*% xx)) %*% (t(xx) %*% y)) ^ 2), silent = T)
    })

    solution <- data.frame(cbind(solution, candidate[ , which.min(rss)]))
    colnames(solution)[step] <- colnames(candidate)[which.min(rss)]
    candidate <- candidate[ , -which.min(rss)]
    if(dim(solution)[2] > 1) {
      interactions <- data.frame(model.matrix(as.formula(paste("~0+(", paste(colnames(solution)[colnames(solution) %in% colnames(x)], collapse= "+"),")*(", paste(colnames(x),collapse = "+"),")")),x))
      interactions <- interactions[ , setdiff(colnames(interactions), colnames(solution)), drop = F]
      interactions <- interactions[ , setdiff(colnames(interactions), colnames(candidate)), drop = F]
      candidate <- cbind(candidate, interactions)
    }

    bic[step] <- log(as.numeric(min(rss))/n) + (dim(solution)[2]) * (log(n) + 2 * log(p))/n
    if((dim(solution)[2]) > (n/log(n, 2)))break
    step <- step + 1

  }

  model <- data.frame(solution[ , 1 : which.min(bic), drop = FALSE])
  summary(lm(y ~ . + 0 , data = model))

}





 # Developing weak heredity version
  # scenario 1: Weak iForm two-way interaction
  # Two-interaction needs to be included in model for three-way interactions

iForm.weak1.three <- function(data, response){

  y <- data[ , response]
  x <- data[ , -which(colnames(data) == response)]
  p <- dim(x)[2]
  n <- dim(x)[1]
  candidate <- data.frame(x)
  solution <- NULL
  model <- NULL
  step <- 1
  bic <- NULL

  repeat{

    rss <- sapply(colnames(candidate), function(i){
      xx <- as.matrix(cbind(solution, candidate[ , i]))
      try(sum((y - xx %*% (solve(t(xx) %*% xx)) %*% (t(xx) %*% y)) ^ 2),silent = T)
    })

    solution <- data.frame(cbind(solution, candidate[ , which.min(rss)]))
    colnames(solution)[step] <- colnames(candidate)[which.min(rss)]
    candidate <- candidate[ , -which.min(rss)]

    if(dim(solution)[2] > 1) {

      formula <- paste("(", paste(colnames(solution)[colnames(solution) %in% colnames(x)], collapse= "+"),")*(", paste(colnames(x),collapse = "+"),")")

      if(sum(!{colnames(solution) %in% colnames(x)}) > 0) {
        tmp <- Filter(function(x) {length(x) == 2}, strsplit(colnames(solution), "[.]|[:]"))
        tmp <- Map(function(x) {paste0(x, collapse=":")}, Filter(function(x){length(x)==2}, tmp))
        if(length(tmp) > 0) {formula <- paste(formula,paste("(",paste(unlist(tmp),collapse="+"),")*(",paste(colnames(x), collapse = "+"),")"), sep="+")}
      }

      interactions <- data.frame(model.matrix(as.formula(paste("~0+", formula)), x))
      interactions <- interactions[ , setdiff(colnames(interactions), colnames(solution)), drop = F]
      interactions <- interactions[ , setdiff(colnames(interactions), colnames(candidate)), drop = F]
      candidate <- cbind(candidate, interactions)

    }

    bic[step] <- log(as.numeric(min(rss))/n) + (dim(solution)[2]) * (log(n) + 2 * log(p))/n
    if((dim(solution)[2]) > (n/log(n, 2)))break
    step <- step + 1

  }

  model <- data.frame(solution[ , 1 : which.min(bic), drop = FALSE])
  summary(lm(y ~ . + 0 , data = model))

}





# Developing weak heredity version
  # Scenario 2: Weak heredity for all interactions (very weak heredity)


iForm.weak2.three <- function(data, response){

  y <- data[ , response]
  x <- data[ , -which(colnames(data) == response)]
  p <- dim(x)[2]
  n <- dim(x)[1]
  candidate <- data.frame(x)
  solution <- NULL
  model <- NULL
  step <- 1
  bic <- NULL

  repeat{

    rss <- sapply(colnames(candidate), function(i){
      xx <- as.matrix(cbind(solution, candidate[ , i]))
      try(sum((y - xx %*% (solve(t(xx) %*% xx)) %*% (t(xx) %*% y)) ^ 2),silent = T)
    })

    solution <- data.frame(cbind(solution, candidate[ , which.min(rss)]))
    colnames(solution)[step] <- colnames(candidate)[which.min(rss)]
    candidate <- candidate[ , -which.min(rss)]

    if(dim(solution)[2] > 1) {

      formula <- paste("(", paste(colnames(solution)[colnames(solution) %in% colnames(x)], collapse= "+"),")*(", paste(colnames(x),collapse = "+"),")")
      formula <- paste(formula, "+", formula, "*(", paste(colnames(x), collapse = "+"),")")

      interactions <- data.frame(model.matrix(as.formula(paste0("~0+", formula)), x))
      interactions <- interactions[ , setdiff(colnames(interactions), colnames(solution)), drop = F]
      interactions <- interactions[ , setdiff(colnames(interactions), colnames(candidate)), drop = F]
      candidate <- cbind(candidate, interactions)

    }

    bic[step] <- log(as.numeric(min(rss))/n) + (dim(solution)[2]) * (log(n) + 2 * log(p))/n
    if((dim(solution)[2]) > (n/log(n, 2)))break
    step <- step + 1

  }

  model <- data.frame(solution[ , 1 : which.min(bic), drop = FALSE])
  summary(lm(y ~ . + 0 , data = model))

}



check1 <- iForm.weak1.three(data,"response")
check2 <- iForm.weak2.three(data,"response")







# Developing weak heredity version
# scenario 1: Weak iForm two-way interaction with intercept
# Two-interaction needs to be included in model for three-way interactions

iForm.weak1.three.intercept <- function(data, response){

  y <- data[ , response]
  x <- data[ , -which(colnames(data) == response)]
  p <- dim(x)[2]
  n <- dim(x)[1]
  candidate <- data.frame(x)
  solution <- NULL
  model <- NULL
  step <- 1
  bic <- NULL

  repeat{

    rss <- sapply(colnames(candidate), function(i){
      xx <- as.matrix(cbind(solution, candidate[ , i]))
      try(sum((y - xx %*% (solve(t(xx) %*% xx)) %*% (t(xx) %*% y)) ^ 2),silent = T)
    })

    solution <- data.frame(cbind(solution, candidate[ , which.min(rss)]))
    colnames(solution)[step] <- colnames(candidate)[which.min(rss)]
    candidate <- candidate[ , -which.min(rss)]

    if(dim(solution)[2] > 1) {

      formula <- paste("(", paste(colnames(solution)[colnames(solution) %in% colnames(x)], collapse= "+"),")*(", paste(colnames(x),collapse = "+"),")")

      if(sum(!{colnames(solution) %in% colnames(x)}) > 0) {
        tmp <- Filter(function(x) {length(x) == 2}, strsplit(colnames(solution), "[.]|[:]"))
        tmp <- Map(function(x) {paste0(x, collapse=":")}, Filter(function(x){length(x)==2}, tmp))
        if(length(tmp) > 0) {formula <- paste(formula,paste("(",paste(unlist(tmp),collapse="+"),")*(",paste(colnames(x), collapse = "+"),")"), sep="+")}
      }

      interactions <- data.frame(model.matrix(as.formula(paste("~", formula)), x))
      interactions <- interactions[ , setdiff(colnames(interactions), colnames(solution)), drop = F]
      interactions <- interactions[ , setdiff(colnames(interactions), colnames(candidate)), drop = F]
      candidate <- cbind(candidate, interactions)

    }

    bic[step] <- log(as.numeric(min(rss))/n) + (dim(solution)[2]) * (log(n) + 2 * log(p))/n
    if((dim(solution)[2]) > (n/log(n, 2)))break
    step <- step + 1

  }

  model <- data.frame(solution[ , 1 : which.min(bic), drop = FALSE])
  summary(lm(y ~ . , data = model))

}


# ff <- log(Volume) ~ log(Height) + log(Girth)
# utils::str(m <- model.frame(ff, trees))
# mat <- model.matrix(ff, m)
#
# dd <- data.frame(a = gl(3,4), b = gl(4,1,12)) # balanced 2-way
# options("contrasts")
# model.matrix(~ a + b, dd)
# model.matrix(~ a + b, dd, contrasts = list(a = "contr.sum"))
# model.matrix(~ a + b, dd, contrasts = list(a = "contr.sum", b = "contr.poly"))
# m.orth <- model.matrix(~a+b, dd, contrasts = list(a = "contr.helmert"))
# crossprod(m.orth) # m.orth is  ALMOST  orthogonal



# lm function printed
function(formula, data, subset, weights, na.action, method = "qr",
          model = TRUE, x = FALSE, y = FALSE, qr = TRUE, singular.ok = TRUE,
          contrasts = NULL, offset, ...){
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
    z <- list(coefficients = if (is.matrix(y)) matrix(, 0,
                                                      3) else numeric(), residuals = y, fitted.values = 0 *
                y, weights = w, rank = 0L, df.residual = if (!is.null(w)) sum(w !=
                                                                                0) else if (is.matrix(y)) nrow(y) else length(y))
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



