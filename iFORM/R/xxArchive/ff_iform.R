library(foreach)
library(doParallel)
library(iterators)
library(ff)
library(data.table)
library(microbenchmark)
rm(list=ls())


solution <- matrix(rbinom(1000, 1, 0.5), nrow = 200)
candidate2 <- candidate <- matrix(rbinom(1000000, 1, 0.5), nrow = 200)
candidate2 <- as.ff(candidate2)
y <- rnorm(200, 0, 1)

system.time({
  rss_sapply <- sapply(1 : 1 : ncol(candidate), function(i){
    xx <- as.matrix(cbind(solution, candidate[ , i]))
    tryCatch({
      sum((y - xx %*% (solve(t(xx) %*% xx)) %*% (t(xx) %*% y)) ^ 2)
    }, error = function(e) NULL)
  })
})


cl <- makeCluster(4)
registerDoParallel(cl)

system.time({
  rss_foreach <- foreach(i = iter(candidate),
                         .combine = "c",
                         .errorhandling = "pass") %dopar% {
    xx <- as.matrix(cbind(solution, i))
    sum((y - xx %*% (solve(t(xx) %*% xx)) %*% (t(xx) %*% y)) ^ 2)
  }
})

stopCluster(cl)

setequal(rss_sapply, rss_foreach)


rss_sapply_func <- function(solution, candidate, y) {
  sapply(1 : 1 : ncol(candidate), function(i){
    xx <- as.matrix(cbind(solution, candidate[ , i]))
    tryCatch({
      sum((y - xx %*% (solve(t(xx) %*% xx)) %*% (t(xx) %*% y)) ^ 2)
    }, error = function(e) Inf)
  })
}

rss_par_sapply_func <- function(solution, candidate, y, no_cores = 1) {
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)
  clusterExport(cl, c("solution", "candidate", "y"))

  parSapply(cl, 1 : 1 : ncol(candidate), function(i){
    xx <- as.matrix(cbind(solution, candidate[ , i]))
    tryCatch({
      sum((y - xx %*% (solve(t(xx) %*% xx)) %*% (t(xx) %*% y)) ^ 2)
    }, error = function(e) NULL)
  })

  stopCluster(cl)
}


rss_foreach_func <- function(solution, candidate, y, no_cores = 1) {
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)

  foreach(i = 1 : ncol(candidate),
          .combine = "c",
          .errorhandling = "pass") %dopar% {
            xx <- as.matrix(cbind(solution, candidate[ , i]))
            sum((y - xx %*% (solve(t(xx) %*% xx)) %*% (t(xx) %*% y)) ^ 2)
          }
  stopCluster(cl)
  gc(reset = TRUE)
}


system.time({
  rss_sapply_1 <- rss_sapply_func(solution = solution, candidate = candidate, y = y)
})

system.time({
  rss_sapply_2 <- rss_sapply_func(solution = solution, candidate = candidate2, y = y)
})

system.time({
  rss_foreach_1 <- rss_foreach_func(solution = solution, candidate = candidate, y = y)
})

system.time({
  rss_foreach_2 <- rss_foreach_func(solution = solution, candidate = candidate2, y = y)
})


system.time({
  rss_foreach_1_par <- rss_foreach_func(solution = solution, candidate = candidate, y = y, no_cores = 6)
})

system.time({
  rss_foreach_2_par <- rss_foreach_func(solution = solution, candidate = candidate2, y = y, no_cores = 6)
})

microbenchmark(list(rss_sapply_func(solution, candidate, y), rss_sapply_func(solution, candidate2,y)))


dim(candidate2) <- c(200, 1000000)
for( i in 1 : 100 ) {
  candidate2[, (10000 * (i - 1) + 1) : (10000 * i)] <- matrix(rbinom(2000000, 1, 0.5), nrow = 200)
}

system.time({
  sapply_output <- rss_sapply_func(solution = solution, candidate = candidate2, y = y)
})

system.time({
  par_sapply_output <- rss_par_sapply_func(solution = solution, candidate = candidate2, y = y, no_cores = 4)
})

system.time({
  foreach_output <- rss_foreach_func(solution = solution, candidate = candidate2, y = y, no_cores = 6)
})

## iForm w/ ff #####


iForm <- function(formula, data, strong = TRUE, higher.order = FALSE){

  dat <- model.frame(formula, data)
  y <- dat[ , 1]
  x <- dat[ , -1]
  p <- ncol(x)
  n <- nrow(x)
  candidate <- as.ff(as.matrix(x), colnames = colnames(x), overwrite = TRUE)
  solution <- NULL
  model <- NULL
  step <- 1
  bic <- NULL

  fit <- iformselect(x, y, p, n, candidate, solution, model, bic, step, strong, higher.order)
  y <- fit$y
  solution <- fit$solution
  model <- fit$model
  bic <- fit$bic


  model <- data.frame(solution[ , 1 : which.min(bic), drop = FALSE])
  lm(y ~ . + 0 , data = model)

}



iformselect <- function(x, y, p, n, candidate, solution, model, bic, step, strong, higher.order){



  repeat{

    rss <- sapply(colnames(candidate), function(i){
      xx <- as.matrix(cbind(solution, candidate[ , i]))
      tryCatch({
        sum((y - xx %*% (solve(t(xx) %*% xx)) %*% (t(xx) %*% y)) ^ 2)
      }, error = function(e) NULL)
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
        interaction.formula <- paste("~0+(", paste(colnames(solution)[colnames(solution) %in% colnames(x)], collapse= "+"),")^2")
      }else{
        interaction.formula <- paste("~0+(", paste(colnames(solution)[colnames(solution) %in% colnames(x)], collapse= "+"),")*(", paste(colnames(x),collapse = "+"),")")
      }

      if(higher.order){
        if(sum(!{colnames(solution) %in% colnames(x)}) > 2) {

          if(strong){
            tmp <- Filter(function(x) {length(x) == 2}, strsplit(colnames(solution), "[.]|[:]"))
            tmp <- combn(length(tmp), 3, function(x) Reduce(intersect, combn(x, 2, function(y) Reduce(union, tmp[y]), simplify=F)), simplify=F)
            tmp <- Map(function(x) {paste0(x, collapse=":")}, Filter(function(x){length(x)==3}, tmp))
            if(length(tmp) > 0) {interaction.formula <- paste(interaction.formula, paste0(unlist(tmp), collapse = "+"), sep = "+")}

          }else{
            tmp <- Filter(function(x) {length(x) == 2}, strsplit(colnames(solution), "[.]|[:]"))
            tmp <- Map(function(x) {paste0(x, collapse=":")}, Filter(function(x){length(x)==2}, tmp))
            if(length(tmp) > 0) {interaction.formula <- paste(interaction.formula, paste("(",paste(unlist(tmp), collapse="+"), ")*(", paste(colnames(x), collapse = "+"), ")"), sep = "+")}
          }

        }
      }

      interaction.formula <- as.formula(interaction.formula)
      interactions <- model.matrix(interaction.formula, x)
      interactions <- interactions[ , setdiff(colnames(interactions), colnames(solution)), drop = F]
      interactions <- interactions[ , setdiff(colnames(interactions), colnames(candidate)), drop = F]

      if( ncol(interactions) > 0 ){
        newdimstart <- ncol(candidate) + 1
        newdimend <- ncol(candidate) + ncol(interactions)
        dim(candidate) <- c(n, newdimend)
        candidate[, newdimstart : newdimend] <- interactions  # breaking point
        colnames(candidate) <- c(cnames, colnames(interactions))
      }

    }

    bic[step] <- log(as.numeric(min(rss))/n) + (dim(solution)[2]) * (log(n) + 2 * log(p))/n
    if((dim(solution)[2]) > (n / log(n, 2))) break
    step <- step + 1

  }

  list(y = y, solution = solution, model = model, bic = bic)
}



for(input in inputs) {
  tryCatch(print(paste("log of", input, "=", log(input))),
    warning = function(w) {print(paste("negative argument", input));
   log(-input)},
  error = function(e) {print(paste("non-numeric argument", input));
  NaN})
}



