## iForm with ff package #########################

mat <- matrix(rbinom(10000,1,.5), nrow = 200)
y <- mat[,c(1,2,4,5)] %*% rep(3,4) + 
  2*mat[,1]*mat[,2] +
  2*mat[,1]*mat[,4] +
  2*mat[,2]*mat[,5] +
  2*mat[,2]*mat[,4] +
  2*mat[,1]*mat[,2]*mat[,4] + 
  rnorm(200,0,1)

data <- data.frame(y, mat)
formula <- as.formula("y ~ 0 + .")


require(ff)


matff <- as.ff(mat)
dim(matff) <- c(nrow(mat), (ncol(mat) + ncol(newpreds)))
matff[, (ncol(mat) + 1) : (ncol(mat) + ncol(newpreds))] <- newpreds



iForm <- function(formula, data, strong = TRUE, higher.order = FALSE){
  
  dat <- model.frame(formula, data)
  y <- dat[ , 1]
  x <- dat[ , -1]
  p <- ncol(x)
  n <- nrow(x)
  candidate <- as.ff(as.matrix(x))
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





# doesn't work
  # colnames/dimnames are causing an issue

iformselect <- function(x, y, p, n, candidate, solution, model, bic, step, strong, higher.order){
  
  repeat{
    
    rss <- sapply(colnames(candidate), function(i){
      xx <- as.matrix(cbind(solution, candidate[ , i]))
      try(sum((y - xx %*% (solve(t(xx) %*% xx)) %*% (t(xx) %*% y)) ^ 2), silent = T)
    })
    
    solution <- cbind(solution, candidate[ , which.min(rss)])
    colnames(solution)[step] <- colnames(candidate)[which.min(rss)]
    
    # can't have negative subscripts in ff matrices
      # look into update.ff
    candidate <- update(candidate, from = candidate[, which.min(rss)], delete = TRUE)
    # candidate <- as.ff(candidate[ , -which.min(rss)])
    
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
      
      # this part breaks it
      dim(candidate) <- c(n, (ncol(candidate) + ncol(interactions)))
      candidate[, (ncol(candidate) + 1) : (ncol(candidate) + ncol(interactions))] <- interactions
      
      # 
      # candidate <- cbind(candidate, interactions)
    }
    
    bic[step] <- log(as.numeric(min(rss))/n) + (dim(solution)[2]) * (log(n) + 2 * log(p))/n
    if((dim(solution)[2]) > (n / log(n, 2))) break
    step <- step + 1
    
  }
  
  list(y = y, solution = solution, model = model, bic = bic)
}







## ff and foreach together #############################

# load package ####
library(ff)
library(doParallel)
library(foreach)
library(iterators)

mat <- matrix(rbinom(10000,1,.5), nrow = 200)
y <- mat[,c(1,2,4,5)] %*% rep(3,4) + 
  2*mat[,1]*mat[,2] +
  2*mat[,1]*mat[,4] +
  2*mat[,2]*mat[,5] +
  2*mat[,2]*mat[,4] +
  2*mat[,1]*mat[,2]*mat[,4] + 
  rnorm(200,0,1)

data <- data.frame(y, mat)
formula <- as.formula("y ~ 0 + .")

solution <- mat[,samp <- sample(ncol(mat), 10)]
candidate <- mat[,-samp]


  # just foeach
itx <- iter(candidate, by = "col")
rss <- foreach(i = itx, .combine = "c") %do% 
{
  xx <- cbind(solution, i)
  try(sum((y - xx %*% (solve(t(xx) %*% xx)) %*% (t(xx) %*% y)) ^ 2), silent = T)
}


  # foreach in parallel
cl <- makeCluster(2)
registerDoParallel(cl)
itx <- iter(candidate, by = "col")
rss <- foreach(i = itx, .combine = "c") %dopar% 
{
  xx <- cbind(solution, i)
  try(sum((y - xx %*% (solve(t(xx) %*% xx)) %*% (t(xx) %*% y)) ^ 2), silent = T)
}
stopCluster(cl)



  # foreach in parallel with ff package

candidate <- as.ff(candidate)

cl <- makeCluster(2)
registerDoParallel(cl)
# itx <- iter(candidate, by = "col")
rss <- foreach(i = 1:ncol(candidate), .combine = "c", .packages = "ff") %dopar% 
{
  xx <- cbind(solution, candidate[,i])
  try(sum((y - xx %*% (solve(t(xx) %*% xx)) %*% (t(xx) %*% y)) ^ 2), silent = T)
}
stopCluster(cl)

