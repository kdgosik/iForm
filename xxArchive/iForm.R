## iForm Using Model.Matrix Code ###########################

iForm.model.matrix <- function(data, response){
  
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
      interactions <- data.frame(model.matrix(as.formula(paste("~0+(", paste(colnames(solution)[colnames(solution) %in% colnames(x)], collapse= "+"),")^2")), x))
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




## iForm (weak heredity) Using Model.Matrix Code ###########################

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



## Three-way interaction with Model.Matrix #############################


iForm.model.matrix.three <- function(data, response){
  
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
      
      formula <- paste0("(", paste(colnames(solution)[colnames(solution) %in% colnames(x)], collapse= "+"), ")^2")
      
      if(sum(!{colnames(solution) %in% colnames(x)}) > 2) {
        
        tmp <- Filter(function(x) {length(x) == 2}, strsplit(colnames(solution), "[.]|[:]"))
        tmp <- combn(length(tmp), 3, function(x) Reduce(intersect, combn(x, 2, function(y) Reduce(union, tmp[y]), simplify=F)), simplify=F)
        tmp <- Map(function(x) {paste0(x, collapse=":")}, Filter(function(x){length(x)==3}, tmp))
        if(length(tmp) > 0) {formula <- paste(formula,paste0(unlist(tmp),collapse="+"), sep="+")}

      }
      
      interactions <- data.frame(model.matrix(as.formula(paste0("~0+", formula)), solution))
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


## iForm with choice of including three-way interactions #############################


iForm.model.matrix.flex <- function(data, response, three = FALSE){
  
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
      
      formula <- paste0("(", paste(colnames(solution)[colnames(solution) %in% colnames(x)], collapse= "+"), ")^2")
      
      if(three){
        if(sum(!{colnames(solution) %in% colnames(x)}) > 2) {
        
        tmp <- Filter(function(x) {length(x) == 2}, strsplit(colnames(solution), "[.]|[:]"))
        tmp <- combn(length(tmp), 3, function(x) Reduce(intersect, combn(x, 2, function(y) Reduce(union, tmp[y]), simplify=F)), simplify=F)
        tmp <- Map(function(x) {paste0(x, collapse=":")}, Filter(function(x){length(x)==3}, tmp))
        if(length(tmp) > 0) {formula <- paste(formula,paste0(unlist(tmp),collapse="+"), sep="+")}
        
        }
      }
      
      interactions <- data.frame(model.matrix(as.formula(paste0("~0+", formula)), solution))
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



## iFort Using Model.Matrix Code ###########################

iFort.model.matrix <- function(data, response){
  
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
    
    bic[step] <- log(as.numeric(min(rss))/n) + (dim(solution)[2]) * (log(n) + 2 * log(p))/n
    if((dim(solution)[2]) > (n/log(n, 2)))break
    step <- step + 1
    
  }
  
  p2 <- which.min(bic)
  solution <- solution[,1:p2]
  candidate <- model.matrix(as.formula(paste0("~0+", "(", paste0(colnames(solution), collapse="+"), ")^2" )), solution)
  candidate <- candidate[ , setdiff(colnames(candidate), colnames(solution)), drop = F]
  step <- p2 + 1
  
  repeat{
    
    rss <- sapply(colnames(candidate), function(i){
      xx <- as.matrix(cbind(solution, candidate[ , i]))
      try(sum((y - xx %*% (solve(t(xx) %*% xx)) %*% (t(xx) %*% y)) ^ 2), silent = T)
    })
    
    solution <- data.frame(cbind(solution, candidate[ , which.min(rss)]))
    colnames(solution)[step] <- colnames(candidate)[which.min(rss)]
    candidate <- candidate[ , -which.min(rss)]
    
    bic[step] <- log(as.numeric(min(rss))/n) + (dim(solution)[2]) * (log(n) + 2 * log(p))/n
    if((dim(solution)[2]) > (n/log(n, 2)))break
    step <- step + 1
    
  }
  
  model <- data.frame(solution[, 1 : which.min(bic), drop = FALSE])
  summary(lm(y ~ . + 0 , data = model))
  
}




## iFort Using Model.Matrix Code ###########################

iFort.model.matrix.flex <- function(data, response, three=F){
  
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
    candidate <- candidate[ , -which.min(rss), drop = F]
    
    bic[step] <- log(as.numeric(min(rss))/n) + (dim(solution)[2]) * (log(n) + 2 * log(p))/n
    if((dim(solution)[2]) > (n/log(n, 2)))break
    step <- step + 1
    
  }
  
  p2 <- which.min(bic)
  solution <- solution[ , 1 : p2]
  candidate <- model.matrix(as.formula(paste0("~0+", "(", paste0(colnames(solution), collapse="+"), ")^2" )), solution)
  candidate <- candidate[ , setdiff(colnames(candidate), colnames(solution)), drop = F]
  step <- p2 + 1
  
  repeat{
    
    rss <- sapply(colnames(candidate), function(i){
      xx <- as.matrix(cbind(solution, candidate[ , i]))
      try(sum((y - xx %*% (solve(t(xx) %*% xx)) %*% (t(xx) %*% y)) ^ 2), silent = T)
    })
    
    solution <- data.frame(cbind(solution, candidate[ , which.min(rss)]))
    colnames(solution)[step] <- colnames(candidate)[which.min(rss)]
    candidate <- candidate[ , -which.min(rss), drop = F]
    
    bic[step] <- log(as.numeric(min(rss))/n) + (dim(solution)[2]) * (log(n) + 2 * log(p))/n
    if((dim(solution)[2]) > (n/log(n, 2)))break
    step <- step + 1
    
  }

  if(three){
    
  p3 <- which.min(bic)
  solution <- solution[ , 1 : p3]
  tmp <- Filter(function(x) {length(x) == 2}, strsplit(colnames(solution), "[.]|[:]"))
  tmp <- combn(length(tmp), 3, function(x) Reduce(intersect, combn(x, 2, function(y) Reduce(union, tmp[y]), simplify=F)), simplify=F)
  tmp <- Map(function(x) {paste0(x, collapse=":")}, Filter(function(x){length(x)==3}, tmp))
  if(length(tmp) > 0) {formula <- paste0(unlist(tmp), collapse = "+")}
  
  candidate <- model.matrix(as.formula(paste0("~0+", formula)), solution)
  step <- p3 + 1
  
  repeat{
    
    rss <- sapply(colnames(candidate), function(i){
      xx <- as.matrix(cbind(solution, candidate[ , i]))
      try(sum((y - xx %*% (solve(t(xx) %*% xx)) %*% (t(xx) %*% y)) ^ 2), silent = T)
    })
    
    solution <- data.frame(cbind(solution, candidate[ , which.min(rss)]))
    colnames(solution)[step] <- colnames(candidate)[which.min(rss)]
    candidate <- candidate[ , -which.min(rss), drop = F]
    
    bic[step] <- log(as.numeric(min(rss))/n) + (dim(solution)[2]) * (log(n) + 2 * log(p))/n
    if((dim(solution)[2]) > (n/log(n, 2)))break
    if(dim(candidate)[2] == 0)break
    step <- step + 1
    
    }
  }
  
  model <- data.frame(solution[, 1 : which.min(bic), drop = FALSE])
  summary(lm(y ~ . + 0 , data = model))
  
}







