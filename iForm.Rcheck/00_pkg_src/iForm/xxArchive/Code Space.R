data <- data.frame(matrix(rbinom(10000, 1, 0.5), nrow = 200) - 0.5)
names(data) <- paste0("X", 1 : 50) 
response <- 3 * data[ , 1] + 
            3 * data[ , 4] + 
            3 * data[ , 8] + 
            3 * data[ , 9] + 
            3 * data[ , 11] + 
            3 * data[ , 1] * data[ , 4] + 
            3 * data[ , 1] * data[ , 8] + 
            3 * data[ , 4] * data[ , 8] + 
            3 * data[ , 1] * data[ , 4] * data[ , 8] + rnorm(100, 0, 1)
data <- data.frame(response, data)


## Practice Identyfing and Selecting Interacions #################

# selection%in%predictors
# !{selection%in%predictors}
# 
# selection[selection%in%predictors]
# selection[!{selection%in%predictors}]
# 
# unlist(strsplit(selection[!{selection%in%predictors}], ":"))
# unlist(strsplit(selection[!{selection%in%predictors}], ":")) %in% selection[selection%in%predictors]

# c("X6", "X7", "X8" "X6:X7", "X6:X8", "X7:X8")
# cbind(do.call(rbind,strsplit(selection,":")),selection)



predictors <- paste0("X",1:50)

interactions <- NULL
for(i in 1:length(predictors)){
  for(j in i:length(predictors)){
    interactions <- append(interactions, paste0(predictors[i], ":", predictors[j]))
  }
}

selection <- unique(c("X6", "X7", "X8", "X6:X7", "X6:X8", "X7:X8", 
                      sample(c(predictors,interactions),700)))

# This works
  # maybe try to use sapply function instead of a for loop?


threeway.int <- NULL
for(i in 1:length(selection[selection %in% predictors])){
  tmp <- selection[selection %in% predictors][i:(i+2)]
  need <- c(paste0(tmp[1],":",tmp[2]), paste0(tmp[1],":",tmp[3]), paste0(tmp[2],":",tmp[3]))
  if(sum(need%in%selection) == 3) 
  {threeway.int[i] <- paste(tmp[1],tmp[2],tmp[3],sep=":")}
}


threeway.int <- NULL
for(i in 1:length(selection[selection %in% predictors])){
  tmp <- selection[selection %in% predictors][i:(i+2)]
  need <- c(paste0(tmp[1],":",tmp[2]), paste0(tmp[1],":",tmp[3]), paste0(tmp[2],":",tmp[3]))
  if(sum(need%in%selection) == 3) 
  {threeway.int[i] <- paste(tmp[1],tmp[2],tmp[3],sep=":")}
}




x <- read.table(textConnection('
   V1 V2 V3 V4
1  9   25   18
2  5   20   10
3  4   30   12
4  4   34   16'
), header=TRUE)

model.matrix(~V1^2,x)

model.matrix(~(V1+V2+V3+V4)^3,x)

attempt <- list()
for(i in 2:4){
  attempt[[i]] <- model.matrix( ~ (names(x)[1:i])^2,x)
}


xnam <- paste0("x", 1:25)
(fmla <- as.formula(paste(" ~ 0 +", paste(xnam, collapse= "+"))))

xnam <- paste0("x", 1:25)
(fmla <- as.formula(paste(xnam, collapse= "+")))


## iForm Using Model.Matrix Code ###########################

####

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
    
    rss <- sapply(colnames(candidate),function(i){
      xx <- as.matrix(cbind(solution, candidate[ , i]))
      try(sum((y - xx %*% (solve(t(xx) %*% xx)) %*% (t(xx) %*% y)) ^ 2), silent=T)
    })
    
    solution <- data.frame(cbind(solution, candidate[ , which.min(rss)]))
    colnames(solution)[step] <- colnames(candidate)[which.min(rss)]
    candidate <- candidate[ , -which.min(rss)]
    if(dim(solution)[2] > 1) {
      interactions <- data.frame(model.matrix(as.formula(paste("~0+", paste("(", paste(colnames(solution)[colnames(solution) %in% colnames(x)], collapse= "+"),")^2"))),solution))
      interactions <- interactions[ , setdiff(colnames(interactions), colnames(solution)), drop = F]
      interactions <- interactions[ , setdiff(colnames(interactions), colnames(candidate)), drop = F]      
      candidate <- cbind(candidate, interactions)
    }
    
    bic[step] <- log(as.numeric(min(rss))/n) + (dim(solution)[2]) * (log(n) + 2 * log(p))/n
    if((dim(solution)[2]) > (n/log(n, 2)))break
    step <- step + 1
  
  }
  
  model <- data.frame(solution[,1:which.min(bic)])
  summary(lm(y ~ . + 0 , data = model))
  
}


## Three-way interaction with Model.Matrix #############################


iForm.model.matrix.three<-function(data, response){
  
  y <- data[ , "response"]
  x <- data[ , -which(colnames(data) == "response")]   
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
  
  model <- data.frame(solution[,1:which.min(bic)])
  summary(lm(y ~ . + 0 , data = model))
  
}




solution <- c("X1", "X2", "X3", "X4", "X5", "X1.X2", "X3.X1", "X2.X3", "X2.X4", "X4.X5")

solution.tmp<-strsplit(solution[grepl("[[:punct:]]",solution)],"[[:punct:]]")

combn(5, 3, function(x) Reduce(intersect, combn(x, 2, function(y) Reduce(union, solution.tmp[y]), simplify=F)), simplify=F)


tmp<-combn(5, 3, function(x) Reduce(intersect, combn(x, 2, function(y) Reduce(union, solution.tmp[y]), simplify=F)), simplify=F)


paste0(Filter(function(x){length(x)==3},tmp),collpase=":")
Map(function(x) {paste0(x, collapse=":")}, Filter(function(x){length(x)==3},tmp))


tmp <- Filter(function(x) {length(x) == 2}, strsplit(colnames(solution), "[.]|[:]"))
tmp <- combn(length(tmp), 3, function(x) Reduce(intersect, combn(x, 2, function(y) Reduce(union, tmp[y]), simplify=F)), simplify=F)
tmp <- Map(function(x) {paste0(x, collapse=":")}, Filter(function(x){length(x)==3},tmp))


solut <- c("X1", "X2", "X3", "X4", "X5", "X1.X2", "X3.X1", "X2.X3", "X2.X4", "X4.X5")
solut <- c("X1", "X2", "X3", "X4", "X5", "X1.X2", "X2.X3", "X2.X4", "X4.X5")
tmp <- Filter(function(x) {length(x) == 2}, strsplit(solut, "[.]|[:]"))
tmp <- combn(length(tmp), 3, function(x) Reduce(intersect, combn(x, 2, function(y) Reduce(union, tmp[y]), simplify=F)), simplify=F)
tmp <- Map(function(x) {paste0(x, collapse=":")}, Filter(function(x){length(x)==3},tmp))
if(length(tmp) > 0) {formula <- paste(formula,paste0(unlist(tmp),collapse="+"), sep="+")}
