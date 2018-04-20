library(data.table)
library(parallel)

sim.data <-lapply(1:100, function(i){
data <- data.frame(matrix(rbinom(100000, 1, 0.5), nrow = 200) - 0.5)
names(data) <- paste0("X", 1 : 500) 
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
data
})

ptm <- proc.time()
iform.output <- Map(function(x) {iForm.model.matrix.flex(x, "response", three=T)}, sim.data)
proc.time() - ptm


iform.output[[1]]$"coefficients"
iform.output[[1]]$"adj.r.squared"
iform.output[[1]]$"df"[1] # model size
iform.output[[1]]$"sigma" # RSE

sum(rownames(iform.output[[1]]$"coefficients")%like%"X1.X8.X4|X1.X4.X8|X4.X1.X8|X4.X8.X1|X8.X1.X4|X8.X4.X1")

sum(rownames(iform.output[[1]]$"coefficients")%like%"X1|X4|X8|X11")

sum(rownames(iform.output[[1]]$"coefficients")%like%"X1.X8|X1.X4|X4.X1|X4.X8|X8.X1|X8.X4")

Filter(function(x) length(x)==1, strsplit(rownames(iform.output[[1]]$"coefficients"),"[.]"))%like%"X1|X4|X8|X9|X11"

iform.coef.summary <- Map(function(x){
c(sum(rownames(x$"coefficients")%in%"X1|X8|X4"),
  sum(rownames(x$"coefficients")%like%"X1.X8|X1.X4|X4.X1|X4.X8|X8.X1|X8.X4"),
  sum(rownames(x$"coefficients")%like%"X1.X8.X4|X1.X4.X8|X4.X1.X8|X4.X8.X1|X8.X1.X4|X8.X4.X1"))

iform.summary <- Map(function(x){
  c(length(Filter(function(y) length(y)==1, strsplit(rownames(x$"coefficients"),"[.]"))),
    sum(Filter(function(y) length(y)==1, strsplit(rownames(x$"coefficients"),"[.]"))%like%"X1|X4|X8|X9|X11"),
    length(Filter(function(y) length(y)==2, strsplit(rownames(x$"coefficients"),"[.]"))),
    sum(Filter(function(y) length(y)==2, strsplit(rownames(x$"coefficients"),"[.]"))%like%"X1.X8|X1.X4|X4.X1|X4.X8|X8.X1|X8.X4"),    
    length(Filter(function(y) length(y)==3, strsplit(rownames(x$"coefficients"),"[.]"))),
    sum(rownames(x$"coefficients")%like%"X1.X8.X4|X1.X4.X8|X4.X1.X8|X4.X8.X1|X8.X1.X4|X8.X4.X1"),     
  x$"adj.r.squared",
  x$"df"[1], # model size
  x$"sigma") # RSE  
}, iform.output)

iform.summary.df <- do.call(rbind,iform.summary)
colnames(iform.summary.df) <- c("Adj.R.Squared", "Size", "RSE")






## Running 3-way Interaction in Parallel ################################

setwd("~/School/PhD/Dissertation")

source("iForm.R")
library(snow)

gene.expr<-t(read.table("Data Files/M_quantile.txt", header=TRUE))
gene.expr <- as.data.frame(gene.expr)
gene.type<-read.table("Data Files/genotype.txt", header=TRUE)

n=dim(gene.type)[1]
candidate<-data.frame(gene.type)
j<-1
repeat{
  prop.alike<-NULL
  for(i in 1:(dim(candidate)[2])){
    prop.alike[i]<-sum(candidate[,j]==candidate[,i])/n
  }
  same<-which(prop.alike==1)
  same<-same[same!=j]
  if(length(same>1))candidate<-candidate[,-same]
  j<-j+1
  if(j>=(dim(candidate)[2]))break
}

gene.type.nodups<-candidate
gene.type.nodups[gene.type.nodups!=0]<-2
gene.type.nodups[gene.type.nodups==0]<-1
gene.type.nodups <- gene.type.nodups - 1.5


  # Calculate the number of cores
  no_cores <- 6
  
  # Initiate cluster
  cl <- makeCluster(no_cores, type="FORK")

ptm <- proc.time()

output.list <- parLapply(cl, 3001:5000, function(i){
  data <- data.frame(response = gene.expr[,i], gene.type.nodups)
iForm.model.matrix.flex(data, "response", three = TRUE)
})

proc.time() - ptm

stopCluster(cl)

saveRDS(output.list, "celegans.iform.3way.rds")


output.list <- c(output.list, output.list2)







## Tree Marker Data ##############################################


# Developing weak heredity version
# scenario 1: Weak iForm two-way interaction with intercept
# Two-interaction needs to be included in model for three-way interactions

iForm.weak.three.intercept <- function(data, response){
  
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
  lm(y ~ . , data = model)
  
}


# checking to see if intercept is included
fit_iform <- iForm.weak.three.intercept(sim.data[[1]],"response")

## iForm w/ ff and intercept ########

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
  lm(y ~ ., data = model)
  
}



iformselect <- function(x, y, p, n, candidate, solution, model, bic, step, strong, higher.order){
  
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


## working with tree data ################


library(data.table)
library(softImpute)
library(mice) 

dia <- fread("L-abr-1_DIA.csv")
setnames(dia, names(dia)[1], "Tree")

  # adding growth curves
for( t in 1 : 16 ) {
  eval(parse(text = paste0("dia[, g_", t, ":= a/(1+b*exp(-r * " , t, "))]")))
}

ht <- fread("L-abr-1_HT.csv")
setnames(ht, names(ht)[1], "Tree")

  # adding growth curves
for( t in 1 : 16 ) {
  eval(parse(text = paste0("ht[, g_", t, ":= a/(1+b*exp(-r * " , t, "))]")))
}


impute_two_levels <- fread("markers_impute_two_levels.csv")
impute_two_levels <- impute_two_levels[, -1, with = F]
impute_three_levels <- fread("markers_impute_three_levels.csv")
impute_three_levels <- impute_three_levels[, -1, with = F]
markers_impute <- cbind(impute_two_levels, impute_three_levels)

markers.t <- fread("L-marker.csv")
markers <- t(markers.t)
colnames(markers) <- markers[3,]
markers <- as.data.table(markers[4:207,])
setnames(markers, names(markers)[1], "Tree")
# markers[,paste(names(markers)[-1]) := lapply(.SD, function(x){
#   ifelse(x == "--", NA, x)
# }), .SDcols = names(markers)[-1]]
# 
# markers[,paste(names(markers)[-1]) := lapply(.SD, function(x) factor(x)), 
#         .SDcols = names(markers)[-1]]

markers <- data.frame(Tree = markers$Tree, markers_impute)
# removing duplicate markers

data.t <- t(markers)
data.t <- unique(data.t)
markers <- as.data.table(t(data.t))

markers[,paste(names(markers)[-1]) := lapply(.SD, function(x) factor(x)), 
        .SDcols = names(markers)[-1]]

markers[,paste(names(markers)[-1]) := lapply(.SD, function(x) as.numeric(x)), 
        .SDcols = names(markers)[-1]]


  # by subtracting 2 and multiplying by negative 1 we get the coding for the addtive effects
markers[,paste(names(markers)[-1]) := lapply(.SD, function(x) (x - 2)*(-1)),
        .SDcols = names(markers)[-1]]

# na.rows <- apply(markers, 1, function(x) mean(is.na(x)))
# na.rows <- which(na.rows > 0.15)
# 
# markers <-markers[-na.rows,]
# 
#   # checking for NA columns less than 5%
# na.check <- markers[,lapply(.SD, function(x) mean(is.na(x))),.SDcols=names(markers)[-1]]
# low.na <- which(unlist(na.check) < 0.05)
# 
# markers <- markers[,.SD,.SDcols = c("Tree", names(low.na))]

  # Search for 3 levels and then create dominant marker

three.cols <- which(unlist(markers[,lapply(.SD, function(x) length(unique(na.omit(x)))==3),.SDcols = names(markers)])==TRUE)

# need to code additive only markers as (1, 0)
# need to code addtive-dominant markers as additive: (1, 0 , -1) 
# and dominant (0,1,0)

dom.effects <- markers[,three.cols, with = FALSE]
dom.effects[,paste(names(dom.effects)) := lapply(.SD, function(x) (abs(x) - 1)*(-1)), .SDcols = names(dom.effects)]
setnames(dom.effects, names(dom.effects), paste0(names(dom.effects), "_d"))

setnames(markers, names(markers)[-1], paste0(names(markers)[-1], "_a"))

markers <- data.frame(markers, dom.effects)
dia[, Tree := as.character(Tree)]

data.mrg <- merge(dia, markers, by = "Tree")
setnames(data.mrg, c("a", "b", "r"), c("dia_a", "dia_b", "dia_r"))

ht[, Tree := as.character(Tree)]
data.mrg <- merge(ht, data.mrg, by = "Tree")
setnames(data.mrg, c("a", "b", "r"), c("ht_a", "ht_b", "ht_r"))


ptm <- proc.time()

data1 <- as.data.frame(data.mrg[,.SD, .SDcols = c("dia_a", colnames(markers)[-1])])
dia.a.output <- iForm(dia_a ~ ., data = data1, strong = FALSE, higher.order = TRUE)

proc.time() - ptm



rm(data1); gc(reset = TRUE)
data1 <- as.data.frame(data.mrg[,.SD, .SDcols = c("dia_b", names(markers)[-1])])
dia.b.output <- iForm(dia_b ~ ., data = data1, strong = FALSE, higher.order = TRUE)

rm(data1); gc(reset = TRUE)
data1 <- as.data.frame(data.mrg[,.SD, .SDcols = c("dia_r", names(markers)[-1])])
dia.r.output <- iForm(dia_r ~ ., data = data1, strong = FALSE, higher.order = TRUE)

rm(data1); gc(reset = TRUE)
data1 <- as.data.frame(data.mrg[,.SD, .SDcols = c("ht_a", names(markers)[-1])])
ht.a.output <- iForm(ht_a ~ ., data = data1, strong = FALSE, higher.order = TRUE)

rm(data1); gc(reset = TRUE)
data1 <- as.data.frame(data.mrg[,.SD, .SDcols = c("ht_b", names(markers)[-1])])
ht.b.output <- iForm(ht_b ~ ., data = data1, strong = FALSE, higher.order = TRUE)

rm(data1); gc(reset = TRUE)
data1 <- as.data.frame(data.mrg[,.SD, .SDcols = c("ht_r", names(markers)[-1])])
ht.r.output <- iForm(ht_r ~ ., data = data1, strong = FALSE, higher.order = TRUE)



output.list <- list(dia.a.output, dia.b.output, dia.r.output,
                    ht.a.output, ht.b.output, ht.r.output)

saveRDS(output.list, "Mei_Tree.iFormThreeWayOutput_2.rds")


ht.r.output.2way <- iForm.weak(data, "ht_r")



iForm.weak1.three.intercept(data, "dia.a")


output <- lapply(c("dia.a", "dia.b", "dia.r", "ht.a", "ht.b", "ht.r"), function(x){
 data <- as.data.frame(data.mrg[,.SD, .SDcols = c(x, names(markers)[-1])])
 iForm.weak1.three.intercept(data, x)
})



mei <- readRDS("Mei_Tree.iFormThreeWayOutput.rds")
celegans <- readRDS("iform.celegans.output.3way.rds")




data.mrg <- fread("MeiTrees_imputed_merged.csv")
data.mrg[,lapply(.SD, mean), .SD = c("ht_a", "ht_b", "ht_r")]

data.mrg[,lapply(.SD, mean), .SD = c("ht_a", "ht_b", "ht_r"), by = AATTC_nn_np_2815_a]
data.mrg[,lapply(.SD, mean), .SD = c("ht_a", "ht_b", "ht_r"), by = AATTC_lm_ll_3034_a]
data.mrg[,lapply(.SD, mean), .SD = c("ht_a", "ht_b", "ht_r"), by = AATTC_nn_np_1615_a]

data.mrg[,lapply(.SD, mean), .SD = c("ht_a", "ht_b", "ht_r"), by = .(AATTC_nn_np_2815_a,
                                                                     AATTC_lm_ll_3034_a,
                                                                     AATTC_nn_np_1615_a)]


## Imputation Attempt #####

 ## Helper Functions ####
impute_na <- function(data, g_dist, obs) {
  na.locf(data[data$`Genetic_Distances(cM)` == g_dist, obs], na.rm = FALSE, fromLast = TRUE) %>%
    na.locf(., na.rm = FALSE, fromLast = FALSE)
}


na_impute_data <- function( data ) {
  non_imputed_data <- data
  imputed_data <- data
  
  for( i in unique(data$`Genetic_Distances(cM)`) ) {
    for( j in 1 : 207) {
      imputed_data[imputed_data$`Genetic_Distances(cM)` == i, j] <- impute_na(imputed_data, i, j)
    }
  }
  
  change <- which(apply(imputed_data, 1, function(x) {
    sum(!{unique(na.omit(x[4:207])) %in% unlist(strsplit(x[3], "_"))[2:3]})
  }) > 0)
  
  imputed_data[change, ] <- non_imputed_data[change, ]
  
  list(non_imputed_data, imputed_data)
}


library(data.table)
library(linkim)
library(zoo)
library(magrittr)

markers.t <- fread("L-marker.csv")
setnames(markers.t, names(markers.t), as.character(markers.t[1,]))
markers.t <- markers.t[-1, ]
markers.t[,paste(names(markers.t)[-1]) := lapply(.SD, function(x){
  ifelse(x == "--", NA, x)
}), .SDcols = names(markers.t)[-1]]

markers.t_df <- as.data.frame(markers.t)
markers.t <- na_impute_data(markers.t_df)

markers <- t(markers.t[[2]]) # the imputed dataset
colnames(markers) <- markers[3,]
r <- markers[2, ]
markers <- as.data.table(markers[4:207,], keep.rownames = TRUE)
setnames(markers, names(markers)[1], "Tree")
markers[,paste(names(markers)[-1]) := lapply(.SD, function(x){
  ifelse(x == "--", NA, x)
}), .SDcols = names(markers)[-1]]

markers[,paste(names(markers)[-1]) := lapply(.SD, function(x) factor(x)), 
        .SDcols = names(markers)[-1]]

data.t <- t(markers)
data.t <- unique(data.t)
markers <- as.data.table(t(data.t))

markers[,paste(names(markers)[-1]) := lapply(.SD, function(x) factor(x)), 
        .SDcols = names(markers)[-1]]

two_levels <- which(sapply(markers,function(x) length(levels(x)))==2)
three_levels <- which(sapply(markers,function(x) length(levels(x)))==3)

markers_two <- as.data.frame(markers[,lapply(.SD,as.numeric),.SDcols = two_levels]) - 1
markers_three <- as.data.frame(markers[,.SD,.SDcols = three_levels])

markers_new_two <- link.im(markers_two, as.numeric(r[two_levels]))


# removing duplicate markers


markers[,paste(names(markers)[-1]) := lapply(.SD, function(x) factor(x)), 
        .SDcols = names(markers)[-1]]




## example
library(linkim) 

data(barley)
dat <- barley[,-1]
r <- as.numeric(dat[1,])
data <- dat[-1,]
new.data <- link.im(data,r)



## reading imputed data


