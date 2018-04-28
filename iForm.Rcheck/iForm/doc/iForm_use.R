## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
library(iForm)
data("mtcars")
names(mtcars)  # Variable names
dim(mtcars)  # dimensions of the dataset
sum(is.na(mtcars$hp)) # checking number of NA values
help("mtcars") # more info on mtcars

  # Running iForm on R dataset mtcars using hourse power as response (hp)
iForm.fit1 <- iForm(hp ~ ., mtcars)
iForm.fit1

## ------------------------------------------------------------------------
  # fitting the iForm procedure with the weak heredity assumption
iForm.weak.fit1 <- iForm(hp ~ ., mtcars, heredity = "weak")
iForm.weak.fit1

## ------------------------------------------------------------------------
library(ISLR)
data("Hitters")
names(Hitters)  # Variable names
dim(Hitters)  # dimensions of the dataset
sum(is.na(Hitters$Salary))  # number of NA values in the Salary variable
Hitters <- na.omit(Hitters)  # removing all na values
dim(Hitters)  # rechecking dimension size after removal of NA's
sum(is.na(Hitters))   # rechecking number of NA values

## ------------------------------------------------------------------------
help("Hitters") # more information on Hitters dataset

iForm.fit2 <- iForm(Salary ~ ., Hitters)
iForm.fit2

## ------------------------------------------------------------------------
help("Hitters") # more information on Hitters dataset

iForm.fit2_weak <- iForm(Salary ~ ., Hitters, heredity = "weak")
iForm.fit2_weak

