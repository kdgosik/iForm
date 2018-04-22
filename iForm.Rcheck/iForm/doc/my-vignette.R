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
iForm.weak.fit1 <- iForm(hp ~ ., mtcars, strong = FALSE)
iForm.weak.fit1

## ------------------------------------------------------------------------
  # Using a slightly larger dataset from the package ISLR

  # If package isn't installed this code will install it

list.of.packages <- c("ISLR")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

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

