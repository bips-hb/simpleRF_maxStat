
##' @title Data class
##' @description Wrapper class for \code{data.frame}. 
##' Allows to access a \code{data.frame} by reference.
Data <- setRefClass("Data", 
  fields = list(
    data = "data.frame",
    glmdata = "data.frame",
    glmformula = "formula",
    ncol = "integer", 
    nrow = "integer",
    nconf = "integer",
    names = "character",
    confounders = "character"),  
  methods = list(
    
    initialize = function(...) {
      callSuper(...)
      ncol <<- ncol(data)
      nconf <<- as.integer(ncol(glmdata)-1)
      nrow <<- nrow(data)
      names <<- colnames(data)
      confounders <<- colnames(glmdata)[-1]
    },
    
    column = function(col) {
      return(data[, col])   
    }, 
    
    subset = function(row, col) {
      return(data[row, col])
    },
  
    glmcolumn = function(col) {
      return(glmdata[, col])   
    }, 
  
    glmsubset = function(row, col) {
     return(glmdata[row, col, drop=FALSE])
    })
  
  
)