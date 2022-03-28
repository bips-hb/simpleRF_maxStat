
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
      nrow <<- nrow(data)
      names <<- colnames(data)
    },
    
    column = function(col) {
      return(data[, col])   
    }, 
    
    subset = function(row, col) {
      return(data[row, col])
    })
)