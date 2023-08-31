
##' @title Forest class
##' @description Virtual class for Random forest. 
##' Contains all fields and methods used in all Forest subclasses.
##' @importFrom parallel mclapply
##' @import methods
Forest <- setRefClass("Forest", 
  fields = list(
    num_trees = "integer", 
    mtry = "integer", 
    minsplit = "integer", 
    minbucket = "integer",
    splitrule = "character",
    always_split_varIDs = "integer",
    unordered_factors = "character",
    data = "Data",
    predict_data = "Data",
    formula = "formula",
    trees = "list",
    treetype = "character",
    replace = "logical", 
    predleaf = "character",
    maxstat = "logical",
    minprop = "numeric",
    alpha = "numeric",
    pmethod = "character",
    covariate_levels = "list"),
  methods = list(
    
    grow = function(num_threads) { 
      
      ## Init trees
      temp <- lapply(trees, function(x) {
        x$mtry <- mtry
        x$minsplit <- minsplit
        x$minbucket <- minbucket
        x$splitrule <- splitrule
        x$always_split_varIDs <- always_split_varIDs
        x$unordered_factors <- unordered_factors
        x$data <- data
        x$predleaf <- predleaf
        x$maxstat <- maxstat
        x$minprop <- minprop
        x$alpha <- alpha
        x$pmethod <- pmethod
      })
      
      ## Grow trees
      trees <<- mclapply(trees, function(x) {
        x$grow(replace)
        x
      }, mc.cores = num_threads)
    }, 
    
    predict = function(newdata) {
      model.data <- model.frame(formula, newdata)
      glm.data <- model.frame(data$glmformula, newdata)
      
      #remove confounders from model.data if necessary
      model.data <- model.data[,which(colnames(model.data) %in% data$names)]

      ## Recode factors if forest grown 'order_once' mode
      if (unordered_factors == "order_once" & length(covariate_levels) > 0) {
        model.data[, -1] <- mapply(function(x, y) {
          if(is.null(y)) {
            x
          } else {
            new.levels <- setdiff(levels(x), y)
            factor(x, levels = c(y, new.levels), ordered = TRUE)
          }
        }, model.data[, -1], covariate_levels, SIMPLIFY = FALSE)
      }

      ## Save prediction data in model
      predict_data <<- Data$new(data = model.data, glmdata=glm.data)
      
      ## Predict in trees
      predictions <- simplify2array(lapply(trees, function(x) {
        x$predict(predict_data)
      }))
      
      ## Aggregate predictions
      return(aggregatePredictions(predictions))
    }, 
    
    aggregatePredictions = function(predictions) {
      ## Empty virtual function
    }, 
    
    predictionError = function() {
      ## Empty virtual function
    },
    
    variableImportance = function(type="permutation", start=2, num_threads=1) {
      ## Calculate tree VIM
      vim_trees <- mclapply(trees, function(x) {
        x$variableImportance(type, start)
      }, mc.cores = num_threads)
      
      ## Aggregate over trees
      rowMeans(simplify2array(vim_trees))
      #apply(simplify2array(vim_trees), MARGIN = 1, median)
      #simplify2array(vim_trees)
    },
    
    show = function() {
      cat("simpleRF Forest\n")
      cat("Type:                               ", treetype, "\n")
      cat("Splitrule:                          ", splitrule, "\n")
      cat("Maxstat splitting:                  ", maxstat, "\n")
      cat("Confounders:                        ", data$confounders, "\n")
      cat("Number of trees:                    ", num_trees, "\n")
      cat("Sample size:                        ", data$nrow, "\n")
      cat("Number of independent variables:    ", data$ncol-1, "\n")
      cat("Mtry:                               ", mtry, "\n")
      cat("Minimal node size for splitting:    ", minsplit, "\n")
      cat("Minimal terminal node size:         ", minbucket, "\n")
      cat("Minprop:                            ", minprop, "\n")
      cat("Alpha:                              ", alpha, "\n")
      cat("Replace:                            ", replace, "\n")
      cat("Predictions in leaves:              ", predleaf, "\n")
      cat("Always split variables (IDs):       ", always_split_varIDs, "\n")
      cat("Unordered factor handling:          ", unordered_factors, "\n")
      cat("OOB prediction error:               ", predictionError(), "\n")
    }, 
    
    print = function() {
      show()
    })
)


