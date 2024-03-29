##' Implements Random Forests (Breiman 2001) with emphasis on simplicity. 
##' Uses reference classes and only plain \code{R}. 
##' Not optimized for computation speed. 
##' Allows rapid prototyping of RF-type algorithms.
##' 
##' Unordered factor variables can be handled in different ways. 
##' Use "ignore" to treat them as ordered in the order of the factor levels. 
##' With "order_once" and "order_split" they are ordered by their response values. For "order_once" this is done once before the analysis, for "order_split" this is done in each split.
##' With "partition" all 2-partitions of the factor levels are considered for splitting.
##' 
##' @title simpleRF
##' @param formula Object of class \code{formula} or \code{character} describing the model to fit.
##' @param data Training data of class \code{data.frame}.
##' @param num_trees Number of trees.
##' @param mtry Number of variables to possibly split at in each node.
##' @param minsplit Minimal node size for a node to be considered for splitting. Default 1 for classification, 5 for regression, 3 for survival and 10 for probability estimation.
##' @param minbucket Minimal node size for a terminal node. Default 0, i.e. no restrictions on the terminal node size.
##' @param replace Sample with replacement. Default TRUE.
##' @param probability Grow a probability forest. Default FALSE.
##' @param splitrule Splitrule to use in trees. Default "Gini" for classification and probability forests, "Variance" for regression forests and "Logrank" for survival forests.
##' @param always_split_variables Character vector with variable names to be always selected in addition to the mtry variables tried for splitting. 
##' @param unordered_factors How to handle unordered factor variables. One of "ignore", "order_once", "order_split" and "partition" with default "ignore".
##' @param leafpred Prediction in the leaves. Default "Meanresid" for splitrule="Residuals" & "Meanout" for splitrule="Variance".
##' @param maxstat Use maximally selected statistics for splitting. Default FALSE.
##' @param minprop Lower quantile of covariate distribtuion to be considered for splitting.
##' @param alpha Significance threshold to allow splitting.
##' @param pmethod P-value approximation method. Can either be one of the generic methods "approximation", "permutation" and "asymptotic" or a specific method of \code{\link{maxstat_chisq}} for classification or \code{\link{maxstat}} for regression and survival. See the respective function for details. 
##' @param num_threads Number of threads used for mclapply, set to 1 for debugging.
##' @examples 
##' \donttest{
##' library(simpleRF)
##' 
##' # Classification
##' simpleRF(Species ~ ., iris)
##' 
##' # Prediction
##' train_idx <- sample(nrow(iris), 2/3 * nrow(iris))
##' iris_train <- iris[train_idx, ]
##' iris_test <- iris[-train_idx, ]
##' rf_iris <- simpleRF(Species ~ ., data = iris_train)
##' pred_iris <- rf_iris$predict(iris_test)
##' table(iris_test$Species, pred_iris)
##' }
##' 
##' @author Marvin N. Wright
##' @references
##' Breiman, L. (2001). Random forests. Mach Learn, 45(1), 5-32. \cr
##' @import stats
##' @export
simpleRF <- function(formula, data, num_trees = 50, mtry = NULL, 
                     minsplit = NULL, minbucket = 1, replace = TRUE, probability = FALSE, 
                     splitrule = NULL, always_split_variables = NULL, unordered_factors = "ignore", 
                     predleaf = NULL, maxstat = FALSE, minprop =  0.1, alpha=0.05, 
                     pmethod = "approximation", num_threads = 1) {
  
  #distinguish between RF with and without confounders
  if (grepl("|", deparse1(formula), fixed = TRUE)){
    #formula includes confounders (indicated by | on the right side of the formula)
    response <- attr(terms(formula), "variables")[[2]]
    confounders <- as.character(attr(terms(formula), "variables")[[3]])[2]
    predictors <- as.character(attr(terms(formula), "variables")[[3]])[3]
    glmformula <- as.formula(paste(response, confounders, sep=" ~ "))
    #glmformula is needed as an extra component of the data object to fit the glm at each split
    #otherwise it is not possible to consider interactions
    splitformula <- as.formula(paste(response, predictors, sep=" ~ "))
    glm.data <- model.frame(glmformula, data)
    model.data <- model.frame(splitformula, data)
    if (predictors == "."){
      #default: if no splitting variables are specified then confounders (main effects) are not
      #considered as splitting variables
      mainconf <- attr(terms(glmformula), "term.labels")[which(attr(terms(glmformula), "order")==1)]
      model.data <- model.data[,-which(colnames(model.data) %in% mainconf)]
    }
    
  } else {
    #formula does not include confounders
    splitformula <- formula
    model.data <- model.frame(splitformula, data)
    #glmformula <- as.formula(NULL)
    #glm.data <- as.data.frame(NULL)
    glmformula <- as.formula(paste(as.character(formula)[[2]], "1", sep=" ~ "))
    glm.data <- model.frame(glmformula, data)
  }
  
  always_split_varIDs <- which(colnames(model.data) %in% always_split_variables)
  
  if (class(model.data[, 1]) == "factor") {
    if (probability) {
      treetype <- "Probability" 
    } else {
      treetype <- "Classification"
    }
  } else if (class(model.data[, 1]) == "numeric") {
    treetype <- "Regression"
  } else if (class(model.data[, 1]) == "Surv") {
    treetype <- "Survival"
  } else {
    stop("Unkown response type.")
  }
  
  ## Check parameters
  if (is.null(mtry)) {
    mtry <- sqrt(ncol(model.data)-length(always_split_variables)-1)
  } else if (mtry + length(always_split_variables) > ncol(model.data)-1) {
    stop("Mtry plus the number of always_split_variables cannot be larger than number of independent variables.")
  }
  if (is.null(minsplit)) {
    if (treetype == "Classification") {
      minsplit <- 1
    } else if (treetype == "Probability") {
      minsplit <- 10
    } else if (treetype == "Regression") {
      minsplit <- 5
    } else if (treetype == "Survival") {
      minsplit <- 3
    }
  }
  
  ## Splitrule
  if (is.null(splitrule)) {
    if (treetype == "Classification") {
      splitrule <- "Gini"
    } else if (treetype == "Probability") {
      splitrule <- "Gini"
    } else if (treetype == "Regression") {
      splitrule <- "Variance"
    } else if (treetype == "Survival") {
      splitrule <- "Logrank"
    }
  }
  
  ## Predictions in leaves
  if (is.null(predleaf)) {
    if (splitrule == "Residuals") {
      predleaf <- "Meanresid"
    } else {
      predleaf <- "Meanout"
    }
  }
  
  ## Unordered factors
  if (!(unordered_factors %in% c("ignore", "order_once", "order_split", "partition"))) {
    stop("Unknown value for unordered_factors.")
  }
  covariate_levels <- list()
  
  if (unordered_factors == "order_once") {
    ## Reorder factor columns depending on response type
    model.data <- reorder.factor.columns(model.data)
    
    ## Save levels
    covariate_levels <- lapply(model.data[, -1], levels)
  }
  else if (unordered_factors == "ignore") {
    ## Just set to ordered if "ignore"
    character.idx <- sapply(model.data[, -1], is.character)
    ordered.idx <- sapply(model.data[, -1], is.ordered)
    factor.idx <- sapply(model.data[, -1], is.factor)
    recode.idx <- character.idx | (factor.idx & !ordered.idx)
    model.data[, -1][, recode.idx] <- lapply(model.data[, -1][, recode.idx, drop=FALSE], as.ordered)
    
    ## Save levels
    covariate_levels <- lapply(model.data[, -1], levels)
  }
  
  ## Select pmethod
  if (treetype == "Classification") {
    if (pmethod == "approximation") {
      pmethod <- "betensky"
    } else if (pmethod == "asymptotic") {
      stop("Currently no asymptotic method for classification implemented.")
    }
  } else if (treetype == "Regression" | treetype == "Survival") {
    if (pmethod == "approximation") {
      pmethod <- "minLau92Lau94"
    } else if (pmethod == "permutation") {
      pmethod <- "condMC"
    } else if (pmethod == "asymptotic") {
      pmethod <- "exactGauss"
    }
  } else {
    stop("Unkown tree type.")
  }
  
    
  ## Create forest object
  if (treetype == "Classification") {
    forest <- ForestClassification$new(num_trees = as.integer(num_trees), mtry = as.integer(mtry), 
                                       minsplit = as.integer(minsplit), minbucket = as.integer(minbucket),
                                       replace = replace, splitrule = splitrule, predleaf=predleaf,
                                       data = Data$new(data = model.data, glmdata=glm.data, glmformula=glmformula), 
                                       formula = splitformula, always_split_varIDs = always_split_varIDs,
                                       unordered_factors = unordered_factors, 
                                       covariate_levels = covariate_levels,
                                       response_levels = levels(model.data[, 1]),
                                       maxstat = maxstat, minprop = minprop, alpha = alpha, pmethod = pmethod)
  } else if (treetype == "Probability") {
    forest <- ForestProbability$new(num_trees = as.integer(num_trees), mtry = as.integer(mtry), 
                                    minsplit = as.integer(minsplit), minbucket = as.integer(minbucket),
                                   replace = replace, splitrule = splitrule, predleaf=predleaf,
                                   data = Data$new(data = model.data, glmdata=glm.data, glmformula=glmformula), 
                                   formula = splitformula, always_split_varIDs = always_split_varIDs,
                                   unordered_factors = unordered_factors,
                                   covariate_levels = covariate_levels,
                                   response_levels = levels(model.data[, 1]),
                                   maxstat = maxstat, minprop = minprop, alpha = alpha, pmethod = pmethod)
  } else if (treetype == "Regression") {
    forest <- ForestRegression$new(num_trees = as.integer(num_trees), mtry = as.integer(mtry), 
                                   minsplit = as.integer(minsplit), minbucket = as.integer(minbucket),
                                   replace = replace, splitrule = splitrule, predleaf=predleaf,
                                   data = Data$new(data = model.data, glmdata=glm.data, glmformula=glmformula), 
                                   formula = splitformula, always_split_varIDs = always_split_varIDs,
                                   unordered_factors = unordered_factors, 
                                   covariate_levels = covariate_levels,
                                   maxstat = maxstat, minprop = minprop, alpha = alpha, pmethod = pmethod)
  } else if (treetype == "Survival") {
    idx.death <- model.data[, 1][, 2] == 1
    timepoints <- sort(unique(model.data[idx.death, 1][, 1]))
    forest <- ForestSurvival$new(num_trees = as.integer(num_trees), mtry = as.integer(mtry), 
                                 minsplit = as.integer(minsplit), minbucket = as.integer(minbucket),
                                 replace = replace, splitrule = splitrule, predleaf=predleaf,
                                 data = Data$new(data = model.data, glmdata=glm.data, glmformula=glmformula), 
                                 formula = splitformula, always_split_varIDs = always_split_varIDs,
                                 unordered_factors = unordered_factors, 
                                 covariate_levels = covariate_levels,
                                 timepoints = timepoints,
                                 maxstat = maxstat, minprop = minprop, alpha = alpha, pmethod = pmethod)
  } else {
    stop("Unkown tree type.")
  }

  ## Grow forest
  forest$grow(num_threads = num_threads)

  ## Return forest
  return(forest) 
}