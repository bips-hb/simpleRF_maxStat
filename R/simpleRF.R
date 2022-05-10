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
##' @param min_node_size Minimal node size. Default 1 for classification, 5 for regression, 3 for survival and 10 for probability estimation.
##' @param min_daughter Use the minimal node size for the daughter nodes. Default FALSE.
##' @param replace Sample with replacement. Default TRUE.
##' @param probability Grow a probability forest. Default FALSE.
##' @param splitrule Splitrule to use in trees. Default "Gini" for classification and probability forests, "Variance" for regression forests and "Logrank" for survival forests.
##' @param unordered_factors How to handle unordered factor variables. One of "ignore", "order_once", "order_split" and "partition" with default "ignore".
##' @param glmleaf Fit glm to leaves. Default FALSE.
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
                     min_node_size = NULL, min_daughter=FALSE, replace = TRUE, probability = FALSE, 
                     splitrule = NULL, unordered_factors = "ignore", 
                     glmleaf = FALSE, num_threads = 1) {
  
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
    mtry <- sqrt(ncol(model.data)-1)
  } else if (mtry > ncol(model.data)-1) {
    stop("Mtry cannot be larger than number of independent variables.")
  }
  if (is.null(min_node_size)) {
    if (treetype == "Classification") {
      min_node_size <- 1
    } else if (treetype == "Probability") {
      min_node_size <- 10
    } else if (treetype == "Regression") {
      min_node_size <- 5
    } else if (treetype == "Survival") {
      min_node_size <- 3
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
    
  ## Create forest object
  if (treetype == "Classification") {
    forest <- ForestClassification$new(num_trees = as.integer(num_trees), mtry = as.integer(mtry), 
                                       min_node_size = as.integer(min_node_size), min_daughter=min_daughter,
                                       replace = replace, splitrule = splitrule, glmleaf=glmleaf,
                                       data = Data$new(data = model.data, glmdata=glm.data, glmformula=glmformula), 
                                       formula = splitformula, unordered_factors = unordered_factors, 
                                       covariate_levels = covariate_levels,
                                       response_levels = levels(model.data[, 1]))
  } else if (treetype == "Probability") {
    forest <- ForestProbability$new(num_trees = as.integer(num_trees), mtry = as.integer(mtry), 
                                   min_node_size = as.integer(min_node_size), min_daughter=min_daughter,
                                   replace = replace, splitrule = splitrule, glmleaf=glmleaf,
                                   data = Data$new(data = model.data, glmdata=glm.data, glmformula=glmformula), 
                                   formula = splitformula, unordered_factors = unordered_factors,
                                   covariate_levels = covariate_levels,
                                   response_levels = levels(model.data[, 1]))
  } else if (treetype == "Regression") {
    forest <- ForestRegression$new(num_trees = as.integer(num_trees), mtry = as.integer(mtry), 
                                   min_node_size = as.integer(min_node_size), min_daughter=min_daughter,
                                   replace = replace, splitrule = splitrule, glmleaf=glmleaf,
                                   data = Data$new(data = model.data, glmdata=glm.data, glmformula=glmformula), 
                                   formula = splitformula, unordered_factors = unordered_factors, 
                                   covariate_levels = covariate_levels)
  } else if (treetype == "Survival") {
    idx.death <- model.data[, 1][, 2] == 1
    timepoints <- sort(unique(model.data[idx.death, 1][, 1]))
    forest <- ForestSurvival$new(num_trees = as.integer(num_trees), mtry = as.integer(mtry), 
                                 min_node_size = as.integer(min_node_size), min_daughter=min_daughter,
                                 replace = replace, splitrule = splitrule, glmleaf=glmleaf,
                                 data = Data$new(data = model.data, glmdata=glm.data, glmformula=glmformula), 
                                 formula = splitformula, unordered_factors = unordered_factors, 
                                 covariate_levels = covariate_levels,
                                 timepoints = timepoints)
  } else {
    stop("Unkown tree type.")
  }

  ## Grow forest
  forest$grow(num_threads = num_threads)

  ## Return forest
  return(forest) 
}