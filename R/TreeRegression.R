
##' @title Regression tree class
##' @description Subclass for regression tree.
##' Contains all fields and methods used special for regression trees.
TreeRegression <- setRefClass("TreeRegression", 
  contains = "Tree",
  fields = list(),
  methods = list(
 
    splitNodeInternal = function(nodeID, possible_split_varIDs) {
      ## Check node size, stop if maximum reached
      if (length(sampleIDs[[nodeID]]) <= minsplit) {
        return(NULL)
      }
      
      ## Stop if node is pure      
      unique_response <- unique(data$subset(sampleIDs[[nodeID]], 1))
      if (length(unique_response) == 1) {
        return(NULL)
      }
      
      ## Find best split, stop if no decrease of impurity
      return(findBestSplit(nodeID, possible_split_varIDs))
    }, 
    
    findBestSplit = function(nodeID, possible_split_varIDs) {
      ## Initialize
      best_split <- NULL
      if (maxstat) {
        #decrease corresponds to the pvalue, thus smaller values are preferable
        best_split$decrease <- 2
      } else {
        best_split$decrease <- -1
      }
      best_split$varID <- -1
      best_split$value <- -1
      best_split$meanresid_left <- NULL
      best_split$meanresid_right <- NULL
      glm <- NULL
      #to perform the Benjamini-Hochberg adjustment a list of all p-values is required
      best_split$pvalues <- cbind(possible_split_varIDs,numeric(length(possible_split_varIDs)))
      
      ## Get response
      response <- data$subset(sampleIDs[[nodeID]], 1)
      
      ## Fit lm to the parent node & obtain residuals
      resid <- NULL
      if (splitrule == "Residuals") {
        #GLM can only be fitted for factor/ character predictors with more than one class in the node
        #test if glm includes factor/ character predictors & exclude factor/ character predictors with
        #only a single class in the node
        char_fact_ind <- which(attr(terms(data$glmdata), "dataClasses") %in% c("factor","character","ordered"))
        glmsubset <- data$glmsubset(sampleIDs[[nodeID]],)
        if (length(char_fact_ind)==0){
          #no characters or factors included in the formula
          glm <- lm(data$glmformula, glmsubset)
          resid <- glm$residuals
        } else{
          #formula includes characters or factors
          nunique <- sapply(glmsubset[,char_fact_ind, drop=FALSE], function(x){length(unique(x))})
          singleclass <- colnames(glmsubset)[char_fact_ind[which(nunique==1)]]
          if (length(singleclass)==0){
            #no need for changes, model can be fitted as usual
            glm <- lm(data$glmformula, glmsubset)
            resid <- glm$residuals
          } else{
            #factors/ characters with only one level must be excluded from glm
            #this also holds for interactions including them
            confounders_mod <- attr(terms(data$glmdata), "term.labels")[which(!grepl(singleclass, attr(terms(data$glmdata), "term.labels"),fixed=TRUE))]
            response <- as.character(attr(terms(data$glmdata), "variables"))[[2]][1]
            glmformula_mod <- as.formula(paste(response, paste(confounders_mod, collapse =" + "), sep=" ~ "))
            glm <- lm(glmformula_mod, glmsubset)
            resid <- glm$residuals
          }
        }
      }
      
      ## For all possible variables
      for (i in 1:length(possible_split_varIDs)) {
        split_varID <- possible_split_varIDs[i]
        data_values <- data$subset(sampleIDs[[nodeID]], split_varID)
        
        ## Handle ordered factors
        if (!is.numeric(data_values) & !is.ordered(data_values) & unordered_factors == "order_split") {
          ## Order factor levels
          num.response <- as.numeric(response)
          means <- sapply(levels(data_values), function(x) {
            mean(num.response[data_values == x])
          })
          levels.ordered <- as.character(levels(data_values)[order(means)])
          
          ## Get all levels not in node
          levels.missing <- setdiff(levels(data_values), levels.ordered)
          levels.ordered <- c(levels.missing, levels.ordered)
          
          ## Return reordered factor
          data_values <- factor(data_values, levels = levels.ordered, ordered = TRUE)
        }
        
        ## If still not ordered, use partition splitting
        if (!is.numeric(data_values) & !is.ordered(data_values)) {
          best_split = findBestSplitValuePartition(split_varID, data_values, best_split, response, resid)
          
          ## Set split levels left
          if (best_split$varID == split_varID) {
            split_levels_left[[nodeID]] <<- best_split$values_left
          }
        } else {
          best_split = findBestSplitValueOrdered(split_varID, data_values, best_split, response, resid)
          
          ## Set split levels left (empty if ordered splitting)
          if (unordered_factors == "order_split") {
            if (best_split$varID == split_varID) {
              split_levels_left[[nodeID]] <<- unique(data_values[data_values <= best_split$value])

              if (is.factor(data_values)) {
                ## Use same splits as in partition
                ints <- as.integer(factor(split_levels_left[[nodeID]], levels = levels(data$subset(sampleIDs[[nodeID]], split_varID))))
                if (sum(2^(ints-1)) >= 2^(max(as.numeric(data$subset(sampleIDs[[nodeID]], split_varID))) - 1)) {
                  split_levels_left[[nodeID]] <<- unique(data_values[data_values > best_split$value])
                }
              }
            }
          } else {
            if (best_split$varID == split_varID) {
              split_levels_left[[nodeID]] <<- list()
            }
          }
          
        }
      }
      
      if (best_split$varID < 0) {
        ## Stop if no good split found
        return(NULL)
      } else if (maxstat) {
        ## Adjust p values with Benjamini/Hochberg, use smallest
        adjusted_pvalue <- min(p.adjust(best_split$pvalues[,2], "fdr"), na.rm = TRUE)
        if (adjusted_pvalue > alpha) {
          ## Stop if no good split found
          return(NULL)
        } else {
          ## Return best split
          result <- NULL
          result$varID <- as.integer(best_split$varID)
          result$value <- best_split$value
          return(result)
        }
        
      } else {
        ## Return best split
        result <- NULL
        result$varID <- as.integer(best_split$varID)
        result$value <- best_split$value
        result$parentglm <- glm
        result$meanresid_left <- best_split$meanresid_left
        result$meanresid_right <- best_split$meanresid_right
        return(result)
      }      
    }, 
    
    findBestSplitValueOrdered = function(split_varID, data_values, best_split, response, resid) {
      ## maxstat=TRUE
      if (maxstat) {#test to determine split point
        if (splitrule == "Variance") {
          ## Use response for splitting
          maxstat_result <- maxstat_wrapper(y = response, x = data_values, smethod = "Wilcoxon", 
                          pmethod = pmethod, minprop = minprop, maxprop = 1-minprop)
          
        } else if (splitrule == "Residuals") { 
          ## Use residuals for splitting
          maxstat_result <- maxstat_wrapper(y = resid, x = data_values, smethod = "Wilcoxon", 
                          pmethod = pmethod, minprop = minprop, maxprop = 1-minprop)
        } else {
          stop("Unknown splitrule.")
        } 
        split_value <- maxstat_result$estimate
        #pvalue saved in decrease to save storage
        decrease <- maxstat_result$p.value
        
        if (!is.null(decrease) & !is.na(decrease)) {
          best_split$pvalues[which(best_split$pvalues[,1]==split_varID),2] <- decrease
          ## Use this split if better than before
          #careful: decrease is the pvalue, thus smaller values are preferable
          if (decrease < best_split$decrease) {
            best_split$value <- split_value
            best_split$varID <- split_varID
            best_split$decrease <- decrease
          }
        } else {
          best_split$pvalues[which(best_split$pvalues[,1]==split_varID),2] <- NA
        }
        
        return(best_split)
        
      } else {#exhaustive search
        ## For all possible splits
        possible_split_values <- unique(data_values)
        for (j in 1:length(possible_split_values)) {
          split_value <- possible_split_values[j]
          
          ## Sum responses & residuals in childs
          idx <- data_values <= split_value
          response_left <- response[idx]
          response_right <- response[!idx]
          resid_left <- resid[idx]
          resid_right <- resid[!idx]
          
          ## Skip if one child smaller than minbucket or if child empty if minbucket=1
          if (length(response_left) < minbucket | length(response_right) < minbucket) {
            next
          }
          
          if (splitrule == "Variance") {
            ## Decrease of impurity
            decrease <- sum(response_left)^2/length(response_left) + 
              sum(response_right)^2/length(response_right)
          } else if (splitrule == "Residuals") { 
            ## Decrease of impurity
            decrease <- sum(resid_left)^2/length(resid_left) + 
              sum(resid_right)^2/length(resid_right)
          } else {
            stop("Unknown splitrule.")
          }
          
          ## Use this split if better than before
          if (decrease > best_split$decrease) {
            best_split$value <- split_value
            best_split$varID <- split_varID
            best_split$decrease <- decrease
            if (splitrule == "Residuals"){
              best_split$meanresid_left <- mean(resid_left)
              best_split$meanresid_right <- mean(resid_right)
            }
          }
        }
        return(best_split)
      }
    },
    
    findBestSplitValuePartition = function(split_varID, data_values, best_split, response, resid) {
      ## For all possible splits
      possible_split_values <- sort(unique(data_values))

      ## For all 2^(n-1)-1 2-partitions
      num_partitions <- 2^(length(possible_split_values) - 1) - 1
      for (j in 1:num_partitions) {
        ## Convert number to logic vector
        left_idx <- as.bitvect(j, length = length(possible_split_values))
        values_left <- possible_split_values[left_idx]
        
        ## Sum responses & residuals in childs
        idx <- data_values %in% values_left
        response_left <- response[idx]
        response_right <- response[!idx]
        resid_left <- resid[idx]
        resid_right <- resid[!idx]
        
        ## Skip if one child smaller than minbucket or if child empty if minbucket=1
        if (length(response_left) < minbucket | length(response_right) < minbucket) {
          next
        }
        
        if (splitrule == "Variance") {
          ## Decrease of impurity
          decrease <- sum(response_left)^2/length(response_left) + 
            sum(response_right)^2/length(response_right)
        } else if (splitrule == "Residuals") { 
          ## Decrease of impurity
          decrease <- sum(resid_left)^2/length(resid_left) + 
            sum(resid_right)^2/length(resid_right)
        } else {
          stop("Unknown splitrule.")
        }
        
        ## Use this split if better than before
        if (decrease > best_split$decrease) {
          best_split$values_left <- values_left
          best_split$varID <- split_varID
          best_split$decrease <- decrease
        }
      }
      return(best_split)
    },
    
    estimate = function(nodeID) {      
      if (splitrule == "Residuals"){
        ## Return mean of the residuals
        return(mean_resid[[nodeID]])
      } else {
        ## Return mean response
        return(mean(data$subset(sampleIDs[[nodeID]], 1)))
      }
    }, 
    
#    fit_glm = function(nodeID) {      
#      ## Return glm
#      #GLM can only be fitted for factor/ character predictors with more than one class in the node
#      #test if glm includes factor/ character predictors & exclude factor/ character predictors with
#      #only a single class in the node
#      char_fact_ind <- which(attr(terms(data$glmdata), "dataClasses") %in% c("factor","character","ordered"))
#      glmsubset <- data$glmsubset(sampleIDs[[nodeID]],)
#      if (length(char_fact_ind)==0){
#        #no characters or factors included in the formula
#        glm <- lm(data$glmformula, glmsubset)
#      } else{
#        #formula includes characters or factors
#        nunique <- sapply(glmsubset[,char_fact_ind, drop=FALSE], function(x){length(unique(x))})
#        singleclass <- colnames(glmsubset)[char_fact_ind[which(nunique==1)]]
#        if (length(singleclass)==0){
#          #no need for changes, model can be fitted as usual
#          glm <- lm(data$glmformula, glmsubset)
#        } else{
#          #factors/ characters with only one level must be excluded from glm
#          #this also holds for interactions including them
#          confounders_mod <- attr(terms(data$glmdata), "term.labels")[which(!grepl(singleclass, attr(terms(data$glmdata), "term.labels"),fixed=TRUE))]
#          response <- as.character(attr(terms(data$glmdata), "variables"))[[2]][1]
#          glmformula_mod <- as.formula(paste(response, paste(confounders_mod, collapse =" + "), sep=" ~ "))
#          glm <- lm(glmformula_mod, glmsubset)
#        }
#      }
#      return(glm)
#    },
    
    getNodePrediction = function(nodeID, predictobs) {
      #if (splitrule=="Residuals"){
      #  pred <- predict.lm(parentglm[[nodeID]], newdata=predictobs) + split_values[nodeID]
      #} else{
        pred <- split_values[nodeID]
      #}
      return(pred)
    }, 
    
    predictionError = function(pred = NULL) {
      if (is.null(pred)) {
        pred <- predictOOB()
      }
      sum((pred - data$subset(oob_sampleIDs, 1))^2, na.rm = TRUE) / length(oob_sampleIDs)
    })
    
)