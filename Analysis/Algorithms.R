#### Algorithms.R #####
#### Wrapper functions for different versions of the RF algorithm 

#0) Help functions ####
coeff_effectx1 = function(tree){
  #if tree splits at x1:
  #get effect estimates for the confounder z1 & distinguish between glms affected & unaffected by 
  #splits at x1
  if (2 %in% tree$split_varIDs) {
    glm_x1 <- NULL
    remove_glm_x1 <- NULL
    
    for (i in 1:length(tree$child_nodeIDs)){
      if (!is.na(tree$split_varIDs[i])){
        #no terminal node
        if (tree$split_varIDs[i]==2 | i %in% glm_x1){
          glm_x1 <- c(glm_x1,tree$child_nodeIDs[[i]])
        }
        if (tree$split_varIDs[i]==2 & !(i %in% glm_x1)){
          #only remove daughters split at x1 that are NOT affected by previous splits at x1
          remove_glm_x1 <- c(remove_glm_x1,tree$child_nodeIDs[[i]])
        }
      }
    }
    
    if (tree$predleaf=="Meanresid"){
      #the leaves include the glm from the parent node which is not affected by the split
      #at x1 for the first daughter generation
      glm_x1 <- glm_x1[-which(glm_x1 %in% remove_glm_x1)]
    }
    glm_nox1 <- seq(1,length(tree$split_varIDs))
    coef_x1 <- NULL
    coef_nox1 <- NULL
    leaves <- which(is.na(tree$split_varIDs))
    glm_x1 <- glm_x1[which(glm_x1 %in% leaves)]
    glm_nox1 <- glm_nox1[which(glm_nox1 %in% leaves)]
    
    if (length(glm_x1)>0){
      glm_nox1 <- glm_nox1[-which(glm_nox1 %in% glm_x1)]
      #get estimate for regression parameter
      coef_x1 <- sapply(glm_x1, function(x){tree$nodeglm[[x]]$coefficients[2]})
    }
    
    #get estimate for regression parameter
    if (length(glm_nox1)>0){
      coef_nox1 <- sapply(glm_nox1, function(x){tree$nodeglm[[x]]$coefficients[2]})
    }
    
    betaz1 <- data.frame(x1split=c(rep("no",times=length(glm_nox1)),rep("yes",times=length(glm_x1))),
                         betaz1=c(coef_nox1,coef_x1))
    return(betaz1)
  } else {
    return(NULL)
  } 
}

#The following three functions are required for the runtime comparison with microbenchmark:
vim_rf_var_real = function(data, ...){
  #confounder z1 & z2 are included as splitting variables
  rf <- simpleRF(bmi ~ ., data=data, splitrule="Variance", replace=FALSE, predleaf="Meanout",
                 minsplit=50, unordered_factors="order_once", num_trees=200, num_threads=1)
  #return permutation vim
  rf$variableImportance(start=6, num_threads=1)
}

vim_rf_zhao_real = function(data, ...){
  #calculate residuals for outcome and predictors
  resdata <- as.data.frame(apply(data[,-1:-4], 2, function(x){lm(x ~ sex + age + country + isced, data=data)$residuals}))
  rf <- simpleRF(bmi ~ ., data=resdata, splitrule="Variance", replace=FALSE, predleaf="Meanout",
                 minsplit=50, unordered_factors="order_once", num_trees=200, num_threads=1)
  #return permutation vim
  rf$variableImportance(num_threads=1)
}

vim_rf_maxstat_real = function(data, ...){
  rf <- simpleRF(bmi ~ sex + age + country + isced |., data=data, splitrule="Residuals", maxstat=TRUE, replace=FALSE,
                 predleaf="Meanresid", minsplit=50, alpha=0.05, unordered_factors="order_once",
                 num_trees=200, num_threads=1)
  #return permutation vim
  rf$variableImportance(num_threads=1)
}

#1) Simple Simulation ####
#*VIM ####
rf_var_wrapper = function(data, job, instance, ...){
  #confounder z1 is included as splitting variable
  library(maxstat)
  library(devtools)
  load_all("~/Masterarbeit/simpleRF_maxStat")
  rf <- simpleRF(y ~ ., data=instance, splitrule="Variance", replace=FALSE, predleaf="Meanout", ...)
  #return permutation vim
  rf$variableImportance()
}

rf_zhao_wrapper = function(data, job, instance, ...){
  library(maxstat)
  library(devtools)
  load_all("~/Masterarbeit/simpleRF_maxStat")
  #calculate residuals for outcome and predictors
  resdata <- as.data.frame(apply(instance[,-2], 2, function(x){lm(x ~ z1, data=instance)$residuals}))
  rf <- simpleRF(y ~ ., data=resdata, splitrule="Variance", replace=FALSE, predleaf="Meanout", ...)
  #return permutation vim
  rf$variableImportance()
}

rf_res_wrapper = function(data, job, instance, ...){
  library(maxstat)
  library(devtools)
  load_all("~/Masterarbeit/simpleRF_maxStat")
  rf <- simpleRF(y ~ z1 |., data=instance, splitrule="Residuals", replace=FALSE, ...)
  #return permutation vim
  rf$variableImportance()
}

rf_maxstat_wrapper = function(data, job, instance, ...){
  library(maxstat)
  library(devtools)
  load_all("~/Masterarbeit/simpleRF_maxStat")
  rf <- simpleRF(y ~ z1 |., data=instance, splitrule="Residuals", maxstat=TRUE, replace=FALSE, ...)
  #return permutation vim
  rf$variableImportance()
}

#*Prediction ####
pred_rf_base_wrapper = function(data, job, instance, ...){
  #only the confounder z1 is included as splitting variable
  library(maxstat)
  library(devtools)
  load_all("~/Masterarbeit/simpleRF_maxStat")
  rf <- simpleRF(y ~ ., data=instance$train[,1:2], splitrule="Variance", replace=FALSE, predleaf="Meanout", ...)
  pred <- rf$predict(instance$test)
  #return RMSE
  rmse <- sqrt(mean((pred - instance$test$y)^2))
  names(rmse) <- "error"
  rmse
}

pred_rf_var_wrapper = function(data, job, instance, ...){
  #confounder z1 is included as splitting variable
  library(maxstat)
  library(devtools)
  load_all("~/Masterarbeit/simpleRF_maxStat")
  rf <- simpleRF(y ~ ., data=instance$train, splitrule="Variance", replace=FALSE, predleaf="Meanout", ...)
  pred <- rf$predict(instance$test)
  #return RMSE
  rmse <- sqrt(mean((pred - instance$test$y)^2))
  names(rmse) <- "error"
  rmse
}

pred_rf_zhao_wrapper = function(data, job, instance, accuracy="Outcome", ...){
  #allows to  measure the prediction accuracy based on the outcome or residuals
  library(maxstat)
  library(devtools)
  load_all("~/Masterarbeit/simpleRF_maxStat")
  #fit linear model for outcome and predictors
  lm <- apply(instance$train[,-2], 2, function(x){lm(x ~ z1, data=instance$train)})
  restrain <- as.data.frame(sapply(lm, residuals.lm))
  rf <- simpleRF(y ~ ., data=restrain, splitrule="Variance", replace=FALSE, predleaf="Meanout", ...)
  #predictions and residuals for test data
  predtest <- as.data.frame(sapply(lm, function(x){predict.lm(x, newdata=instance$test)}))
  restest <- cbind(instance$test[,-2] - predtest, z1=instance$test$z1)
  #predictions from RF
  pred <- rf$predict(restest)
  if (accuracy=="Outcome"){
    rmse <- sqrt(mean((predtest$y + pred - instance$test$y)^2))
  } else {
    rmse <- sqrt(mean((pred - restest$y)^2))
  }
  names(rmse) <- "error"
  rmse
}

pred_rf_res_wrapper = function(data, job, instance, ...){
  library(maxstat)
  library(devtools)
  load_all("~/Masterarbeit/simpleRF_maxStat")
  rf <- simpleRF(y ~ z1 |., data=instance$train, splitrule="Residuals", replace=FALSE, ...)
  pred <- rf$predict(instance$test)
  #return RMSE
  rmse <- sqrt(mean((pred - instance$test$y)^2))
  names(rmse) <- "error"
  rmse
}

pred_rf_maxstat_wrapper = function(data, job, instance, ...){
  library(maxstat)
  library(devtools) 
  load_all("~/Masterarbeit/simpleRF_maxStat")
  rf <- simpleRF(y ~ z1 |., data=instance$train, splitrule="Residuals", maxstat=TRUE, replace=FALSE, ...)
  pred <- rf$predict(instance$test)
  #return RMSE
  rmse <- sqrt(mean((pred - instance$test$y)^2))
  names(rmse) <- "error"
  rmse
}

#*Selection frequency ####
freq_rf_var_wrapper = function(data, job, instance, ...){
  #confounder z1 is included as splitting variable
  library(maxstat)
  library(devtools)
  load_all("~/Masterarbeit/simpleRF_maxStat")
  rf <- simpleRF(y ~ ., data=instance, splitrule="Variance", replace=FALSE, predleaf="Meanout",  ...)
  splitvar <- unlist(sapply(rf$trees, function(x){na.omit(x$split_varIDs)-1}))
  npred <- rf$data$ncol - 1
  freq <- table(factor(splitvar, levels = 1:npred))
  names(freq) <- rf$data$names[-1]
  #return selection frequency
  freq
}

freq_rf_zhao_wrapper = function(data, job, instance, ...){
  library(maxstat)
  library(devtools)
  load_all("~/Masterarbeit/simpleRF_maxStat")
  #calculate residuals for outcome and predictors
  resdata <- as.data.frame(apply(instance[,-2], 2, function(x){lm(x ~ z1, data=instance)$residuals}))
  rf <- simpleRF(y ~ ., data=resdata, splitrule="Variance", replace=FALSE, predleaf="Meanout", ...)
  splitvar <- unlist(sapply(rf$trees, function(x){na.omit(x$split_varIDs)-1}))
  npred <- rf$data$ncol - 1
  freq <- table(factor(splitvar, levels = 1:npred))
  names(freq) <- rf$data$names[-1]
  #return selection frequency
  freq
}

freq_rf_res_wrapper = function(data, job, instance, ...){
  library(maxstat)
  library(devtools)
  load_all("~/Masterarbeit/simpleRF_maxStat")
  rf <- simpleRF(y ~ z1 |., data=instance, splitrule="Residuals", replace=FALSE, ...)
  splitvar <- unlist(sapply(rf$trees, function(x){na.omit(x$split_varIDs)-1}))
  npred <- rf$data$ncol - 1
  freq <- table(factor(splitvar, levels = 1:npred))
  names(freq) <- rf$data$names[-1]
  #return selection frequency
  freq
}

freq_rf_maxstat_wrapper = function(data, job, instance, ...){
  library(maxstat)
  library(devtools)
  load_all("~/Masterarbeit/simpleRF_maxStat")
  rf <- simpleRF(y ~ z1 |., data=instance, splitrule="Residuals", maxstat=TRUE, replace=FALSE, ...)
  splitvar <- unlist(sapply(rf$trees, function(x){na.omit(x$split_varIDs)-1}))
  npred <- rf$data$ncol - 1
  freq <- table(factor(splitvar, levels = 1:npred))
  names(freq) <- rf$data$names[-1]
  #return selection frequency
  freq
}

#*Effect estimate ####
beta_rf_res_wrapper = function(data, job, instance, ...){
  #so far this is only designed for leafpred="Meanresid"
  library(maxstat)
  library(devtools)
  load_all("~/Masterarbeit/simpleRF_maxStat")
  source("Algorithms.R")
  rf <- simpleRF(y ~ z1 |., data=instance, splitrule="Residuals", replace=FALSE, ...)
  betaz1 <- lapply(rf$trees, coeff_effectx1)
  #return effect estimates
  do.call(rbind.data.frame,betaz1)
}

beta_rf_maxstat_wrapper = function(data, job, instance, ...){
  #so far this is only designed for leafpred="Meanresid"
  library(maxstat)
  library(devtools)
  load_all("~/Masterarbeit/simpleRF_maxStat")
  source("Algorithms.R")
  rf <- simpleRF(y ~ z1 |., data=instance, splitrule="Residuals", maxstat=TRUE, replace=FALSE, ...)
  betaz1 <- lapply(rf$trees, coeff_effectx1)
  #return effect estimates
  do.call(rbind.data.frame,betaz1)
}


#2) Genetic Simulation ####
#*VIM ####
rf_var_snp_wrapper = function(data, job, instance, ...){
  #confounder z1 & z2 are included as splitting variables
  library(maxstat)
  library(devtools)
  load_all("~/Masterarbeit/simpleRF_maxStat")
  rf <- simpleRF(y ~ ., data=instance, splitrule="Variance", replace=FALSE, predleaf="Meanout", ...)
  #return permutation vim
  rf$variableImportance()
}

rf_zhao_snp_wrapper = function(data, job, instance, ...){
  library(maxstat)
  library(devtools)
  load_all("~/Masterarbeit/simpleRF_maxStat")
  #calculate residuals for outcome and predictors
  resdata <- as.data.frame(apply(instance[,c(-2,-3)], 2, function(x){lm(x ~ z1+z2, data=instance)$residuals}))
  rf <- simpleRF(y ~ ., data=resdata, splitrule="Variance", replace=FALSE, predleaf="Meanout", ...)
  #return permutation vim
  rf$variableImportance()
}

rf_res_snp_wrapper = function(data, job, instance, ...){
  library(maxstat)
  library(devtools)
  load_all("~/Masterarbeit/simpleRF_maxStat")
  rf <- simpleRF(y ~ z1+z2 |., data=instance, splitrule="Residuals", replace=FALSE, ...)
  #return permutation vim
  rf$variableImportance()
}

rf_maxstat_snp_wrapper = function(data, job, instance, ...){
  library(maxstat)
  library(devtools)
  load_all("~/Masterarbeit/simpleRF_maxStat")
  rf <- simpleRF(y ~ z1+z2 |., data=instance, splitrule="Residuals", maxstat=TRUE, replace=FALSE, ...)
  #return permutation vim
  rf$variableImportance()
}

#*Prediction ####
pred_rf_base_snp_wrapper = function(data, job, instance, ...){
  #confounder z1 & z2 are included as splitting variable
  library(maxstat)
  library(devtools)
  load_all("~/Masterarbeit/simpleRF_maxStat")
  rf <- simpleRF(y ~ ., data=instance$train[,1:3], splitrule="Variance", replace=FALSE, predleaf="Meanout", ...)
  pred <- rf$predict(instance$test)
  #return RMSE
  rmse <- sqrt(mean((pred - instance$test$y)^2))
  names(rmse) <- "error"
  rmse
}

pred_rf_var_snp_wrapper = function(data, job, instance, ...){
  #confounder z1 & z2 are included as splitting variable
  library(maxstat)
  library(devtools)
  load_all("~/Masterarbeit/simpleRF_maxStat")
  rf <- simpleRF(y ~ ., data=instance$train, splitrule="Variance", replace=FALSE, predleaf="Meanout", ...)
  pred <- rf$predict(instance$test)
  #return RMSE
  rmse <- sqrt(mean((pred - instance$test$y)^2))
  names(rmse) <- "error"
  rmse
}

pred_rf_zhao_snp_wrapper = function(data, job, instance, accuracy="Outcome", ...){
  #allows to  measure the prediction accuracy based on the outcome or residuals
  library(maxstat)
  library(devtools)
  load_all("~/Masterarbeit/simpleRF_maxStat")
  #fit linear model for outcome and predictors
  lm <- apply(instance$train[,c(-2,-3)], 2, function(x){lm(x ~ z1+z2, data=instance$train)})
  restrain <- as.data.frame(sapply(lm, residuals.lm))
  rf <- simpleRF(y ~ ., data=restrain, splitrule="Variance", replace=FALSE, predleaf="Meanout", ...)
  #predictions and residuals for test data
  predtest <- as.data.frame(sapply(lm, function(x){predict.lm(x, newdata=instance$test)}))
  restest <- cbind(instance$test[,c(-2,-3)] - predtest, z1=instance$test$z1, z2=instance$test$z2)
  #predictions from RF
  pred <- rf$predict(restest)
  if (accuracy=="Outcome"){
    rmse <- sqrt(mean((predtest$y + pred - instance$test$y)^2))
  } else {
    rmse <- sqrt(mean((pred - restest$y)^2))
  }
  names(rmse) <- "error"
  rmse
}

pred_rf_res_snp_wrapper = function(data, job, instance, ...){
  library(maxstat)
  library(devtools)
  load_all("~/Masterarbeit/simpleRF_maxStat")
  rf <- simpleRF(y ~ z1+z2 |., data=instance$train, splitrule="Residuals", replace=FALSE, ...)
  pred <- rf$predict(instance$test)
  #return RMSE
  rmse <- sqrt(mean((pred - instance$test$y)^2))
  names(rmse) <- "error"
  rmse
}

pred_rf_maxstat_snp_wrapper = function(data, job, instance, ...){
  library(maxstat)
  library(devtools)
  load_all("~/Masterarbeit/simpleRF_maxStat")
  rf <- simpleRF(y ~ z1+z2 |., data=instance$train, splitrule="Residuals", maxstat=TRUE, replace=FALSE, ...)
  pred <- rf$predict(instance$test)
  #return RMSE
  rmse <- sqrt(mean((pred - instance$test$y)^2))
  names(rmse) <- "error"
  rmse
}


#3) Realdata ####
#*VIM ####
rf_var_real_wrapper = function(data, job, instance, ...){
  #confounders are included as splitting variables
  library(maxstat)
  library(devtools)
  load_all("~/Masterarbeit/simpleRF_maxStat")
  rf <- simpleRF(bmi ~ ., data=data, splitrule="Variance", replace=FALSE, predleaf="Meanout",
                 minsplit=50, unordered_factors="order_once", num_trees=200, num_threads=1)
  #return permutation vim
  rf$variableImportance(num_threads=1)
}

rf_zhao_real_wrapper = function(data, job, instance, ...){
  library(maxstat)
  library(devtools)
  load_all("~/Masterarbeit/simpleRF_maxStat")
  #calculate residuals for outcome and predictors
  resdata <- as.data.frame(apply(data[,-1:-4], 2, function(x){lm(x ~ sex + age + country + isced, data=data)$residuals}))
  rf <- simpleRF(bmi ~ ., data=resdata, splitrule="Variance", replace=FALSE, predleaf="Meanout",
                 minsplit=50, unordered_factors="order_once", num_trees=200, num_threads=1)
  #return permutation vim
  rf$variableImportance(num_threads=1)
}

rf_maxstat_real_wrapper = function(data, job, instance, ...){
  library(maxstat)
  library(devtools)
  load_all("~/Masterarbeit/simpleRF_maxStat")
  rf <- simpleRF(bmi ~ sex + age + country + isced |., data=data, splitrule="Residuals", maxstat=TRUE, replace=FALSE,
                 predleaf="Meanresid", minsplit=50, alpha=0.05, unordered_factors="order_once",
                 num_trees=200, num_threads=1)
  #return permutation vim
  rf$variableImportance(num_threads=1)
}

#*Prediction ####
pred_rf_base_real_wrapper = function(data, job, instance, ...){
  #only the confounders are included as splitting variables
  library(maxstat)
  library(devtools)
  load_all("~/Masterarbeit/simpleRF_maxStat")
  rf <- simpleRF(bmi ~ ., data=data[,1:5], splitrule="Variance", replace=FALSE, predleaf="Meanout",
                 minsplit=50, unordered_factors="order_once", num_trees=200, num_threads=1)
  #return OOB error
  rf$predictionError()
}

pred_rf_var_real_wrapper = function(data, job, instance, ...){
  #confounders are included as splitting variables
  library(maxstat)
  library(devtools)
  load_all("~/Masterarbeit/simpleRF_maxStat")
  rf <- simpleRF(bmi ~ ., data=data, splitrule="Variance", replace=FALSE, predleaf="Meanout",
                 minsplit=50, unordered_factors="order_once", num_trees=200, num_threads=1)
  #return OOB error
  rf$predictionError()
}

pred_rf_zhao_real_wrapper = function(data, job, instance, ...){
  library(maxstat)
  library(devtools)
  load_all("~/Masterarbeit/simpleRF_maxStat")
  #calculate residuals for outcome and predictors
  resdata <- as.data.frame(apply(data[,-1:-4], 2, function(x){lm(x ~ sex + age + country + isced, data=data)$residuals}))
  rf <- simpleRF(bmi ~ ., data=resdata, splitrule="Variance", replace=FALSE, predleaf="Meanout",
                 minsplit=50, unordered_factors="order_once", num_trees=200, num_threads=1)
  #return OOB error
  rf$predictionError()
}

pred_rf_maxstat_real_wrapper = function(data, job, instance, ...){
  library(maxstat)
  library(devtools)
  load_all("~/Masterarbeit/simpleRF_maxStat")
  rf <- simpleRF(bmi ~ sex + age + country + isced |., data=data, splitrule="Residuals", maxstat=TRUE, replace=FALSE,
                 predleaf="Meanresid", minsplit=50, alpha=0.05, unordered_factors="order_once",
                 num_trees=200, num_threads=1)
  #return OOB error
  rf$predictionError()
}

#*Runtime ####
run_rf_var_real_wrapper = function(data, job, instance, ...){
  #confounders are included as splitting variables
  library(microbenchmark)
  library(maxstat)
  library(devtools)
  load_all("~/Masterarbeit/simpleRF_maxStat")
  source("Algorithms.R")
  microbenchmark(rf_var = vim_rf_var_real(data=data, minsplit=minsplit, num_trees=ntree, mtry=mtry, num_threads=1),
                 times=1)
}

run_rf_zhao_real_wrapper = function(data, job, instance, ...){
  library(microbenchmark)
  library(maxstat)
  library(devtools) 
  load_all("~/Masterarbeit/simpleRF_maxStat")
  source("Algorithms.R")
  microbenchmark(rf_zhao = vim_rf_zhao_real(data=data, minsplit=minsplit, num_trees=ntree, mtry=mtry, num_threads=1),
                 times=1) 
}

run_rf_maxstat_real_wrapper = function(data, job, instance, ...){
  library(microbenchmark)
  library(maxstat)
  library(devtools)
  load_all("~/Masterarbeit/simpleRF_maxStat")
  source("Algorithms.R")
  microbenchmark(rf_maxstat = vim_rf_maxstat_real(data=data, minsplit=minsplit, num_trees=ntree, mtry=mtry, num_threads=1),
                 times=1) 
}