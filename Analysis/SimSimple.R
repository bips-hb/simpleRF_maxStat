## ***************************************************************************
## SimSimple.R #####
##
## Script name: SimSimple
##
## Purpose of script: Analyse the bias of the RF_VIM for the simple simulation
##
## Author: Annika Swenne
##
## Date Created: 2022-06-20
##
## Version/ Date: V03 / 2022-08-09 
##                V02 / 2022-07-27 
##                V01 / 2022-07-25
##
## ***************************
##
## Modification Notes:
## V02/ 2022-07-27: only include the most relevant simulations (computational issues)
##   
## ***************************
##
## Notes/ to do:
##
## ***************************************************************************

#set working directory for HPC
setwd("~/Masterarbeit")
#setwd("C:/Users/Annika Swenne/Documents/Uni_Bremen/MedicalBiometry/Masterarbeit/Code/Batchtools")

#load the required packages:  (uncomment as required)
library(batchtools)

set.seed(42)

## ***************************************************************************
#*Registry ####
reg <-  makeExperimentRegistry(file.dir="SimSimple", seed = 1)
#reg <- loadRegistry(file.dir="SimSimple", writeable=TRUE)

#*Problems ####
source("Problems.R")
addProblem(name="simsimple_null", fun=SimSimple, seed=43)
addProblem(name="simsimple_power", fun=SimSimplePower, seed=44)
addProblem(name="simsimple_local", fun=SimSimpleLocal, seed=46)
addProblem(name="pred_simsimple_null", fun=Pred_SimSimple, seed=47)
addProblem(name="pred_simsimple_power", fun=Pred_SimSimplePower, seed=48)
addProblem(name="pred_simsimple_local", fun=Pred_SimSimpleLocal, seed=49)

#*Algorithms ####
source("Algorithms.R")
addAlgorithm(name="rf_var", fun=rf_var_wrapper)
addAlgorithm(name="rf_zhao", fun=rf_zhao_wrapper)
addAlgorithm(name="rf_res", fun=rf_res_wrapper)
addAlgorithm(name="rf_maxstat", fun=rf_maxstat_wrapper)
addAlgorithm(name="freq_rf_var", fun=freq_rf_var_wrapper)
addAlgorithm(name="freq_rf_zhao", fun=freq_rf_zhao_wrapper)
addAlgorithm(name="freq_rf_res", fun=freq_rf_res_wrapper)
addAlgorithm(name="freq_rf_maxstat", fun=freq_rf_maxstat_wrapper)
addAlgorithm(name="beta_rf_res", fun=beta_rf_res_wrapper)
addAlgorithm(name="beta_rf_maxstat", fun=beta_rf_maxstat_wrapper)
addAlgorithm(name="pred_rf_var", fun=pred_rf_var_wrapper)
addAlgorithm(name="pred_rf_zhao", fun=pred_rf_zhao_wrapper)
addAlgorithm(name="pred_rf_res", fun=pred_rf_res_wrapper)
addAlgorithm(name="pred_rf_maxstat", fun=pred_rf_maxstat_wrapper)

#*Experiments ####
ntree <- 100
nsim <- 100

#**VIM ####
#Null case
prob_design <- list(simsimple_null=expand.grid(stringsAsFactors = FALSE))
algo_design <- list(rf_var = expand.grid(minsplit=10, num_trees=ntree, stringsAsFactors = FALSE),
                    rf_zhao = expand.grid(minsplit=10, num_trees=ntree, stringsAsFactors = FALSE),
                    rf_res = data.frame(minsplit=c(10,1,10), minbucket=c(1,10,1), num_trees=ntree, 
                                        predleaf=c(rep("Meanout",times=1),rep("GLM",times=1),rep("Meanresid",times=1)), stringsAsFactors = FALSE),
                    rf_maxstat = data.frame(alpha=rep(c(0.05,0.3,0.5,0.7),times=3), num_trees=ntree, minsplit=c(10,10,10,10,20,20,20,20,10,10,10,10), 
                                            minbucket=c(1,1,1,1,10,10,10,10,1,1,1,1),
                                            predleaf=c(rep("Meanout",times=4),rep("GLM",times=4),rep("Meanresid",times=4)), stringsAsFactors = FALSE))
addExperiments(prob.designs=prob_design, algo.designs=algo_design, repls=nsim)
#Power case
prob_design <- list(simsimple_power=expand.grid(betax2=0.5, stringsAsFactors = FALSE))
algo_design <- list(rf_var = expand.grid(minsplit=10, num_trees=ntree, stringsAsFactors = FALSE),
                    rf_zhao = expand.grid(minsplit=10, num_trees=ntree, stringsAsFactors = FALSE),
                    rf_res = data.frame(minsplit=c(10,1,10), minbucket=c(1,10,1), num_trees=ntree, 
                                        predleaf=c(rep("Meanout",times=1),rep("GLM",times=1),rep("Meanresid",times=1)), stringsAsFactors = FALSE),
                    rf_maxstat = data.frame(alpha=rep(c(0.05,0.3,0.5,0.7),times=3), num_trees=ntree, minsplit=c(10,10,10,10,20,20,20,20,10,10,10,10), 
                                            minbucket=c(1,1,1,1,10,10,10,10,1,1,1,1),
                                            predleaf=c(rep("Meanout",times=4),rep("GLM",times=4),rep("Meanresid",times=4)), stringsAsFactors = FALSE))
addExperiments(prob.designs=prob_design, algo.designs=algo_design, repls=nsim)
#Power case sensitivity
prob_design <- list(simsimple_power=expand.grid(betax2=0.15, stringsAsFactors = FALSE))
algo_design <- list(rf_var = expand.grid(minsplit=10, num_trees=ntree, stringsAsFactors = FALSE),
                    rf_zhao = expand.grid(minsplit=10, num_trees=ntree, stringsAsFactors = FALSE),
                    rf_maxstat = data.frame(alpha=c(0.05,0.3,0.5,0.7), num_trees=ntree, minsplit=10, minbucket=1,
                                            predleaf=rep("Meanresid",times=4), stringsAsFactors = FALSE))
addExperiments(prob.designs=prob_design, algo.designs=algo_design, repls=nsim)
#Local
prob_design <- list(simsimple_local=expand.grid(stringsAsFactors = FALSE))
algo_design <- list(rf_var = expand.grid(minsplit=10, num_trees=ntree, mtry=c(3,10), stringsAsFactors = FALSE),
                    rf_zhao = expand.grid(minsplit=10, num_trees=ntree, mtry=c(3,10), stringsAsFactors = FALSE),
                    rf_maxstat = expand.grid(alpha=c(0.05,0.3,0.5,0.7), mtry=c(3,10), 
                                            num_trees=ntree, minsplit=10, minbucket=1,
                                            predleaf="Meanresid", stringsAsFactors = FALSE))
addExperiments(prob.designs=prob_design, algo.designs=algo_design, repls=nsim)
#Local sensitivity
prob_design <- list(simsimple_local=expand.grid(addx2=TRUE, stringsAsFactors = FALSE))
algo_design <- list(rf_var = expand.grid(minsplit=10, num_trees=ntree, mtry=c(3,10), stringsAsFactors = FALSE),
                    rf_zhao = expand.grid(minsplit=10, num_trees=ntree, mtry=c(3,10), stringsAsFactors = FALSE),
                    rf_maxstat = expand.grid(alpha=0.05, mtry=c(3,10), 
                                             num_trees=ntree, minsplit=10, minbucket=1,
                                             predleaf="Meanresid", stringsAsFactors = FALSE))
addExperiments(prob.designs=prob_design, algo.designs=algo_design, repls=nsim)

#**Selection frequency ####
#Null case
prob_design <- list(simsimple_null=expand.grid(stringsAsFactors = FALSE))
algo_design <- list(freq_rf_var = expand.grid(minsplit=10, num_trees=ntree, stringsAsFactors = FALSE),
                    freq_rf_zhao = expand.grid(minsplit=10, num_trees=ntree, stringsAsFactors = FALSE),
                    freq_rf_res = expand.grid(minsplit=10, num_trees=ntree, predleaf="Meanresid",stringsAsFactors = FALSE),
                    freq_rf_maxstat = data.frame(alpha=c(0.05,0.3,0.5,0.7), num_trees=ntree, minsplit=10, minbucket=1,
                                            predleaf="Meanresid", stringsAsFactors = FALSE))
addExperiments(prob.designs=prob_design, algo.designs=algo_design, repls=nsim)

#**Effect estimate ####
#Null case
prob_design <- list(simsimple_null=expand.grid(stringsAsFactors = FALSE))
algo_design <- list(beta_rf_res = data.frame(minsplit=c(1,1,10,1), minbucket=c(10,20,1,10), num_trees=ntree, 
                                             predleaf=c(rep("GLM",times=2),rep("Meanresid",times=2)), stringsAsFactors = FALSE),
                    beta_rf_maxstat = data.frame(alpha=1, minsplit=c(20,40,10,20), minbucket=c(10,20,1,10), num_trees=ntree, 
                                             predleaf=c(rep("GLM",times=2),rep("Meanresid",times=2)), stringsAsFactors = FALSE))
addExperiments(prob.designs=prob_design, algo.designs=algo_design, repls=nsim)

#**Prediction ####
#Null case
prob_design <- list(pred_simsimple_null=expand.grid(stringsAsFactors = FALSE))
algo_design <- list(pred_rf_var = expand.grid(minsplit=10, num_trees=ntree, stringsAsFactors = FALSE),
                    pred_rf_zhao = expand.grid(minsplit=10, num_trees=ntree, stringsAsFactors = FALSE),
                    pred_rf_res = data.frame(minsplit=c(10,1,10), minbucket=c(1,10,1), num_trees=ntree, 
                                        predleaf=c(rep("Meanout",times=1),rep("GLM",times=1),rep("Meanresid",times=1)), stringsAsFactors = FALSE),
                    pred_rf_maxstat = data.frame(alpha=rep(c(0.05,0.3,0.5,0.7),times=3), num_trees=ntree, minsplit=c(10,10,10,10,20,20,20,20,10,10,10,10), 
                                            minbucket=c(1,1,1,1,10,10,10,10,1,1,1,1),
                                            predleaf=c(rep("Meanout",times=4),rep("GLM",times=4),rep("Meanresid",times=4)), stringsAsFactors = FALSE))
addExperiments(prob.designs=prob_design, algo.designs=algo_design, repls=nsim)
#Power case
prob_design <- list(pred_simsimple_power=expand.grid(stringsAsFactors = FALSE))
algo_design <- list(pred_rf_var = expand.grid(minsplit=10, num_trees=ntree, stringsAsFactors = FALSE),
                    pred_rf_zhao = expand.grid(minsplit=10, num_trees=ntree, stringsAsFactors = FALSE),
                    pred_rf_res = data.frame(minsplit=c(10,1,10), minbucket=c(1,10,1), num_trees=ntree, 
                                             predleaf=c(rep("Meanout",times=1),rep("GLM",times=1),rep("Meanresid",times=1)), stringsAsFactors = FALSE),
                    pred_rf_maxstat = data.frame(alpha=rep(c(0.05,0.3,0.5,0.7),times=3), num_trees=ntree, minsplit=c(10,10,10,10,20,20,20,20,10,10,10,10), 
                                            minbucket=c(1,1,1,1,10,10,10,10,1,1,1,1),
                                            predleaf=c(rep("Meanout",times=4),rep("GLM",times=4),rep("Meanresid",times=4)), stringsAsFactors = FALSE))
addExperiments(prob.designs=prob_design, algo.designs=algo_design, repls=nsim)
#Local
prob_design <- list(pred_simsimple_local=expand.grid(stringsAsFactors = FALSE))
algo_design <- list(pred_rf_var = expand.grid(minsplit=10, num_trees=ntree, mtry=c(3,10), stringsAsFactors = FALSE),
                    pred_rf_zhao = expand.grid(minsplit=10, num_trees=ntree, mtry=c(3,10), stringsAsFactors = FALSE),
                    pred_rf_maxstat = expand.grid(alpha=c(0.05,0.3,0.5,0.7), mtry=c(3,10), 
                                                 num_trees=ntree, minsplit=10, minbucket=1,
                                                 predleaf="Meanresid", stringsAsFactors = FALSE))
addExperiments(prob.designs=prob_design, algo.designs=algo_design, repls=nsim)

#*Test jobs ####
#Test one job for each algorithm
#-VIM
#Null case
#testJob(id=head(findExperiments(prob.name = "simsimple_null", algo.name = "rf_var"), 1))
#testJob(id=head(findExperiments(prob.name = "simsimple_null", algo.name = "rf_zhao"), 1))
#testJob(id=head(findExperiments(prob.name = "simsimple_null", algo.name = "rf_res", algo.pars = (minsplit == 20)), 1))
#testJob(id=head(findExperiments(prob.name = "simsimple_null", algo.name = "rf_maxstat", algo.pars = (alpha == 0.05)), 1))
#Power case
#testJob(id=head(findExperiments(prob.name = "simsimple_power", algo.name = "rf_res", algo.pars = (minsplit == 20)), 1))
#testJob(id=head(findExperiments(prob.name = "simsimple_power", algo.name = "rf_var"), 1))
#testJob(id=head(findExperiments(prob.name = "simsimple_power", algo.name = "rf_zhao"), 1))
#testJob(id=head(findExperiments(prob.name = "simsimple_power", algo.name = "rf_maxstat", algo.pars = (alpha == 0.05)), 1))

#-Prediction
#Null case
#testJob(id=head(findExperiments(prob.name = "pred_simsimple_null", algo.name = "pred_rf_res", algo.pars = (minsplit == 20)), 1))
#testJob(id=head(findExperiments(prob.name = "pred_simsimple_null", algo.name = "pred_rf_var"), 1))
#testJob(id=head(findExperiments(prob.name = "pred_simsimple_null", algo.name = "pred_rf_zhao"), 1))
#testJob(id=head(findExperiments(prob.name = "pred_simsimple_null", algo.name = "pred_rf_maxstat", algo.pars = (alpha == 0.05)), 1))
#Power case
#testJob(id=head(findExperiments(prob.name = "pred_simsimple_power", algo.name = "pred_rf_res", algo.pars = (minsplit == 20)), 1))
#testJob(id=head(findExperiments(prob.name = "pred_simsimple_power", algo.name = "pred_rf_var"), 1))
#testJob(id=head(findExperiments(prob.name = "pred_simsimple_power", algo.name = "pred_rf_zhao"), 1))
#testJob(id=head(findExperiments(prob.name = "pred_simsimple_power", algo.name = "pred_rf_maxstat", algo.pars = (alpha == 0.05)), 1))

#*Submit ####
if (grepl("node\\d{2}|bipscluster", system("hostname", intern = TRUE))) {
  ids <- findNotStarted()
  #ids <- findNotDone()
  ids[, chunk := chunk(job.id, chunk.size = 50)]
  submitJobs(ids = ids, # walltime in seconds, 10 days max, memory in MB
             resources = list(name = "SimSimple", chunks.as.arrayjobs = TRUE, 
                              ncpus = 1, memory = 6000, walltime = 10*24*3600, 
                              max.concurrent.jobs = 40))
} else {
  ids <- findNotDone()
  submitJobs(ids)
}
waitForJobs()

#*Get results ####
ids <- findExperiments()
ids_beta1 <- findExperiments(algo.name ="beta_rf_res")
ids_beta2 <- findExperiments(algo.name="beta_rf_maxstat")
ids2 <- ids[-which(ids$job.id %in% c(ids_beta1$job.id,ids_beta2$job.id)),]

#all other algorithms
res <-  flatten(ijoin(getJobPars(), reduceResultsDataTable(ids2)))
save(res,file="Res_simple.RData")

#effect estimates for beta
temp <- reduceResultsDataTable(rbind(ids_beta1,ids_beta2))
nrow <- sapply(temp$result, nrow)
effect <- do.call(rbind.data.frame,temp$result)
res <- cbind(job.id=rep(temp$job.id, times=nrow),
             effect)
res <-  ijoin(getJobPars(), res)
res <- flatten(res)
save(res,file="Res_simple_effect.RData")





