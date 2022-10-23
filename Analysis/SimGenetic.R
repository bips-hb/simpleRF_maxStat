## ***************************************************************************
## SimGenetic.R #####
##
## Script name: SimGenetic
##
## Purpose of script: Analyse the bias of the RF_VIM for the genetic simulation
##
## Author: Annika Swenne
##
## Date Created: 2022-07-26
##
## Version/ Date: V03 / 2022-08-09
##                V02 / 2022-07-27
##                V01 / 2022-07-26
##
## ***************************
##
## Modification Notes:
## V01/ 2022-07-26: only include the most relevant simulations (computational issues)
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
reg <-  makeExperimentRegistry(file.dir="SimSNP", seed = 1)
#reg <- loadRegistry(file.dir="SimSNP", writeable=TRUE)

#*Problems ####
source("Problems.R")
addProblem(name="simsnp_null", fun=SimSNP, seed=43)
addProblem(name="simsnp_power", fun=SimSNPPower, seed=44)
addProblem(name="pred_simsnp_null", fun=Pred_SimSNP, seed=45)
addProblem(name="pred_simsnp_power", fun=Pred_SimSNPPower, seed=46)

#*Algorithms ####
source("Algorithms.R")
addAlgorithm(name="rf_var", fun=rf_var_snp_wrapper)
addAlgorithm(name="rf_zhao", fun=rf_zhao_snp_wrapper)
addAlgorithm(name="rf_res", fun=rf_res_snp_wrapper)
addAlgorithm(name="rf_maxstat", fun=rf_maxstat_snp_wrapper)
addAlgorithm(name="pred_rf_var", fun=pred_rf_var_snp_wrapper)
addAlgorithm(name="pred_rf_zhao", fun=pred_rf_zhao_snp_wrapper)
addAlgorithm(name="pred_rf_res", fun=pred_rf_res_snp_wrapper)
addAlgorithm(name="pred_rf_maxstat", fun=pred_rf_maxstat_snp_wrapper)

#*Experiments ####
ntree <- 200
nsim <- 100

#**VIM ####
#Null case
prob_design <- list(simsnp_null=expand.grid(stringsAsFactors = FALSE))
algo_design <- list(rf_var = expand.grid(minsplit=10, mtry=33, num_trees=ntree, stringsAsFactors = FALSE),
                    rf_zhao = expand.grid(minsplit=10, mtry=33, num_trees=ntree, stringsAsFactors = FALSE),
                    rf_maxstat = data.frame(alpha=c(0.05,0.3,0.5,0.7), num_trees=ntree, minsplit=10, minbucket=1,
                                            predleaf=rep("Meanresid",times=4), mtry=33, stringsAsFactors = FALSE))
addExperiments(prob.designs=prob_design, algo.designs=algo_design, repls=nsim)
#Power case
prob_design <- list(simsnp_power=expand.grid(stringsAsFactors = FALSE, betax2=c(4,20)))
algo_design <- list(rf_var = expand.grid(minsplit=10, mtry=33, num_trees=ntree, stringsAsFactors = FALSE),
                    rf_zhao = expand.grid(minsplit=10, mtry=33, num_trees=ntree, stringsAsFactors = FALSE),
                    rf_maxstat = data.frame(alpha=c(0.05,0.3,0.5,0.7), num_trees=ntree, minsplit=10, minbucket=1,
                                           predleaf=rep("Meanresid",times=4), mtry=33, stringsAsFactors = FALSE))
addExperiments(prob.designs=prob_design, algo.designs=algo_design, repls=nsim)

#**Prediction ####
#Null case
prob_design <- list(pred_simsnp_null=expand.grid(stringsAsFactors = FALSE))
algo_design <- list(pred_rf_var = expand.grid(minsplit=10, mtry=33, num_trees=ntree, stringsAsFactors = FALSE),
                    pred_rf_zhao = expand.grid(minsplit=10, mtry=33, num_trees=ntree, stringsAsFactors = FALSE),
                    pred_rf_maxstat = data.frame(alpha=c(0.05,0.3,0.5,0.7), num_trees=ntree, minsplit=10, minbucket=1,
                                            predleaf=rep("Meanresid",times=4), mtry=33, stringsAsFactors = FALSE))
addExperiments(prob.designs=prob_design, algo.designs=algo_design, repls=nsim)
#Power case
prob_design <- list(pred_simsnp_power=expand.grid(stringsAsFactors = FALSE, betax2=20))
algo_design <- list(pred_rf_var = expand.grid(minsplit=10, mtry=33, num_trees=ntree, stringsAsFactors = FALSE),
                    pred_rf_zhao = expand.grid(minsplit=10, mtry=33, num_trees=ntree, stringsAsFactors = FALSE),
                    pred_rf_maxstat = data.frame(alpha=c(0.05,0.3,0.5,0.7), num_trees=ntree, minsplit=10, minbucket=1,
                                                 predleaf=rep("Meanresid",times=4), mtry=33, stringsAsFactors = FALSE))
addExperiments(prob.designs=prob_design, algo.designs=algo_design, repls=nsim)

#*Test jobs ####
#Test one job for each algorithm
#-VIM
#Null case
#testJob(id=head(findExperiments(prob.name = "simsnp_null", algo.name = "rf_var"), 1))
#testJob(id=head(findExperiments(prob.name = "simsnp_null", algo.name = "rf_zhao"), 1))
#testJob(id=head(findExperiments(prob.name = "simsnp_null", algo.name = "rf_maxstat", algo.pars = (alpha == 0.05)), 1))
#Power case
#testJob(id=head(findExperiments(prob.name = "simsnp_power", algo.name = "rf_var"), 1))
#testJob(id=head(findExperiments(prob.name = "simsnp_power", algo.name = "rf_zhao"), 1))
#testJob(id=head(findExperiments(prob.name = "simsnp_power", algo.name = "rf_maxstat", algo.pars = (alpha == 0.05)), 1))

#-Prediction
#Null case
#testJob(id=head(findExperiments(prob.name = "pred_simsnp_null", algo.name = "pred_rf_var"), 1))
#testJob(id=head(findExperiments(prob.name = "pred_simsnp_null", algo.name = "pred_rf_zhao"), 1))
#testJob(id=head(findExperiments(prob.name = "pred_simsnp_null", algo.name = "pred_rf_maxstat", algo.pars = (alpha == 0.05)), 1))
#Power case
#testJob(id=head(findExperiments(prob.name = "pred_simsnp_power", algo.name = "pred_rf_var"), 1))
#testJob(id=head(findExperiments(prob.name = "pred_simsnp_power", algo.name = "pred_rf_zhao"), 1))
#testJob(id=head(findExperiments(prob.name = "pred_simsnp_power", algo.name = "pred_rf_maxstat", algo.pars = (alpha == 0.05)), 1))

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
res <-  flatten(ijoin(getJobPars(), reduceResultsDataTable(ids)))
save(res,file="Res_snp.RData")





