## ***************************************************************************
## Realdata.R #####
##
## Script name: Realdata
##
## Purpose of script: Analyse the GAS for BMI for the IDEFICS/I.Family data &
##                    compute the runtime
##
## Author: Annika Swenne
##
## Date Created: 2022-08-25
##
## Version/ Date: V01 / 2022-08-25
##
## ***************************
##
## Modification Notes:
## V01/ 2022-08-23: 
##   
## ***************************
##
## Notes/ to do:
##
## ***************************************************************************

#set working directory for HPC
setwd("~/Masterarbeit")
#setwd("C:/Users/Annika Swenne/Documents/Uni_Bremen/MedicalBiometry/Masterarbeit/Code/Realdata")

#load the required packages:  (uncomment as required)
library(batchtools)

#read the data
data <- readRDS("~/Masterarbeit/selected_snps_cc.Rds")
data$isced <- as.factor(data$isced)

set.seed(42)

## ***************************************************************************
#1) Analysis ####
#*Registry ####
reg <-  makeExperimentRegistry(file.dir="Realdata", seed = 1)
#reg <- loadRegistry(file.dir="Realdata", writeable=TRUE)

#*Problems ####
source("Problems.R")
addProblem(name="realdata", data=data, seed=43)

#*Algorithms ####
source("Algorithms.R")
addAlgorithm(name="rf_var", fun=rf_var_real_wrapper)
addAlgorithm(name="rf_zhao", fun=rf_zhao_real_wrapper)
addAlgorithm(name="rf_maxstat", fun=rf_maxstat_real_wrapper)
addAlgorithm(name="pred_rf_base", fun=pred_rf_base_real_wrapper)
addAlgorithm(name="pred_rf_var", fun=pred_rf_var_real_wrapper)
addAlgorithm(name="pred_rf_zhao", fun=pred_rf_zhao_real_wrapper)
addAlgorithm(name="pred_rf_maxstat", fun=pred_rf_maxstat_real_wrapper)

#*Experiments ####
ntree <- 200
nsim <- 10
num_threads <- 12

##**Tune alpha ####
prob_design <- list(realdata=expand.grid(stringsAsFactors = FALSE))
algo_design <- list(pred_rf_maxstat = data.frame(alpha=c(0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), mtry=33, num_trees=ntree,
                                            minsplit=50, minbucket=1, predleaf="Meanresid",
                                            num_threads=num_threads, stringsAsFactors = FALSE))
addExperiments(prob.designs=prob_design, algo.designs=algo_design, repls=1)

#**VIM ####
prob_design <- list(realdata=expand.grid(stringsAsFactors = FALSE))
algo_design <- list(rf_var = expand.grid(minsplit=50, mtry=33, num_trees=ntree, num_threads=num_threads, stringsAsFactors = FALSE),
                    rf_zhao = expand.grid(minsplit=50, mtry=33, num_trees=ntree, num_threads=num_threads, stringsAsFactors = FALSE),
                    rf_maxstat = data.frame(alpha=0.2, mtry=33, num_trees=ntree, minsplit=50, minbucket=1,
                                            predleaf="Meanresid", num_threads=num_threads, stringsAsFactors = FALSE))
addExperiments(prob.designs=prob_design, algo.designs=algo_design, repls=nsim)

#**Prediction ####
prob_design <- list(realdata=expand.grid(stringsAsFactors = FALSE))
algo_design <- list(pred_rf_base = expand.grid(minsplit=50, mtry=1, num_trees=ntree, num_threads=num_threads, stringsAsFactors = FALSE),
		      pred_rf_var = expand.grid(minsplit=50, mtry=33, num_trees=ntree, num_threads=num_threads, stringsAsFactors = FALSE),
                    pred_rf_zhao = expand.grid(minsplit=50, mtry=33, num_trees=ntree, num_threads=num_threads, stringsAsFactors = FALSE),
                    pred_rf_maxstat = data.frame(alpha=0.2, mtry=33, num_trees=ntree, minsplit=50, minbucket=1,
                                            predleaf="Meanresid", num_threads=num_threads, stringsAsFactors = FALSE))
addExperiments(prob.designs=prob_design, algo.designs=algo_design, repls=nsim)

#*Test jobs ####
#Test one job for each algorithm
#-VIM
#testJob(id=head(findExperiments(prob.name = "realdata", algo.name = "rf_var"), 1))
#testJob(id=head(findExperiments(prob.name = "realdata", algo.name = "rf_zhao"), 1))
#testJob(id=head(findExperiments(prob.name = "realdata", algo.name = "rf_maxstat", algo.pars = (alpha == 0.05)), 1))
#-Prediction
#testJob(id=head(findExperiments(prob.name = "realdata", algo.name = "pred_rf_var"), 1))
#testJob(id=head(findExperiments(prob.name = "realdata", algo.name = "pred_rf_zhao"), 1))
#testJob(id=head(findExperiments(prob.name = "realdata", algo.name = "pred_rf_maxstat", algo.pars = (alpha == 0.05)), 1))

if (grepl("node\\d{2}|bipscluster", system("hostname", intern = TRUE))) {
  ids <- findNotStarted()
  #ids <- findNotDone()
  ids[, chunk := chunk(job.id, chunk.size = 5)]
  submitJobs(ids = ids, # walltime in seconds, 10 days max, memory in MB
             resources = list(name = "SimSimple", chunks.as.arrayjobs = TRUE, 
                              ncpus = 12, memory = 6000, walltime = 10*24*3600, 
                              max.concurrent.jobs = 10))
} else {
  ids <- findNotDone()
  submitJobs(ids)
}
waitForJobs()

#*Get results ####
#the results are automatically saved in the registry directory
ids <- findExperiments()
res <-  flatten(ijoin(getJobPars(), reduceResultsDataTable(ids)))
save(res,file="Res_real.RData")


## ***************************************************************************
#2) Runtime ####
#*Registry ####
reg <-  makeExperimentRegistry(file.dir="Runtime", seed = 1)
#reg <- loadRegistry(file.dir="Runtime", writeable=TRUE)

#*Problems ####
source("Problems.R")
addProblem(name="realdata", data=data, seed=43)

#*Algorithms ####
source("Algorithms.R")
addAlgorithm(name="run_rf_var", fun=run_rf_var_real_wrapper)
addAlgorithm(name="run_rf_zhao", fun=run_rf_zhao_real_wrapper)
addAlgorithm(name="run_rf_maxstat", fun=run_rf_maxstat_real_wrapper)

#*Experiments ####
ntree <- 200
nsim <- 10
num_threads <- 1
minsplit <- 50
alpha <- 0.2
mtry <- 33

#**Runtime ####
prob_design <- list(realdata=expand.grid(stringsAsFactors = FALSE))
algo_design <- list(run_rf_var = expand.grid(minsplit=minsplit, mtry=mtry, num_trees=ntree, num_threads=num_threads, stringsAsFactors = FALSE),
                    run_rf_zhao = expand.grid(minsplit=minsplit, mtry=mtry, num_trees=ntree, num_threads=num_threads, stringsAsFactors = FALSE),
                    run_rf_maxstat = data.frame(alpha=alpha, mtry=mtry, num_trees=ntree, minsplit=minsplit, minbucket=1,
                                                predleaf="Meanresid", num_threads=num_threads, stringsAsFactors = FALSE))
addExperiments(prob.designs=prob_design, algo.designs=algo_design, repls=nsim)

#*Test jobs ####
#Test one job for each algorithm
#-Runtime
#testJob(id=head(findExperiments(prob.name = "realdata", algo.name = "run_rf_var"), 1))
#testJob(id=head(findExperiments(prob.name = "realdata", algo.name = "run_rf_var"), 1))
#testJob(id=head(findExperiments(prob.name = "realdata", algo.name = "run_rf_var"), 1))

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
#the results are automatically saved in the registry directory
ids <- findExperiments()
res <-  flatten(ijoin(getJobPars(), reduceResultsDataTable(ids)))
save(res,file="Res_runtime.RData")






