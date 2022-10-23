## ***************************************************************************
## SimGenetic_Plots.R #####
##
## Script name: SimGenetic_Plots
##
## Purpose of script: Plot the results from the Genetic simulation
##
## Author: Annika Swenne
##
## Date Created: 2022-06-20
##
## Version/ Date: V01 / 2022-06-20
##                V02 / 2022-09-02
##
## ***************************
##
## Modification Notes:
## V02/ 2022-09-02: Different plots
##   
## ***************************
##
## Notes/ to do:
##
## ***************************************************************************

#load the required packages:  (uncomment as required)
library(batchtools)
library(ggplot2)
library(dplyr)

set.seed(42)

## ***************************************************************************
load(file="C:/Users/Annika Swenne/Documents/Uni_Bremen/MedicalBiometry/Masterarbeit/Code/Batchtools/Res_snp.RData")

npred <- 100
predictors <- paste(rep("x", times=npred), seq(1,npred), sep="")

#change the structure of the results for the vim to long format to visualize them with ggplot2
reslong <- res[which(!is.na(res$x1)),]
if ("error" %in% colnames(res)){
  reslong <- reshape(reslong, varying=predictors, direction="long", v.names="vim", times=predictors, drop=c("z1","error","accuracy"))
} else {
  reslong <- reshape(reslong, varying=predictors, direction="long", v.names="vim", times=predictors)
}
reslong <- reslong[order(reslong$job.id),]
reslong <- subset(reslong, select=-id)
colnames(reslong)[names(reslong)=="time"] <- "predictor"

reslong$predictor <- factor(reslong$predictor, levels=unique(reslong$predictor))

#assign ranks to the VIM (highest VIM gets rank 1)
ranks <- res
temp <- t(apply(-ranks[,predictors,with=FALSE],1,FUN = rank, ties.method="random"))
colnames(temp) <- predictors
for(i in 1:npred){
  ind  <- predictors[i]
  ranks[,which(colnames(ranks)==ind)] <- temp[,i]
}

#change the structure of the results to long format to visualize them with ggplot2
rankslong <- ranks[which(!is.na(ranks$x1)),]
if ("error" %in% colnames(ranks)){
  rankslong <- reshape(rankslong, varying=predictors, direction="long", v.names="rank", times=predictors, drop=c("z1","error","accuracy"))
} else {
  rankslong <- reshape(rankslong, varying=predictors, direction="long", v.names="rank", times=predictors)
}
rankslong <- rankslong[order(rankslong$job.id),]
rankslong <- subset(rankslong, select=-id)
colnames(rankslong)[names(rankslong)=="time"] <- "predictor"

rankslong$predictor <- factor(rankslong$predictor, levels=unique(rankslong$predictor))

#create dataset for the prediction accuracy
respred <- res[which(!is.na(res$error)),-which(colnames(res) %in% predictors),with=FALSE]
#rename the algorithms
respred$rf <- NA
respred$rf[which(respred$algorithm=="pred_rf_var")] <- "Variance \nsplitting"
respred$rf[which(respred$algorithm=="pred_rf_zhao")] <- "Zhao et al. \n(2012)"

respred$rf[which(respred$algorithm=="pred_rf_maxstat" & respred$predleaf=="Meanresid" & respred$alpha==0.05)] <- "Maxstat \nsplitting \nalpha=0.05"
respred$rf[which(respred$algorithm=="pred_rf_maxstat" & respred$predleaf=="Meanresid" & respred$alpha==0.3)] <- "Maxstat \nsplitting \nalpha=0.3"
respred$rf[which(respred$algorithm=="pred_rf_maxstat" & respred$predleaf=="Meanresid" & respred$alpha==0.5)] <- "Maxstat \nsplitting \nalpha=0.5"
respred$rf[which(respred$algorithm=="pred_rf_maxstat" & respred$predleaf=="Meanresid" & respred$alpha==0.7)] <- "Maxstat \nsplitting \nalpha=0.7"

respred$rf <- factor(respred$rf, levels=unique(respred$rf))

setwd("C:/Users/Annika Swenne/Documents/Uni_Bremen/MedicalBiometry/Masterarbeit/Latex_MA")

#*21) Ranks Null case ----
id_zhao <- which(rankslong$problem=="simsnp_null" & !(rankslong$predictor %in% c("x1","x2","x3","x4","x5"))
                 & rankslong$algorithm=="rf_zhao")
id_var <- which(rankslong$problem=="simsnp_null" & !(rankslong$predictor %in% c("x1","x2","x3","x4","x5"))
                & rankslong$algorithm=="rf_var")
id_max <- which(rankslong$problem=="simsnp_null" & !(rankslong$predictor %in% c("x1","x2","x3","x4","x5"))
                & rankslong$algorithm=="rf_maxstat" & rankslong$alpha==0.05)
zhao <- data.frame(algorithm=rep("rf_zhao",times=npred),
                   prop.table(table(factor(rankslong$rank[id_zhao],levels=1:npred))))
var <- data.frame(algorithm=rep("rf_var",times=npred),
                  prop.table(table(factor(rankslong$rank[id_var],levels=1:npred))))
max <- data.frame(algorithm=rep("rf_max",times=npred),
                  prop.table(table(factor(rankslong$rank[id_max],levels=1:npred))))
rankdata <- rbind(var,zhao,max)
rankdata$SNP <- "Irrelevant SNPs"
id_zhao <- which(rankslong$problem=="simsnp_null" & rankslong$predictor %in% c("x1","x2","x3","x4","x5")
                 & rankslong$algorithm=="rf_zhao")
id_var <- which(rankslong$problem=="simsnp_null"& rankslong$predictor %in% c("x1","x2","x3","x4","x5")
                & rankslong$algorithm=="rf_var")
id_max <- which(rankslong$problem=="simsnp_null" & rankslong$predictor %in% c("x1","x2","x3","x4","x5") 
                & rankslong$algorithm=="rf_maxstat" & rankslong$alpha==0.05)
zhao <- data.frame(algorithm=rep("rf_zhao",times=npred),
                   prop.table(table(factor(rankslong$rank[id_zhao],levels=1:npred))))
var <- data.frame(algorithm=rep("rf_var",times=npred),
                  prop.table(table(factor(rankslong$rank[id_var],levels=1:npred))))
max <- data.frame(algorithm=rep("rf_max",times=npred),
                  prop.table(table(factor(rankslong$rank[id_max],levels=1:npred))))
temp <- rbind(var,zhao,max)
temp$SNP <- "Confounded SNPs"
rankdata <- rbind(rankdata,temp)
colnames(rankdata) <- c("algorithm","rank","proportion","SNP")
rankdata$rank <- as.numeric(rankdata$rank)

#create names for the algorithms
rankdata$rf <- NA
rankdata$rf[which(rankdata$algorithm=="rf_var")] <- "Variance splitting"
rankdata$rf[which(rankdata$algorithm=="rf_zhao")] <- "Zhao et al. (2012)"
rankdata$rf[which(rankdata$algorithm=="rf_max")] <- "Maximally selected \nresidual rank statistic"
rankdata$rf <- factor(rankdata$rf, levels=unique(rankdata$rf))

ggplot(rankdata, aes(x=rank, y=proportion)) + 
  geom_bar(stat="identity", position=position_dodge(), width=0.5) + 
  facet_grid(cols=vars(rf), rows=vars(SNP), scales="free") +
  theme_bw() + 
  labs(x="Rank", y="Proportion")
ggsave("img/Figure21_null.png", height=4, width=7)

#*22) Ranks Null case: Alpha ----
#**a) Ranks ----
#id_alpha005 <- which(rankslong$problem=="simsnp_null" & !(rankslong$predictor %in% c("x1","x2","x3","x4","x5"))
#                & rankslong$algorithm=="rf_maxstat" & rankslong$alpha==0.05)
id_alpha03 <- which(rankslong$problem=="simsnp_null" & !(rankslong$predictor %in% c("x1","x2","x3","x4","x5"))
                      & rankslong$algorithm=="rf_maxstat" & rankslong$alpha==0.3)
id_alpha05 <- which(rankslong$problem=="simsnp_null" & !(rankslong$predictor %in% c("x1","x2","x3","x4","x5"))
                      & rankslong$algorithm=="rf_maxstat" & rankslong$alpha==0.5)
id_alpha07 <- which(rankslong$problem=="simsnp_null" & !(rankslong$predictor %in% c("x1","x2","x3","x4","x5"))
                      & rankslong$algorithm=="rf_maxstat" & rankslong$alpha==0.7)
#alpha005 <- data.frame(alpha=rep("alpha: 0.05", times=npred),
#                   prop.table(table(factor(rankslong$rank[id_alpha005],levels=1:npred))))
alpha03 <- data.frame(alpha=rep("alpha: 0.3", times=npred),
                       prop.table(table(factor(rankslong$rank[id_alpha03],levels=1:npred))))
alpha05 <- data.frame(alpha=rep("alpha: 0.5", times=npred),
                      prop.table(table(factor(rankslong$rank[id_alpha05],levels=1:npred))))
alpha07 <- data.frame(alpha=rep("alpha: 0.7", times=npred),
                       prop.table(table(factor(rankslong$rank[id_alpha07],levels=1:npred))))
#rankdata <- rbind(alpha005,alpha03,alpha05,alpha07)
rankdata <- rbind(alpha03,alpha05,alpha07)
rankdata$SNP <- "Irrelevant SNPs"
#id_alpha005 <- which(rankslong$problem=="simsnp_null" & rankslong$predictor %in% c("x1","x2","x3","x4","x5")
#                     & rankslong$algorithm=="rf_maxstat" & rankslong$alpha==0.05)
id_alpha03 <- which(rankslong$problem=="simsnp_null" & rankslong$predictor %in% c("x1","x2","x3","x4","x5")
                    & rankslong$algorithm=="rf_maxstat" & rankslong$alpha==0.3)
id_alpha05 <- which(rankslong$problem=="simsnp_null" & rankslong$predictor %in% c("x1","x2","x3","x4","x5")
                    & rankslong$algorithm=="rf_maxstat" & rankslong$alpha==0.5)
id_alpha07 <- which(rankslong$problem=="simsnp_null" & rankslong$predictor %in% c("x1","x2","x3","x4","x5")
                    & rankslong$algorithm=="rf_maxstat" & rankslong$alpha==0.7)
#alpha005 <- data.frame(alpha=rep("alpha: 0.05", times=npred),
#                       prop.table(table(factor(rankslong$rank[id_alpha005],levels=1:npred))))
alpha03 <- data.frame(alpha=rep("alpha: 0.3", times=npred),
                      prop.table(table(factor(rankslong$rank[id_alpha03],levels=1:npred))))
alpha05 <- data.frame(alpha=rep("alpha: 0.5", times=npred),
                      prop.table(table(factor(rankslong$rank[id_alpha05],levels=1:npred))))
alpha07 <- data.frame(alpha=rep("alpha: 0.7", times=npred),
                      prop.table(table(factor(rankslong$rank[id_alpha07],levels=1:npred))))
#temp <- rbind(alpha005,alpha03,alpha05,alpha07)
temp <- rbind(alpha03,alpha05,alpha07)
temp$SNP <- "Confounded SNPs"
rankdata <- rbind(rankdata,temp)
colnames(rankdata) <- c("alpha","rank","proportion","SNP")
rankdata$rank <- as.numeric(rankdata$rank)

ggplot(rankdata, aes(x=rank, y=proportion)) + 
  geom_bar(stat="identity", position=position_dodge(), width=0.5) + 
  facet_grid(cols=vars(alpha), rows=vars(SNP), scales="free") +
  theme_bw() + 
  labs(x="Rank", y="Proportion")
ggsave("img/Figure22_null_alpha.png", height=4, width=7)

#**b) VIM ----
#randomly select 5 irrelevant SNPs
set.seed(11)
sample(seq(11,100),size=5) #44 66 35 26 47
id <- which(reslong$problem=="simsnp_null" & reslong$algorithm=="rf_maxstat"
             & reslong$predictor %in% c("x1","x2","x3","x4","x5","x44","x66","x35","x26","x47"))
ggplot(reslong[id,], aes(x=predictor, y=vim)) + 
stat_boxplot(geom = "errorbar", width = 0.3) +
  geom_boxplot(fill="lightgrey", notch=FALSE) + 
  geom_hline(yintercept = 0, linetype="dotted") +
  facet_wrap(~alpha, ncol=1, scales="free", labeller=label_both) +
  theme_bw() + 
  labs(x="Predictor", y="Permutation importance")
  #labs(title="Nullcase: Zhao et al. (nsim=100)", 
  #   x="Predictor", y="Permutation importance")
ggsave("img/Figure22_null_vim.png", height=8, width=7)


#*23) MSE Null case ----
id <- which(respred$problem=="pred_simsnp_null")
ggplot(respred[id,], aes(x=rf, y=error)) + 
  stat_boxplot(geom = "errorbar", width = 0.3) +
  geom_boxplot(fill="lightgrey", notch=FALSE) + 
  theme_bw() + 
  labs(x="Algorithm", y="MSE")
ggsave("img/Figure23_null_mse.png", height=3, width=7)


#*24) Ranks Power case ----
betax2 <- 20
id_zhao <- which(rankslong$problem=="simsnp_power" & !(rankslong$predictor %in% c("x1","x2","x3","x4","x5","x6","x7","x8","x9","x10"))
                 & rankslong$algorithm=="rf_zhao" & rankslong$betax2==betax2)
id_var <- which(rankslong$problem=="simsnp_power" & !(rankslong$predictor %in% c("x1","x2","x3","x4","x5","x6","x7","x8","x9","x10"))
                & rankslong$algorithm=="rf_var" & rankslong$betax2==betax2)
id_max <- which(rankslong$problem=="simsnp_power" & !(rankslong$predictor %in% c("x1","x2","x3","x4","x5","x6","x7","x8","x9","x10"))
                & rankslong$algorithm=="rf_maxstat" & rankslong$alpha==0.05 & rankslong$betax2==betax2)
zhao <- data.frame(algorithm=rep("rf_zhao",times=npred),
                   prop.table(table(factor(rankslong$rank[id_zhao],levels=1:npred))))
var <- data.frame(algorithm=rep("rf_var",times=npred),
                  prop.table(table(factor(rankslong$rank[id_var],levels=1:npred))))
max <- data.frame(algorithm=rep("rf_max",times=npred),
                  prop.table(table(factor(rankslong$rank[id_max],levels=1:npred))))
rankdata <- rbind(zhao,var,max)
rankdata$SNP <- "Irrelevant SNPs"
id_zhao <- which(rankslong$problem=="simsnp_power" & rankslong$predictor %in% c("x1","x2","x3","x4","x5")
                 & rankslong$algorithm=="rf_zhao" & rankslong$betax2==betax2)
id_var <- which(rankslong$problem=="simsnp_power"& rankslong$predictor %in% c("x1","x2","x3","x4","x5")
                & rankslong$algorithm=="rf_var" & rankslong$betax2==betax2)
id_max <- which(rankslong$problem=="simsnp_power" & rankslong$predictor %in% c("x1","x2","x3","x4","x5") 
                & rankslong$algorithm=="rf_maxstat" & rankslong$alpha==0.05 & rankslong$betax2==betax2)
zhao <- data.frame(algorithm=rep("rf_zhao",times=npred),
                   prop.table(table(factor(rankslong$rank[id_zhao],levels=1:npred))))
var <- data.frame(algorithm=rep("rf_var",times=npred),
                  prop.table(table(factor(rankslong$rank[id_var],levels=1:npred))))
max <- data.frame(algorithm=rep("rf_max",times=npred),
                  prop.table(table(factor(rankslong$rank[id_max],levels=1:npred))))
temp <- rbind(zhao,var,max)
temp$SNP <- "Confounded SNPs"
rankdata <- rbind(rankdata,temp)
id_zhao <- which(rankslong$problem=="simsnp_power" & rankslong$predictor %in% c("x6","x7","x8","x9","x10")
                 & rankslong$algorithm=="rf_zhao" & rankslong$betax2==betax2)
id_var <- which(rankslong$problem=="simsnp_power"& rankslong$predictor %in% c("x6","x7","x8","x9","x10")
                & rankslong$algorithm=="rf_var" & rankslong$betax2==betax2)
id_max <- which(rankslong$problem=="simsnp_power" & rankslong$predictor %in% c("x6","x7","x8","x9","x10") 
                & rankslong$algorithm=="rf_maxstat" & rankslong$alpha==0.05 & rankslong$betax2==betax2)
zhao <- data.frame(algorithm=rep("rf_zhao",times=npred),
                   prop.table(table(factor(rankslong$rank[id_zhao],levels=1:npred))))
var <- data.frame(algorithm=rep("rf_var",times=npred),
                  prop.table(table(factor(rankslong$rank[id_var],levels=1:npred))))
max <- data.frame(algorithm=rep("rf_max",times=npred),
                  prop.table(table(factor(rankslong$rank[id_max],levels=1:npred))))
temp <- rbind(zhao,var,max)
temp$SNP <- "Causal SNPs"
rankdata <- rbind(rankdata,temp)
colnames(rankdata) <- c("algorithm","rank","proportion","SNP")
rankdata$rank <- as.numeric(rankdata$rank)

#create names for the algorithms
rankdata$rf <- NA
rankdata$rf[which(rankdata$algorithm=="rf_var")] <- "Variance splitting"
rankdata$rf[which(rankdata$algorithm=="rf_zhao")] <- "Zhao et al. (2012)"
rankdata$rf[which(rankdata$algorithm=="rf_max")] <- "Maximally selected \nresidual rank statistic"
rankdata$rf <- factor(rankdata$rf, levels=unique(rankdata$rf))

ggplot(rankdata, aes(x=rank, y=proportion)) + 
  geom_bar(stat="identity", position=position_dodge(), width=0.5) + 
  facet_grid(cols=vars(rf), rows=vars(SNP), scales="free") +
  theme_bw() + 
  labs(x="Rank", y="Proportion")
ggsave("img/Figure24_power.png", height=6, width=7)


#*25) Ranks Power case: Alpha ----
betax2 <- 20
#**a) Ranks ----
#id_alpha005 <- which(rankslong$problem=="simsnp_power" & !(rankslong$predictor %in% c("x1","x2","x3","x4","x5","x6","x7","x8","x9","x10"))
#                     & rankslong$algorithm=="rf_maxstat" & rankslong$betax2==betax2 & rankslong$alpha==0.05)
id_alpha03 <- which(rankslong$problem=="simsnp_power" & !(rankslong$predictor %in% c("x1","x2","x3","x4","x5","x6","x7","x8","x9","x10"))
                    & rankslong$algorithm=="rf_maxstat" & rankslong$betax2==betax2 & rankslong$alpha==0.3)
id_alpha05 <- which(rankslong$problem=="simsnp_power" & !(rankslong$predictor %in% c("x1","x2","x3","x4","x5","x6","x7","x8","x9","x10"))
                    & rankslong$algorithm=="rf_maxstat" & rankslong$betax2==betax2 & rankslong$alpha==0.5)
id_alpha07 <- which(rankslong$problem=="simsnp_power" & !(rankslong$predictor %in% c("x1","x2","x3","x4","x5","x6","x7","x8","x9","x10"))
                    & rankslong$algorithm=="rf_maxstat" & rankslong$betax2==betax2 & rankslong$alpha==0.7)
#alpha005 <- data.frame(alpha=rep("alpha: 0.05", times=npred),
#                       prop.table(table(factor(rankslong$rank[id_alpha005],levels=1:npred))))
alpha03 <- data.frame(alpha=rep("alpha: 0.3", times=npred),
                      prop.table(table(factor(rankslong$rank[id_alpha03],levels=1:npred))))
alpha05 <- data.frame(alpha=rep("alpha: 0.5", times=npred),
                      prop.table(table(factor(rankslong$rank[id_alpha05],levels=1:npred))))
alpha07 <- data.frame(alpha=rep("alpha: 0.7", times=npred),
                      prop.table(table(factor(rankslong$rank[id_alpha07],levels=1:npred))))
rankdata <- rbind(alpha03,alpha05,alpha07)
#rankdata <- rbind(alpha005,alpha03,alpha05,alpha07)
rankdata$SNP <- "Irrelevant SNPs"
#id_alpha005 <- which(rankslong$problem=="simsnp_power" & rankslong$predictor %in% c("x1","x2","x3","x4","x5")
#                     & rankslong$algorithm=="rf_maxstat" & rankslong$betax2==betax2 & rankslong$alpha==0.05)
id_alpha03 <- which(rankslong$problem=="simsnp_power" & rankslong$predictor %in% c("x1","x2","x3","x4","x5")
                    & rankslong$algorithm=="rf_maxstat" & rankslong$betax2==betax2 & rankslong$alpha==0.3)
id_alpha05 <- which(rankslong$problem=="simsnp_power" & rankslong$predictor %in% c("x1","x2","x3","x4","x5")
                    & rankslong$algorithm=="rf_maxstat" & rankslong$betax2==betax2 & rankslong$alpha==0.5)
id_alpha07 <- which(rankslong$problem=="simsnp_power" & rankslong$predictor %in% c("x1","x2","x3","x4","x5")
                    & rankslong$algorithm=="rf_maxstat" & rankslong$betax2==betax2 & rankslong$alpha==0.7)
#alpha005 <- data.frame(alpha=rep("alpha: 0.05", times=npred),
#                       prop.table(table(factor(rankslong$rank[id_alpha005],levels=1:npred))))
alpha03 <- data.frame(alpha=rep("alpha: 0.3", times=npred),
                      prop.table(table(factor(rankslong$rank[id_alpha03],levels=1:npred))))
alpha05 <- data.frame(alpha=rep("alpha: 0.5", times=npred),
                      prop.table(table(factor(rankslong$rank[id_alpha05],levels=1:npred))))
alpha07 <- data.frame(alpha=rep("alpha: 0.7", times=npred),
                      prop.table(table(factor(rankslong$rank[id_alpha07],levels=1:npred))))
#temp <- rbind(alpha005,alpha03,alpha05,alpha07)
temp <- rbind(alpha03,alpha05,alpha07)
temp$SNP <- "Confounded SNPs"
rankdata <- rbind(rankdata,temp)
#id_alpha005 <- which(rankslong$problem=="simsnp_power" & rankslong$predictor %in% c("x6","x7","x8","x9","x10")
#                    & rankslong$algorithm=="rf_maxstat" & rankslong$betax2==betax2 & rankslong$alpha==0.05)
id_alpha03 <- which(rankslong$problem=="simsnp_power" & rankslong$predictor %in% c("x6","x7","x8","x9","x10")
                    & rankslong$algorithm=="rf_maxstat" & rankslong$betax2==betax2 & rankslong$alpha==0.3)
id_alpha05 <- which(rankslong$problem=="simsnp_power" & rankslong$predictor %in% c("x6","x7","x8","x9","x10")
                    & rankslong$algorithm=="rf_maxstat" & rankslong$betax2==betax2 & rankslong$alpha==0.5)
id_alpha07 <- which(rankslong$problem=="simsnp_power" & rankslong$predictor %in% c("x6","x7","x8","x9","x10")
                    & rankslong$algorithm=="rf_maxstat" & rankslong$betax2==betax2 & rankslong$alpha==0.7)
#alpha005 <- data.frame(alpha=rep("alpha: 0.05", times=npred),
#                       prop.table(table(factor(rankslong$rank[id_alpha005],levels=1:npred))))
alpha03 <- data.frame(alpha=rep("alpha: 0.3", times=npred),
                      prop.table(table(factor(rankslong$rank[id_alpha03],levels=1:npred))))
alpha05 <- data.frame(alpha=rep("alpha: 0.5", times=npred),
                      prop.table(table(factor(rankslong$rank[id_alpha05],levels=1:npred))))
alpha07 <- data.frame(alpha=rep("alpha: 0.7", times=npred),
                      prop.table(table(factor(rankslong$rank[id_alpha07],levels=1:npred))))
temp <- rbind(alpha03,alpha05,alpha07)
#temp <- rbind(alpha005,alpha03,alpha05,alpha07)
temp$SNP <- "Causal SNPs"
rankdata <- rbind(rankdata,temp)
colnames(rankdata) <- c("alpha","rank","proportion","SNP")
rankdata$rank <- as.numeric(rankdata$rank)

ggplot(rankdata, aes(x=rank, y=proportion)) + 
  geom_bar(stat="identity", position=position_dodge(), width=0.5) + 
  facet_grid(cols=vars(alpha), rows=vars(SNP), scales="free") +
  theme_bw() + 
  labs(x="Rank", y="Proportion")
ggsave("img/Figure25_power_alpha.png", height=6, width=7)

##**b) VIM ----
#randomly select 5 irrelevant SNPs
set.seed(11)
sample(seq(11,100),size=5) #44 66 35 26 47
id <- which(reslong$problem=="simsnp_power" & reslong$algorithm=="rf_maxstat" & reslong$betax2==20
            & reslong$predictor %in% c("x1","x2","x3","x4","x5","x6","x7","x8","x9","x10","x44","x66","x35","x26","x47"))
ggplot(reslong[id,], aes(x=predictor, y=vim)) + 
  stat_boxplot(geom = "errorbar", width = 0.3) +
  geom_boxplot(fill="lightgrey", notch=FALSE) + 
  geom_hline(yintercept = 0, linetype="dotted") +
  facet_wrap(~alpha, ncol=1, scales="free", labeller=label_both) +
  theme_bw() + 
  labs(x="Predictor", y="Permutation importance")
ggsave("img/Figure25_power_vim.png", height=8, width=7)

#close-up
set.seed(11)
sample(seq(11,100),size=5) #44 66 35 26 47
idmax <- which(reslong$problem=="simsnp_power" & reslong$algorithm=="rf_maxstat" & reslong$betax2==20 & reslong$alpha==0.05
               & reslong$predictor %in% c("x1","x2","x3","x4","x5","x6","x7","x8","x9","x10","x44","x66","x35","x26","x47"))
ggplot(reslong[idmax,], aes(x=predictor, y=vim)) + 
  stat_boxplot(geom = "errorbar", width = 0.3) +
  geom_boxplot(fill="lightgrey", notch=FALSE) + 
  geom_hline(yintercept = 0, linetype="dotted") +
  facet_wrap(~alpha, ncol=1, scales="free", labeller=label_both) +
  theme_bw() + 
  labs(x="Predictor", y="Permutation importance")
ggsave("img/Figure25_power_vim_closeup.png", height=4, width=7)


#*26) Ranks Power case: Sensitivity ----
betax2 <- 4
id_zhao <- which(rankslong$problem=="simsnp_power" & !(rankslong$predictor %in% c("x1","x2","x3","x4","x5","x6","x7","x8","x9","x10"))
                 & rankslong$algorithm=="rf_zhao" & rankslong$betax2==betax2)
id_var <- which(rankslong$problem=="simsnp_power" & !(rankslong$predictor %in% c("x1","x2","x3","x4","x5","x6","x7","x8","x9","x10"))
                & rankslong$algorithm=="rf_var" & rankslong$betax2==betax2)
id_max <- which(rankslong$problem=="simsnp_power" & !(rankslong$predictor %in% c("x1","x2","x3","x4","x5","x6","x7","x8","x9","x10"))
                & rankslong$algorithm=="rf_maxstat" & rankslong$alpha==0.05 & rankslong$betax2==betax2)
zhao <- data.frame(algorithm=rep("rf_zhao",times=npred),
                   prop.table(table(factor(rankslong$rank[id_zhao],levels=1:npred))))
var <- data.frame(algorithm=rep("rf_var",times=npred),
                  prop.table(table(factor(rankslong$rank[id_var],levels=1:npred))))
max <- data.frame(algorithm=rep("rf_max",times=npred),
                  prop.table(table(factor(rankslong$rank[id_max],levels=1:npred))))
rankdata <- rbind(zhao,var,max)
rankdata$SNP <- "Irrelevant SNPs"
id_zhao <- which(rankslong$problem=="simsnp_power" & rankslong$predictor %in% c("x1","x2","x3","x4","x5")
                 & rankslong$algorithm=="rf_zhao" & rankslong$betax2==betax2)
id_var <- which(rankslong$problem=="simsnp_power"& rankslong$predictor %in% c("x1","x2","x3","x4","x5")
                & rankslong$algorithm=="rf_var" & rankslong$betax2==betax2)
id_max <- which(rankslong$problem=="simsnp_power" & rankslong$predictor %in% c("x1","x2","x3","x4","x5") 
                & rankslong$algorithm=="rf_maxstat" & rankslong$alpha==0.05 & rankslong$betax2==betax2)
zhao <- data.frame(algorithm=rep("rf_zhao",times=npred),
                   prop.table(table(factor(rankslong$rank[id_zhao],levels=1:npred))))
var <- data.frame(algorithm=rep("rf_var",times=npred),
                  prop.table(table(factor(rankslong$rank[id_var],levels=1:npred))))
max <- data.frame(algorithm=rep("rf_max",times=npred),
                  prop.table(table(factor(rankslong$rank[id_max],levels=1:npred))))
temp <- rbind(zhao,var,max)
temp$SNP <- "Confounded SNPs"
rankdata <- rbind(rankdata,temp)
id_zhao <- which(rankslong$problem=="simsnp_power" & rankslong$predictor %in% c("x6","x7","x8","x9","x10")
                 & rankslong$algorithm=="rf_zhao" & rankslong$betax2==betax2)
id_var <- which(rankslong$problem=="simsnp_power"& rankslong$predictor %in% c("x6","x7","x8","x9","x10")
                & rankslong$algorithm=="rf_var" & rankslong$betax2==betax2)
id_max <- which(rankslong$problem=="simsnp_power" & rankslong$predictor %in% c("x6","x7","x8","x9","x10") 
                & rankslong$algorithm=="rf_maxstat" & rankslong$alpha==0.05 & rankslong$betax2==betax2)
zhao <- data.frame(algorithm=rep("rf_zhao",times=npred),
                   prop.table(table(factor(rankslong$rank[id_zhao],levels=1:npred))))
var <- data.frame(algorithm=rep("rf_var",times=npred),
                  prop.table(table(factor(rankslong$rank[id_var],levels=1:npred))))
max <- data.frame(algorithm=rep("rf_max",times=npred),
                  prop.table(table(factor(rankslong$rank[id_max],levels=1:npred))))
temp <- rbind(zhao,var,max)
temp$SNP <- "Causal SNPs"
rankdata <- rbind(rankdata,temp)
colnames(rankdata) <- c("algorithm","rank","proportion","SNP")
rankdata$rank <- as.numeric(rankdata$rank)

#create names for the algorithms
rankdata$rf <- NA
rankdata$rf[which(rankdata$algorithm=="rf_var")] <- "Variance splitting"
rankdata$rf[which(rankdata$algorithm=="rf_zhao")] <- "Zhao et al. (2012)"
rankdata$rf[which(rankdata$algorithm=="rf_max")] <- "Maximally selected \nresidual rank statistic"
rankdata$rf <- factor(rankdata$rf, levels=unique(rankdata$rf))

ggplot(rankdata, aes(x=rank, y=proportion)) + 
  geom_bar(stat="identity", position=position_dodge(), width=0.5) + 
  facet_grid(cols=vars(rf), rows=vars(SNP), scales="free") +
  theme_bw() + 
  labs(x="Rank", y="Proportion")
ggsave("img/Figure26_power_sensitivity.png", height=6, width=7)


#*27) Ranks Power case: Sensitivity Alpha----
betax2 <- 4
#id_alpha005 <- which(rankslong$problem=="simsnp_power" & !(rankslong$predictor %in% c("x1","x2","x3","x4","x5","x6","x7","x8","x9","x10"))
#                     & rankslong$algorithm=="rf_maxstat" & rankslong$betax2==betax2 & rankslong$alpha==0.05)
id_alpha03 <- which(rankslong$problem=="simsnp_power" & !(rankslong$predictor %in% c("x1","x2","x3","x4","x5","x6","x7","x8","x9","x10"))
                    & rankslong$algorithm=="rf_maxstat" & rankslong$betax2==betax2 & rankslong$alpha==0.3)
id_alpha05 <- which(rankslong$problem=="simsnp_power" & !(rankslong$predictor %in% c("x1","x2","x3","x4","x5","x6","x7","x8","x9","x10"))
                    & rankslong$algorithm=="rf_maxstat" & rankslong$betax2==betax2 & rankslong$alpha==0.5)
id_alpha07 <- which(rankslong$problem=="simsnp_power" & !(rankslong$predictor %in% c("x1","x2","x3","x4","x5","x6","x7","x8","x9","x10"))
                    & rankslong$algorithm=="rf_maxstat" & rankslong$betax2==betax2 & rankslong$alpha==0.7)
#alpha005 <- data.frame(alpha=rep("alpha: 0.05", times=npred),
#                       prop.table(table(factor(rankslong$rank[id_alpha005],levels=1:npred))))
alpha03 <- data.frame(alpha=rep("alpha: 0.3", times=npred),
                      prop.table(table(factor(rankslong$rank[id_alpha03],levels=1:npred))))
alpha05 <- data.frame(alpha=rep("alpha: 0.5", times=npred),
                      prop.table(table(factor(rankslong$rank[id_alpha05],levels=1:npred))))
alpha07 <- data.frame(alpha=rep("alpha: 0.7", times=npred),
                      prop.table(table(factor(rankslong$rank[id_alpha07],levels=1:npred))))
rankdata <- rbind(alpha03,alpha05,alpha07)
#rankdata <- rbind(alpha005,alpha03,alpha05,alpha07)
rankdata$SNP <- "Irrelevant SNPs"
#id_alpha005 <- which(rankslong$problem=="simsnp_power" & rankslong$predictor %in% c("x1","x2","x3","x4","x5")
#                     & rankslong$algorithm=="rf_maxstat" & rankslong$betax2==betax2 & rankslong$alpha==0.05)
id_alpha03 <- which(rankslong$problem=="simsnp_power" & rankslong$predictor %in% c("x1","x2","x3","x4","x5")
                    & rankslong$algorithm=="rf_maxstat" & rankslong$betax2==betax2 & rankslong$alpha==0.3)
id_alpha05 <- which(rankslong$problem=="simsnp_power" & rankslong$predictor %in% c("x1","x2","x3","x4","x5")
                    & rankslong$algorithm=="rf_maxstat" & rankslong$betax2==betax2 & rankslong$alpha==0.5)
id_alpha07 <- which(rankslong$problem=="simsnp_power" & rankslong$predictor %in% c("x1","x2","x3","x4","x5")
                    & rankslong$algorithm=="rf_maxstat" & rankslong$betax2==betax2 & rankslong$alpha==0.7)
#alpha005 <- data.frame(alpha=rep("alpha: 0.05", times=npred),
#                       prop.table(table(factor(rankslong$rank[id_alpha005],levels=1:npred))))
alpha03 <- data.frame(alpha=rep("alpha: 0.3", times=npred),
                      prop.table(table(factor(rankslong$rank[id_alpha03],levels=1:npred))))
alpha05 <- data.frame(alpha=rep("alpha: 0.5", times=npred),
                      prop.table(table(factor(rankslong$rank[id_alpha05],levels=1:npred))))
alpha07 <- data.frame(alpha=rep("alpha: 0.7", times=npred),
                      prop.table(table(factor(rankslong$rank[id_alpha07],levels=1:npred))))
temp <- rbind(alpha03,alpha05,alpha07)
#temp <- rbind(alpha005,alpha03,alpha05,alpha07)
temp$SNP <- "Confounded SNPs"
rankdata <- rbind(rankdata,temp)
#id_alpha005 <- which(rankslong$problem=="simsnp_power" & rankslong$predictor %in% c("x6","x7","x8","x9","x10")
#                     & rankslong$algorithm=="rf_maxstat" & rankslong$betax2==betax2 & rankslong$alpha==0.05)
id_alpha03 <- which(rankslong$problem=="simsnp_power" & rankslong$predictor %in% c("x6","x7","x8","x9","x10")
                    & rankslong$algorithm=="rf_maxstat" & rankslong$betax2==betax2 & rankslong$alpha==0.3)
id_alpha05 <- which(rankslong$problem=="simsnp_power" & rankslong$predictor %in% c("x6","x7","x8","x9","x10")
                    & rankslong$algorithm=="rf_maxstat" & rankslong$betax2==betax2 & rankslong$alpha==0.5)
id_alpha07 <- which(rankslong$problem=="simsnp_power" & rankslong$predictor %in% c("x6","x7","x8","x9","x10")
                    & rankslong$algorithm=="rf_maxstat" & rankslong$betax2==betax2 & rankslong$alpha==0.7)
#alpha005 <- data.frame(alpha=rep("alpha: 0.05", times=npred),
#                       prop.table(table(factor(rankslong$rank[id_alpha005],levels=1:npred))))
alpha03 <- data.frame(alpha=rep("alpha: 0.3", times=npred),
                      prop.table(table(factor(rankslong$rank[id_alpha03],levels=1:npred))))
alpha05 <- data.frame(alpha=rep("alpha: 0.5", times=npred),
                      prop.table(table(factor(rankslong$rank[id_alpha05],levels=1:npred))))
alpha07 <- data.frame(alpha=rep("alpha: 0.7", times=npred),
                      prop.table(table(factor(rankslong$rank[id_alpha07],levels=1:npred))))
#temp <- rbind(alpha005,alpha03,alpha05,alpha07)
temp <- rbind(alpha03,alpha05,alpha07)
temp$SNP <- "Causal SNPs"
rankdata <- rbind(rankdata,temp)
colnames(rankdata) <- c("alpha","rank","proportion","SNP")
rankdata$rank <- as.numeric(rankdata$rank)

ggplot(rankdata, aes(x=rank, y=proportion)) + 
  geom_bar(stat="identity", position=position_dodge(), width=0.5) + 
  facet_grid(cols=vars(alpha), rows=vars(SNP), scales="free") +
  theme_bw() + 
  labs(x="Rank", y="Proportion")
ggsave("img/Figure27_power_sensitivity_alpha.png", height=6, width=7)


#*28) MSE Power case ----
id <- which(respred$problem=="pred_simsnp_power")
ggplot(respred[id,], aes(x=rf, y=error)) + 
  stat_boxplot(geom = "errorbar", width = 0.3) +
  geom_boxplot(fill="lightgrey", notch=FALSE) + 
  theme_bw() + 
  labs(x="Algorithm", y="MSE")
ggsave("img/Figure28_power_mse.png", height=3, width=7)