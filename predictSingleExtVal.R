#####################################################################
############################# Libraries #############################
#####################################################################

require(stringr)
require(dplyr)
require(IlluminaHumanMethylation450kanno.ilmn12.hg19)
require(pROC)
require(doParallel)
require(e1071)
require(glmnet)
require(PRROC)
require(ROCR)
require(randomForest)
require(argparse)
require(caret)

source('utils.R')
set.seed(123)

## set up parser ##
parser <- ArgumentParser()
parser$add_argument("-o","--outdir", action="store")
parser$add_argument("-i","--id", action="store")

outdir <- args$outdir

#####################################################################
############################ Main Script ############################
#####################################################################

## Read in data ## 
# data_test <- readRDS(paste0(outdir,'rds/',id,'TestSet.rds'))
# This can be replaced later in the pipeline so that this part starts once removing confounders is done
data_test <- readRDS('example_data.rda')

## Read in features ## 
features <- read.csv('features.txt',sep='\t')
features$probe <- as.character(features$probe)
features$gene <- as.character(features$gene)
## Aggregate probes  ## 
data_test <- aggregate_probes(data_test,features)
## Scale data ## 
data_test <- scale_df(data_test,features$gene)

## Predict on new data ## 
xgboost_results <- pred_cancer_xgboost_test(data_test,features$gene)

## Calibrate results ## 
calibrated_results <- platt_scaling(xgboost_results)
## Generate prediction metrics ## 
ROCobj_test <- ROCInfo_atcutoff(calibrated_results,other_title = 'Predictions')

## Save results object ## 
saveRDS(ROCobj_test,paste0(outdir,'rds/',id,'_AgeOfOnsetResults.rds'))
