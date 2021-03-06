library(caret)
library(tidyverse)
source('helpers.R')

# get model and validation data for probability calibration below (otherwise valid not used)
mod <- readRDS('valid_model.rds')
valid <- mod[[2]]
mod <- mod[[1]]

# get features
feat_list <- mod[[11]]$xNames

# get our results 
test <- readRDS('NoobCorrected_after_covs_beta_ProjPC2Adj_lfs_5UTR_scaled_canceratdraw_systreat_xgboost_ageofonset_ROCInfoTest.rds')
test <- test[[1]]
test$test_label <- ifelse(test$test_label==1, 'Yes', 'No')

# calibrate probabilities
temp_calibrated <- calibrate_probs_valli(train_results = valid, test_results = test)

# validation and test data
valid <- temp_calibrated[[1]]
test <- temp_calibrated[[2]]

# remove 
rm(temp_calibrated, valid)


