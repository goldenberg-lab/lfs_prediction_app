# LFS Prediction App #

Li-Fraumeni syndrome (LFS) is an autosomal dominant cancer predisposition syndrome. Approximately 80% of LFS patients harbour a germline TP53 mutation rendering them susceptible to a wide spectrum of early onset malignancies. A comprehensive surveillance regimen termed the ‘Toronto Protocol’, has recently been adopted for early tumor detection, demonstrating significant improvement in survival among these patients. However, the protocol’s “one-size-fits-all” approach fails to consider an individual patient's risk of cancer. We have built a machine learning model that predicts early onset of primary tumors in LFS patients by estimating the probability of cancer onset before the age of six, leveraging a patient's peripheral blood leukocyte methylation profile. 

## LFS Age of Onset Pipeline ##

This directory contains the following scripts which are used to preprocess methylation data and subsequently use that methylation data to predict early tumor onset in LFS patients. 

~~~
    utils.R
~~~

This script contains all the helper functions to perform preprocess, predict and plot. 

~~~
    remove_confounders.R
~~~ 
This script runs the preprocessing which involves three steps:

1. Outlier removal
2. Removal of batch confounder
3. Removal of array confounder

~~~
    predictSingleExtVal.R
~~~

This script runs the feature selection, model fitting and outputs the test results.




