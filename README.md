---
title: "Case Based Reasoning with Confidence"
author: "Christopher Bartlett"
date: "11/24/2020"
output: 
  html_document:
  theme: journal
  toc: true
  toc_float:
    collapsed: false
    smooth_scroll: true
  preserve_yaml: true
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Case Based Reasoning with Confidence

Case Based Reasoning with Confidence (CBR-CONF) is an R package for survival prediction that uses a novel confidence weight. Confidence is 1/2 the Euclidean distance between training samples to a prototypical case of its own risk group (low, medium or high risk of having the event) plus the average distance of the training sample to all testing samples. This weight is normalized between 0 and 1. During retrieval, training samples are retrieved for a test sample until confidence is greater than 1.0. The predicted survival time for a test sample is then the mean survival time of the retrieved cases. 

The prognostic score for determining risk is computed by performing a multivariate cox regression and multiplying the beta coefficients by the feature value (such as a gene expression value) to determine the value of each feature and then summing these values to determine the prognostic score of a sample. Risk groups are then computed by taking a tertile ranking. 

During feature selection, the dataset is split into n training and testing sets and a multivariate cox regression is performed on each training set. Significant features with a log rank < 0.01 recurring in a specified percentage of the training sets are kept at the end of the runs. The prognostic score is then calculated from the retained features. 

### Requirements
Only a dataset is required that has "time" as the first column and "status" as the second column. These column identifiers are case sensitive. "time" is a numeric variable that depicts the time to the event, and "status" is a binary value pertaining to whether or not the event occurred. The rest of the columns are features, and each row is a sample. 

A sample dataset is provided.

### Loading CBR-CONF
```{r settingup, echo=FALSE}
suppressMessages(require(caret))
suppressMessages(require(RegParallel))
suppressMessages(require(survival))
suppressMessages(require(survcomp))
suppressMessages(require(stats))
suppressMessages(require(textshape))
suppressMessages(require(gbm))
suppressMessages(require(rms))
suppressMessages(require(ranger))
suppressMessages(require(randomForestSRC))
```

```{r loading, echo = FALSE}
library(CBRCONF)
```

### Included data
A sample dataset of a subset of TCGA-BRCA gene expression data has been provided. This dataset, called GeneSub, has 500 random cancer tissue samples and 500 random gene expression features. The first column are overall survival times, and the second column is the survival status. 1 means that the sample was deceased and 0 means that the sample was living at the time of the last follow up. The features were found to be differentially expressed between the cancer and normal tissue samples and are after batch effect analysis but no other preprocessing has been performed. 

```{r retrieve}
data(GeneSub)
```

### Dimensions of the included dataset
```{r preddim, echo = FALSE}
dim(GeneSub)
```
### CBRCONF object
The CBRCONF object is simply an R6 class that stores objects that will be used throughout the various functions, as well as results. R6 classes have accessors to access their information, used with \$. An example is CBRCONF\$testing_results.

The accessors are:

* dataset
  + The dataset where time and status are the first and second columns and subsequent columns are features
* sig_feature_list
  + List of all significant features found in each iteration of feature selection
* sig_feature_data
  + Significant features that are retained after all iterations of feature selection
* feature_scores
  + A dataframe to contain all of the individual feature scores. This is the feature's expression value multiplied by the mean beta coefficient
* prognostic_scores
  + Dataframe of prognostic scores for all samples
* results_list
  + Holds the results for feature selection
* risk_groups
  + Holds the risk groups, which is assigned by the prognostic risk scores
* group_means
  + The features of each sample within each risk group is averaged to create one prototypical representation of each group. This stores that information
* low_risk_features
  + A subset of the original dataset that is just the low risk group, and their corresponding features
* medium_risk_features
  + A subset of the original dataset that is just the medium risk group, and their corresponding feature
* high_risk_features
  + A subset of the original dataset that is just the high risk group, and their corresponding feature
* low_risk_prototype
  + A vector of the averaged features of the low risk group
* medium_risk_prototype
  + A vector of the averaged features of the medium risk group
* high_risk_prototype
  + A vector of the averaged features of the high risk group
* confidence_weights_low
  + The confidence weights for just the low risk group
* confidence_weights_medium
  + The confidence weights for just the medium risk group
* confidence_weights_high
  + The confidence weights for just the high risk group
* confidence_weights
  + A dataframe that combines the confidence weights of the three risk groups
* data.folds
  + Each of the training and testing folds
* training_folds
  + Each of the training folds
* testing_folds
  + Each of the testing folds
* case_retrieval_information
  + Which cases were retrieved by which case
* testing_number_retrieved
  + How many cases each testing sample retrieved
* testing_concordance
  + The concordance values for each testing fold
* testing_sample_information
  + The risk group, predicted survival times, status, and testing_results




### Creating a CBRCONF object
To get started, you'll need to create a CBRCONF object to hold all the data. It can be named anything. 


```{r create}
CBRCONF <- create_CBRCONF(GeneSub)
```

### Cox Module
This module performs the feature selection step, in which it runs for Num iterations. At each iteration, a training set
and testing set will be created and every feature will run through a multivariate cox regression. After Num runs, features found to be significant in the percent of training sets specified with Percent are kept (log ranking < 0.01 is used for significance). Cores and blocksize are used for a parallel process (uses RegParallel). Cores is the number of computing cores assigned to the task, and blocksize is the number of features to run at once. Blocksize should not exceed the number of features in your dataset. 

```{r cox}
CBRCONF <- cox_module(CBRCONF, Num = 10, Percent = 0.5, Cores = 1)
```


### Risk Group Module
This module constructs 3 prognostic risk groups, taken from the ile function. In other words, each risk group is derived from a tertile of the prognostic score. A KM plot can be constructed at this point from the actual survival data.

```{r risk groups}
CBRCONF <- risk_groups_module(CBRCONF, KM.Plot = FALSE)
```

### Prototype Module
This module creates 3 prototypical representations of the risk groups. This is performed by using all samples within a risk
group and taking an average of each feature. Specifically, a vector is constructed for low, medium and high risk that has
averaged values for each feature.

```{r prototype}
CBRCONF <- prototype_module(CBRCONF)
```

### Prototype Distance Moduke
This module builds the distance matrices that determine how similar the samples are. Euclidean distance is always used.
Confidence weights will start to be built in this module (at this stage, only the distance from the prototype is calculated)

```{r prototype distance}
CBRCONF <- prototype_distance_module(CBRCONF)
```

### CBR Module
This function performs n fold cross validation using the CBR-CONF framework. Confidence weights are used to determine
when retrieval ends for each test sample. Upon conclusion, the concordance index will be displayed.

```{r CBR}
CBRCONF <- CBR_module(CBRCONF, n.folds = 10, KM.Plot = FALSE)
```
### CBRCONF All in One
You can also perform all of these options in one comprehensive function. Please note that you won't have access to all of the results that you would have if you ran the modules. 

```{r CBR All in One}
AllinOne <- CBRCONF_AllinOne(dataset = GeneSub, Num = 10, Percent = 0.5, Cores = 1, blocksize = 200, n.folds = 10, KM.Plot = FALSE)
```


