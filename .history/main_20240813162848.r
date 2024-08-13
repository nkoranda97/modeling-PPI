#modeling
library(caret)
#all tidyr libraries
library(tidyverse)
#ggformula specific graphing
library(ggformula)
#ROC graphs
library(pROC)
#compute models in parallel by using all CPU cores
library(doParallel)
num_cores <- detectCores()

PPI <-read_csv("protein_analysis_results.csv", show_col_types = FALSE) %>%
  mutate(pred_hel_diff = pred_hel_1 - pred_hel_2) #forgot to add this in python script