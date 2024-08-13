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

PPI_T <- PPI %>% 
    mutate(
      interact = as.factor(interact),
      across(
        .cols = where(is.numeric) & !matches("(_1|_2)$"), #select all difference measurements
    .fns = ~ abs(.))) %>% #take absolute value, difference values created confusing interpretations
  select(
    !nearZeroVar(.)) %>% 
  drop_na() 