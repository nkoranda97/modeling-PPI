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

cl <- makeCluster(num_cores) #identify cpu cores
registerDoParallel(cl) #tell r to run in parallel

#xgb model, testing range or roungs and depth
set.seed(99)
ctrl = trainControl(method = "cv", number = 5)
fit_xgb = train(interact ~ .,
           data = PPI_T,
           method = "xgbTree",
                     tuneGrid = expand.grid(nrounds = c(80,90,100),
                                            max_depth = 8:10,
                                            eta = 0.3,
                                            gamma = 0.3,
                                            colsample_bytree = 0.9,
                                            min_child_weight = 1,
                                            subsample = 1),
           preProc = c("center", "scale"),
                    trControl = ctrl)

#ann model, testing number hidden nodes and penalization
set.seed(99)

ctrl = trainControl(method = "cv", "number" = 5)

fit_ann = train(interact ~ .,
                    data = PPI_T,
                    method = "nnet",
                    tuneGrid = expand.grid(size = 1:2,
                                           decay = c(0.5, 0.55, 0.6)),
                    trace = FALSE,
                    maxit = 5000,
                    preProc = c("center", "scale"),
                    trControl = ctrl)

