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

PPI = read_csv("protein_analysis_results.csv", show_col_types = FALSE) %>%
  mutate(pred_hel_diff = pred_hel_1 - pred_hel_2) #forgot to add this in python script

PPI_T = PPI %>% 
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

#swap all _1 with _2
PPI_T_reversed = PPI_T %>%
  rename_with(~ gsub("_1$", "_temp", .), ends_with("_1")) %>% # switch all _1 variables with _2
  rename_with(~ gsub("_2$", "_1", .), ends_with("_2")) %>%
  rename_with(~ gsub("_temp$", "_2", .), ends_with("_temp")) 

#create predictions
actual = PPI_T$interact
predictions = predict(fit_xgb, newdata = PPI_T, type = "prob")[ ,2]
roc_1 = roc(actual, predictions)
predictions_reversed = predict(fit_xgb, newdata = PPI_T_reversed, type = "prob")[ ,2]
roc_rev = roc(actual, predictions_reversed)

#put ROC data into graph
roc_1 = data.frame(
  Sensetivity = roc_1$sensitivities,
  Specificity = roc_1$specificities,
  Data = "Normal Data"
)

roc_rev = data.frame(
  Sensetivity = roc_rev$sensitivities,
  Specificity = roc_rev$specificities,
  Data = "Switched Order"
)

bind_rows(roc_1, roc_rev) %>% 
  gf_line(Sensetivity ~ Specificity, color = ~ Data, size = 1) %>% 
  gf_abline(linetype = "dashed", slope = 1, intercept = 1, color = "grey") %>%
  gf_theme(theme_minimal()) %>%
  gf_refine(scale_x_reverse()) %>% 
  gf_labs(title = "Model is Inaccurate When Switching Order of Proteins")

#add origninal and reversed rows
PPI_T_2 = bind_rows(PPI_T, PPI_T_reversed)

#repeat xgb on full data
set.seed(99)
ctrl = trainControl(method = "cv", number = 5)
fit_xgb_2 = train(interact ~ .,
           data = PPI_T_2,
           method = "xgbTree",
                     tuneGrid = expand.grid(nrounds = c(80, 90, 100),
                                            max_depth = 8:10,
                                            eta = 0.3,
                                            gamma = 0.3,
                                            colsample_bytree = 0.9,
                                            min_child_weight = 1,
                                            subsample = 1),
           preProc = c("center", "scale"),
                    trControl = ctrl)

