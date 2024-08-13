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

#repeat ann on full data
set.seed(99)

ctrl = trainControl(method = "cv", "number" = 5)

fit_ann_2 = train(interact ~ .,
                    data = PPI_T_2,
                    method = "nnet",
                    tuneGrid = expand.grid(size = 1:2,
                                           decay = c(0.5, 0.55, 0.6)),
                    trace = FALSE,
                    maxit = 5000,
                    preProc = c("center", "scale"),
                    trControl = ctrl)

#calculate new predictions
actual_full = PPI_T_2$interact
predictions_full = predict(fit_xgb_2, type = "prob")[ ,2]
predictions_2 = predict(fit_xgb_2, newdata = PPI_T, type = "prob")[ ,2]
predictions_reversed_2 = predict(fit_xgb_2, newdata = PPI_T_reversed, type = "prob")[ ,2]

#roc data for new graphs
roc_full = roc(actual_full, predictions_full)
roc_2 = roc(actual, predictions_2)
roc_rev_2 = roc(actual, predictions_reversed_2)

roc_full = data.frame(
  Sensetivity = roc_full$sensitivities,
  Specificity = roc_full$specificities,
  Data = "Full Data"
)

roc_2 = data.frame(
  Sensetivity = roc_2$sensitivities,
  Specificity = roc_2$specificities,
  Data = "Normal Data"
)

roc_rev_2 = data.frame(
  Sensetivity = roc_rev_2$sensitivities,
  Specificity = roc_rev_2$specificities,
  Data = "Reversed Data"
)

bind_rows(roc_full, bind_rows(roc_2, roc_rev_2)) %>% 
  gf_line(Sensetivity ~ Specificity, color = ~ Data, size = 1) %>% 
  gf_abline(linetype = "dashed", slope = 1, intercept = 1, color = "grey") %>%
  gf_theme(theme_minimal()) %>%
  gf_refine(scale_x_reverse(),
            theme(legend.position = "none")) %>% 
  gf_labs(title = "Training Model With Original and Reversed Data Set Imporoves Robustness") %>% 
  gf_facet_wrap( ~ Data)

#extract varImp data from xgb model for nice graph
importance = varImp(fit_xgb_2)
importance_df = as.data.frame(importance$importance) #turn varImp into a dataframe
importance_df$Variable = rownames(importance_df)

importance_df %>% 
  filter(Overall > 10) %>% #filter for most important variables
  mutate(Variable = fct_rev(fct_inorder(Variable))) %>% #order by importance
  mutate(base_name = sapply(str_split(Variable, "_"), `[`, 1)) %>% #extract base name
  gf_col(Variable ~ Overall, fill = ~ base_name) %>% 
  gf_refine(guides(fill = "none"),
            theme = theme_minimal()) %>% 
  gf_labs(x = "Importance",
          title = "Variable Importance")

#set every predictor variable at median except pred_hel_diff, predict response and then graph
pred_hel_diff = seq(0, 15, by = 0.05)

median_df = PPI_T_2 %>% 
  select(-c(interact, pred_hel_diff)) %>% 
  summarise_all(median) %>% 
  dplyr::slice(rep(1,length(pred_hel_diff))) %>% 
  mutate(pred_hel_diff = pred_hel_diff)

hel_pred = predict(fit_xgb_2, newdata = median_df, type = "prob")[ ,2]

gf_line(hel_pred ~ pred_hel_diff, size = 1) %>% 
  gf_theme(theme_minimal()) %>% 
  gf_labs(y = "Predicted Probability of Interaction",
          x = "Magnitude of Predicted Helicies",
          title = "Magnitude of Predicted Helicies Effect on PPI")

#set every predictor variable at median except pairwise alignment score, predict response and then graph

pred_align = seq(80, 4000, by = 10)

median_df = PPI_T_2 %>% 
  select(-c(interact, pairwise_alignment)) %>% 
  summarise_all(median) %>% 
  dplyr::slice(rep(1,length(pred_align))) %>% 
  mutate(pairwise_alignment = pred_align)

hel_align= predict(fit_xgb_2, newdata = median_df, type = "prob")[ ,2]

gf_line(hel_align ~ pred_align, size = 1) %>% 
  gf_theme(theme_minimal()) %>% 
  gf_labs(y = "Predicted Probability of Interaction",
          x = "Pairwise Alignment score",
          title = "Pairwise Alignment Score Effect on PPI")

# double cv validation
set.seed(99)

dataused = PPI_T_2
n = nrow(dataused)
nfolds = 5

cv_pred = factor(vector(length = n), levels(PPI_T$interact))

groups = rep(1:nfolds, length = n)
cv_groups = sample(groups, n)

best_fold_model = list(rep(NA, nfolds))

ctrl = trainControl(method = "cv", number = 5)

for(i in 1:nfolds){
  in_train = (cv_groups != i)
  in_test = (cv_groups == i)
  
  traindata = dataused[in_train, ]
  testdata = dataused[in_test, ]
  
  fit_xgb_inner= train(interact ~ .,
             data = traindata,
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
  
  fit_ann_inner = train(interact ~ .,
                    data = traindata,
                    method = "nnet",
                    tuneGrid = expand.grid(size = 1:2,
                                           decay = c(0.5, 0.55, 0.6)),
                    trace = FALSE,
                    maxit = 5000,
                    preProc = c("center", "scale"),
                    trControl = ctrl)
  
  accuaracies = list(max(fit_xgb_inner$results$Accuracy),
                     max(fit_ann_inner$results$Accuracy))
  model_list = list(fit_xgb_inner, fit_ann_inner)
  
  best_fold_model[[i]] = model_list[[which.max(accuaracies)]]
  
  cv_pred[in_test] = predict(best_fold_model[[i]], newdata = testdata)
  
  }

stopCluster(cl) #stop parallel computing

sum(diag(table(cv_pred, actual_full )))/nrow(PPI_T_2)
