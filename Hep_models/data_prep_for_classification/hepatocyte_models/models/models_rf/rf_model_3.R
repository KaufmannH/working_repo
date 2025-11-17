
# RF moodel 3: tidymodels RF model with cross validation and tuning 

# split
set.seed(123)

df_matched <- read_tsv('data_prep_for_classification/hepatocyte_models/results/matched_df.tsv')
head(df_matched)
# select cols i need
data <- df_matched %>% 
  select(var_non_amb, TATA_box, Inr,CCAAT_box, GC_box) %>%  
  mutate(var_non_amb = factor(var_non_amb)) %>%
  drop_na()

data_split <- initial_split(data, prop = 0.7, strata = var_non_amb)
train_data <- training(data_split)
test_data  <- testing(data_split)


cv_folds <- vfold_cv(train_data, v = 5, strata = var_non_amb) 

# specify model: RF
tree_model <- decision_tree(
  cost_complexity = tune(),
  tree_depth = tune(),
  min_n = tune()) %>%
  set_engine("rpart") %>%
  set_mode("classification")

# define workflow
tree_wf <- workflow() %>%
  add_model(tree_model) %>%
  add_formula(var_non_amb ~ TATA_box + GC_box + Inr + CCAAT_box)

# define tuning grid
  tree_grid <- grid_regular(
              cost_complexity(),
              tree_depth(),
              min_n(),
              levels = 5 ) # 5x5x5 = 125 Kombinationen)

# tune model
tuned_results <- tune_grid(
                tree_wf,
                resamples = cv_folds,
                grid = tree_grid,
                metrics = metric_set(accuracy, roc_auc),
                control = control_grid(save_pred = TRUE))

# evaluate
best_params <- select_best(tuned_results, metric = "roc_auc")
final_wf <- finalize_workflow(tree_wf, best_params)

# retrain 
final_fit <- final_wf %>%
  last_fit(split = data_split)

tuned_metrics <- collect_metrics(final_fit)

final_predictions <- collect_predictions(final_fit)
metrics_selected <- metric_set(accuracy, sens, yardstick::spec)
tuned_metrics <- metrics_selected(final_predictions, truth = var_non_amb, estimate = .pred_class)
tuned_metrics

saveRDS(tuned_metrics, 'data_prep_for_classification/hepatocyte_models/models/models_rf/rf_model_3')

