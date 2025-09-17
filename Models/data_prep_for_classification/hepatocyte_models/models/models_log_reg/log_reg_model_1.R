library(tidymodels)
library(dplyr)
library(pROC)
library(readr)


df_matched <- read_tsv('data_prep_for_classification/hepatocyte_models/results/matched_df.tsv')
head(df_matched)
# select cols i need
data <- df_matched %>% 
  select(var_non_amb, TATA_box, Inr,CCAAT_box, GC_box) %>%  
  mutate(var_non_amb = factor(var_non_amb)) %>%
  drop_na()

# I. basic model

# split
set.seed(123)
data_split <- initial_split(data, prop = 0.7, strata = var_non_amb)
train_data <- training(data_split)
test_data  <- testing(data_split)

# specify model: RF
log_model <- logistic_reg() %>% 
    set_engine('glm') %>% 
    set_mode('classification')

# normal fit
log_fit <- log_model %>% 
	  fit(var_non_amb ~ TATA_box + GC_box + Inr + CCAAT_box,
	      data = train_data)

prediction_class <-  predict(log_fit, new_data = test_data, type = 'class') 
prediction_prob <-  predict(log_fit, new_data = test_data, type = 'prob') 

results <- test_data |>
    select(var_non_amb) |>
    bind_cols(prediction_class, prediction_prob)


# speed variant
last_fit <- last_fit(log_model, var_non_amb ~ TATA_box + GC_box + Inr + CCAAT_box, split = data_split)
last_fit_results <- collect_predictions(last_fit)
last_fit_results

# evaluation
metrics_selected <- metric_set(sens, spec)
metrics <- metrics_selected(results, truth = var_non_amb, estimate = .pred_class)
confusiton_matrix <- conf_mat(results, truth = var_non_amb, estimate = .pred_class)


# II. Model with cross validation and tuning


# split
set.seed(123)
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

collect_metrics(final_fit)







# plot
#conf_plot <- autoplot(confusiton_matrix, type = 'mosaic') # mosaic type too

roc_data <- roc_curve(results, truth = var_non_amb, .pred_LVG)
roc_auc_val <- roc_auc(results, truth = var_non_amb, .pred_LVG)
roc_auc_val

roc_plot <- ggplot(roc_data, aes(x = 1 - specificity, y = sensitivity)) +
  geom_line(color = "lightblue", size = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey") +
  coord_equal() +
  labs(
    title = sprintf("HVG ident. by promoter motifs: AUC = %.3f", roc_auc_val),
    x = "Specificity (TP)",
    y = "Sensitivity (FP)") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(size = 14, hjust = 0.5))

ggsave("data_prep_for_classification/hepatocyte_models/results/predictions/roc_log_reg_1.png", roc_plot, width = 8, height = 8)

# confusion matrix
cm_df <- as.data.frame(confusiton_matrix$table)

conf_plot <- ggplot(cm_df, aes(x = Prediction, y = Truth, fill = Freq)) +
  geom_tile(colour = "black") +
  geom_text(aes(label = Freq), size = 6) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(
    title = "Confusion Matrix",
    x = "Predicted Class",
    y = "True Class",
    fill = "Count"
  ) +
  theme_classic()
ggsave("data_prep_for_classification/hepatocyte_models/results/predictions/conf_log_reg_1.png", conf_plot, width = 8, height = 8)
