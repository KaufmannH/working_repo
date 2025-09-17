library(tidymodels)
library(dplyr)
#library(pROC)
library(readr)


# model: basic tidymodel rf model (rf model 2)
# predicts: LVG/HVG across all ages (just based on gmean)
# input data: only gmean (QC)

df_matched <- readRDS('data_prep_for_classification/hepatocyte_models/data/matched_tf_binding_df.rds')  
head(df_matched)
# select cols i need
data <- df_matched %>% 
  select(var_non_amb, gmean) %>%  
  mutate(var_non_amb = factor(var_non_amb)) %>%
  drop_na()
data

# split
set.seed(123)
data_split <- initial_split(data, prop = 0.7, strata = var_non_amb)
train_data <- training(data_split)
test_data  <- testing(data_split)

# specify model: RF
rf_model <- decision_tree() %>% 
    set_engine('rpart') %>% 
    set_mode('classification')

# normal fit
rf_fit <- rf_model %>% 
          fit(var_non_amb ~ gmean,
          data = train_data)

prediction_class <-  predict(rf_fit, new_data = test_data, type = 'class') 
prediction_prob <-  predict(rf_fit, new_data = test_data, type = 'prob') 

results <- test_data |>
    select(var_non_amb) |>
    bind_cols(prediction_class, prediction_prob)
results


# speed variant
last_fit <- last_fit(rf_model, var_non_amb ~ gmean, split = data_split)
last_fit_results <- collect_predictions(last_fit)
last_fit_results

# evaluation
metrics_selected <- metric_set(accuracy, sens, yardstick::spec)
basic_metrics <- metrics_selected(results, truth = var_non_amb, estimate = .pred_class)
roc_result <- roc_auc(results, truth = var_non_amb, .pred_HVG)
metrics_combined <- bind_rows( basic_metrics,roc_result)

saveRDS(metrics_combined, 'data_prep_for_classification/hepatocyte_models/models/models_rf/rf_model_2_3')



# ROC curve
#conf_plot <- autoplot(confusiton_matrix, type = 'mosaic') # mosaic type too

roc_data <- roc_curve(results, truth = var_non_amb, .pred_HVG)
roc_auc_val <- roc_auc(results, truth = var_non_amb, .pred_HVG)$.estimate
roc_auc_val

roc_plot <- ggplot(roc_data, aes(x = 1 - specificity, y = sensitivity)) +
  geom_line(color = "#39a379", size = 1) +
  geom_ribbon(aes(ymin = 0, ymax = sensitivity), fill = "#39a379", alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey") +
  labs(
    title = sprintf("HVG ident. by gmean: AUC = %.3f", roc_auc_val),
    x = "Specificity (TP)",
    y = "Sensitivity (FP)") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(size = 14, hjust = 0.5))

ggsave("data_prep_for_classification/hepatocyte_models/results/predictions/roc_rf_M2_3.png", roc_plot, width = 5, height = 5)



# plot feature importance
importance <- rf_fit$fit$variable.importance
importance_df <- enframe(importance, name = "feature", value = "importance")


familial_profile_metadata <- readRDS("data_prep_for_classification/familial_binding_profile_predictions/jaspar/tf_names_df_250_50.rds")
importance_annotated_df <- familial_profile_metadata |>
  rename(feature = familial_profile) |>
  right_join(importance_df, by = 'feature')

feat_plot <- ggplot(importance_annotated_df, aes(x = reorder(feature, importance), y = importance)) +
  geom_col(fill = "#39a379") +
  coord_flip() +
  labs(
    title = "Feature Importance (Decision Tree)",
    x = "Feature",
    y = "Importance") +
  theme_classic()

ggsave("data_prep_for_classification/hepatocyte_models/results/predictions/feat_rf_M2_2.png", feat_plot, width = 5, height = 5)

important <- importance_annotated_df |>
  filter(feature == 'Familial_profile_122')
important  
