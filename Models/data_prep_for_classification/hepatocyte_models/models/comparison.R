library(tidymodels)
library(dplyr)
library(pROC)
library(readr)


# load model results
rf2_results <- readRDS('data_prep_for_classification/hepatocyte_models/models/models_rf/rf_model_2')
rf3_results <- readRDS('data_prep_for_classification/hepatocyte_models/models/models_rf/rf_model_3')

rf2_results$model <- "rf 2 (basic)"
rf3_results$model <- "rf 3(tuned)"
rf3_results

all_metrics <- bind_rows(
  rf2_results %>% select(.metric, .estimate, model),
  rf3_results %>% select(.metric, .estimate, model))

all_metrics_arranged <- all_metrics %>%
  arrange(.metric, desc(.estimate))


# bar plot

comp_plot <- ggplot(all_metrics_arranged, aes(x = model, y = .estimate, fill = model)) +
  geom_col(position = "dodge", width = 0.7) +
  facet_wrap(~ .metric) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = "Model",
    y = "Metric value",
    title = "Model comparison across metrics"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
ggsave("data_prep_for_classification/hepatocyte_models/models/models_rf/comparison_1.png", comp_plot, width = 8, height = 8)


# ROC curve
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

confusiton_matrix <- conf_mat(results, truth = var_non_amb, estimate = .pred_class)
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
