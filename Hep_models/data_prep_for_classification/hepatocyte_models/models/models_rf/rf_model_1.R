library(randomForest)
library(dplyr)
library(pROC)
library(readr)
library(ggplot2)

# RF model 1: first try of basic rf model (randomForest package)
df_main <- read_tsv('data_prep_for_classification/hepatocyte_models/results/hvg_lvg_predictors_counts.tsv') # og eszter data
df_main <- read_tsv('data_prep_for_classification/hepatocyte_models/results/matched_df.tsv') # new 16thjune


df <- df_main %>%
  filter(var_non_amb %in% c("HVG", "LVG")) %>%       
  mutate(class = factor(var_non_amb))   

# split: 70/30
set.seed(123) 
train_indices <- sample(nrow(df), 0.7 * nrow(df))
train_data <- df[train_indices, ]
test_data  <- df[-train_indices, ]

# fit rf model
rf_model <- randomForest(
  class ~ TATA_box + GC_box + Inr + CCAAT_box,
  data = train_data,
  importance = TRUE,
  ntree = 500)

# predict
prob <- predict(rf_model, newdata = test_data, type = "prob")[, "HVG"]

# plot
roc_obj <- roc(response = test_data$var_non_amb, predictor = prob, levels = c("LVG", "HVG"))
auc_val <- auc(roc_obj)

roc_df <- data.frame(
  specificity  = roc_obj$specificities,
  sensitivity  = roc_obj$sensitivities)

roc_plot <- ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity)) +
  geom_line(color = "lightblue", size = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey") +
  coord_equal() +
  labs(
    title = sprintf("HVG ident. by promoter motifs: AUC = %.3f", auc_val),
    x = "Specificity (TP)",
    y = "Sensitivity (FP)") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(size = 14, hjust = 0.5))

ggsave("data_prep_for_classification/hepatocyte_models/results/predictions/roc_rf_1.png", roc_plot, width = 8, height = 8)

