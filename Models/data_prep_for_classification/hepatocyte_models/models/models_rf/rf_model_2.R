library(tidymodels)
library(dplyr)
library(pROC)
library(readr)


# basic tidymodel rf model (rf model 2)

df_matched <- read_tsv('data_prep_for_classification/hepatocyte_models/results/matched_df.tsv')
head(df_matched)
# select cols i need
data <- df_matched %>% 
  select(var_non_amb, TATA_box, Inr,CCAAT_box, GC_box) %>%  
  mutate(var_non_amb = factor(var_non_amb)) %>%
  drop_na()


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
	  fit(var_non_amb ~ TATA_box + GC_box + Inr + CCAAT_box,
	      data = train_data)

prediction_class <-  predict(rf_fit, new_data = test_data, type = 'class') 
prediction_prob <-  predict(rf_fit, new_data = test_data, type = 'prob') 

results <- test_data |>
    select(var_non_amb) |>
    bind_cols(prediction_class, prediction_prob)
results

# speed variant
last_fit <- last_fit(rf_model, var_non_amb ~ TATA_box + GC_box + Inr + CCAAT_box, split = data_split)
last_fit_results <- collect_predictions(last_fit)
last_fit_results

# evaluation
metrics_selected <- metric_set(accuracy, sens, yardstick::spec)
basic_metrics <- metrics_selected(results, truth = var_non_amb, estimate = .pred_class)
roc_result <- roc_auc(results, truth = var_non_amb, .pred_HVG)
metrics_combined <- bind_rows( basic_metrics,roc_result)

saveRDS(metrics_combined, 'data_prep_for_classification/hepatocyte_models/models/models_rf/rf_model_2')




