


```{r}
library(survival)
set.seed(123)
mgus_df <- complete(mice(mgus, m=1, maxit = 50, method = 'pmm', seed = 500),1)
mgus_split <- mgus_df %>% 
  select(-futime,-death,-id,-ptime, -alb) %>% 
  initial_split(., strata = pstat)
mgus_train <- training(mgus_split)
mgus_test <- testing(mgus_split)
```

```{r}
xgb_spec <- boost_tree(
  trees = 1000,
  tree_depth = tune(), min_n = tune(), loss_reduction = tune(),
  sample_size = tune(), mtry = tune(), learn_rate = tune()
) %>%
  set_engine("xgboost",scale_pos_weight = 1000000) %>%
  set_mode('regression')

xgb_spec
```

```{r}
xgb_grid <- grid_latin_hypercube(
  tree_depth(),
  min_n(),
  loss_reduction(),
  sample_size = sample_prop(),
  finalize(mtry(), train),
  learn_rate(),
  size = 30
)

xgb_grid
```

```{r}
xgb_wf <- workflow() %>%
  add_formula(pstat ~ .) %>%
  add_model(xgb_spec)
```

```{r}
set.seed(123)
mgus_folds <- vfold_cv(mgus_train, strata = pstat)
mgus_folds
```

```{r}
doParallel::registerDoParallel()

set.seed(234)
xgb_res <- tune_grid(
  xgb_wf,
  resamples = mgus_folds,
  grid = xgb_grid, 
  control = control_grid(save_pred = TRUE)
)
```

```{r}
xgb_res %>%
  collect_metrics() %>%
  filter(.metric == "rmse") %>%
  select(mean, mtry:sample_size) %>%
  pivot_longer(mtry:sample_size,
               values_to = "value",
               names_to = "parameter"
  ) %>%
  ggplot(aes(value, mean, color = parameter)) +
  geom_point(alpha = 0.8, show.legend = FALSE) +
  facet_wrap(~parameter, scales = "free_x") +
  labs(x = NULL, y = "RMSE")

show_best(xgb_res, "rmse")

best_rmse <- select_best(xgb_res, "rmse")
best_rmse

final_xgb <- finalize_workflow(
  xgb_wf,
  best_rmse
)

final_xgb
```

```{r}
library(vip)
final_xgb %>%
  fit(data = mgus_train) %>%
  pull_workflow_fit() %>%
  vip(geom = "point")

final_res <- last_fit(final_xgb, mgus_split)

collect_metrics(final_res)

final_res %>%
  collect_predictions() %>%
  mutate(pstat = as.factor(pstat))      %>%
  roc_curve(pstat, .pred) %>%
  ggplot(aes(x = 1 - specificity, y = sensitivity)) +
  geom_line(size = 1.5, color = "midnightblue") +
  geom_abline(
    lty = 2, alpha = 0.5,
    color = "gray50",
    size = 1.2
  )

#Plot metrics
library(runway)
final_res %>%
  collect_predictions() %>% 
  threshperf_plot(.,
                  outcome = 'pstat',
                  prediction = '.pred')

final_res %>%
  collect_predictions() %>% 
  cal_plot(.,
           outcome = 'pstat', 
           prediction = '.pred',
           show_loess = TRUE)

final_res %>%
  collect_predictions() %>% 
  roc_plot(., 
           outcome = 'pstat', 
           prediction = '.pred',
           ci = TRUE, 
           plot_title = 'Single ROC curve w/CI ribbon')
```