fit_and_predict <- function(rsplit, wflow, parameters) {
  model <- finalize_workflow(wflow, parameters)
  model <- fit(model, training(rsplit))
  
  preprocessor <- extract_preprocessor(model)
  outcome <- preprocessor$term_info |> 
    filter(role == "outcome") |> 
    pull(variable)
  
  holdout_pred <-
    predict(model, testing(rsplit)) %>%
    bind_cols(testing(rsplit) |> select(!!outcome))
  
  holdout_pred
}

fit_and_score_inner <-
  function(resamples, wflow, grid, metrics, metric, control) {
    res <- tune_grid(
      wflow,
      resamples = resamples,
      grid = grid,
      metrics = metrics,
      control = control
    )
    select_best(res, metric = metric)
  }

#' Perform nested cross-validation
#'
#' @param wflow a workflow object with a preprocessor and model.
#' @param resamples a nested_cv resampling object created using the
#'   `rsample::nested_cv` function.
#' @param grid Hyperparameter tuning grid. Default is 25L. If an integer is
#'   supplied, then an automatic tuning grid is generated for the workflow based
#'   on `tune::latin_hypercube`.
#' @param metrics a metric_set object
#' @param metric character specifying the name of the metric in `metrics` to use
#'   for model selection as part of the hyperparameter tuning.
#' @param control fitting control object. Default is `tune::control_grid`.
#'
#' @return a nested_cv object with the additional columns 'parameters' to store
#'   the best hyperparameters for each inner resample, 'predictions' to store
#'   the out-of-fold outer predictions, and 'metrics' to store the outer
#'   resampling scoring results.
#' @export
cross_validate <-
  function(wflow,
           resamples,
           grid = 25L,
           metrics = metric_set(rmse, mae),
           metric = "rmse",
           control = control_grid()) {
    # fit inner resamples
    resamples$parameters <-
      map(
        resamples$inner_resamples,
        fit_and_score_inner,
        wflow = wflow,
        grid = grid,
        metrics = metrics,
        metric = metric,
        control = control
      )
    
    # fit outer resamples
    resamples$predictions <-
      map2(resamples$splits,
           resamples$parameters,
           fit_and_predict,
           wflow = wflow)
    
    # score outer resamples
    resamples$metrics <-
      map(
        resamples$predictions,
        ~ metrics(
          data = .x,
          truth = bedrock_dep,
          estimate = .pred
        )
      )
    
    return(resamples)
  }
