# From https://www.rpubs.com/kaz_yos/alasso
adaptive_lasso <- function(x_cont, y_cont, gamma = 1, type.measure = "mse", nfolds = 10, intercept){
  ridge1_cv <- cv.glmnet(x = x_cont, y = y_cont, type.measure = type.measure, nfolds = nfolds, alpha = 0)
  best_ridge_coef <- as.numeric(coef(ridge1_cv, s = ridge1_cv$lambda.min))[-1]
  
  cv.glmnet(x = x_cont, y = y_cont, type.measure = type.measure, nfolds = nfolds, alpha = 1,
            penalty.factor = 1 / (abs(best_ridge_coef)^gamma), intercept = intercept)
}

lasso <- function(x_cont, y_cont, gamma = 1, type.measure = "mse", nfolds = 10, intercept){
  cv.glmnet(x = x_cont, y = y_cont, type.measure = type.measure, nfolds = nfolds, alpha = 1, intercept = intercept)
}


regression_wrapper <- function(my_ev, regression_method, column, shuffle, mat_sf_expression, quants, intercept){
  
  # get y data
  y <- quants[quants$event_id == my_ev, c("sample_id", column)] |>
    column_to_rownames("sample_id") |>
    filter(!is.na(.data[[column]])) |>
    as.matrix()
  
  # get x data
  x <- mat_sf_expression[rownames(y),]
  
  if(shuffle){
    x <- x[sample(nrow(x)),]
  }
  
  n <- nrow(x)
  train <- sample(n, round(.7*n))
  
  
  
  fit <- get(regression_method)(x[train,], y[train], nfolds = 20, intercept = intercept)
  
  # Estimate on test data
  prediction_on_test <- predict(fit, newx = x[-train,], s = "lambda.1se") |>
    as.data.frame() |>
    as_tibble(rownames = "sample_id") |>
    rename(predicted = lambda.1se) |>
    add_column(measured = y[-train])
  
  rsquare <- summary(lm(predicted ~ measured, data = prediction_on_test))$adj.r.squared
  
  coefs_sf <- coef(fit, s = "lambda.1se") |>
    as.matrix() |>
    as_tibble(rownames = "transcript_id")
  
  
  
  list(rsquare = rsquare, nb_coefs = sum(coefs_sf$s1 != 0),
       prediction_on_test = prediction_on_test,
       coefs_sf = coefs_sf, fit = fit)
}


do_regression <- function(my_ev, regression_method, column, shuffle = FALSE, mat_sf_expression, quants, intercept){
  possibly(regression_wrapper,
           otherwise = list(NA))(my_ev, regression_method, column, shuffle, mat_sf_expression, quants, intercept)
}
