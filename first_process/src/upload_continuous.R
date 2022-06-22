# Upload evidence score results to a central place
library(mrbrt001, lib.loc = "/ihme/code/mscm/R/packages/")
globals <- new.env()
source("config.R", local = globals)
ARCHIVE <- globals$ARCHIVE


# get exposure names
get_exp_names <- function(
    signal_model
) {
    cov_model <- signal_model$sub_models[[1]]$cov_models[[1]]
    ref_exp_names <- cov_model$ref_cov
    alt_exp_names <- cov_model$alt_cov
    return(list(ref_exp_names, alt_exp_names))
}


# function get signal
get_signal <- function(
    signal_model,
    ref_val,
    alt_val
) {
    exp_names <- get_exp_names(signal_model)
    
    df <- data.frame(
        c(sapply(exp_names[[1]], function(x) ref_val, simplify = FALSE, USE.NAMES = TRUE),
          sapply(exp_names[[2]], function(x) alt_val, simplify = FALSE, USE.NAMES = TRUE))
    )

    data <- MRData()
    data$load_df(df, col_covs = unlist(exp_names))

    signal_model$predict(data)
}


get_exp <- function(signal_model) {
    df <- signal_model$data$to_df()
    exp_names <- get_exp_names(signal_model)
    ref_exp <- df[exp_names[[1]]]
    alt_exp <- df[exp_names[[2]]]

    list(ref_exp = ref_exp, alt_exp = alt_exp)
}


get_exp_mid <- function(signal_model) {
    exp <- get_exp(signal_model)
    list(ref_mid = rowMeans(exp[[1]]), alt_mid = rowMeans(exp[[2]]))
}


get_exp_limits <- function(signal_model) {
    exp <- get_exp(signal_model)
    c(min(unlist(exp)), max(unlist(exp)))
}


get_exp_bounds <- function(signal_model) {
    exp_mid <- get_exp_mid(signal_model)
    c(quantile(exp_mid[[1]], 0.15), quantile(exp_mid[[2]], 0.85))
}


# function that get data
get_data <- function(
    rei_id,
    cause_id,
    signal_model,
    linear_model,
    normalize_to_tmrel = FALSE
) {
    # extract data frame
    df <- signal_model$data$to_df()
    
    # get exposure names
    exp_names <- get_exp_names(signal_model)

    # get midpoints for data exposures
    exp_mid <- get_exp_mid(signal_model)
    df$ref_risk <- exp_mid[[1]]
    df$alt_risk <- exp_mid[[2]]

    # get outlier indices
    df$is_outlier <- as.integer(signal_model$get_w_soln() < 0.1)

    # get the exposure limit
    exp_limits <- get_exp_limits(signal_model)

    # get signal for reference and alternative midpoints
    df$log_ref_cause <- get_signal(signal_model,
                                   ref_val = exp_limits[1],
                                   alt_val = df$ref_risk)
    if (normalize_to_tmrel) {
        pred <- get_pred(signal_model)
        df$log_ref_cause <- df$log_ref_cause - min(pred[[2]])
    }
    df$log_alt_cause <- df$log_ref_cause + df$obs

    # get upper and lower bounds for the log rr
    df$log_alt_cause_lower <- df$log_alt_cause - 1.96*df$obs_se
    df$log_alt_cause_upper <- df$log_alt_cause + 1.96*df$obs_se

    # rename columns
    names(df)[names(df) == "obs"] <- "log_rr"
    names(df)[names(df) == "obs_se"] <- "log_rr_se"
    names(df)[names(df) == exp_names[[1]][1]] <- "ref_risk_lower"
    names(df)[names(df) == exp_names[[2]][1]] <- "alt_risk_lower"
    if (length(exp_names[[1]]) == 2) {
        names(df)[names(df) == exp_names[[1]][2]] <- "ref_risk_upper"
    } else {
        df$ref_risk_upper <- df$ref_risk_lower
    }
    if (length(exp_names[[2]]) == 2) {
        names(df)[names(df) == exp_names[[2]][2]] <- "alt_risk_upper"
    } else {
        df$alt_risk_upper <- df$alt_risk_lower
    }

    # adding new columns
    df$rei_id <- rei_id
    df$cause_id <- cause_id
    
    df$linear_ref_cause <- exp(df$log_ref_cause)
    df$linear_alt_cause_lower <- exp(df$log_alt_cause_lower)
    df$linear_alt_cause <- exp(df$log_alt_cause)
    df$linear_alt_cause_upper <- exp(df$log_alt_cause_upper)

    df[c(
        "rei_id",
        "cause_id",
        "study_id",
        "log_rr",
        "log_rr_se",
        "ref_risk_lower",
        "ref_risk",
        "ref_risk_upper",
        "alt_risk_lower",
        "alt_risk",
        "alt_risk_upper",
        "log_ref_cause",
        "linear_ref_cause",
        "log_alt_cause_lower",
        "log_alt_cause",
        "log_alt_cause_upper",
        "linear_alt_cause_lower",
        "linear_alt_cause",
        "linear_alt_cause_upper",
        "is_outlier"
    )]
}


get_beta <- function(linear_model) {
    beta_soln <- linear_model$beta_soln
    names(beta_soln) <- linear_model$cov_names
    beta <- beta_soln[["signal"]]

    model_specs <- mrbrt001::core$other_sampling$extract_simple_lme_specs(linear_model)
    beta_soln_hessian <- mrbrt001::core$other_sampling$extract_simple_lme_hessian(model_specs)
    beta_soln_sd <- sqrt(diag(solve(beta_soln_hessian)))
    names(beta_soln_sd) <- linear_model$cov_names
    beta_sd <- beta_soln_sd[["signal"]]
    return(c(beta, beta_sd))
}


get_gamma <- function(linear_model) {
    gamma_soln <- linear_model$gamma_soln
    gamma_soln_fisher <- linear_model$lt$get_gamma_fisher(gamma_soln)
    gamma_soln_sd <- sqrt(diag(solve(gamma_soln_fisher)))

    gamma <- gamma_soln[[1]]
    gamma_sd <- gamma_soln_sd[[1]]
    return(c(gamma, gamma_sd))
}


get_pred <- function(signal_model, linear_model, pred_exp = NULL) {
    exp_limits <- get_exp_limits(signal_model)
    if (is.null(pred_exp)) {
      pred_exp <- seq(exp_limits[1], exp_limits[2], length.out = 100)
    }
    pred_lrr <- get_signal(signal_model, exp_limits[1], pred_exp)
    list(pred_exp = pred_exp, pred_lrr = pred_lrr)
}


get_beta_bounds <- function(linear_model) {
    beta_info <- get_beta(linear_model)
    gamma_info <- get_gamma(linear_model)
    
    beta <- beta_info[1]
    inner_beta_sd <- beta_info[2]
    outer_beta_sd <- sqrt(beta_info[2]**2 + gamma_info[1] + 2.0*gamma_info[2])
    
    c(qnorm(0.05, mean = beta, sd = outer_beta_sd),
      qnorm(0.05, mean = beta, sd = inner_beta_sd),
      qnorm(0.95, mean = beta, sd = inner_beta_sd),
      qnorm(0.95, mean = beta, sd = outer_beta_sd))
}


get_pred_bounds <- function(signal_model, linear_model, pred_exp = NULL) {
    beta_bounds <- get_beta_bounds(linear_model)
    pred <- get_pred(signal_model, linear_model, pred_exp = pred_exp)
    outer(beta_bounds, pred[[2]])
}

get_score <- function(signal_model, linear_model, normalize_to_tmrel = FALSE) {
    pred <- get_pred(signal_model, linear_model)
    pred_bounds <- get_pred_bounds(signal_model, linear_model)
    if (normalize_to_tmrel) {
        index <- which.min(pred[[2]])
        pred[[2]] <- pred[[2]] - pred[[2]][index]
        pred_bounds <- pred_bounds - pred_bounds[, index]
    }

    exp_bounds <- get_exp_bounds(signal_model)
    
    index <- (pred[[1]] >= exp_bounds[1]) & (pred[[1]] <= exp_bounds[2])
    ssign <- sign(mean(pred[[2]][index]))
    min(rowMeans(pred_bounds[c(1, 4),index])*ssign)
}


get_pval <- function(mean, sd, one_sided = FALSE) {
    zscore <- abs(mean/sd)
    if (one_sided) {
        pval <- 1 - pnorm(zscore)
    } else {
        pval <- 2*(1 - pnorm(zscore))
    }
    pval
}


egger_regression <- function(residual, residual_sd, one_sided = TRUE) {
    weighted_residual <- residual/residual_sd
    r_mean <- mean(weighted_residual)
    r_sd <- 1/sqrt(length(weighted_residual))
    r_pval <- get_pval(r_mean, r_sd, one_sided = one_sided)
    list(mean = r_mean, sd = r_sd, pval = r_pval)
}


get_pub_bias <- function(signal_model, linear_model, one_sided = FALSE) {
    signal <- signal_model$predict(signal_model$data)
    gamma_info <- get_gamma(linear_model)
    residual <- signal_model$data$obs - signal
    residual_sd <- sqrt(signal_model$data$obs_se**2 + signal**2*gamma_info[1])
    index <- signal_model$get_w_soln() > 0.1
    residual <- residual[index]
    residual_sd <- residual_sd[index]
    egger_result <- egger_regression(residual, residual_sd, one_sided = one_sided)
    return(as.integer(egger_result[[3]] < 0.05))
}


# function that get meta data
get_meta <- function(
    rei_id,
    cause_id,
    risk_unit,
    signal_model,
    linear_model,
    normalize_to_tmrel = FALSE
) {
    df <- data.frame(
        rei_id = rei_id,
        cause_id = cause_id,
        risk_unit = risk_unit,
        risk_type = "continuous"
    )

    # get midpoints for data exposures
    exp_bounds <- get_exp_bounds(signal_model)
    df$risk_lower <- exp_bounds[[1]]
    df$risk_upper <- exp_bounds[[2]]

    # get scores
    df$score <- get_score(signal_model, linear_model,
                          normalize_to_tmrel = normalize_to_tmrel)

    # publication bias
    df$pub_bias <- get_pub_bias(signal_model, linear_model)
    
    df[c(
        "rei_id",
        "cause_id",
        "risk_lower",
        "risk_upper",
        "risk_unit",
        "risk_type",
        "score",
        "pub_bias"
    )]
}


get_draw <- function(
    rei_id,
    cause_id,
    age_group_id,
    sex_id,
    location_id,
    signal_model,
    linear_model,
    normalize_to_tmrel = FALSE,
    pred_exp = NULL
) {
    pred <- get_pred(signal_model, linear_model, pred_exp = pred_exp)
    pred_bounds <- get_pred_bounds(signal_model, linear_model, pred_exp = pred_exp)
    if (normalize_to_tmrel) {
        index <- which.min(pred[[2]])
        pred[[2]] <- pred[[2]] - pred[[2]][index]
        pred_bounds <- pred_bounds - pred_bounds[, index]
    }
    
    data.frame(
        rei_id = rei_id,
        cause_id = cause_id,
        age_group_id = age_group_id,
        sex_id = sex_id,
        location_id = location_id,
        risk = pred[[1]],
        outer_log_cause_lower = pred_bounds[1,],
        inner_log_cause_lower = pred_bounds[2,],
        log_cause = pred[[2]],
        inner_log_cause_upper = pred_bounds[3,],
        outer_log_cause_upper = pred_bounds[4,],
        outer_linear_cause_lower = exp(pred_bounds[1,]),
        inner_linear_cause_lower = exp(pred_bounds[2,]),
        linear_cause = exp(pred[[2]]),
        inner_linear_cause_upper = exp(pred_bounds[3,]),
        outer_linear_cause_upper = exp(pred_bounds[4,])
    )
}


upload_results <- function(
    signal_model_path,
    linear_model_path,
    results_folder,
    rei_id,
    cause_id,
    risk_unit,
    age_group_id = 22,
    sex_id = 3,
    location_id = 1,
    normalize_to_tmrel = FALSE,
    pred_exp = NULL
) {
    signal_model <- py_load_object(filename = signal_model_path, pickle = "dill")
    linear_model <- py_load_object(filename = linear_model_path, pickle = "dill")
    
    py_save_object(object = signal_model,
                   filename = file.path(results_folder, "signal_model.pkl"),
                   pickle = "dill")
    py_save_object(object = linear_model,
                   filename = file.path(results_folder, "linear_model.pkl"),
                   pickle = "dill")
    
    dataframes <- list(
        study_data = get_data(rei_id, cause_id,
                              signal_model, linear_model,
                              normalize_to_tmrel = normalize_to_tmrel),
        risk_cause_metadata = get_meta(rei_id, cause_id, risk_unit,
                                       signal_model, linear_model,
                                       normalize_to_tmrel = normalize_to_tmrel),
        output_data = get_draw(rei_id, cause_id, age_group_id, sex_id, location_id,
                               signal_model, linear_model,
                               normalize_to_tmrel = normalize_to_tmrel,
                               pred_exp = pred_exp)
    )
    
    for (name in names(dataframes)) {
        write.csv(dataframes[[name]],
                  file.path(results_folder, paste0(name, ".csv")),
                  row.names = FALSE)
    }
    return(dataframes)
}
