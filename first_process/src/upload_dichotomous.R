# Upload evidence score results to a central place
library(mrbrt001, lib.loc = "/ihme/code/mscm/R/packages/")
globals <- new.env()
source("config.R", local = globals)
ARCHIVE <- globals$ARCHIVE

# function that get data
get_data <- function(
    rei_id,
    cause_id,
    model
) {
    # extract data frame
    df <- model$data$to_df()

    # rename columns
    names(df)[names(df) == "obs"] <- "log_rr"
    names(df)[names(df) == "obs_se"] <- "log_rr_se"

    # fill in other columns with NA
    df$rei_id <- rei_id
    df$cause_id <- cause_id

    df$ref_risk_lower <- NA_real_
    df$ref_risk <- NA_real_
    df$ref_risk_upper <- NA_real_

    df$alt_risk_lower <- NA_real_
    df$alt_risk <- NA_real_
    df$alt_risk_upper <- NA_real_

    df$log_ref_cause <- NA_real_
    df$linear_ref_cause <- NA_real_

    df$log_alt_cause_lower <- NA_real_
    df$log_alt_cause <- NA_real_
    df$log_alt_cause_upper <- NA_real_

    df$linear_alt_cause_lower <- NA_real_
    df$linear_alt_cause <- NA_real_
    df$linear_alt_cause_upper <- NA_real_

    df$is_outlier <- as.integer(model$w_soln < 0.1)

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


get_meta <- function(
    rei_id,
    cause_id,
    model
) {
    df <- data.frame(
        rei_id = rei_id,
        cause_id = cause_id,
        risk_unit = NA_character_,
        risk_type = "dichotomous"
    )

    # get score
    df$score <- get_score(model)

    # publication bias
    df$pub_bias <- get_pub_bias(model)

    # fill other column na
    df$risk_lower <- NA_real_
    df$risk_upper <- NA_real_

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
    model
) {
    beta_info <- get_beta(model)
    beta_bounds <- get_beta_bounds(model)
    data.frame(
        rei_id = rei_id,
        cause_id = cause_id,
        age_group_id = age_group_id,
        sex_id = sex_id,
        location_id = location_id,
        risk = NA_real_,
        outer_log_cause_lower = beta_bounds[1],
        inner_log_cause_lower = beta_bounds[2],
        log_cause = beta_info[1],
        inner_log_cause_upper = beta_bounds[3],
        outer_log_cause_upper = beta_bounds[4],
        outer_linear_cause_lower = NA_real_,
        inner_linear_cause_lower = NA_real_,
        linear_cause = NA_real_,
        inner_linear_cause_upper = NA_real_,
        outer_linear_cause_upper = NA_real_
    )
}


upload_results <- function(
    model_path,
    results_folder,
    rei_id,
    cause_id,
    age_group_id = 22,
    sex_id = 3,
    location_id = 1
) {
    model <- py_load_object(filename = model_path, pickle = "dill")
    
    py_save_object(object = model,
                   filename = file.path(results_folder, "model.pkl"),
                   pickle = "dill")
    
    dataframes <- list(
        study_data = get_data(rei_id, cause_id, model),
        risk_cause_metadata = get_meta(rei_id, cause_id, model),
        output_data = get_draw(rei_id, cause_id, age_group_id, sex_id, location_id, model)
    )
    
    for (name in names(dataframes)) {
        write.csv(dataframes[[name]],
                  file.path(results_folder, paste0(name, ".csv")),
                  row.names = FALSE)
    }
    return(dataframes)
}


get_beta <- function(model) {
    beta_soln <- model$beta_soln
    names(beta_soln) <- model$cov_names
    beta <- beta_soln["intercept"]

    model_specs <- mrbrt001::core$other_sampling$extract_simple_lme_specs(model)
    beta_soln_hessian <- mrbrt001::core$other_sampling$extract_simple_lme_hessian(model_specs)
    beta_soln_sd <- sqrt(diag(solve(beta_soln_hessian)))
    names(beta_soln_sd) <- model$cov_names
    beta_sd <- beta_soln_sd[["intercept"]]
    return(c(beta, beta_sd))
}

get_gamma <- function(model) {
    gamma_soln <- model$gamma_soln
    gamma_soln_fisher <- model$lt$get_gamma_fisher(gamma_soln)
    gamma_soln_sd <- sqrt(diag(solve(gamma_soln_fisher)))

    gamma_names <- vector("character", length=length(gamma_soln))
    i <- 1L
    for (cov_model in model$cov_models) {
        if (cov_model$use_re) {
            gamma_names[i] <- cov_model$alt_cov[[1]]
            i <- i + 1L
        }
    }
    gamma <- gamma_soln[gamma_names == "intercept"]
    gamma_sd <- gamma_soln_sd[gamma_names == "intercept"]
    return(c(gamma, gamma_sd))
}

get_beta_bounds <- function(model) {
    beta_info <- get_beta(model)
    gamma_info <- get_gamma(model)

    beta <- beta_info[1]
    inner_beta_sd <- beta_info[2]
    outer_beta_sd <- sqrt(beta_info[2]**2 + gamma_info[1] + 2.0*gamma_info[2])
    
    c(qnorm(0.05, mean = beta, sd = outer_beta_sd),
      qnorm(0.05, mean = beta, sd = inner_beta_sd),
      qnorm(0.95, mean = beta, sd = inner_beta_sd),
      qnorm(0.95, mean = beta, sd = outer_beta_sd))
}

get_score <- function(model) {
    beta_info <- get_beta(model)
    beta_bounds <- get_beta_bounds(model)
    ssign <- sign(beta_info[1])
    0.5*min(beta_bounds[c(1, 4)]*ssign)
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


get_pub_bias <- function(model, one_sided = FALSE) {
    beta_info <- get_beta(model)
    gamma_info <- get_gamma(model)
    residual <- model$data$obs - beta_info[1]
    residual_sd <- sqrt(model$data$obs_se**2 + gamma_info[1])
    index <- model$w_soln > 0.1
    residual <- residual[index]
    residual_sd <- residual_sd[index]
    egger_result <- egger_regression(residual, residual_sd, one_sided = one_sided)
    return(as.integer(egger_result[[3]] < 0.05))
}
