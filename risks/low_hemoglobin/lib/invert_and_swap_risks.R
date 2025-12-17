#' Invert ln_rr, swap exposure bounds, and relabel risk type where applicable
invert_and_swap_risks <- function(dat) {
  
  # Define conditions for inversion and swapping
  condition_low_hb <- dat$alt_risk_upper > 130 & dat$hb_risk_type == "low_hb"
  condition_high_hb <- dat$alt_risk_lower < 130 & dat$hb_risk_type == "high_hb"
  
  # Invert ln_rr
  dat <- dplyr::mutate(
    dat,
    ln_rr = dplyr::case_when(
      condition_low_hb | condition_high_hb ~ -ln_rr,
      TRUE ~ ln_rr
    )
  )
  
  # Swap exposure bounds and update hb_risk_type
  dat <- dplyr::mutate(
    dat,
    original_ref_risk_upper = ref_risk_upper,
    original_ref_risk_lower = ref_risk_lower
  )
  
  dat <- dplyr::mutate(
    dat,
    ref_risk_upper = dplyr::if_else(condition_low_hb | condition_high_hb, alt_risk_upper, original_ref_risk_upper),
    alt_risk_upper = dplyr::if_else(condition_low_hb | condition_high_hb, original_ref_risk_upper, alt_risk_upper),
    
    ref_risk_lower = dplyr::if_else(condition_low_hb | condition_high_hb, alt_risk_lower, original_ref_risk_lower),
    alt_risk_lower = dplyr::if_else(condition_low_hb | condition_high_hb, original_ref_risk_lower, alt_risk_lower),
    
    hb_risk_type = dplyr::case_when(
      condition_low_hb ~ "high_hb",
      condition_high_hb ~ "low_hb",
      TRUE ~ hb_risk_type
    )
  )
  
  return(dat)
}
