# Convert reported HR with 95% confidence intervals and single exposure range 
# to log RR with associated log SE
# When we are given a per-unit HR, calculate the log effect per unit increase
# Then compute average log RR across the reported exposure range
# This means predicting the log RR at the exposure midpoint using the correct slope in log space
# This average log RR value is what will be used for the reported log RR value in the pipeline
# Values extracted from study

##  example 1:  https://pubmed.ncbi.nlm.nih.gov/31883872/. Fig 1 ARIC Males 

## ARIC males 

hr = 1.02;         # Reported HR
unit = 355;        # Scale of unit increase that the HR is reported for (e.g., a HR of 0.85 per 0.5 unit)
ci_low = 0.86;     # Lower bound of 95% CI
ci_hi = 1.21;      # Upper bound of 95% CI
exp_low = 0;     # Smallest reported exposure level: not reported but we can assume zero
#exp_hi = unknow ;      # Largest reported exposure level: Unknown
# Step 1: calculate log effect (log RR) per unit increase
# Calculating the slopes of the lines for the log RR
hr_m = log(hr)/unit;
# Step 2: calculate average log RR value across reported exposure range
# Predict log RR at the midpoint of the exposure range (slope m = hr_m)
# Use a line that goes through log RR = 0 at the lowest exposure

midpt = 127  ## reported in the study
#midpt = (exp_low + exp_hi)/2;    # calculate midpoint of exposure range
b = -hr_m*exp_low;               # y = mx + b => logRR = 0 at lowest exposure => 0 = hr_m*exp_low + b
log_RR = hr_m*midpt + b;     # evaluate actual line y = mx + b at the midpoint to get the average log RR
# Step 3: transform 95% CIs the same way and calculate standard error
# Calculate slopes
ci_low_m = log(ci_low)/unit;
ci_hi_m = log(ci_hi)/unit;
# Calculate intercepts of line so the 95% CIs also pass through 0 at lowest exposure ranges
b1 = -ci_low_m*exp_low;
b2 = -ci_hi_m*exp_low;
# Evaluate average upper and lower 95% CI boundaries
log_ci_low = ci_low_m*midpt + b1;
log_ci_hi = ci_hi_m*midpt + b2;
# Calculate log standard error: (log(CI_upper) - log(CI_lower))/(2*1.96)
log_SE = (log_ci_hi - log_ci_low)/(2*1.96)
##########

### ````````````````````````````````````````````````````````````````````
### Male: midpoint exposure = 127
## log_RR =0.00708 = exp(Log_RR ) = 1.007, 
## log_ci_low = -0.05395636 = exp(log_ci_low)= 0.9474735,
## log_ci_high =  0.06819376 = exp(log_ci_hi) = 1.070573

## newly calculated HR UI is narrower compared to the orginal 

## newHR = 1.007(0.947, 1.071)

## orginalHR = 1.02(0.86, 1.21)


#````````````````````````````````````````````````````````````````````````




#### ARIC female https://pubmed.ncbi.nlm.nih.gov/31883872/. Fig 1 ARIC females  

hr = 1.06;         # Reported HR
unit = 355;        # Scale of unit increase that the HR is reported for (e.g., a HR of 0.85 per 0.5 unit)
ci_low = 0.82;     # Lower bound of 95% CI
ci_hi = 1.37;      # Upper bound of 95% CI
exp_low = 0;     # Smallest reported exposure level: We can assume zero
#exp_hi = unknown;      # Largest reported exposure level: Unknown
# Step 1: calculate log effect (log RR) per unit increase
# Calculating the slopes of the lines for the log RR
hr_m = log(hr)/unit;
# Step 2: calculate average log RR value across reported exposure range
# Predict log RR at the midpoint of the exposure range (slope m = hr_m)
# Use a line that goes through log RR = 0 at the lowest exposure

midpt = 41
#midpt = (exp_low + exp_hi)/2;    # calculate midpoint of exposure range
b = -hr_m*exp_low;               # y = mx + b => logRR = 0 at lowest exposure => 0 = hr_m*exp_low + b
log_RR = hr_m*midpt + b;     # evaluate actual line y = mx + b at the midpoint to get the average log RR
# Step 3: transform 95% CIs the same way and calculate standard error
# Calculate slopes
ci_low_m = log(ci_low)/unit;
ci_hi_m = log(ci_hi)/unit;
# Calculate intercepts of line so the 95% CIs also pass through 0 at lowest exposure ranges
b1 = -ci_low_m*exp_low;
b2 = -ci_hi_m*exp_low;
# Evaluate average upper and lower 95% CI boundaries
log_ci_low = ci_low_m*midpt + b1;
log_ci_hi = ci_hi_m*midpt + b2;
# Calculate log standard error: (log(CI_upper) - log(CI_lower))/(2*1.96)
log_SE = (log_ci_hi - log_ci_low)/(2*1.96)


##### --------------------------------------------------

## Female: midpoint exposure: 41 ml from the study
## log_RR =0.006729649 = exp(Log_RR ) = 1.00675,
## log_ci_low = -0.02291969= exp(log_ci_low)= 0.977341,
## log_ci_high =  0.06819376 = exp(log_ci_hi) = 1.037027

## newHR = 1.007(0.947, 1.071)

## orginalHR  = 1.06(0.82, 1.37)

#``````````````````````````````````````````


#### ATBC study ------------


##  males 

hr = 1.08;         # Reported HR
unit = 355;        # Scale of unit increase that the HR is reported for (e.g., a HR of 0.85 per 0.5 unit)
ci_low = 0.96;     # Lower bound of 95% CI
ci_hi = 1.21;      # Upper bound of 95% CI
exp_low = 0;     # Smallest reported exposure level
#exp_hi = 370;      # Largest reported exposure level: Unknown
# Step 1: calculate log effect (log RR) per unit increase
# Calculating the slopes of the lines for the log RR
hr_m = log(hr)/unit;
# Step 2: calculate average log RR value across reported exposure range
# Predict log RR at the midpoint of the exposure range (slope m = hr_m)
# Use a line that goes through log RR = 0 at the lowest exposure

midpt = 47
#midpt = (exp_low + exp_hi)/2;    # calculate midpoint of exposure range
b = -hr_m*exp_low;               # y = mx + b => logRR = 0 at lowest exposure => 0 = hr_m*exp_low + b
log_RR = hr_m*midpt + b;     # evaluate actual line y = mx + b at the midpoint to get the average log RR
# Step 3: transform 95% CIs the same way and calculate standard error
# Calculate slopes
ci_low_m = log(ci_low)/unit;
ci_hi_m = log(ci_hi)/unit;
# Calculate intercepts of line so the 95% CIs also pass through 0 at lowest exposure ranges
b1 = -ci_low_m*exp_low;
b2 = -ci_hi_m*exp_low;
# Evaluate average upper and lower 95% CI boundaries
log_ci_low = ci_low_m*midpt + b1;
log_ci_hi = ci_hi_m*midpt + b2;
# Calculate log standard error: (log(CI_upper) - log(CI_lower))/(2*1.96)
log_SE = (log_ci_hi - log_ci_low)/(2*1.96)
##########

exp(log_RR)
exp(log_ci_low)
exp(log_ci_hi)





##### Iowa Women Health study---------- 


#### female 

hr = 0.88;         # Reported HR
unit = 355;        # Scale of unit increase that the HR is reported for (e.g., a HR of 0.85 per 0.5 unit)
ci_low = 0.65;     # Lower bound of 95% CI
ci_hi = 1.18;      # Upper bound of 95% CI
exp_low = 0;     # Smallest reported exposure level
#exp_hi = 370;      # Largest reported exposure level: Unknown
# Step 1: calculate log effect (log RR) per unit increase
# Calculating the slopes of the lines for the log RR
hr_m = log(hr)/unit;
# Step 2: calculate average log RR value across reported exposure range
# Predict log RR at the midpoint of the exposure range (slope m = hr_m)
# Use a line that goes through log RR = 0 at the lowest exposure

midpt = 29.8
#midpt = (exp_low + exp_hi)/2;    # calculate midpoint of exposure range
b = -hr_m*exp_low;               # y = mx + b => logRR = 0 at lowest exposure => 0 = hr_m*exp_low + b
log_RR = hr_m*midpt + b;     # evaluate actual line y = mx + b at the midpoint to get the average log RR
# Step 3: transform 95% CIs the same way and calculate standard error
# Calculate slopes
ci_low_m = log(ci_low)/unit;
ci_hi_m = log(ci_hi)/unit;
# Calculate intercepts of line so the 95% CIs also pass through 0 at lowest exposure ranges
b1 = -ci_low_m*exp_low;
b2 = -ci_hi_m*exp_low;
# Evaluate average upper and lower 95% CI boundaries
log_ci_low = ci_low_m*midpt + b1;
log_ci_hi = ci_hi_m*midpt + b2;
# Calculate log standard error: (log(CI_upper) - log(CI_lower))/(2*1.96)
log_SE = (log_ci_hi - log_ci_low)/(2*1.96)


exp(log_RR)
exp(log_ci_low)
exp(log_ci_hi)



