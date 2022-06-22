# Post Analysis

Here, we document the processes after the model fitting, including

* detect and adjust for publication-bias
* compute evidence score
* output diagnostic figures and data

## Detect Publication-Bias

Publication-bias analysis is an important part of systematic review.
As a metric of evaluating the evidence in the dataset, evidence score needs to take
publication-bias into account.

To detect publication-bias, we use a data-driven approach known as [Egger's Regression](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2127453/).
The idea is simple. We want to detect if there is a significant correlation between the
residuals and their standard errors. Egger's regression function can be found [here](https://stash.ihme.washington.edu/users/jiaweihe/repos/evidence_score_pipeline/browse/src/utils/continuous_functions.R#106-112).

Interestingly, we find out that the trimming algorithm helps against the publication-bias.
In the process, we apply Egger's Regression on both untrimmed and trimmed data.
Examples can be found [here](https://stash.ihme.washington.edu/users/jiaweihe/repos/evidence_score_pipeline/browse/src/05_evidence_score_continuous.R#39-42).

## Adjust for Publication-Bias

If publication-bias has been detected, to adjust for it, we use an algorithm called [The Trim-and-Fill Method](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6571372/).
The method is based on the assumption that residual should be distributed symmetrically, and if there is asymmetry, it means that there are certain studies missing.
The algorithm involve two major steps, iteratively trim data to get accurate mean estimation, and fill in the "missing" data based on the mean estimation and get final result.

Our cases are slightly different than the ones considered in "The Trim-and-Fill Method".
Nonetheless, we could modify and applied to our problem.
One major change we need to make is that we remove the "trim" step, since we have our own trimming and we trust our mean estimation. And our adjustment process involves

* use the rank statistics and the residual to compute the number of points need to be filled
* fill the data and re-fit the model

Create filled data function can be found [here](https://stash.ihme.washington.edu/users/jiaweihe/repos/evidence_score_pipeline/browse/src/utils/continuous_functions.R#90-104).
And re-fit the model step can be found [here](https://stash.ihme.washington.edu/users/jiaweihe/repos/evidence_score_pipeline/browse/src/05_evidence_score_continuous.R#53-70).

## Get Scores and Diagnostics

Finally, we need to get [the evidence scores](https://stash.ihme.washington.edu/users/jiaweihe/repos/evidence_score_pipeline/browse/src/05_evidence_score_continuous.R#73-78) and diagnostics including [risk function plot](https://stash.ihme.washington.edu/users/jiaweihe/repos/evidence_score_pipeline/browse/src/05_evidence_score_continuous.R#88-93), [residual funnel plot](https://stash.ihme.washington.edu/users/jiaweihe/repos/evidence_score_pipeline/browse/src/05_evidence_score_continuous.R#86) and [summary dataframe](https://stash.ihme.washington.edu/users/jiaweihe/repos/evidence_score_pipeline/browse/src/05_evidence_score_continuous.R#96-103).
Notice that for dichotomous outcomes there is no plot residual step, because it will be exactly the same with the model plot.

## File Structure

These processes is organized in `src/05_evidence_score_*.R` and their corresponding functions is in `src/utils/*_functions.R`. Currently we have the scripts for continuous and dichotomous outcomes.
The old evidence score step is saved as `src/05_evidence_score_legacy.R`.

