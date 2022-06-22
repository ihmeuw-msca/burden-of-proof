# Evidence score pipeline

**Note: The pipeline is still under development.**

## Configuration

In `config.R`, the variables are global variables that will be used in scripts in the pipeline. Users can change the value of variables based on their tasks.

Users need to change some settings before running the pipeline, including `WORK_DIR`, `PROJ`, `VERSION_ID`, `OUT_DIR`. Other settings are task-specific, including, but not limited to, `INPUT_DATA_DIR`, `OBS_VAR`  and an array of settings for models in each stage. Users are advised to read through the `config.R` file to make sure the settings match their needs.

## Run script

`run_pipeline.R` is the main script to run the pipeline. Users need to change `WORK_DIR` to specify the location where the scripts are saved.

The pipeline consists of five stages:
1. Ensemble model with exposure only to get signal; no random effects
2. Log-linear model to get slope prior for covariate selection
3. Covariate selection model
4. Final mixed effects model that combines the signal and selected covariates
5. Get evidence score with signal model and final model; plots

Each stage is run sequentially but paralleled for risk-outcome pairs. Note that Step 5 cannot be run in parallel on cluster. It need to be run mannully after typing "repl_python()" and "esc". Details can be found in script.

## Results

The plots of evidence score will be saved in `05_evidence_score` under the `OUT_DIR` specified.

## Useful links
[GBD 2020 Guidance on Evidence Score](https://docs.google.com/document/d/1gP7-T6cxah2rLfjaTxWZaO0ejRw7Wk4eVFDGoszetj4/edit)

[Introduction to MR-BRT](https://rpubs.com/rsoren/mrbrt_gbd2020)

[MR-BRT examples](https://rpubs.com/rsoren/mrbrt_examples_gbd2020)