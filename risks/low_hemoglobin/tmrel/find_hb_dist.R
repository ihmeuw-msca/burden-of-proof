source('~/repos/bophf/R/hb_model_ensemble_draws.R')
library(ggplot2)

draws_mean <- ihme::get_draws(
  gbd_id_type = 'modelable_entity_id',
  source = 'epi',
  gbd_id = 28869,
  year_id = 2022,
  release_id = 16,
  num_workers = parallelly::availableCores()
)

draws_sd <- ihme::get_draws(
  gbd_id_type = 'modelable_entity_id',
  source = 'epi',
  gbd_id = 28870,
  year_id = 2022,
  release_id = 16,
  num_workers = parallelly::availableCores()
)

merge_draws <- data.table::merge.data.table(
  x = nch::pivot_draws_longer(draws_mean),
  y = nch::pivot_draws_longer(draws_sd),
  by = c('age_group_id', 'year_id', 'sex_id', 'location_id', 'draw_id'),
  suffixes = c('.mean', '.sd')
)

NUM_HB_DRAWS <- 1000

future::plan(
  future::multisession, 
  workers = parallelly::availableCores()
)

hb_dist_df <- furrr::future_map2(
  .x = merge_draws$value.mean,
  .y = merge_draws$value.sd,
  \(mn, stdev) {
    hb_draws <- withr::with_seed(
      seed = 123,
      code = hb_me_ensemble_draws(
        num_draws = NUM_HB_DRAWS,
        mn = mn,
        vr = stdev ^ 2
      )
    ) 
    
    data.table::data.table(
      high_hb = quantile(hb_draws, pnorm(2)),
      low_hb = quantile(hb_draws, pnorm(-2))
    )
  }
) |> data.table::rbindlist()

fst::write.fst(
  x = hb_dist_df,
  path = '~/anemia/bop_pipeline_notes/tmrel/hb_dist.fst'
)

hist(hb_dist_df$high_hb)

low_hb_bound <- hb_dist_df |>
  purrr::chuck('low_hb') |>
  quantile(0.025)

hb_dist_df |>
  purrr::chuck('high_hb') |>
  quantile(0.975)

ggplot(hb_dist_df[low_hb > 0], aes(x = low_hb)) +
  geom_histogram(fill = 'darkgreen', colour = 'black') +
  geom_vline(xintercept = low_hb_bound, colour = 'red', linetype = 'dashed') +
  labs(
    title = 'Distribution of 2.5th percentile of hemoglobin in pregnant WRA for all \nage/sex/loctations in 2022',
    x = 'Hemoglobin (g/L)'
  ) +
  theme_minimal()

ax <- 1:100
y <- x ^ 2

x <- c(-100, x)
y <- c(-2, y)

some_function <- approxfun(x = x, y = y)
some_function(-50)

plot(
  x = -100:100,
  y = some_function(-100:100)
)
