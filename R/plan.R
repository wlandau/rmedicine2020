plan <- drake_plan(
  sim = target(
    seq_len(10),
    hpc = FALSE
  ),
  patients = target(
    simulate_trial(
      mean_control = 20,
      mean_treatment = mean_treatment,
      patients_per_arm = 100,
      censor = 30
    ),
    dynamic = map(sim),
    transform = map(mean_treatment = c(10, 15, 20)),
    format = "fst_tbl"
  ),
  models = target(
    model_hazard(patients, 50),
    dynamic = map(patients),
    transform = map(patients, .id = mean_treatment),
    format = "fst_tbl"
  ),
  summaries = target(
    summarize_models(models),
    transform = map(models, .id = mean_treatment),
    format = "fst_tbl"
  ),
  results = target(
    bind_rows(summaries),
    transform = combine(summaries),
    format = "fst_tbl",
    hpc = FALSE
  )
)
