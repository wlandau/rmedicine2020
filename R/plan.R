plan <- drake_plan(
  sim = target(
    seq_len(10),
    hpc = FALSE
  ),
  patients = target(
    simulate_trial(
      mean_control = 15,
      mean_treatment = 10,
      patients_per_arm = patients_per_arm,
      censor = 30
    ),
    dynamic = map(sim),
    transform = map(patients_per_arm = c(100, 200, 300, 400)),
    format = "fst_tbl"
  ),
  models = target(
    model_hazard(patients, 100),
    dynamic = map(patients),
    transform = map(patients, .id = patients_per_arm),
    format = "fst_tbl"
  ),
  summaries = target(
    summarize_models(models),
    transform = map(models, .id = patients_per_arm),
    format = "fst_tbl"
  ),
  results = target(
    bind_rows(summaries),
    transform = combine(summaries),
    format = "fst_tbl",
    hpc = FALSE
  )
)
