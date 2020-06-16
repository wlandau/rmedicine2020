simulate_patients <- function(mean_control, mean_treatment, n, censor) {
  bind_rows(
    simulate_arm(mean_control, censor, n, "control"),
    simulate_arm(mean_control, censor, n, "treatment")
  ) %>%
    mutate(n = n)
}

simulate_arm <- function(mean, censor, n, arm) {
  tibble(days = rexp(n = n, rate = 1 / mean)) %>%
    mutate(censored = days >= censor) %>%
    mutate(days = pmin(days, censor)) %>%
    mutate(arm = arm)
}

model_hazard <- function(patients) {
  summarize_fit(map(seq_len(4), ~run_chain(patients)), n)
}

run_chain <- function(patients) {
  indeptCoxph(
    formula = survival::Surv(days, !censored) ~ arm,
    data = patients,
    mcmc = list(nburn = 1e3, nsave = 1e3, ndisplay = 2e3, nskip = 0)
  )
}

summarize_fit <- function(fit, n) {
  hr_list <- map(fit, ~as.mcmc(t(exp(.x$beta))))
  hr <- unlist(hr_list)
  tibble(
    prob_success = mean(hr > 1.25),
    median = median(hr),
    psrf = gelman.diag(hr_list, multivariate = FALSE)$psrf[, 1],
    n = n
  )
}
