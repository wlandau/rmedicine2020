#' @title Simulate a single iteration of a clinical trial.
#' @description Simulated patient-level data of time to hospital discharge.
#' @return Data frame with one row per patient and the following fields:
#' * `days`: Observed days until hospital discharged, right-censored at `censor`
#' * `censored`: Logical, whether the patient's `days` value was right-censored.
#' * `arm`: Study arm of the patient.
#' * `patients_per_arm`: Number of patients per arm. Here, there is one
#'   control arm and one treatment arm, with 1:1 allocation.
#' * `mean_control`: Mean time until hospital discharge for the control arm.
#' * `mean_treatment`: Same for the treatment arm.
#' @param mean_control Mean time until hospital discharge for the control arm.
#' @param mean_treatment Mean time until hospital discharge for the treatment arm.
#' @inheritParams simulate_arm
#' @examples
#' simulate_trial(patients_per_arm = 200)
simulate_trial <- function(
  mean_control = 15,
  mean_treatment = 10,
  patients_per_arm = 100,
  censor = 30
) {
  bind_rows(
    simulate_arm(mean_control, censor, patients_per_arm, "control"),
    simulate_arm(mean_treatment, censor, patients_per_arm, "treatment")
  ) %>%
    mutate(
      patients_per_arm = patients_per_arm,
      mean_control = mean_control,
      mean_treatment = mean_treatment
    )
}

#' @title Simulate a single arm of a single iteration of a clinical trial.
#' @description Simulated patient-level data of time to hospital discharge.
#' @return Data frame with one row per patient and the following fields:
#' * `days`: Observed days until hospital discharged, right-censored at `censor`
#' * `censored`: Logical, whether the patient's `days` value was right-censored.
#' * `arm`: Study arm of the patient.
#' @param patients_per_arm Number of patients per arm. Here, there is one
#'   control arm and one treatment arm, with 1:1 allocation.
#' @param censor Number of days at which to right-censor the data.
#' @examples
#' simulate_arm()
simulate_arm <- function(
  mean = 15,
  censor = 30,
  patients_per_arm = 100,
  arm = "control"
) {
  days <- rtruncnorm(
    n = patients_per_arm,
    a = 1,
    mean = mean,
    sd = 10
  )
  tibble(days = days) %>%
    mutate(censored = days > censor) %>%
    mutate(days = pmin(days, censor)) %>%
    mutate(arm = arm)
}

#' @title Analyze simulated trial data with a Bayesian proportional hazards model.
#' @description Uses the `spBayesSurv` package by Zhou, Hanson, and Zhang (2020),
#'   <https://www.jstatsoft.org/article/view/v092i09>.
#' @return Same value as from [summarize_samples()].
#' @inheritParams run_chain
#' @examples
#' patients <- simulate_trial()
#' model_hazard(patients, iterations = 100)
model_hazard <- function(patients, iterations) {
  samples <- map(seq_len(4), ~run_chain(patients, iterations))
  summarize_samples(samples, patients)
}

#' @title Run a single chain of the Bayesian proportional hazards model.
#' @description Uses the `spBayesSurv` package by Zhou, Hanson, and Zhang (2020),
#'   <https://www.jstatsoft.org/article/view/v092i09>.
#' @return A list of output from one call to `spBayesSurv::indeptCoxph()`.
#' @param patients Data frame, a single iteration of a simulated trial
#'   from [simulate_trial()]
#' @param iterations Number of non-warmup MCMC iterations.
#' @examples
#' patients <- simulate_trial()
#' run_chain(patients, iterations = 100)
run_chain <- function(patients, iterations) {
  mcmc <- list(
    nburn = iterations / 2,
    nsave = iterations,
    ndisplay = 2 * iterations,
    nskip = 0
  )
  indeptCoxph(
    formula = survival::Surv(days, !censored) ~ arm,
    data = patients,
    mcmc = mcmc
  )
}

#' @title Summarize a batch of MCMC chains on a single trial.
#' @description Summarize a batch of MCMC chains on a single trial.
#' @return A single data frame with the following fields:
#' * `prob_effect`: Posterior probability that the hazard ratio of
#'   discharge from the hospital is greater than 1.
#' * `median`: Posterior median hazard ratio.
#' * `psrf`: Potential scale reduction factor, an MCMC convergence diagnostic.
#' * `patients_per_arm`: Number of patients per arm.
#' * `mean_control`: True mean time to discharge for the control arm.
#' * `mean_treatment`: True mean time to discharge for the treatment arm.
#' @param samples A list of MCMC chain objects from calls to [run_chain()].
#' @param patients Data frame, a single iteration of a simulated trial
#'   from [simulate_trial()]
#' @examples
#' patients <- simulate_trial()
#' samples <- map(seq_len(4), ~run_chain(patients, 100))
#' summarize_samples(samples, patients)
summarize_samples <- function(samples, patients) {
  hazard_ratio_list <- map(samples, ~as.mcmc(t(exp(.x$beta))))
  hazard_ratio <- unlist(hazard_ratio_list)
  tibble(
    prob_effect = mean(hazard_ratio > 1.5),
    median = median(hazard_ratio),
    psrf = gelman.diag(hazard_ratio_list, multivariate = FALSE)$psrf[, 1],
    patients_per_arm = patients$patients_per_arm[1],
    mean_control = patients$mean_control[1],
    mean_treatment = patients$mean_treatment[1]
  )
}

#' @title Summarize multiple trial simulations.
#' @description Summarize multiple trial simulations from
#'   independent calls to [model_hazard()].
#' @return A single data frame with the following fields:
#' * `prob_success`: Prob(Prob(hazard ratio > 1 | data) > 0.6).
#' * `median`: Median of the posterior median hazard ratios.
#' * `max_psrf`: Max potential scale reduction factor across all simulations.
#' * `patients_per_arm`: Number of patients per arm.
#' * `mean_control`: True mean time to discharge for the control arm.
#' * `mean_treatment`: True mean time to discharge for the treatment arm.
#' @param models A data frame of row-binded output from independent
#'   calls to [model_hazard()].
#' @examples
#' models <- bind_rows(
#'   model_hazard(simulate_trial(), 100),
#'   model_hazard(simulate_trial(), 100)
#' )
#' models
#' summarize_models(models)
summarize_models <- function(models) {
  summarize(
    models,
    prob_success = mean(prob_effect > 0.6),
    median = median(median),
    max_psrf = max(psrf),
    patients_per_arm = patients_per_arm[1],
    mean_control = mean_control[1],
    mean_treatment = mean_treatment[1]
  )
}
