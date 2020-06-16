source("R/packages.R")
source("R/functions.R")
source("R/plan.R")
options(clustermq.scheduler = "sge", clustermq.template = "sge.tmpl")
drake_config(
  plan,
  parallelism = "clustermq",
  jobs = 1000,
  caching = "worker",
  recover = TRUE,
  history = FALSE
)
