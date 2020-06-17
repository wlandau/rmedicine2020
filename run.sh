rm -f .RData
module load R/3.6.3
nohup R CMD BATCH run.R &
