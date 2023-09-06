# sbatch slurm/run_TEMREG.job
# sbatch --account=short4hr slurm/run_TEMREG.job
# Job Wall-clock time: 00:38:54
# Memory Utilized: 28.83 GB 
# requires extra memory; specify on JASMIN:
#SBATCH --mem=60G

here::i_am("./run.R")
library(targets)

source("./_targets.R")

tar_outdated() # what is out of date?

# serial
# system.time(tar_make())
# parallel
system.time(tar_make_future(workers = 6L))

quit(save = "no")
