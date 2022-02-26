if (!require("sleuth")) {
  install.packages("devtools")
  library(devtools)
  install_github("pachterlab/sleuth")
}
ags<- commandArgs(trailingOnly = TRUE)
if (length(ags)==0){
  cat("Usage: runSleuth sample_name1 sample_name2 sample_name3 ...\n
      Will run sleuth on the kallisto results in current directory and output diff genes to SleuthResults.csv\n
      Make sure you have experiment_info.txt in the directory and it has the following format:\n
      sample\tcondition\n
      ctl_1\tctl\n
      ctl_2\tctl\n
      ctl_3\tctl\n
      pull_1\tpull\n
      pull_2\tpull\n
      pull_3\tpull\n")
}
sample_id <- ags
kal_dirs <- file.path(".", sample_id)
s2c <- read.table(file.path(".", "metadata", "experiment_info.txt"), header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::select(s2c, sample = run_accession, condition)
so <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE)
so <- sleuth_fit(so, ~condition, 'full')
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
head(sleuth_significant,25)
write.csv(x = sleuth_significant,file = "SleuthResult.csv",quote = NA,row.names = TRUE,col.names = TRUE)
