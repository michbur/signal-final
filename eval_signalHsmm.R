library(seqinr)
library(signalHsmm)
library(dplyr)
library(pbapply)

files_path <- "/home/michal/Dropbox/signal_final_dat/"

res <- pblapply(list.files(paste0(files_path, "data/")), function(ith_file) {
  dat <- read.fasta(paste0(files_path, "data/", ith_file), seqtype = "AA")
  try({run_signalHsmm(dat) %>% 
    pred2df %>% 
    mutate(file_name = ith_file, seq_name = rownames(.)) %>% 
    select(file_name, seq_name, sp.probability)
  }, silent = TRUE)
})

save(res, file = paste0(files_path, "results/res.RData"))
