library(seqinr)
library(signalHsmm)
library(dplyr)
library(pbapply)

files_path <- "/home/michal/Dropbox/signal_final_dat/"

res <- pblapply(list.files(paste0(files_path, "data/")), function(ith_file) {

  dat <- read.fasta(paste0(files_path, "data/", ith_file), seqtype = "AA")
  res <- unlist(lapply(dat, function(ith_seq) try(run_signalHsmm(ith_seq), silent = TRUE)), recursive = FALSE)
  
  lapply(res, function(ith_res) {
    sp.probability  <- if("sp_probability" %in% names(ith_res)) {
      ith_res[["sp_probability"]]
    } else {
      NA
    }
    data.frame(name = ith_res[["name"]], signalHsmm = sp.probability)
  }) %>% 
    do.call(rbind, .) %>% 
    mutate(id = 1L:nrow(.),
           file = ith_file) %>% 
    select(file, id, name, signalHsmm)
}) %>% 
  do.call(rbind, .)

write.csv(x = res, file = paste0(files_path, "results/res.csv"), row.names = FALSE)
