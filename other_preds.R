library(dplyr)
library(pbapply)

files_path <- "/home/michal/Dropbox/signal_final_dat/"

res <- pblapply(list.files(paste0(files_path, "data/")), function(ith_file) {
  
  deepSig_command <- paste0(files_path, "deepsig-1.0/runDeepSig.sh ", files_path, ith_file, " euk tmp.txt")
  system(deepSig_command)
  
  res <- read.table("tmp.txt") %>% 
    mutate(deepSig = ifelse(V2 == "SignalPeptide", V3, NA),
           id = 1L:nrow(.),
           file = ith_file) %>% 
    select(file, id, name = V1, deepSig)

  file.remove("tmp.txt")
  res
}) %>% 
  do.call(rbind, .)

write.csv(x = res, file = paste0(files_path, "results/res.csv"), row.names = FALSE)
