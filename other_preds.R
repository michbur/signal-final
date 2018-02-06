library(dplyr)
library(pbapply)

files_path <- "/home/michal/Dropbox/signal_final_dat/"

res <- pblapply(list.files(paste0(files_path, "data/")), function(ith_file) {
  
  browser()
  signalP_command <- paste0(files_path, "signalp-4.1/signalp -t euk -f short ", 
                            files_path, "data/", ith_file, " > signalP.txt")
  deepSig_command <- paste0(files_path, "deepsig-1.0/runDeepSig.sh ", 
                            files_path, "data/", ith_file, " euk deepSig.txt")
  system(deepSig_command, ignore.stdout = TRUE)
  system(signalP_command, ignore.stdout = TRUE)
  #V1            V2 V3 V4
  #1 XP_001613005.1; SignalPeptide  1 25
  #2 XP_001613005.1; SignalPeptide  1 23
  res_deepSig <- read.table("deepSig.txt") %>% 
    mutate(deepSig = ifelse(V2 == "SignalPeptide", V3, NA),
           file = ith_file) %>% 
    filter(!duplicated(.)) %>% 
    select(file, name = V1, deepSig)

  res_signalP <- read.table("signalP.txt", skip = 1) %>% 
    mutate(signalP = V9,
           id = 1L:nrow(.),
           file = ith_file) %>% 
    select(file, id, name = V1, signalP)
  
  setdiff(as.character(res_deepSig[["name"]]),
          as.character(res_signalP[["name"]]))
  filter(res_deepSig, name == "XP_001612976.1;")
  
  read.table("deepSig.txt")
  
  file.remove("deepSig.txt")
  res
}) %>% 
  do.call(rbind, .)

write.csv(x = res, file = paste0(files_path, "results/res.csv"), row.names = FALSE)
