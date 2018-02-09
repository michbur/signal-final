library(seqinr)
library(signalHsmm)
library(dplyr)
library(pbapply)

files_path <- "/home/michal/Dropbox/signal_final_dat/"


all_seqs_res <- pblapply(list.files(paste0(files_path, "data/")), function(ith_file) {
  
  dat <- read.fasta(paste0(files_path, "data/", ith_file), seqtype = "AA")
  res <- lapply(dat, function(ith_seq) try({
    signalHsmm_pred <- run_signalHsmm(ith_seq)
    
    write.fasta(ith_seq, name = attr(ith_seq, "name"), file.out = "tmp_seq.fasta")
    
    signalP_command <- paste0(files_path, "signalp-4.1/signalp -t euk -f short tmp_seq.fasta > signalP.txt")
    deepSig_command <- paste0(files_path, "deepsig-1.0/runDeepSig.sh tmp_seq.fasta euk deepSig.txt")
    
    system(deepSig_command, ignore.stdout = TRUE, ignore.stderr = TRUE)
    system(signalP_command)
    
    res_deepSig <- read.table("deepSig.txt") %>% 
      mutate(deepSig = ifelse(V2 == "SignalPeptide", V3, 1 - V3)) %>% 
      select(deepSig) %>% 
      unlist
    
    res_signalP <- read.table("signalP.txt", skip = 1) %>% 
      mutate(signalP = V9) %>% 
      select(signalP) %>% 
      unlist
    
    file.remove(c("tmp_seq.fasta", "deepSig.txt", "signalP.txt"))
    
    data.frame(file = ith_file, 
               name = attr(ith_seq, "name"),
               signalHsmm = signalHsmm_pred[[1]][["sp_probability"]],
               deepSig = res_deepSig,
               signalP = res_signalP,
               stringsAsFactors = FALSE)
  }, silent = TRUE)) 
  
  nonproblematic <- sapply(res, function(i) class(i) != "try-error")
  do.call(rbind, res[nonproblematic]) %>% 
    mutate(id = cumsum(nonproblematic)) %>% 
    select(file, id, name, signalHsmm, signalP, deepSig)
}) 

save(all_seqs_res, file = "plasmodium_signals.RData")
#write.csv(x = all_seqs_res, file = paste0(files_path, "results/full_res.csv"), row.names = FALSE)
