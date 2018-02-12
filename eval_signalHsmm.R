library(seqinr)
library(signalHsmm)
library(dplyr)
library(pbapply)

files_path <- "/home/michal/Dropbox/signal_final_dat/"

# in i[["name"]] : subscript out of bounds
all_seqs_res <- pblapply(list.files(paste0(files_path, "data/")), function(ith_file) 
  try({
    dat <- read.fasta(paste0(files_path, "data/", ith_file), seqtype = "AA")
    
    res <- run_signalHsmm(dat[lengths(dat) > 1])
    
    res_signalHsmm <- pred2df(res) %>% 
      mutate(name = rownames(.), 
             id = 1L:nrow(.),
             file = ith_file,
             signalHsmm = sp.probability) %>% 
      select(file, id, name, signalHsmm)
    
    signalP_command <- paste0(files_path, "signalp-4.1/signalp -t euk -f short ", 
                              paste0(files_path, "data/", ith_file), " > signalP.txt")
    system(signalP_command)
    res_signalP <- read.table("signalP.txt", skip = 1, stringsAsFactors = FALSE) %>% 
      mutate(signalP = V9) %>% 
      select(name = V1, signalP)
    
    deepSig_command <- paste0(files_path, "deepsig-1.0/runDeepSig.sh ", 
                              paste0(files_path, "data/", ith_file), " euk deepSig.txt")
    system(deepSig_command, ignore.stdout = TRUE, ignore.stderr = TRUE)
    res_deepSig <- read.table("deepSig.txt", stringsAsFactors = FALSE) %>% 
      mutate(deepSig = ifelse(V2 == "SignalPeptide", V3, 1 - V3)) %>% 
      select(name = V1, deepSig)
    
    file.remove(c("deepSig.txt", "signalP.txt"))
    
    res_signalHsmm %>% 
      inner_join(res_signalP, by = c("name" = "name")) %>% 
      inner_join(res_deepSig, by = c("name" = "name"))
  }, silent = TRUE)
)

#write.csv(x = all_seqs_res, file = paste0(files_path, "results/full_res.csv"), row.names = FALSE)
save(all_seqs_res, file = paste0(files_path, "results/full_res.RData"))

load(paste0(files_path, "results/full_res.RData"))

write.csv(x = do.call(rbind, all_seqs_res), 
          file = paste0(files_path, "results/full_res.csv"), row.names = FALSE)
