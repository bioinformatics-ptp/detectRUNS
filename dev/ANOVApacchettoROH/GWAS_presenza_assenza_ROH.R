#### ROH ####
library(dplyr) ; library(detectRUNS) ; library(tidyr) ; library(purrr) ; library(broom)
setwd("/Users/alberto/Library/CloudStorage/OneDrive-UniversitàdegliStudidiSassari/LAVORI/CULUCCIA/")
rm(list=ls()) ; gc() ; if(!is.null(dev.list())) {dev.off()} ; cat("\14")

pheno=data.table::fread("pheno.txt")
roh=consecutiveRUNS.run("BIRU_puliti.ped","BIRU_puliti.map",minSNP=103,minLengthBps=2e6) %>%
  mutate(mb=lengthBps/1e6) %>% mutate(nome=paste0("chr",chrom,"_",from,"_",to))

nIDS=nlevels(as.factor(roh$id))

# tenere le ROH trovate in almeno x% di animali
n_id=n_distinct(roh$ID)
roh_valide=roh %>% distinct(id, nome) %>% count(nome) %>% filter(n / n_id >= 0.10) %>% pull(nome)
roh_condivise=roh %>% 
  filter(nome %in% roh_valide) %>% select(id,nome) %>% 
  mutate(val = 1) %>% distinct(id, nome, .keep_all = T) %>%
  pivot_wider(names_from=nome,values_from=val,values_fill = 0)

# testo effetto della presenza/assenza delle ROH sul fenotipo
ok=pheno %>% inner_join(roh_condivise,"id")
cols_test=names(ok)[(which(names(ok) == "pheno") + 1):ncol(ok)]
risultati=map_dfr(cols_test, function(col) {
  roh_trovate=sum(ok[[col]] == 1, na.rm = T)
  fit=lm(pheno ~ factor(ok[[col]]),data = ok)
  tidy(fit) %>% filter(term == "factor(ok[[col]])1") %>%
    transmute(roh=col,roh_trovate=roh_trovate,beta=estimate,pvalue=p.value)})




