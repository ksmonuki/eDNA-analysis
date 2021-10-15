#Thesis Fix Hashes

library(tidyverse)
library(stringr)
#Hash Keys
input_hash_path_c19 <- "/Users/zackgold/Documents/UCLA_phd/Projects/California/Keira_honors_thesis/anacapa/august2020/hashes_keira.txt"

Hash.key_c19 <- read.table(input_hash_path_c19, header = 1, sep = "\t", stringsAsFactors = F)

Hash.key_c19 %>% View()
head(Hash.key_c19)
dim(Hash.key_c19)
tail(Hash.key_c19)

Hash.key_c19 %>% 
  mutate(.,match_species = if_else( sum.taxonomy_c19 == sum.taxonomy_fishcard ,TRUE,FALSE)) %>% 
  group_by(match_species) %>% 
  tally()

Hash.key_c19 %>% 
  mutate(.,match_species = if_else(sum.taxonomy_c19 == sum.taxonomy_fishcard ,TRUE,FALSE)) %>%
  mutate(., mammal = str_detect(sum.taxonomy_c19,"Mammalia")) %>% 
  mutate(., empty_fish = str_detect(sum.taxonomy_fishcard,"")) %>% 
  mutate(., bird = str_detect(sum.taxonomy_c19,"Aves;")) %>%
  mutate(., actino = str_detect(sum.taxonomy_fishcard,"Eukaryota;Chordata;Actinopteri;;;;")) %>% 
  mutate(., Sebastes = str_detect(sum.taxonomy_c19,"Sebastes")) %>% 
  mutate(., Serranidae = str_detect(sum.taxonomy_c19,"Serranidae")) %>% 
  mutate(., Syngnathiformes = str_detect(sum.taxonomy_c19,"Syngnathiformes")) %>% 
  mutate(., Myliobatiformes = str_detect(sum.taxonomy_c19,"Myliobatiformes")) %>% 
  mutate(., Pleuronectiformes = str_detect(sum.taxonomy_c19,"Pleuronectiformes")) %>% 
  mutate(., Pomacentridae = str_detect(sum.taxonomy_c19,"Pomacentridae")) %>% 
  mutate(., Acanthuriformes = str_detect(sum.taxonomy_c19,"Acanthuriformes")) %>% 
  mutate(., Cottidae = str_detect(sum.taxonomy_c19,"Cottidae")) %>% 
  mutate(., Labriformes = str_detect(sum.taxonomy_c19,"Labriformes")) %>% 
  mutate(., Stromateidae = str_detect(sum.taxonomy_c19,"Stromateidae")) %>% 
  mutate(., Kyphosidae = str_detect(sum.taxonomy_c19,"Kyphosidae")) %>% 
  mutate(., Cichlidae = str_detect(sum.taxonomy_c19,"Cichlidae")) %>% 
  mutate(., Gobiiformes = str_detect(sum.taxonomy_c19,"Gobiiformes")) %>% 
  mutate(., surf_c19 = str_detect(sum.taxonomy_c19,"Eukaryota;Chordata;Actinopteri;;Embiotocidae;;")) %>% 
  mutate(., surf_fish = str_detect(sum.taxonomy_fishcard,"Eukaryota;Chordata;Actinopteri;;Embiotocidae;;")) %>% 
  mutate(., grunt_c19 = str_detect(sum.taxonomy_c19,"Eukaryota;Chordata;Actinopteri;;Sciaenidae;;")) %>% 
  mutate(., grunt_fish = str_detect(sum.taxonomy_fishcard,"Eukaryota;Chordata;Actinopteri;;Sciaenidae;;")) %>% 
  filter(., match_species==FALSE) %>% 
  filter(., empty_fish==TRUE) %>% 
  filter(., bird ==FALSE) %>%
  filter(., mammal == FALSE) %>%
  filter(., actino == FALSE) %>%
  filter(., Sebastes == FALSE) %>% 
  filter(., Serranidae == FALSE) %>% 
  filter(., Syngnathiformes == FALSE) %>%
  filter(., Myliobatiformes == FALSE) %>% 
  filter(., Pleuronectiformes == FALSE) %>%
  filter(., Pomacentridae == FALSE) %>% 
  filter(., Acanthuriformes == FALSE) %>%
  filter(., Cottidae == FALSE) %>%
  filter(., Labriformes == FALSE) %>%
  filter(., Stromateidae == FALSE) %>% 
  filter(., Kyphosidae == FALSE) %>% 
  filter(., Cichlidae == FALSE) %>%
  filter(., Gobiiformes == FALSE) %>%
  filter(., surf_c19 == FALSE) %>% 
  filter(., surf_fish == FALSE) %>% 
  filter(., grunt_c19 == FALSE) %>%
  filter(., grunt_fish == FALSE) %>% View()

Hash.key_c19 %>% 
  mutate(.,match_species = if_else(sum.taxonomy_c19 == sum.taxonomy_fishcard ,TRUE,FALSE)) %>%
  mutate(., mammal = str_detect(sum.taxonomy_c19,"Mammalia")) %>% 
  mutate(., empty_fish = str_detect(sum.taxonomy_fishcard,"")) %>% 
  mutate(., bird = str_detect(sum.taxonomy_c19,"Aves;")) %>%
  mutate(., actino = str_detect(sum.taxonomy_fishcard,"Eukaryota;Chordata;Actinopteri;;;;")) %>% 
  mutate(., Sebastes = str_detect(sum.taxonomy_c19,"Sebastes")) %>% 
  mutate(., Serranidae = str_detect(sum.taxonomy_c19,"Serranidae")) %>% 
  mutate(., Myliobatiformes = str_detect(sum.taxonomy_c19,"Myliobatiformes")) %>% 
  mutate(., Pomacentridae = str_detect(sum.taxonomy_c19,"Pomacentridae")) %>% 
  mutate(., Acanthuriformes = str_detect(sum.taxonomy_c19,"Acanthuriformes")) %>% 
  mutate(., Cottidae = str_detect(sum.taxonomy_c19,"Cottidae")) %>% 
  mutate(., Labriformes = str_detect(sum.taxonomy_c19,"Labriformes")) %>% 
  mutate(., Stromateidae = str_detect(sum.taxonomy_c19,"Stromateidae")) %>% 
  mutate(., Kyphosidae = str_detect(sum.taxonomy_c19,"Kyphosidae")) %>% 
  mutate(., Cichlidae = str_detect(sum.taxonomy_c19,"Cichlidae")) %>% 
  mutate(., Gobiiformes = str_detect(sum.taxonomy_c19,"Gobiiformes")) %>% 
  mutate(., surf_c19 = str_detect(sum.taxonomy_c19,"Eukaryota;Chordata;Actinopteri;;Embiotocidae;;")) %>% 
  mutate(., surf_fish = str_detect(sum.taxonomy_fishcard,"Eukaryota;Chordata;Actinopteri;;Embiotocidae;;")) %>% 
  mutate(., grunt_c19 = str_detect(sum.taxonomy_c19,"Eukaryota;Chordata;Actinopteri;;Sciaenidae;;")) %>% 
  mutate(., grunt_fish = str_detect(sum.taxonomy_fishcard,"Eukaryota;Chordata;Actinopteri;;Sciaenidae;;")) %>% 
  mutate(., to_keep = case_when(match_species==TRUE ~ sum.taxonomy_c19,
                                mammal == TRUE ~sum.taxonomy_c19,
                                bird == TRUE ~sum.taxonomy_c19,
                                actino == TRUE ~sum.taxonomy_c19,
                                Sebastes == TRUE ~sum.taxonomy_c19,
                                Serranidae == TRUE ~sum.taxonomy_c19,
                                Myliobatiformes == TRUE ~sum.taxonomy_c19,
                                Pomacentridae == TRUE ~sum.taxonomy_c19,
                                Cottidae ==TRUE ~sum.taxonomy_c19,
                                Labriformes ==TRUE ~sum.taxonomy_c19,
                                Stromateidae ==TRUE ~sum.taxonomy_c19,
                                Kyphosidae ==TRUE ~sum.taxonomy_fishcard,
                                Cichlidae ==TRUE ~sum.taxonomy_c19,
                                Gobiiformes ==TRUE ~sum.taxonomy_c19,
                                surf_c19 ==TRUE ~sum.taxonomy_c19,
                                surf_fish ==TRUE ~sum.taxonomy_fishcard,
                                grunt_c19 ==TRUE ~sum.taxonomy_fishcard,
                                grunt_fish ==TRUE ~sum.taxonomy_c19
                                )) %>%
  dplyr::select(seq_number=seq_number1,to_keep) -> best_taxa


saveRDS(best_taxa, file="/Users/zackgold/Documents/UCLA_phd/Projects/California/Keira_honors_thesis/anacapa/august2020/best_hashes_chosen_08052020")
