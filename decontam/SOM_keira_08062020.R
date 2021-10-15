#This code is adapted from Ryan Kelly and Ramon Gallego by Zachary Gold
#4/13/20

######################
library(tidyverse)
library(rstan)
library(shinystan)
library(bayesplot)
library(broom)
library(vegan)
library(proxy)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
setwd("/Users/zackgold/Documents/UCLA_phd/Projects/California/Keira_honors_thesis/anacapa/august2020/")
######################

#Load in Data
ASV.nested <- read_rds("/Users/zackgold/Documents/UCLA_phd/Projects/California/Keira_honors_thesis/anacapa/august2020/decontam/Cleaning.before.Occ.model")

input_meta_path <- "/Users/zackgold/Documents/UCLA_phd/Projects/California/Keira_honors_thesis/analysis/decontamination/KM_metadata_4.2.19_edited.txt"
metadata <- read.table(input_meta_path, header = 1, sep = "\t", stringsAsFactors = F)

Hash.key <- readRDS(file="/Users/zackgold/Documents/UCLA_phd/Projects/California/Keira_honors_thesis/anacapa/august2020/best_hashes_chosen_08052020")
Hash.key %>% 
  dplyr::select(seq_number, sum.taxonomy= to_keep) -> Hash.key


#Ryan's STAN Model
#set up data with different hierarchical levels, in case you want to play with these later
ASV.nested %>% 
  dplyr::select(Step2.tibble) %>% 
  unnest(Step2.tibble) %>% # This stage is a long table
  ungroup() %>% 
  mutate(., seq_number=Hash) %>%
  left_join(Hash.key) %>%
  mutate(nReads = 1) %>% 
  left_join(metadata, by=c("sample"="New_name")) %>% 
  dplyr::select(sum.taxonomy, Collected_Date, Location,Replicate, nReads) %>% distinct() -> distinct_species_reps

sum_NA <- function(x) {if (all(is.na(x))) x[NA_integer_] else sum(x, na.rm = TRUE)}

distinct_species_reps %>%  
  group_by(Collected_Date,Location) %>% 
  mutate(tot_rep = n_distinct(Replicate)) %>%
  ungroup() %>% 
  pivot_wider(values_from = nReads, names_from = Replicate, values_fill = list (nReads = 0)) %>% 
  mutate(ndetections = `A`+`B`+`C`) %>% 
  unite(c(`A`,`B`,`C`), col="pattern_presence") %>%
  unite(sum.taxonomy,Collected_Date,Location, col="unique_sample_sum.taxonomy",sep = ":") %>% ungroup() -> Pattern.of.presence

distinct_species_reps %>% 
  unite(sum.taxonomy,Collected_Date,Location, col="unique_sample_sum.taxonomy",sep = ":", remove = FALSE) %>% 
  unite(Collected_Date,Location, col="unique_sample_name",sep = ";", remove = FALSE) %>% 
  left_join(Pattern.of.presence,by ="unique_sample_sum.taxonomy") -> distinct_species_reps


distinct_species_reps %>%  
  pivot_wider(values_from = nReads, names_from = Replicate, values_fill = list (nReads = 0)) %>% 
  rename(K= tot_rep) -> data_tech_summarized

data_tech_summarized %>% 
  group_by(ndetections, K, pattern_presence) %>% 
  slice(1) -> unique_patterns


unique_patterns -> unique_data

# Problem is there are different levels of replication now that we have dropped samples
unique_data %>% View()
#create unique identifier for species; for use in hierarchical modeling
#apparently Ryan's STAN model requires species names to be numbers.
SS_species <- unite(data = unique_data,
                    col = SS_species,
                    c("sum.taxonomy")
) %>% pull(SS_species)
unique_data$Species <- match(SS_species, unique(SS_species)) #index for unique site-species combinations

#create unique identifier for combinations of site-species; for use in hierarchical modeling
SS <- unite(data = unique_data,
            col = SS,
            c("Collected_Date","Location", "sum.taxonomy")
) %>% pull(SS)
unique_data$SiteSpecies <- match(SS, unique(SS)) #index for unique site-species combinations

#Fix to match Ryan's STAN model format
unique_data %>%
  ungroup() %>% 
  mutate(., N= ndetections) %>%
  dplyr::select(-ndetections)-> unique_data

unique_data$N <- as.numeric(unique_data$N)

##Stan Model
sink("Stan_SOM_hierarchical_keira.stan")
cat(
  "data{/////////////////////////////////////////////////////////////////////
  int<lower=1> S;    // number of samples (nrow)
  int<lower=1> Species[S];    // index of species, each of which will have a different value for p11 and p10
  int<lower=1> Nspecies;    // number of species, each of which will have a different value for p11 and p10
  int<lower=1> L[S];   // index of locations or species/site combinations, each of which will have a different value psi
  int<lower=1> Nloc;   // number of locations or species/site combinations, each of which will have a different value psi
  int<lower=1> K[S];   // number of replicates per site (ncol)
  int<lower=0> N[S]; // number of detections among these replicates
  int z[S];   // integer flag to help estimate psi parameter
  }
  parameters{/////////////////////////////////////////////////////////////////////
  real<lower=0,upper=1> psi[Nloc];  //commonness parameter
  real<lower=0,upper=1> p11[Nspecies]; //true positive detection rate
  real<lower=0,upper=1> p10[Nspecies]; //false positive detection rate
  }
  transformed parameters{/////////////////////////////////////////////////////////////////////
  }
  model{/////////////////////////////////////////////////////////////////////
  real p[S];
  
  for (i in 1:S){
  z[i] ~ bernoulli(psi[L[i]]);
  p[i] = z[i]*p11[Species[i]] + (1-z[i])*p10[Species[i]];
  N[i] ~ binomial(K[i], p[i]);
  }; 
  
  //priors
  psi ~ beta(2,2); 
  p11 ~ beta(2,2); 
  p10 ~ beta(1,10);
  }
  generated quantities{
  real<lower=0,upper=1> Occupancy_prob[S];    //after inferring parameters above, now calculate occupancy probability for each observation. Equation from Lahoz-Monfort et al. 2015
  
  for (i in 1:S){
  Occupancy_prob[i]  = (psi[L[i]]*(p11[Species[i]]^N[i])*(1-p11[Species[i]])^(K[i]-N[i])) 
  / ((psi[L[i]]*(p11[Species[i]]^N[i])*(1-p11[Species[i]])^(K[i]-N[i])) 
  + (((1-psi[L[i]])*(p10[Species[i]]^N[i]))*((1-p10[Species[i]])^(K[i]-N[i])))
  );
  }
  }
  
  ",
  fill=TRUE)
sink()

#####################
#run Stan model
#note this will take a while the first time you run a particular model, because it needs to compile from C++
#####################      
myHierarchicalModel <- stan(file = "Stan_SOM_hierarchical_keira.stan", 
                            data = list(
                              S = nrow(unique_data),
                              Species = unique_data$Species,
                              Nspecies = length(unique(unique_data$Species)),
                              L = unique_data$SiteSpecies,
                              Nloc = length(unique(unique_data$SiteSpecies)),
                              K = unique_data$K,
                              N = unique_data$N,
                              z = ifelse(unique_data$N > 0, 1, 0)
                            ), 
                            chains = 10,   #number of chains
                            iter = 10000   #number of iterations per chain
)
unique_data %>%  View()
myHierarchicalStanResults <- tidy(as.data.frame(myHierarchicalModel))   
#launch_shinystan(myHierarchicalModel)
saveRDS(myHierarchicalStanResults, file="keira_myHierarchicalStanResults.rds")
myHierarchicalStanResults <- readRDS(, file="keira_myHierarchicalStanResults.rds")

myHierarchicalStanResults %>% View()

#see histogram of occupancy probabilities
myHierarchicalStanResults %>% 
  filter(grepl("Occupancy_prob", column)) %>% 
  pull(mean) %>% 
  hist()
head(myHierarchicalStanResults)
##NOTE -- below here doesn't work now, because the Occupancy Probs aren't numbered in the same way as the actual parameter estimates. 

#Now Re-attach SOM to all Hash-Site combinations
myHierarchicalStanResults %>% 
  filter(grepl("Occupancy_prob", column)) %>%
  separate(column, into=c("column","SiteSpecies"), sep="([\\[\\]])") %>% 
  pivot_wider(., names_from=column, values_from = c(mean, sd))-> myHierarchicalStanResults_wide
myHierarchicalStanResults_wide$SiteSpecies <- as.numeric(myHierarchicalStanResults_wide$SiteSpecies)

unique_data %>% 
  left_join(myHierarchicalStanResults_wide) %>%  
  dplyr::select(K,pattern_presence, mean_Occupancy_prob,sd_Occupancy_prob) %>% 
  unite(K, pattern_presence, sep=":", col="identif") -> som_results

som_results %>%  View()

data_tech_summarized %>% 
  unite(K, pattern_presence, sep=":", col="identif") %>% 
  left_join(som_results, by="identif") %>%  separate(identif, sep = ":", into = c("K","pattern_presence")) -> som_results_merged


som_results_merged %>% 
  pull(mean_Occupancy_prob) %>% 
  hist()

som_results_merged %>%
  dplyr::select(K,pattern_presence, mean_Occupancy_prob) %>% 
  group_by(pattern_presence,K) %>% 
  summarise(max_Occupancy_probmax=max(mean_Occupancy_prob)) %>% arrange(desc(max_Occupancy_probmax))%>% View()

#Observations
#0.9 is 2 single detections in 3 and 2 tech reps

som_results_merged %>%
  dplyr::select(sum.taxonomy,Collected_Date,Location,K,pattern_presence, mean_Occupancy_prob) %>% 
  group_by(sum.taxonomy,Collected_Date,Location) %>% 
  summarise(max_Occupancy_prob=max(mean_Occupancy_prob)) %>% arrange(desc(max_Occupancy_prob)) -> som_results_merged_summarized



#Summarize
som_results_merged_summarized %>% 
  ungroup() %>% 
  mutate(., to_keep = case_when(max_Occupancy_prob < 0.75 ~ "Toss",
                                max_Occupancy_prob > 0.90 ~ "Keep",
                                TRUE ~ "Borderline")) %>% 
  count(to_keep) %>% 
  mutate(., total=sum(n)) %>% 
  mutate(per=paste0(round(100*n/total,2),'%')) %>% 
  dplyr::select(-total)

#Save Data
som_results_merged_summarized -> hashsite_som_all

saveRDS(hashsite_som_all,"keira_occupancy_results_all.RDS")

hashsite_som_all %>% 
  write.csv(paste0("keira_Data_wOccupancy_all_",Sys.Date(),".csv"), row.names = F)

hashsite_som_all %>% 
  filter(., max_Occupancy_prob > 0.9)-> hashsite_som_to_keep

saveRDS(hashsite_som_to_keep,"keira_occupancy_results_9.RDS")

hashsite_som_to_keep %>% 
  write.csv(paste0("keira_Data_wOccupancy_9_",Sys.Date(),".csv"), row.names = F)
