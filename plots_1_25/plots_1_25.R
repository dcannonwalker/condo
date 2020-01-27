#packages 
library(dplyr)
library(ggplot2)
library(kdensity)
library(car)
library(tidyr)

#setting up data
dat <- read.csv("core microbiome/plots_1_25/asv112.csv")
head(dat)
colnames(dat) <- c(paste0("BTx623_", 1:8), 
                    paste0("Grassl_", 1:8), 
                    paste0("PI655972_", 1:8), 
                    paste0("Chinese_Amber_", 1:8),
                    paste0("Rio_",1:8))
dat <- dat %>% mutate(propzero_all = apply(dat, 1, function(s){sum(s==0)/length(s)}), 
                        asv = rownames(dat))  %>% filter(propzero_all < 1)

dat_longer <- dat %>% 
  pivot_longer(cols = (-c("asv", "propzero_all")),
               names_to = c("geno", "sample"),
               names_pattern = "(.*)_(.)",
               values_to = "count") %>% 
  group_by(geno, sample) %>% 
  mutate(sampleSum = sum(count), relAbund = count/sampleSum) %>% 
  group_by(geno, sample, asv) %>%
  mutate(samplePresence = count > 0) %>%
  ungroup()

save(dat_longer, file = 'core microbiome/plots_1_25/asv112_longer.RData')

dat_summaries <- dat_longer %>% 
  group_by(geno, asv) %>% 
  summarise(genoPresence = sum(samplePresence),
            cumRelAbund = sum(relAbund), 
            totalCount = sum(count))

save(dat_summaries, file = 'core microbiome/plots_1_25/asv112_summaries.RData')


