#Data analysis for Story Lake Indiana
#Nora Schlenker
#October 2023

#rm(list = ls())
# Load Packages --------
library(tidyverse)
library(readxl)
library(Bchron)
library(neotoma2)
library(topicmodels)
library(bcp)


# Load data ----------
# read in all files from the excel sheet STOR19_final and put them all in a list 
sheet_names <- excel_sheets("./data_raw/STOR19_final.xlsx")
STOR19_data <- list()
for (i in 1:length(sheet_names)) {
  st <- sheet_names[i]
  temp <- read_xlsx("./data_raw/STOR19_final.xlsx", sheet = st, col_names = TRUE)
  STOR19_data[[i]] <- temp
}

names(STOR19_data) <- sheet_names

rm(i, st, temp, sheet_names)

#NOTE: data can also be downloaded from the Neotoma Paleoecology Database, for example on how to transform Neotoma data see the Spicer and/or Pretty Lake portion of this code

# Age-depth Model---------
# Story Lake age-depth model
# run age depth model with Bchron using 15 AMS dates, the coretop, and Ambroisa Rise
story_ams <- STOR19_data[["AMS"]]
story_ams_ages <- as.numeric(unlist(story_ams[,2]))
story_ams_error <- as.numeric(unlist(story_ams[,3]))
story_ams_depth <- as.numeric(unlist(story_ams[,4]))
story_ams_calcurve <- c(rep("normal", 4), rep("intcal20", 13))

story_cal_ages <- BchronCalibrate(ages=story_ams_ages,
                ageSds=story_ams_error,
                calCurves=story_ams_calcurve)
apply(sampleAges(story_cal_ages), 2, quantile, prob=c(0.025,0.5, 0.975))

# Calculate chronology with Bchron using default settings 
story_chronology <- Bchronology(ages = story_ams_ages, 
                                ageSds =  story_ams_error, 
                                positions = story_ams_depth, 
                                calCurves = story_ams_calcurve)
plot(story_chronology) #+ scale_x_reverse(limits = c(13500, -500))

# Use chronology to predict midpoint ages every centimeter 
story_ages_predict <- predict(story_chronology, 
                                 newPositions = seq(from = 0.5, to = 1200.5, by = 1))
story_ages_midpoint <- apply(story_ages_predict, 2, mean)
story_ages_midpoint <- data.frame(depth =  seq(from = 0.5, to = 1200.5, by = 1), ages = story_ages_midpoint)

# Calculate the 95th confidence of age estimates
story_ages_95conf <- apply(story_ages_predict, 2, function(x){
  quantile(x, probs = c(0.025, 0.975))
  })
story_ages_95conf <- t(story_ages_95conf)
story_ages_combo <- data.frame(story_ages_midpoint, "lower" = story_ages_95conf[,1], "upper" = story_ages_95conf[,2])

# Calculate sedimentation rate 
story_sedrate <- summary(story_chronology, type = 'sed_rate', useExisting = FALSE,
                         probs=c(0.25,0.5,0.75))

# remove temporary elements
rm(story_ams, story_ams_ages, story_ams_calcurve, story_ams_depth, story_ams_error)

# Pollen data wrangling -------------
ecol_groups <- read.csv("./data_raw/ecolgroup_taxa.csv") #NOTE: if downloaded from Neotoma ecological groups are already included

story_pollen_raw <- STOR19_data[["Pollen"]]
story_control <-  dplyr::select(story_pollen_raw, c(Depth, Control, Control_conc)) #extract control
story_pollen_raw <- dplyr::select(story_pollen_raw, c(-Control, -Control_conc))#remove control

#convert to long format, remove taxa that have no counts, remove aquatic taxa
story_pollen_long <- tibble(story_pollen_raw[, -which(colSums(story_pollen_raw, na.rm = T) == 0)]) %>%  
  pivot_longer(cols = -Depth, names_to = "variablename") %>%
  left_join(distinct(ecol_groups), by = c("variablename" = "taxonname")) %>%
  dplyr::filter(ecolgroupid != "AQVP")
  
# group certain taxa into genus gropus and group all unknowns together
story_pollen_grouped_counts_long <- story_pollen_long %>%
  mutate(variablename = replace(variablename, stringr::str_detect(variablename, "Pinus*"), "Pinus"),
         variablename = replace(variablename, stringr::str_detect(variablename, "Picea*"), "Picea"),
         variablename = replace(variablename, stringr::str_detect(variablename, "Acer*"), "Acer"),
         variablename = replace(variablename, stringr::str_detect(variablename, "Juglans*"), "Juglans"),
         variablename = replace(variablename, stringr::str_detect(variablename, "Fraxinus*"), "Fraxinus"),
         variablename = replace(variablename, stringr::str_detect(ecolgroupid, "UNID"), "Unknown")
  ) %>%
  group_by(variablename, Depth) %>%
  summarise(value = sum(value), .groups = 'keep') %>%
  ungroup()
#convert pollen counts to proportions and save as wide format and add ages
story_pollen_grouped_prop_wide <- story_pollen_grouped_counts_long %>%
  group_by(Depth) %>%
  mutate(count = sum(value, na.rm = TRUE)) %>%
  group_by(variablename) %>% 
  mutate(prop = value / count) %>% 
  arrange(desc(Depth)) %>%
  tidyr::pivot_wider(id_cols = Depth,
                     names_from = variablename, 
                     values_from = prop,
                     values_fill = 0) %>%
  left_join(story_ages_midpoint, by = c("Depth" = "depth")) 

#convert grouped pollen counts to wide format and add ages
story_pollen_grouped_counts_wide <- story_pollen_grouped_counts_long %>%
  pivot_wider(id_cols = Depth, names_from = variablename, values_from = value, values_fill = 0) %>%
  left_join(story_ages_midpoint, by = c("Depth" = "depth"))

#extract fagus and quercus pollen with depth and age
story_fagus <- dplyr::select(story_pollen_grouped_prop_wide, Depth, ages, "Fagus grandifolia")
colnames(story_fagus) <- c("depth", "ages", "fagus")
story_quercus <- dplyr::select(story_pollen_grouped_prop_wide, Depth, ages, "Quercus")
colnames(story_quercus) <- c("depth", "ages", "quercus")

story_pollen_age_conf <- left_join(story_fagus[, c(1,3)], story_ages_combo, by = c("depth" = "depth"))

## Pollen accumulation rate ------------------
story_control$age <- story_pollen_grouped_counts_wide$ages
story_control$total_pollen <- rowSums(story_pollen_grouped_counts_wide[,c(-1, -length(story_pollen_grouped_counts_wide))]) 
#story_accrate <- rbind(c(0.5, 1),story_accrate)
story_control <- left_join(story_control, story_accrate, by = c("Depth" = "age"))
story_control <- mutate(story_control, accu_rate = (((total_pollen*Control_conc)/Control)*acc_rate))
pollen_accurate_func <- function(pollen, control_conc, control, accu_rate) {
  (((pollen*control_conc)/control)*accu_rate)
}
story_pollen_grouped_counts_wide_only <- story_pollen_grouped_counts_wide[,c(-1, -length(story_pollen_grouped_counts_wide))]
story_pollen_acc <- ((story_pollen_grouped_counts_wide_only*story_control$Control_conc)/story_control$Control)*story_control$acc_rate

pollen_sample_ageres <- story_pollen_grouped_prop_wide$ages[1:(length(story_pollen_grouped_prop_wide$ages)-1)]-story_pollen_grouped_prop_wide$ages[2:length(story_pollen_grouped_prop_wide$ages)]
mean(pollen_sample_ageres)
sd(pollen_sample_ageres)

A <- story_pollen_grouped_counts_wide[,c(2:4,8:13,17:19,21:28,32,33,35,36,38:41)]
NonA <- story_pollen_grouped_counts_wide[,c(5:7,14:16,20,29:31,34,37,43,44)]

A_NAP <- data.frame(ages = story_pollen_grouped_counts_wide$ages, "AP/NAP" = rowSums(A)/rowSums(NonA))

# Topic Analysis - CTM ----------
#create defaults for multiple numbers of topics 
nlist <- 2:8 #set number of topics
reptimes <- length(nlist) #set number or replications for below
#Used controls suggested for CTM from Grün and Hornik 2011
myctrl <-
  list(
    estimate.beta = TRUE,
    verbose = 0,
    prefix = tempfile(),
    save = 0,
    keep = 0,
    seed = as.integer(Sys.time()),
    nstart = 1L,
    best = TRUE,
    var = list(iter.max = 500, tol = 10 ^ -6),
    em = list(iter.max = 1000, tol = 10 ^ -4),
    initialize = "random",
    cg = list(iter.max = 500, tol = 10 ^ -5)
  )

#make repeated lists of pollen counts and controls to input into mapply 
poll_only <- dplyr::select(story_pollen_grouped_counts_wide, -Depth, -ages)
poll_reps <- do.call("list", replicate(reptimes, as.matrix(poll_only), simplify = FALSE))
ctrls <- do.call("list", replicate(reptimes, myctrl, simplify = FALSE))

# Run the CTM topic model
ctm_mods <- mapply(CTM, k=nlist, x=poll_reps, control = ctrls)

# Evaluate the CTM models
aicsctm <- do.call(rbind,lapply(ctm_mods, AIC))
bicsctm <- do.call(rbind,lapply(ctm_mods, BIC))
eval_ctm <- data.frame(k = nlist, aic = aicsctm, bic = bicsctm)
eval_ctm_plot <- ggplot(eval_ctm) + geom_point(aes(k, aicsctm), col = "blue") + geom_point(aes(k, bicsctm), col = "orange") + xlim(c(3,8)) + xlab("Number of Topics") + ylab("Model Fit (AIC and BIC)") + theme_minimal() 
eval_ctm_plot
ggsave("./figures/topic_plots/eval_ctm.jpeg", eval_ctm_plot)

#function to extract terms and topics and make basic plots for each number of groups
termtopic_func <- function(lda_result, age_vector = story_pollen_grouped_counts_wide$ages){
  postr <- lapply(lda_result, posterior)
  
  ## Plot "terms" 
  terms_dfs <- lapply(postr, function(x){
    x <- x[["terms"]]
    x <- data.frame(x)
    x$ngroups <- c(rep(nrow(x), nrow(x)))
    x$topic <- factor(1:nrow(x))
    pivot_longer(x, c(-ngroups, -topic), names_to = "taxa", values_to = "beta")
  })
  
  term_plot <- lapply(terms_dfs, function(x){
    x<- x %>%
      group_by(ngroups, topic) %>%
      top_n(5, beta) %>%
      ungroup() %>%
      mutate(taxa = reorder_within(taxa, beta, topic))
    ggplot(x, aes(x = reorder(taxa, beta), beta, fill = topic)) + geom_col() +
      facet_wrap(~topic, scales = "free", drop = TRUE, ncol = 2) +
      scale_x_reordered() +
      coord_flip()
  })
  
  ## Plot "topics"
  topic_dfs <- lapply(postr, function(x){
    x <- x[["topics"]]
    x <- data.frame(x)
    colnames(x) <- factor(1:ncol(x))
    x$sample <- 1:nrow(x)
    #$depth <- story_depth
    x$age <- age_vector
    x<- pivot_longer(x, c(-sample, -age
                          #, -depth
    ), names_to = "topic", values_to = "gamma")
    x
  })
  
  
  topic_plot <- lapply(topic_dfs, function(x){
    ggplot(x, aes(x = age, y = gamma, col = topic)) + geom_line() +
      scale_x_reverse()+
      theme_minimal()
  })
  list(terms_dfs, term_plot, topic_dfs, topic_plot)
}

ctm_termtopic <- termtopic_func(ctm_mods)

ctm_termtopic[[4]][[2]]
ctm_termtopic[[2]][[2]]

# Save all graphs for sup figs
# save term graphs
for (i in 1:length(nlist)) {
  n <- i + (min(nlist)-1)
  p <- ctm_termtopic[[2]][[i]]
  jpeg(paste0("./figures/topic_plots/terms_N", n, ".jpeg"))
  print(p)
  dev.off()
}
# save topic graphs
for (i in 1:length(nlist)) {
  n <- i + (min(nlist)-1)
  p <- ctm_termtopic[[4]][[i]]
  jpeg(paste0("./figures/topic_plots/topics_N", n, ".jpeg"))
  print(p)
  dev.off()
}

rm(p, n, i, nlist, reptimes, ctrls, aicsctm, bicsctm, myctrl, poll_only, poll_reps)
#Will use 3 topics moving forward with main analysis
ctm_topics3 <- ctm_termtopic[[3]][[2]]
ctm_topics3 <- ctm_topics3 %>%
  pivot_wider(names_from = topic, values_from = gamma)
ctm_terms3 <- ctm_termtopic[[1]][[2]]

# BCP Analysis -----------
# Bayesian change point detection for story lake topic analysis
ctm_topics3_bcp <- apply(ctm_topics3[,3:5], MARGIN = 2, bcp, burnin = 100, mcmc = 1000, p0 = 0.12)

ctm_topics3_postprob <- data.frame(age = story_pollen_grouped_counts_wide$ages, Topic_1 = ctm_topics3_bcp[[1]]$posterior.prob, 
                                        Topic_2 = ctm_topics3_bcp[[2]]$posterior.prob,
                                        Topic_3 = ctm_topics3_bcp[[3]]$posterior.prob)

write.csv(ctm_topics3_postprob, file = "./output/final/bcp3_posterior.csv")

# Charcoal Analysis ---------
char_raw <- STOR19_data[["CHAR"]]

char_depths <- seq(from = 1.5, to = max(char_raw$CmTop)+0.5, by = 1)
#convert sedimentation rate into sediment accumulation rate
story_accrate <- data.frame(age = story_sedrate$position_grid, acc_rate = 1/story_sedrate$`50%`)
#compile charcoal data and calculate the charcoal accumulation rate (CHAR)
story_char <- data.frame(depth = char_depths, char_count = char_raw$Total_Count[2:length(char_raw$Total_Count)])
story_char$CHAR <- story_accrate$acc_rate[1:length(story_char$char_count)] * story_char$char_count
#Add ages
story_char <- left_join(story_char, story_ages_midpoint, by = c("depth" = "depth"))

#Put fagus pollen portions into 50 year bins
story_fagus$bin50 <- cut(story_fagus$ages, breaks = seq(from = -100, to = 8000, by = 50), labels = seq(from = -50, to = 8000, by = 50))
story_fagus_binned50 <- story_fagus %>%
  group_by(bin50) %>%
  summarise(mean(fagus))
colnames(story_fagus_binned50) <- c("bin50", "fagus")
#Put quercus pollen portions into 50 year bins
story_quercus$bin50 <- cut(story_quercus$ages, breaks = seq(from = -100, to = 8000, by = 50), labels = seq(from = -50, to = 8000, by = 50))
story_quercus_binned50 <- story_quercus %>%
  group_by(bin50) %>%
  summarise(mean(quercus))
colnames(story_quercus_binned50) <-  c("bin50", "quercus")
#Put CHAR into 50 year bins
story_char$bin50 <- cut(story_char$ages, breaks = seq(from = -100, to = 8000, by = 50), labels = seq(from = -50, to = 8000, by = 50))
story_char_binned <- story_char %>%
  group_by(bin50) %>%
  summarise(mean(CHAR, na.rm = TRUE))
colnames(story_char_binned) <- c("bin50", "char")

char_binned <- left_join(story_fagus_binned50, story_quercus_binned50, by = "bin50")
char_binned <- left_join(char_binned, story_char_binned, by = "bin50")

cor(char_binned$fagus, char_binned$char)
summary(lm(char_binned$fagus ~ char_binned$char))

cor(char_binned$quercus, char_binned$char)
summary(lm(char_binned$quercus ~ char_binned$char))

# Isotopic Analysis -------
story_iso_all <- STOR19_data[["ISO"]]
story_iso <- story_iso_all[which(story_iso_all$code != 3),1:2]
story_fagus_iso <- inner_join(story_iso, story_fagus[,1:3], by = c("depth" = "depth"))
colnames(story_fagus_iso) <- c("depth", "d13c", "age", "fagus")

cor(story_fagus_iso$fagus, story_fagus_iso$d13c)
summary(lm(story_fagus_iso$d13c ~ story_fagus_iso$fagus))

# Add average fagus to iso points without direct beech comparison
story_iso_all <- inner_join(story_iso_all, story_ages_midpoint, by = c("depth" = "depth"))
story_fagus_iso_full <- left_join(story_iso_all, story_fagus[,c(1,3)], by = c("depth" = "depth"))
f216 <- mean(c(story_fagus$fagus[which(story_fagus$depth == 212.5 | story_fagus$depth == 220.5)]))
f336 <- mean(c(story_fagus$fagus[which(story_fagus$depth == 332.5 | story_fagus$depth == 340.5)]))
f344 <- mean(c(story_fagus$fagus[which(story_fagus$depth == 340.5 | story_fagus$depth == 356.5)]))
f440 <- mean(c(story_fagus$fagus[which(story_fagus$depth == 436.5 | story_fagus$depth == 444.5)]))
f448 <- mean(c(story_fagus$fagus[which(story_fagus$depth == 444.5 | story_fagus$depth == 456.5)]))
f476 <- mean(c(story_fagus$fagus[which(story_fagus$depth == 473.5 | story_fagus$depth == 480.5)]))

story_fagus_iso_full$fagus[which(story_fagus_iso_full$depth == 216.5)] <- f216
story_fagus_iso_full$fagus[which(story_fagus_iso_full$depth == 336.5)] <- f336
story_fagus_iso_full$fagus[which(story_fagus_iso_full$depth == 344.5)] <- f344
story_fagus_iso_full$fagus[which(story_fagus_iso_full$depth == 440.5)] <- f440
story_fagus_iso_full$fagus[which(story_fagus_iso_full$depth == 448.5)] <- f448
story_fagus_iso_full$fagus[which(story_fagus_iso_full$depth == 476.5)] <- f476

summary(lm(story_fagus_iso_full$fagus~story_fagus_iso_full$d13c))
cor(story_fagus_iso_full$d13c, story_fagus_iso_full$fagus)

# Compare different time periods
story_iso_all[which(story_iso_all$ages < 3000),]

cor(story_fagus_iso[which(story_fagus_iso$age < 3000),]$fagus, story_fagus_iso[which(story_fagus_iso$age < 3000),]$d13c)
plot(story_fagus_iso[which(story_fagus_iso$age < 3000),]$fagus, story_fagus_iso[which(story_fagus_iso$age < 3000),]$d13c)
abline(lm(story_fagus_iso[which(story_fagus_iso$age < 3000),]$d13c ~ story_fagus_iso[which(story_fagus_iso$age < 3000),]$fagus))
summary(lm(story_fagus_iso[which(story_fagus_iso$age < 3000),]$d13c ~ story_fagus_iso[which(story_fagus_iso$age < 3000),]$fagus))

cor(story_fagus_iso[which(story_fagus_iso$age > 3000),]$fagus, story_fagus_iso[which(story_fagus_iso$age > 3000),]$d13c)
plot(story_fagus_iso[which(story_fagus_iso$age > 3000),]$fagus, story_fagus_iso[which(story_fagus_iso$age > 3000),]$d13c)
abline(lm(story_fagus_iso[which(story_fagus_iso$age > 3000),]$d13c ~ story_fagus_iso[which(story_fagus_iso$age > 3000),]$fagus))
summary(lm(story_fagus_iso[which(story_fagus_iso$age > 3000),]$d13c ~ story_fagus_iso[which(story_fagus_iso$age > 3000),]$fagus))


# Comparison to regional lakes ---------

## Spicer Lake data wrangling ----------
# Download Spicer data from Neotoma
spicer_site <- neotoma2::get_sites(sitename="%Spicer%")
spicer_dataset <- neotoma2::get_datasets(spicer_site, all_data = TRUE)
spicer_download <- neotoma2::get_downloads(spicer_dataset, all_data = TRUE)
spicer_samples <- neotoma2::samples(spicer_download)

# Recalculate age-depth model
spicer_ams <- chroncontrols(spicer_download)
spicer_ams <- arrange(spicer_ams, depth)

spicer_chronology <- Bchron::Bchronology(ages = spicer_ams$chroncontrolage,
                                     ageSds = abs(spicer_ams$agelimityounger - 
                                                    spicer_ams$chroncontrolage),
                                     calCurves = c("normal", rep("intcal20", (nrow(spicer_ams)-1))),
                                     positionThicknesses = spicer_ams$thickness,
                                     positions = spicer_ams$depth,
                                     allowOutside = TRUE,
                                     ids = spicer_ams$chroncontrolid)
plot(spicer_chronology)
jpeg("./figures/regional_lakes/spicer_chronology.jpeg", height = 5, width = 7, units = "in", res = 300)
plot(spicer_chronology) + scale_x_reverse() + ylab("Depth (cm below surface)") + xlab("Age (calendar years BP)") + coord_flip()
dev.off()

spicer_ages_predict <- predict(spicer_chronology, 
                                newPositions = seq(from = 0.5, to = 1100.5, by = 1))
spicer_ages_midpoint <- apply(spicer_ages_predict, 2, mean)
spicer_ages_midpoint <- data.frame(depth =  seq(from = 0.5, to = 1100.5, by = 1), ages = spicer_ages_midpoint)
#calculate confidence interval for age depth model
spicer_ages_95conf <- apply(spicer_ages_predict, 2, function(x){
  quantile(x, probs = c(0.025, 0.975))
})

spicer_ages_95conf_low <- apply(spicer_ages_predict, 2, function(x){
  mean(x)-(2*sd(x))
})
spicer_ages_95conf_high <- apply(spicer_ages_predict, 2, function(x){
  mean(x)+(2*sd(x))
})
spicer_ages_95conf_sd <- data.frame(spicer_ages_midpoint, "lower" = spicer_ages_95conf_low, "upper" = spicer_ages_95conf_high)
spicer_ages_95conf <- t(spicer_ages_95conf)
spicer_ages_combo <- data.frame(spicer_ages_midpoint, "lower" = spicer_ages_95conf[,1], "upper" = spicer_ages_95conf[,2])

# Format Pollen data
colnames(spicer_ages_midpoint) <- c("depth", "new_age")

#add new ages to spicer samples
spicer_samples <- spicer_samples %>%
  left_join(spicer_ages_midpoint, by = c("depth" = "depth"))

#extract terrestrial pollen samples from Neotoma data, convert counts to proportions, and convert to wide format
spicer_pollen <- spicer_samples %>%
  dplyr::filter(elementtype == "pollen", ecologicalgroup != "AQVP") %>% 
  dplyr::select(depth, new_age, variablename, value) %>% 
  group_by(new_age) %>%
  mutate(count = sum(value, na.rm = TRUE)) %>%
  group_by(variablename) %>% 
  mutate(prop = value / count) %>% 
  arrange(desc(new_age)) %>%
  tidyr::pivot_wider(id_cols = c(new_age, depth),
                     names_from = variablename, 
                     values_from = prop,
                     values_fill = 0)
#extract fagus proportion 
spicer_fagus <- spicer_pollen %>%
  dplyr::filter(new_age < max(story_fagus$ages))  %>%
  dplyr::select(new_age, "Fagus grandifolia")

colnames(spicer_fagus) <- c("ages", "fagus")

#Remove neotoma specific elements
rm(spicer_site, spicer_dataset, spicer_download, spicer_ams)

## Pretty data wrangling ---------
# download Pretty Lake data from neotoma
pretty_site <- neotoma2::get_sites(sitename="%pretty%", siteid = 2954)
pretty_dataset <- neotoma2::get_datasets(pretty_site, all_data = TRUE)
pretty_download <- neotoma2::get_downloads(pretty_dataset, all_data = TRUE)
pretty_samples <- neotoma2::samples(pretty_download)

# Recalculate age-depth model
pretty_ams <- chroncontrols(pretty_download)
pretty_ams <- arrange(pretty_ams, depth)
pretty_ams <- dplyr::filter(pretty_ams, chronologyid == 1456) #filter to only include one instance of the age controls

pretty_chronology <- Bchron::Bchronology(ages = pretty_ams$chroncontrolage,
                                         ageSds = abs(pretty_ams$agelimityounger - 
                                                        pretty_ams$chroncontrolage),
                                         calCurves = c("normal", "normal", rep("intcal20", (nrow(pretty_ams)-2))),
                                         #positionThicknesses = pretty_ams$thickness, #does not have thickness
                                         positions = pretty_ams$depth,
                                         allowOutside = TRUE,
                                         ids = pretty_ams$chroncontrolid)
plot(pretty_chronology)
jpeg("./figures/regional_lakes/pretty_chronology.jpeg", height = 5, width = 7, units = "in", res = 300)
plot(pretty_chronology)+ scale_x_reverse() + ylab("Depth (cm below surface)") + xlab("Age (calendar years BP)") + coord_flip()
dev.off()

pretty_ages_predict <- predict(pretty_chronology, 
                               newPositions = seq(from = 0, to = 520, by = 1))
pretty_ages_midpoint <- apply(pretty_ages_predict, 2, mean)
pretty_ages_midpoint <- data.frame(depth =  seq(from = 0, to = 520, by = 1), ages = pretty_ages_midpoint)

colnames(pretty_ages_midpoint) <- c("depth", "new_age")

#add new ages to pretty samples
pretty_samples <- pretty_samples %>%
  left_join(pretty_ages_midpoint, by = c("depth" = "depth"))

#extract terrestrial pollen samples from Neotoma data, convert counts to proportions, and convert to wide format
pretty_pollen <- pretty_samples %>%
  dplyr::filter(elementtype == "pollen", ecologicalgroup != "AQVP") %>% 
  dplyr::select(depth, new_age, variablename, value) %>% 
  group_by(new_age) %>%
  mutate(count = sum(value, na.rm = TRUE)) %>%
  group_by(variablename) %>% 
  mutate(prop = value / count) %>% 
  arrange(desc(new_age)) %>%
  tidyr::pivot_wider(id_cols = c(new_age, depth),
                     names_from = variablename, 
                     values_from = prop,
                     values_fill = 0)
#extract fagus proportion 
pretty_fagus <- pretty_pollen %>%
  dplyr::filter(new_age < max(story_fagus$ages))  %>%
  dplyr::select(new_age, "Fagus")

colnames(pretty_fagus) <- c("ages", "fagus")

#Remove neotoma specific elements
rm(pretty_site, pretty_dataset, pretty_download)


## Appleman data wrangling ----------
# Appleman age depth model 
# age depth model for the top of core of appleman lake using core top, ambrosia rise and the bottom two dates from Jensen 2019
apple_ams_ages <- c(-55, 90, 7310, 7600)
apple_ams_error <- c(10, 20, 50, 40)
apple_ams_depth <- c(0.5, 40.5, 618.5, 627.5)
apple_ams_calcurve <- c(rep("normal", 2), rep("intcal20", 2))

apple_cal_ages <- BchronCalibrate(ages=apple_ams_ages,
                                  ageSds=apple_ams_error,
                                  calCurves=apple_ams_calcurve)
apply(sampleAges(apple_cal_ages), 2, quantile, prob=c(0.025,0.5, 0.975))

apple_bchron <- Bchronology(ages = apple_ams_ages, 
                            ageSds =  apple_ams_error, 
                            positions = apple_ams_depth, 
                            calCurves = apple_ams_calcurve)
# plot(apple_chronology)

#plot chronology
jpeg("./figures/regional_lakes/apple_chronology.jpeg", height = 5, width = 7, units = "in", res = 300)
plot(apple_bchron)+ scale_x_reverse() + ylab("Depth (cm below surface)") + xlab("Age (calendar years BP)") + coord_flip()
dev.off()

#predict ages for sample centers (0.5 cm)
apple_predict_ages_midpoint <- predict(apple_bchron, 
                                       newPositions = seq(from = 0.5, to = 10000.5, by = 1))
apple_ages_midpoint <- apply(apple_predict_ages_midpoint, 2, mean)
apple_ages_midpoint <- data.frame(depth =  seq(from = 0.5, to = 10000.5, by = 1), ages = apple_ages_midpoint)

# Appleman pollen
ecol_groups <- read.csv("./data_raw/ecolgroup_taxa.csv")
apple_pollen_raw <- read.csv("./data_raw/APPL05_Pollen.csv")
colnames(apple_pollen_raw) <- gsub("[.]", " ", colnames(apple_pollen_raw))
apple_pollen_raw[is.na(apple_pollen_raw)] <- 0

#convert to long format, remove taxa that have no counts, remove aquatic taxa
apple_pollen_long <- tibble(apple_pollen_raw[, which(colSums(apple_pollen_raw) > 0)]) %>%  
  pivot_longer(cols = -Depth, names_to = "variablename") 

# group certain taxa into genus gropus and group all unknowns together
apple_pollen_grouped_counts_long <- apple_pollen_long %>%
  mutate(variablename = replace(variablename, stringr::str_detect(variablename, "Pinus*"), "Pinus"),
         variablename = replace(variablename, stringr::str_detect(variablename, "Picea*"), "Picea"),
         variablename = replace(variablename, stringr::str_detect(variablename, "Acer*"), "Acer"),
         variablename = replace(variablename, stringr::str_detect(variablename, "Juglans*"), "Juglans"),
         variablename = replace(variablename, stringr::str_detect(variablename, "Fraxinus*"), "Fraxinus"),
         variablename = replace(variablename, stringr::str_detect(variablename, "Ambrosia*"), "Ambrosia")
  ) %>%
  group_by(variablename, Depth) %>%
  summarise(value = sum(value), .groups = 'keep') %>%
  ungroup()
#convert pollen counts to proportions and save as wide format and add ages
apple_pollen_grouped_prop_wide <- apple_pollen_grouped_counts_long %>%
  group_by(Depth) %>%
  mutate(count = sum(value, na.rm = TRUE)) %>%
  group_by(variablename) %>% 
  mutate(prop = value / count) %>% 
  arrange(desc(Depth)) %>%
  tidyr::pivot_wider(id_cols = Depth,
                     names_from = variablename, 
                     values_from = prop,
                     values_fill = 0) %>%
  left_join(apple_ages_midpoint, by = c("Depth" = "depth"))

#convert grouped pollen counts to wide format and add ages
apple_pollen_grouped_counts_wide <- apple_pollen_grouped_counts_long %>%
  pivot_wider(id_cols = Depth, names_from = variablename, values_from = value, values_fill = 0) %>%
  left_join(apple_ages_midpoint, by = c("Depth" = "depth"))

#extract fagus and quercus pollen with depth and age
apple_fagus <- dplyr::select(apple_pollen_grouped_prop_wide, Depth, ages, "Fagus")
colnames(apple_fagus) <- c("depth", "ages", "fagus")
apple_quercus <- dplyr::select(apple_pollen_grouped_prop_wide, Depth, ages, "Quercus")
colnames(apple_quercus) <- c("depth", "ages", "quercus")



## Binned Comparison -----
#bin Story and Spicer data into 250 year bins (minimum needed for every bin to have a value for both sites)
spicer_fagus$bin250 <- cut(spicer_fagus$ages, breaks = seq(from = -250, to = 8000, by = 250), labels = seq(from = 0, to = 8000, by = 250))
story_fagus$bin250 <- cut(story_fagus$ages, breaks = seq(from = -250, to = 8000, by = 250), labels = seq(from = 0, to = 8000, by = 250))
apple_fagus$bin250 <- cut(apple_fagus$ages, breaks = seq(from = -250, to = 8000, by = 250), labels = seq(from = 0, to = 8000, by = 250))
pretty_fagus$bin250 <- cut(pretty_fagus$ages, breaks = seq(from = -250, to = 8000, by = 250), labels = seq(from = 0, to = 8000, by = 250))

story_fagus_binned250 <- story_fagus %>%
  group_by(bin250) %>%
  summarise(mean(fagus))
colnames(story_fagus_binned250) <- c("bin250", "story")

spicer_fagus_binned250 <- spicer_fagus %>%
  group_by(bin250) %>%
  summarise(mean(fagus))
colnames(spicer_fagus_binned250) <- c("bin250", "spicer")

apple_fagus_binned250 <- apple_fagus %>%
  group_by(bin250) %>%
  summarise(mean(fagus))
colnames(apple_fagus_binned250) <- c("bin250", "apple")

pretty_fagus_binned250 <- pretty_fagus %>%
  group_by(bin250) %>%
  summarise(mean(fagus))
colnames(pretty_fagus_binned250) <- c("bin250", "pretty")

fagus_combo <- story_fagus_binned250 %>%
  left_join(spicer_fagus_binned250, by = "bin250") %>%
  left_join(apple_fagus_binned250, by = "bin250") %>%
  left_join(pretty_fagus_binned250, by = "bin250")

cor(fagus_combo$story,fagus_combo$spicer)
summary(lm(fagus_combo$story ~ fagus_combo$spicer))
summary(lm(fagus_combo$story ~ fagus_combo$apple))
summary(lm(fagus_combo$story ~ fagus_combo$pretty))
cor(fagus_combo$story,fagus_combo$apple)
cor(fagus_combo$story[which(fagus_combo$pretty > 0)],fagus_combo$pretty[which(fagus_combo$pretty > 0)])

## Age model comparison -----------
spicer_pollen_age_conf <- left_join(spicer_fagus[,1:2], spicer_ages_combo, by = c("ages" = "new_age"))
ss_age_comp <- data.frame(story_ages_predict[,which(round(story_ages_combo$ages) %in% c(2784, 4615, 7061))], spicer_ages_predict[,which(round(spicer_ages_combo$ages) %in% c(1778, 3216, 6782))])
colnames(ss_age_comp) <- c("story_28", "story_46", "story_70", "spicer_18", "spicer_32", "spicer_68")

t.test(ss_age_comp$story_28, ss_age_comp$spicer_18)
ggplot(ss_age_comp) + geom_histogram(aes(x = story_28), fill = "blue", alpha = 0.5) + geom_histogram(aes(x = spicer_18), fill = "red", alpha = 0.5)+ theme_minimal()
ggsave("./figures/regional_lakes/ss_recovery.jpeg")
t.test(ss_age_comp$story_46, ss_age_comp$spicer_32)
ggplot(ss_age_comp) + geom_histogram(aes(x = story_46), fill = "blue", alpha = 0.5) + geom_histogram(aes(x = spicer_32), fill = "red", alpha = 0.5)+ theme_minimal()
ggsave("./figures/regional_lakes/ss_onset.jpeg")
t.test(ss_age_comp$story_70, ss_age_comp$spicer_68)
ggplot(ss_age_comp) + geom_histogram(aes(x = story_70), fill = "blue", alpha = 0.5) + geom_histogram(aes(x = spicer_68), fill = "red", alpha = 0.5)+ theme_minimal()
ggsave("./figures/regional_lakes/ss_arrival.jpeg")

# Save key data -------
save(story_pollen_grouped_counts_long, file = "./output/final/story_pollen_grouped_counts_long.RData")
save(story_pollen_grouped_counts_wide, file = "./output/final/story_pollen_grouped_counts_wide.RData")



