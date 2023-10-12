#Compilation of Publication Figures for Story Lake Indiana
#Nora Schlenker
#October 2023

library(tidyverse)
library(tidytext)
library(ggtext)
library(rioja)
library(riojaPlot)
library(Bchron)
library(vegan) 
library(gridExtra)
library(gtable)
library(grid)
library(quantreg)

#First run analysis code to get all data
#source("./STOR19_Analysis.R")

# Figure 1: Site Map --------
# Not made in R

# Figure 2: Age-depth model --------
jpeg("./figures/story_chronology.jpeg", height = 5, width = 7, units = "in", res = 300)
chron_plot <- plot(story_chronology)
chron_plot + scale_x_reverse(breaks = seq(0, 14000, 1000), limits = c(13500, -500)) + ylab("Depth (cm below surface)") + xlab("Age (calendar years BP)") + coord_flip()
dev.off()

# Figure 3: Pollen Diagram ----------
# Convert proportion to percent
story_perc <- dplyr::select(story_pollen_grouped_prop_wide, -Depth, -ages) #remove depth and age
story_perc <- story_perc*100 #convert to percent

# filter pollen that have greater than 5% abundance
mx <- apply(story_perc, 2, max, na.rm = TRUE)
poll5 <- story_perc[,mx>3.5]

# Define pollen taxa to include in plot starting with looking at the names of taxa with max>5%
colnames(poll5)
poll_incl <- c("Fagus grandifolia", "Quercus", "Ulmus", "Carya", "Ostrya/Carpinus", "Fraxinus", 
               "Juglans", "Platanus", "Salix", "Acer",  "Cupressaceae/Taxaceae",
               "Ambrosia sp.", "Poaceae", "Cyperaceae")
poll_types <- data.frame(names = poll_incl, group = c(rep("Broad-leaf",10), #broad-leaf = dark grey
                                                      rep("Needle-leaf", 1), #needle-leaf = middle grey
                                                      rep("Herbs", 3))) #herbs = light grey)
poll_types$group <- factor(poll_types$group, levels = c("Broad-leaf", "Needle-leaf", "Herbs"))

taxa_labels <- c(expression(italic("Fagus grandifolia")), expression(italic("Quercus")), expression(italic("Ulmus")), expression(italic("Carya")), 
                 expression(italic("Ostrya/Carpinus")), expression(italic("Fraxinus")), expression(italic("Juglans")),
                 expression(italic("Platanus")), expression(italic("Salix")), expression(italic("Acer")), "Cupressaceae/Taxaceae", 
                 expression(italic("Ambrosia")), "Poaceae", "Cyperaceae")

x_poll <- select(story_perc, any_of(poll_incl))

cols <- c("#3253D9", #fagus = blue
                   "#C75702", #quercus = orange
                   rep("#292929",8), #broad-leaf = dark grey
                   rep("#707070", 1), #needle-leaf = middle grey
                   rep("#b8b6b6", 3)) #herbs = light grey

## extract ages
chron <- data.frame(story_pollen_grouped_prop_wide$ages)
colnames(chron) <- c("Age (calendar years BP)")

jpeg("./figures/full_poll.jpeg", width = 12, height = 6, units = "in", res = 600)
poll_rp <- riojaPlot(x_poll, chron, selVars = poll_incl, groups = poll_types
                     , yvar.name = "Age (calendar years BP)"
                     , ymin = 0, ymax = 8000, yinterval = 1000
                     #, plot.groups = TRUE
                     , scale.percent = TRUE
                     #, plot.cumul = TRUE #plot cumulative plot for all types
                     , srt.xlabel = 45
                     #, labels.italicize = TRUE #this is not working....
                     #, do.clust = TRUE
                     #, plot.clust = TRUE
                     #, plot.zones = "auto"
                     , col.poly=cols
                     , plot.exag=TRUE
                     #, xlabels = y.names
                     #, xRight = 0.8
                     , xlabels = taxa_labels
)
dev.off()

# Figure 4: Topic Analysis --------
## Terms plot -------
#extract top 5 taxa for each topic
ctm_terms3_top <- ctm_terms3 %>%
  group_by(ngroups, topic) %>%
  top_n(5, beta) %>%
  ungroup()
#Plot top 5 taxa for each topic
ctm_terms3_plot <-  ctm_terms3_top %>%
   mutate(
##latin names
  #        taxa = replace(taxa, stringr::str_detect(taxa, "Ambrosia*"), "*Ambrosia*"),
  #        taxa = replace(taxa, stringr::str_detect(taxa, "Quercus*"), "<i>Quercus</i>"),
  #        taxa = replace(taxa, stringr::str_detect(taxa, "Salix*"), "<i>Salix</i>"),
  #        taxa = replace(taxa, stringr::str_detect(taxa, "Fagus*"), "<i>Fagus grandifolia</i>"),
  #        taxa = replace(taxa, stringr::str_detect(taxa, "Fraxinus*"), "<i>Fraxinus</i>"),
  #        taxa = replace(taxa, stringr::str_detect(taxa, "Juglans*"), "<i>Juglans</i>"),
  #        taxa = replace(taxa, stringr::str_detect(taxa, "Ulmus*"), "<i>Ulmus</i>"),
  #        taxa = replace(taxa, stringr::str_detect(taxa, "Platanus*"), "<i>Platanus</i>"),
  #        taxa = replace(taxa, stringr::str_detect(taxa, "Ostrya*"), "<i>Ostrya/Carpinus</i>"),
  #        taxa = replace(taxa, stringr::str_detect(taxa, "Carya*"), "<i>Carya</i>"),
  #        taxa = replace(taxa, stringr::str_detect(taxa, "Cupress*"), "Cupressaceae/<br>Taxaceae"),
  #        taxa = replace(taxa, stringr::str_detect(taxa, "Pinus*"), "<i>Pinus</i>"),
## common names
          taxa = replace(taxa, stringr::str_detect(taxa, "Ambrosia*"), "Ragweed"),
          taxa = replace(taxa, stringr::str_detect(taxa, "Quercus*"), "Oak"),
          taxa = replace(taxa, stringr::str_detect(taxa, "Salix*"), "Willow"),
          taxa = replace(taxa, stringr::str_detect(taxa, "Fagus*"), "American<br>beech"),
          taxa = replace(taxa, stringr::str_detect(taxa, "Fraxinus*"), "Ash"),
          taxa = replace(taxa, stringr::str_detect(taxa, "Juglans*"), "Walnut"),
          taxa = replace(taxa, stringr::str_detect(taxa, "Ulmus*"), "Elm"),
          taxa = replace(taxa, stringr::str_detect(taxa, "Platanus*"), "Sycamore"),
          taxa = replace(taxa, stringr::str_detect(taxa, "Ostrya*"), "Hornbeam/<br>hop-hornbeam"),
          taxa = replace(taxa, stringr::str_detect(taxa, "Carya*"), "Hickory"),
          taxa = replace(taxa, stringr::str_detect(taxa, "Cupress*"), "Cypress/<br>Yew"),
          taxa = replace(taxa, stringr::str_detect(taxa, "Pinus*"), "Pine"),
          taxa = replace(taxa, stringr::str_detect(taxa, "Poaceae*"), "Grass"),
         taxa = reorder_within(taxa, beta, topic)) %>%
  ggplot(aes(x = taxa, beta, fill = topic)) + geom_col(show.legend = FALSE) +
  scale_fill_discrete(type = c("#C75702", "#3253D9",  "#FFBF00")) +
  facet_wrap(~factor(topic, levels = c(2,1,3), labels = c("Beech Hardwood Forest", "Oak Woodland/Forest", "Open or Cleared Vegetation"))
             , scales = "free", drop = TRUE, ncol = 1) +
  coord_flip() +
  theme_minimal()+
  scale_x_reordered() +
  scale_y_continuous(expand = c(0,0))+
  xlab(NULL) +
  ylab("Importance (beta)")+
  theme(axis.text.y = element_markdown(size = 6),
        axis.text.x = element_text(size = 6),
        axis.title = element_text(size = 8),
        strip.text = element_text(size = 6))
ctm_terms3_plot
ggsave("./figures/ctm_terms3.jpeg", ctm_terms3_plot, width = 2, height = 5, units = "in", dpi = 600)

## Topics plot ------
ctm_topics3_plot <- ctm_topics3 %>%
  pivot_longer(c( -sample, -age), names_to = "topic", values_to = "gamma") %>%
  ggplot(aes(x = age, y = gamma, col = topic)) + geom_line(linewidth = 1, show.legend = FALSE) +
  scale_x_reverse()+
  scale_color_discrete(type = c("3" = "#FFBF00","2" = "#3253D9", "1" = "#C75702")) +
  xlab("Age (ka BP)") +
  ylab("Community Prevalence (gamma)") +
  theme_minimal()+
  theme(axis.text.y = element_markdown(size = 6),
        axis.text.x = element_text(size = 6),
        axis.title = element_text(size = 10),
        strip.text = element_text(size = 6))
ctm_topics3_plot
ggsave("./figures/ctm_topics3.jpeg", ctm_topics3_plot, width = 7, height = 3, units = "in", dpi = 600)

## BCP plot ----------
ctm_bcp1 <- data.frame(bcp_prob = ctm_topics3_bcp[[1]]$posterior.prob, topic = rep(1, times = length(ctm_topics3_bcp[[1]]$posterior.prob)), age = ctm_topics3$age)
ctm_bcp2 <- data.frame(bcp_prob = ctm_topics3_bcp[[2]]$posterior.prob, topic = rep(2, times = length(ctm_topics3_bcp[[2]]$posterior.prob)), age = ctm_topics3$age)
ctm_bcp3 <- data.frame(bcp_prob = ctm_topics3_bcp[[3]]$posterior.prob, topic = rep(3, times = length(ctm_topics3_bcp[[3]]$posterior.prob)), age = ctm_topics3$age)

ctm_bcp <- rbind(ctm_bcp1, ctm_bcp2, ctm_bcp3)
ctm_bcp$topic <- factor(ctm_bcp$topic, levels = c("1", "2", "3"))
ctm_bcp_plot <- ggplot(ctm_bcp, aes(age, bcp_prob, col = topic)) + geom_line(linetype = "dashed") +
  scale_x_reverse()+
  scale_color_discrete(type = c( "3" = "#FFBF00","2" = "#3253D9", "1" = "#C75702")) +
  theme_minimal()+
  ylab("BCP Posterior Probability") +
  xlab("Age (ka BP)") +
  geom_hline(yintercept = 0.5) +
  theme(axis.text.y = element_markdown(size = 6),
        axis.text.x = element_text(size = 6),
        axis.title = element_text(size = 8), legend.position = "none")
ctm_bcp_plot
ggsave("./figures/ctm4_bcp.jpeg", ctm_bcp_plot, width = 7, height = 2, units = "in", dpi = 600)

# Figure 5: Isotopes ----
iso_cor_plot <- ggplot(story_fagus_iso_full, aes(d13c, fagus)) +
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) + 
  ylim(0, 0.3)+
  xlim(-30, -20)+
  #annotate("text", label = paste("Pearson's Correlation = ", round(cor(story_fagus_iso$d13c, story_fagus_iso$fagus),2)), x = -27, y = 0.19, size = 2, alpha = 0)+
  theme_minimal() +
  theme(axis.text.y=element_text(size=10), axis.title = element_blank(), panel.grid = element_blank(),
        axis.line = element_line(colour = "gray"), axis.ticks = element_line(colour = "gray"))
iso_cor_plot
ggsave("./figures/iso_fagus_cor_plot.jpeg", iso_cor_plot, width = 3, height = 3, units = "in", dpi = 600)

iso_paired_plot <- ggplot() +
  geom_area(data = story_fagus, aes(ages, fagus), col ="#3253D9", fill ="#3253D9") +
  geom_point(data = story_iso_all, aes(ages, ((d13c/-100)*2)-0.3), size = 0.5) +
  geom_line(data = story_iso_all, mapping = aes(ages, ((d13c/-100)*2)-0.3)) +
  scale_y_continuous(name = "Fagus",
                     sec.axis = sec_axis(trans= ~. + 0.3, name="Fagus d13c",  breaks = c(0.6, 0.55, 0.5, 0.45, 0.4), labels = c("-30.0", "-27.5", "-25.0", "-22.5", "-20.0"))) +
  scale_x_reverse()+
  theme_minimal()+
  theme(axis.text.y=element_text(size=10), axis.title = element_blank(), panel.grid = element_blank(),
        axis.line = element_line(colour = "gray"), axis.ticks = element_line(colour = "gray"))
iso_paired_plot
ggsave("./figures/iso_fagus_plot.jpeg", iso_paired_plot, width = 7, height = 3, units = "in", dpi = 600)

summary(lm(story_fagus_iso$fagus~story_fagus_iso$d13c))

# Figure 6: Charcoal Correlation ----
rq_fagus <- rq(formula = char ~ fagus, data = char_binned, tau = 0.95)
summary(rq_fagus, se="ker")
char_fagus_cor <- ggplot(char_binned, aes(fagus, char))+ geom_point() + 
  theme_minimal() + 
  geom_smooth(method = "lm", se = FALSE, col = "#3253D9") +
  geom_abline(slope=rq_fagus$coefficients[2],
              intercept=rq_fagus$coefficients[1],
              color="#3253D9", linetype = "dashed")+
  xlab("American beech abundance <br/> (pollen percent)") +
  ylab("Charcoal Accumultation Rate <br/> (pieces cm<sup>-2</sup> year<sup>-1</sup>)") +
  theme(axis.title.x = element_markdown(),
        axis.title.y = element_markdown(), panel.grid = element_blank(),
        axis.line = element_line(colour = "gray"), axis.ticks = element_line(colour = "gray"))
char_fagus_cor
ggsave("./figures/char_fagus_cor.jpeg", char_fagus_cor, width = 3.5, height = 3.5, units = "in", dpi = 600)

summary(lm(char_binned$char ~ char_binned$fagus))

rq_quercus <- rq(formula = char ~ quercus, data = char_binned, tau = 0.95)
summary(rq_quercus, se="ker")
char_quercus_cor <- ggplot(char_binned, aes(quercus, char))+ geom_point() + 
  theme_minimal() + 
  geom_smooth(method = "lm", se = FALSE, col = "#C75702") +
  geom_abline(slope=rq_quercus$coefficients[2],
              intercept=rq_quercus$coefficients[1],
              color="#C75702", linetype = "dashed")+
  xlab("Oak abundance <br/> (pollen percent)") +
  ylab("Charcoal Accumultation Rate <br/> (pieces cm<sup>-2</sup> year<sup>-1</sup>)") +
  theme(axis.title.x = element_markdown(),
        axis.title.y = element_markdown(), panel.grid = element_blank(),
        axis.line = element_line(colour = "gray"), axis.ticks = element_line(colour = "gray"))
char_quercus_cor
ggsave("./figures/char_quercus_cor.jpeg", char_quercus_cor, width = 3.5, height = 3.5, units = "in", dpi = 600)

summary(lm(char_binned$char ~ char_binned$quercus))

# Figure 7: Proxy combo plot ----
bcp_points <- c(191, 2784, 4615, 6155, 7061)
ctm_plot <- ctm_topics3 %>%
  pivot_longer(c( -sample, -age), names_to = "topic", values_to = "gamma") %>%
  ggplot(aes(age, gamma, col = topic, alpha = topic))+
  geom_line(linewidth = 1, show.legend = FALSE) +
  scale_color_discrete(type = c("3" = "white","2" = "#3253D9", "1" = "#C75702")) +
  scale_alpha_manual(values = c(1,1,0)) +
  geom_vline(xintercept = bcp_points, linetype = "longdash")+
  scale_x_reverse()+
  xlim(8000, 0)+
  theme_minimal()+
  theme(axis.text.y=element_text(size=10), axis.text.x=element_blank(), axis.title = element_blank(), 
        panel.grid = element_blank(), axis.line.y = element_line(colour = "lightgrey"))
ctm_plot

char_plot <- ggplot(story_char, aes(ages, CHAR))+
  geom_vline(xintercept = bcp_points, linetype = "longdash")+
  geom_col(width = 40, fill = "black")+
  scale_x_reverse()+
  xlim(8000, 0)+
  theme_minimal()+
  theme(axis.text.y=element_text(size=10), axis.text.x=element_blank(), axis.title = element_blank(), 
        panel.grid = element_blank(), axis.line.y = element_line(colour = "lightgrey"))
char_plot

iso_plot <- ggplot(story_iso_all, aes(ages, d13c))+
  geom_vline(xintercept = bcp_points, linetype = "longdash")+
  geom_point() +
  geom_line()+
  #geom_smooth(method = "loess", se = FALSE, linetype = "dotted", color = "darkgrey", span = .5)+
  scale_x_reverse()+
  xlim(8000, 0)+
  scale_y_reverse()+
  theme_minimal()+
  theme(axis.text.y=element_text(size=10), axis.text.x=element_blank(), axis.title = element_blank(), 
        panel.grid = element_blank(), axis.line.y = element_line(colour = "lightgrey"))
iso_plot

ll_data <- STOR19_data[["Lavine_LL"]]
ll_plot <- ggplot(ll_data, aes(Age, LakeElev.cm.below.mod))+
  geom_vline(xintercept = bcp_points, linetype = "longdash")+
  geom_line()+
  geom_ribbon(aes(ymin= LowerBound, ymax=UpperBound), fill = "grey", alpha = 0.3) +
  scale_x_reverse(breaks = seq(8000, 0, by = -200), limits = c(8000, 0),
                  labels = c(8000, rep("",4), 7000, rep("",4), 6000, rep("",4), 5000, rep("",4), 4000, rep("",4), 3000, rep("",4), 2000, rep("",4), 1000, rep("",4), 0))+
  scale_y_reverse()+
  #xlim(8000, 0)+
  theme_minimal()+
  theme(axis.text.y=element_text(size=10),  axis.title = element_blank(), 
        panel.grid = element_blank(), axis.line = element_line(colour = "lightgrey"), axis.ticks = element_line(colour = "lightgrey"))
ll_plot

p1 <- ggplotGrob(ctm_plot)
p2 <- ggplotGrob(char_plot)
p3 <- ggplotGrob(iso_plot)
p4 <- ggplotGrob(ll_plot) 

g4 <- rbind(p1,p2,p3,p4, size = "first")
g4$widths <- unit.pmax(p1$widths, p2$widths, p3$widths, p4$widths)
grid.newpage()

jpeg("./figures/combo_plot_reorg.jpeg", width = 8, height = 6, units = "in", res = 600)
grid.draw(g4)
dev.off()

rm(p1, p2, p3, p4, g4)

#Figure 8: Story/Spicer comparison ----
st <- story_fagus[,2:3]
sp <- spicer_fagus[,1:2]
ap <- apple_fagus[,2:3]
py <- pretty_fagus[,1:2]
st <- left_join(st, story_ages_95conf_sd[,2:4], by = c("ages" = "ages"))
sp <- left_join(sp, spicer_ages_95conf_sd[,2:4], by = c("ages" = "ages"))
ap$lower <- NA
ap$upper <- NA
py$lower <- NA
py$upper <- NA
st$lake <- "story"
sp$lake <- "spicer"
ap$lake <- "apple"
py$lake <- "pretty"
fagus_plot <- data.frame(rbind(st, sp, ap, py))
fagus_plot$lake <- factor(fagus_plot$lake, levels = c("story", "apple", "pretty", "spicer"), labels = c("Story Lake", "Appleman Lake", "Pretty Lake", "Spicer Lake"))
spicer_bcp_points <- 
rm(st, sp, ap, py)

comp_pollen <- 
  fagus_plot %>%
  mutate(lower = ifelse(round(ages) %in% c(2784, 4615, 7061, 6782, 3216, 1778), lower, NA)) %>%
  mutate(upper = ifelse(round(ages) %in% c(2784, 4615, 7061, 6782, 3216, 1778), upper, NA)) %>%
  mutate(point = ifelse(round(ages) %in% c(2784, 4615, 7061, 6782, 3216, 1778), ages, NA)) %>%
  ggplot(aes(x = ages, y = fagus)) + geom_line(col = "black") + geom_area(fill = "grey") +
  geom_errorbar(aes(xmin = lower, xmax = upper, y = fagus), col = "darkred") +
  geom_point(aes(x = point, y = fagus), col = "darkred") +
  facet_wrap(~lake, nrow = 1) +
  xlab("Age (Years BP)") +
  ylab("American beech abundance")+
  coord_flip() +
  scale_x_reverse() +
  theme_minimal()+
  xlim(max(story_fagus$ages),min(story_fagus$ages))+
  theme(axis.title = element_text(size = 11),
        axis.title.x = element_markdown(),
        strip.text = element_text(size = 12))
comp_pollen
ggsave("./figures/comp_pollen.jpeg", comp_pollen, width = 8, height = 8, units = "in", dpi = 600)

ss_binned_cor <- ggplot(fagus_combo, aes(story, spicer))+ geom_point(shape = 2) + 
  theme_minimal() + 
  geom_smooth(method = "lm", se = FALSE, col = "black") +
  xlab("Story Lake<br>beech abundance") +
  ylab("Spicer Lake<br>beech abundance") +
  theme(axis.title.x = element_markdown(),
        axis.title.y = element_markdown(),
        title = element_markdown())
ss_binned_cor
ggsave("./figures/ss_correlation.jpeg", width = 4, height = 4, units = "in", dpi = 600)

spa_binned_cor <- fagus_combo %>%
  select(-spicer) %>%
  pivot_longer(cols = c("apple", "pretty")) %>%
  mutate(name_factor = factor(name, levels = c("apple", "pretty"), labels = c("Appleman Lake", "Pretty Lake")))%>%
  ggplot(aes(story, value, shape = name_factor, linetype = name_factor))+ 
  geom_point() + 
  theme_minimal() + 
  geom_smooth(method = "lm", se = FALSE, col = "black") +
  scale_shape_manual(values = c(1,4))+
  scale_linetype_manual(values = c(2,3))+
  ylim(c(0,0.3)) +
  xlab("Story Lake<br>beech abundance") +
  ylab("Appleman and Pretty Lake <br> beech abundance") +
  theme(legend.position = c(.2, .8), 
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white"),
        axis.title.x = element_markdown(),
        axis.title.y = element_markdown(),
        title = element_markdown())
spa_binned_cor
ggsave("./figures/spa_correlation.jpeg", width = 4, height = 4, units = "in", dpi = 600)



#Supplemental Figures -------------
## Pollen accumulation rate -----------
mx <- apply(story_perc, 2, max)
poll5 <- story_perc[,mx>3.5]

# Define pollen taxa to include in plot starting with looking at the names of taxa with max>5%
colnames(poll5)
poll_incl <- c("Fagus grandifolia", "Quercus", "Ulmus", "Carya", "Ostrya/Carpinus", "Fraxinus", 
               "Juglans", "Platanus", "Salix", "Acer",  "Cupressaceae/Taxaceae",
               "Ambrosia sp.", "Poaceae", "Cyperaceae")
poll_types <- data.frame(names = poll_incl, group = c(rep("Broad-leaf",10), #broad-leaf = dark grey
                                                      rep("Needle-leaf", 1), #needle-leaf = middle grey
                                                      rep("Herbs", 3))) #herbs = light grey)
poll_types$group <- factor(poll_types$group, levels = c("Broad-leaf", "Needle-leaf", "Herbs"))

taxa_labels <- c(expression(italic("Fagus grandifolia")), expression(italic("Quercus")), expression(italic("Ulmus")), expression(italic("Carya")), 
                 expression(italic("Ostrya/Carpinus")), expression(italic("Fraxinus")), expression(italic("Juglans")),
                 expression(italic("Platanus")), expression(italic("Salix")), expression(italic("Acer")), "Cupressaceae/Taxaceae", 
                 expression(italic("Ambrosia")), "Poaceae", "Cyperaceae")

cols <- c("#3253D9", #fagus = blue
                   "#C75702", #quercus = orange
                   rep("#292929",8), #broad-leaf = dark grey
                   rep("#707070", 1), #needle-leaf = middle grey
                   rep("#b8b6b6", 3)) #herbs = light grey

## extract ages
chron_a <- data.frame(sort(story_pollen_grouped_prop_wide$ages))
colnames(chron_a) <- c("Age (calendar years BP)")

a_poll <- select(story_pollen_acc, any_of(poll_incl))
jpeg("./figures/poll_accu.jpeg", width = 12, height = 6, units = "in", res = 600)
poll_ac <- riojaPlot(a_poll, chron_a, selVars = poll_incl, groups = poll_types
                     , yvar.name = "Age (calendar years BP)"
                     , ymin = 0, ymax = 8000, yinterval = 1000
                     #, plot.groups = TRUE
                     , scale.percent = TRUE
                     #, plot.cumul = TRUE #plot cumulative plot for all types
                     , srt.xlabel = 45
                     #, labels.italicize = TRUE #this is not working....
                     #, do.clust = TRUE
                     #, plot.clust = TRUE
                     #, plot.zones = "auto"
                     , col.poly=cols
                     , plot.exag=TRUE
                     #, xlabels = y.names
                     #, xRight = 0.8
                     , xlabels = taxa_labels
)
dev.off()

## Appleman Pollen Diagram -----------
apple_perc <- dplyr::select(apple_pollen_grouped_prop_wide, -Depth, -ages) #remove depth and age
apple_perc <- apple_perc*100 #convert to percent

mx_ap <- apply(apple_perc, 2, max, na.rm = TRUE)
poll5_ap <- apple_perc[,mx_ap>3.5]

# Define pollen taxa to include in plot starting with looking at the names of taxa with max>5%
colnames(poll5_ap)
poll_incl_ap <- c("Fagus", "Quercus", "Ulmus", "Carya", "OstryaCarpinus", "Fraxinus", 
               "Juglans", "Platanus", "Acer", "Pinus", "Ambrosia", "Poaceae undiff ")
poll_types_ap <- data.frame(names = poll_incl_ap, group = c(rep("Broad-leaf",9), #broad-leaf = dark grey
                                                      rep("Needle-leaf", 1), #needle-leaf = middle grey
                                                      rep("Herbs", 2))) #herbs = light grey)
poll_types_ap$group <- factor(poll_types_ap$group, levels = c("Broad-leaf", "Needle-leaf", "Herbs"))

taxa_labels_ap <- c(expression(italic("Fagus")), expression(italic("Quercus")), expression(italic("Ulmus")), expression(italic("Carya")), 
                 expression(italic("Ostrya/Carpinus")), expression(italic("Fraxinus")), expression(italic("Juglans")),
                 expression(italic("Platanus")), expression(italic("Acer")), "Pinus", 
                 expression(italic("Ambrosia")), "Poaceae")

ap_poll <- select(apple_perc, any_of(poll_incl_ap))

cols_ap <- c("#3253D9", #fagus = blue
                   "#C75702", #quercus = orange
                   rep("#292929",7), #broad-leaf = dark grey
                   rep("#707070", 1), #needle-leaf = middle grey
                   rep("#b8b6b6", 2)) #herbs = light grey

## extract ages
chron_ap <- data.frame(apple_pollen_grouped_prop_wide$ages)
colnames(chron_ap) <- c("Age (calendar years BP)")

jpeg("./figures/regional_lakes/poll_apple.jpeg", width = 12, height = 6, units = "in", res = 600)
poll_ap <- riojaPlot(ap_poll, chron_ap, selVars = poll_incl_ap, groups = poll_types_ap
                     , yvar.name = "Age (calendar years BP)"
                     , ymin = 0, ymax = 8000, yinterval = 1000
                     #, plot.groups = TRUE
                     , scale.percent = TRUE
                     #, plot.cumul = TRUE #plot cumulative plot for all types
                     , srt.xlabel = 45
                     #, labels.italicize = TRUE #this is not working....
                     #, do.clust = TRUE
                     #, plot.clust = TRUE
                     #, plot.zones = "auto"
                     , col.poly=cols
                     , plot.exag=TRUE
                     #, xlabels = y.names
                     #, xRight = 0.8
                     , xlabels = taxa_labels_ap
)
dev.off()

