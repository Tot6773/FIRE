setwd("C:/Users/thoma/OneDrive/Bureau/fire_analysis")

library(tidyverse)
library(openxlsx)
library(vegan)
library(ade4)
library(gghalves)
library(glmmTMB)
library(sjPlot)
library(emmeans)
library(iNEXT)
library(gridExtra)


data = read.csv("final_data.csv", header = T, sep = ";", dec = ",")
data$qual_freq = factor(data$qual_freq, levels = c("Vierge", "Peu fréquenté", "Très fréquenté"))
data$route = as.factor(data$route)
data <- data[data$orientation != "W", ] # ?????
data = data %>% group_by(route, sp_code) %>% mutate(rel_abu = sum(rel_abu)) %>% 
  group_by(route) %>% distinct(sp_code, rel_abu, .keep_all = T)
data$qual_freq[data$route == "Initiation_4" | data$route == "Initiation_5" | data$route == "Initiation_6"] = "Très fréquenté"

data3 = data %>% group_by(sp_code) %>% mutate(appear = n()) %>%
  filter(!qual_freq == "Peu fréquenté", appear > 1)
data3$qual_freq = factor(data3$qual_freq, levels = c("Vierge", "Très fréquenté"))
data3 = data3 %>% group_by(route) %>% 
  mutate(SR = n_distinct(sp_code))
data3$grade = factor(data3$grade, levels = c("4a","5a","5b","5c","6a","6b","6b+","6c"))

#
#### Structuring datas ####
data = read.csv("data.csv", sep = ";", dec = ",")

sp = read.csv("sp.csv", sep = ";")
sp = sp[,c(1:6)] #be CARREFUL
sp = sp[c(1:60),] #be CARREFUL

routes = read.csv("route.csv", sep = ";", dec = ",")

final_data = left_join(data, sp, by = "species")
final_data = left_join(final_data, routes, by = "route")
final_data = final_data %>%
  mutate(quan_freq = coalesce(quan_freq.x, quan_freq.y)) %>%
  select(-c(quan_freq.x, quan_freq.y))

final_data$qual_freq = ifelse(final_data$quan_freq == 0, "Vierge", final_data$qual_freq)
final_data$qual_freq = ifelse(final_data$quan_freq == 1 | final_data$quan_freq == 2, "Peu fréquenté", final_data$qual_freq)
final_data$qual_freq = ifelse(final_data$quan_freq == 3 | final_data$quan_freq == 4, "Très fréquenté", final_data$qual_freq)
table(final_data$qual_freq, final_data$quan_freq)

final_data$grade = ifelse(final_data$sector == "Initiation", "5a", final_data$grade) #putting grade on initiation

final_data = final_data %>% mutate(rel_abu = round(abundance / size, digits = 2)) %>%
  group_by(route) %>% 
  mutate(SR = n_distinct(sp_code))

final_data$qual_freq = factor(final_data$qual_freq, levels = c("Vierge", "Peu fréquenté", "Très fréquenté"))
final_data$climbing = ifelse(final_data$qual_freq == "Vierge", "No climbing", "Climbing")

write.xlsx(final_data, "final_data.xlsx")
rm(list = ls())

#### Description of datas ####
data %>% group_by(quan_freq) %>% summarise(number_of_routes = n_distinct(route))
data %>% group_by(qual_freq) %>% summarise(number_of_routes = n_distinct(route)) #better use qualitative

data %>% group_by(observer) %>% summarise(number_of_routes = n_distinct(route)) # bravo gab

data %>% group_by(orientation) %>% summarise(number_of_routes = n_distinct(route)) 

data %>% group_by(angle) %>% summarise(number_of_routes = n_distinct(route))

data %>% group_by(grade) %>% summarise(number_of_routes = n_distinct(route))

hist(data$SR)
hist(data$rel_abu)

#### Alpha diversity ####

# Shannon and equitability for each frequency
data = data %>%
  group_by(qual_freq) %>%
  mutate(SR_freq = n_distinct(sp_code),
         shannon = diversity(rel_abu, index = "shannon"), 
         simpson = diversity(rel_abu, index = "simpson"),
         hill = (simpson^(-1))/exp(shannon),
         equitability = shannon / SR_freq)

data %>% group_by(qual_freq) %>% distinct( SR = SR_freq, shannon, simpson, hill, equitability) #some difference but not that much

# Shannon and equitability and SR for each route

data = data %>%
  group_by(route) %>%
  mutate(shannon2 = diversity(rel_abu, index = "shannon"), 
         simpson2 = diversity(rel_abu, index = "simpson"),
         hill2 = (simpson2^(-1))/exp(shannon2),
         equitability2 = shannon2 / SR)

#shannon simpson hill

data1 = data %>% select(route, qual_freq, shannon2) %>% rename(value = shannon2) %>% mutate(indice = "Shannon")
data2 = data %>% select(route, qual_freq, simpson2) %>% rename(value = simpson2) %>% mutate(indice = "Simpson")
data3 = data %>% select(route, qual_freq, hill2) %>% rename(value = hill2) %>% mutate(indice = "Hill number")
dataf = rbind(data1, data2, data3)
dataf$indice = factor(dataf$indice, levels = c("Shannon", "Simpson", "Hill number"))




#plot
ggplot(dataf, aes(x = qual_freq, y = value, fill = qual_freq)) +
  geom_boxplot() +
  facet_wrap(vars(indice), scales = "free_y") +
  labs(y = "Indice values", x = "") +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
        panel.background = element_rect(fill = 'white'),
        panel.grid.major = element_line(color = 'lightgray'),
        strip.background = element_rect(color="black", fill="lightgrey", size=0.5, linetype="solid")) +
  scale_fill_manual(values = c("#FFBF00", "#2274A5","#E83F6F"))


#SHANNON
##
ggplot(data, aes(x = qual_freq, y = shannon2)) + # Shannon
  geom_boxplot(outlier.shape = NA) 
##

data %>% group_by(qual_freq) %>% summarise(mean(shannon2)) # less diversity on Vierge..
summary(aov(data$shannon2~data$qual_freq))
TukeyHSD(aov(data$shannon2~data$qual_freq)) 

#EQUITABILITY
##
ggplot(data, aes(x = qual_freq, y = equitability2)) + # equitability
  geom_boxplot(outlier.shape = NA) 
##

data %>% group_by(qual_freq) %>% summarise(mean(equitability2)) 
summary(aov(data$equitability2~data$qual_freq))
TukeyHSD(aov(data$equitability2~data$qual_freq)) # Test it if needed
# Less equitability on Peu fréquenté

#SPECIES RICHNESS
##
ggplot(data, aes(x = qual_freq, y = SR)) + # equitability
  geom_boxplot(outlier.shape = NA) 
##

data %>% group_by(qual_freq) %>% summarise(mean(SR)) 
summary(aov(data$SR~data$qual_freq))
TukeyHSD(aov(data$SR~data$qual_freq)) # Test it if needed
# More specific richness on Peu fréquenté

#HILL
##
ggplot(data, aes(x = qual_freq, y = hill2)) + # equitability
  geom_boxplot(outlier.shape = NA) 
##

# Some plots
plots <- list()
for(freq in levels(data$qual_freq)) {
  subset_data <- subset(data, qual_freq == freq)
  subset_data = data %>% filter(qual_freq == freq) %>%
    group_by(sp_code) %>%
    mutate(rel_abu = sum(rel_abu)) %>%
    distinct(sp_code, .keep_all = T)
  subset_data$sp_code <- factor(subset_data$sp_code, levels = subset_data$sp_code[order(subset_data$rel_abu, decreasing = TRUE)])
  g = ggplot(subset_data, aes(x = reorder(sp_code, -rel_abu), y = rel_abu)) +
    geom_bar(stat = "identity", fill = "skyblue", width = 0.7) +
    labs(title = paste("Abundance of Species -", freq), x = "Species", y = "Abundance") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  plots[[freq]] <- g
  rm(g,freq, subset_data)
}
grid.arrange(grobs = plots, ncol = 3)

datag = data %>% group_by(qual_freq, sp_code) %>% mutate(ab_per_qual = sum(rel_abu)) %>%
  select(qual_freq, sp_code, ab_per_qual) %>% distinct(qual_freq, sp_code, ab_per_qual)

ggplot(datag, aes(x = reorder(sp_code, -ab_per_qual), y = ab_per_qual, fill = qual_freq)) +
  geom_bar(position = "dodge", stat = "identity", width = 0.75) +
  theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.25),
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        panel.background = element_rect(fill = 'white'),
        panel.grid.major = element_line(color = 'lightgray')) +
  labs(x = "", y = "Relative abundance") +
  scale_fill_manual(values = c("#5ae66a", "#6a8fe6", "#e66369"))+
  facet_wrap(vars(qual_freq), scales = "free_x") 
  
#creuser
##
#### Abundance GLMM ####
ggplot(data, aes(x = qual_freq, y = rel_abu)) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(title = "Relative abundance by frequentation",
       x = "Frequentation",
       y = "Relative abundance")+
  theme_bw() 
  

ggplot(data, aes(x = qual_freq, y = rel_abu, fill = qual_freq)) + 
  geom_half_boxplot(center=TRUE, errorbar.draw=FALSE, width=0.5) +
  geom_half_violin(side="r", nudge  = 0.05) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#FFBF00", "#2274A5", "#E83F6F")) +
  ylim(0,3) +
  labs(y = "Relative abundance",
       x = "") 

data %>% group_by(qual_freq) %>% summarise(mean(rel_abu))

model <- glmmTMB(round(rel_abu) ~ qual_freq + (1|route), family = poisson, data = data)
model1 <- glm(round(log(rel_abu + 1)) ~ qual_freq, family = poisson, data = data)
model2 <- glm(round(log(rel_abu + 1)) ~ qual_freq + orientation, family = poisson, data = data)
model3 <- glmmTMB(round(log(rel_abu + 1)) ~ qual_freq + (1|route), family = poisson, data = data)
model4 <- glmmTMB(round(log(rel_abu + 1)) ~ qual_freq + orientation + (1|route), family = poisson, data = data)

AIC(model, model1, model2, model3, model4)
summary(model1)
tab_model(model2)

emmeans(model2, list(pairwise~qual_freq), adjust = "none") #ok model 2 

#### AFC ####

data = data %>% group_by(sp_code) %>% mutate(appear = n()) %>%
  filter(appear > 1)

#AFC (rom)
data$abu_afc = ifelse(data$rel_abu > 0, 1, 0)
data_afc = data[,c("route", "sp_code", "abu_afc")]
data_afc <- unique(data_afc)
data_afc = data_afc %>% pivot_wider(names_from = "sp_code", values_from = "abu_afc")
data_afc[is.na(data_afc)] <- 0
data_afc <- as.data.frame(data_afc)
rownames(data_afc) <- data_afc[,1]
data_afc <-  data_afc[,-1]

data_freq <- data[,c("route", "qual_freq")]
data_freq <- unique(data_freq)
data_afc$qual_freq <- data_freq$qual_freq

afc = dudi.coa(data_afc, scannf = F, nf = 4)
s.label(afc$co)
s.class(afc$li, data_afc$qual_freq,  csta = 1, cgrid = 0, cpoint = 1.4, pch = 20,
        clab = 0, col = c("#FFBF00", "#2274A5", "#E83F6F"), axesell = F,
        grid = F, xlim = c(0,0.3), ylim = c(-1.5,1.5), cellipse = 2)


# CCA
data_afc2 = subset(data_afc, select = - qual_freq)
cca <- cca(data_afc2, data_freq) #constrained correspondance analysis

category_colors <- c("#5ae66a", "#6a8fe6", "#e66369")
plot(cca, type = "n", xlim = c(-3,1), ylim = c(-2,2))  # type = "n" to create an empty plot
points(cca, display = "sites", pch = 20, col = category_colors[as.factor(data_freq$qual_freq)], cex = 0.75)
#text(cca, display = "sites", col = category_colors[as.factor(data_freq$qual_freq)], cex = 0.5, pos = 3)
ordiellipse(cca, groups = data_freq$qual_freq, kind = "sd", col = c("green", "blue","red"), label = F) #kind = "sd"




#
#### Beta diversity ####

# Extract species abundance matrix
data2 = read.csv("data2_relab_freq.csv")
community_data <- data2[, -1]
rownames(community_data) <- unique(data$qual_freq)

# Calculate dissimilarity matrix using Bray-Curtis index
bc_dissimilarity <- vegdist(community_data, method = "bray")
print(bc_dissimilarity)

# clustering
hclust_result <- hclust(bc_dissimilarity)
plot(hclust_result)

#### Functional role ####
#specialist richness ratio
data$ri_spe = ifelse(data$function. == "specialiste", 1, 0)
data = data %>% group_by(qual_freq, route) %>% 
  mutate(ratio_ri_spe = sum(ri_spe) / SR)

ggplot(data, aes(x = qual_freq, y = ratio_ri_spe, fill = qual_freq)) + 
  geom_half_boxplot(center=TRUE, errorbar.draw=FALSE, width=0.5) +
  geom_half_violin(side="r", nudge  = 0.05) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#5ae66a", "#6a8fe6", "#e66369")) +
  labs(title = "Specialist ratio by frequentation modalities", y = "Specialist richness ratio",
       x = "") 

summary(aov(data$ratio_ri_spe ~ data$qual_freq))
TukeyHSD(aov(data$ratio_ri_spe ~ data$qual_freq))

# specialist abundance ratio
data$ab_spe = ifelse(data$function. == "specialiste", data$rel_abu, 0)
data = data %>% group_by(qual_freq, route) %>% 
  mutate(ratio_ab_spe = sum(ab_spe) / sum(rel_abu))

ggplot(data, aes(x = qual_freq, y = ratio_ab_spe, fill = qual_freq)) + 
  geom_half_boxplot(center=TRUE, errorbar.draw=FALSE, width=0.5) +
  geom_half_violin(side="r", nudge  = 0.05) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#5ae66a", "#6a8fe6", "#e66369")) +
  labs(title = "Specialist ratio by frequentation modalities", y = "Specialist abundance ratio",
       x = "") 

summary(aov(data$ratio_ab_spe ~ data$qual_freq))
TukeyHSD(aov(data$ratio_ab_spe ~ data$qual_freq))
plot(aov(data$ratio_ab_spe ~ data$qual_freq))

kruskal.test(data$ratio_ab_spe ~ data$qual_freq)



#### Get data2's ####
data <- data[data$orientation != "W", ] # ?????

#sp_relabu
data2 = data %>% select(route, sp_code, rel_abu) %>%
  group_by(route, sp_code) %>% mutate(rel_abu = sum(rel_abu)) %>% 
  group_by(route) %>% distinct(sp_code, rel_abu) %>%
  pivot_wider(names_from = "sp_code", values_from = "rel_abu")
data2 = as.data.frame(data2)
rownames(data2) <- data2[,1]
data2 = data2[,-1]
data2[is.na(data2)] <- 0
write.csv(data2, file = "data2_relab.csv", row.names = F)

#sp_relabu par freq
data2 = data %>% select(qual_freq, sp_code, rel_abu) %>%
  group_by(qual_freq, sp_code) %>% mutate(rel_abu = sum(rel_abu)) %>% 
  group_by(qual_freq) %>% distinct(sp_code, rel_abu) %>%
  pivot_wider(names_from = "sp_code", values_from = "rel_abu")
data2 = as.data.frame(data2)
rownames(data2) <- data2[,1]
data2 = data2[,-1]
data2[is.na(data2)] <- 0
write.csv(data2, file = "data2_relab_freq.csv", row.names = F)

#sp_relabu for each modality
abundance_matrix = list()
for(i in unique(data$qual_freq)) {
  data2 = data %>% filter(qual_freq == i) %>% select(route, sp_code, rel_abu) %>%
    group_by(route, sp_code) %>% mutate(rel_abu = sum(rel_abu)) %>% 
    group_by(route) %>% distinct(sp_code, rel_abu) %>%
    pivot_wider(names_from = "sp_code", values_from = "rel_abu")
  data2 = as.data.frame(data2)
  rownames(data2) <- data2[,1]
  data2 = data2[,-1]
  data2[is.na(data2)] <- 0
  
  abundance_matrix[[i]] = data2
  rm(i, data2)
}


# sp_abu
data2 = data %>% select(route, sp_code, abundance) %>%
  group_by(route, sp_code) %>% mutate(abundance = sum(abundance)) %>% 
  group_by(route) %>% distinct(sp_code, abundance) %>%
  pivot_wider(names_from = "sp_code", values_from = "abundance")
data2 = as.data.frame(data2)
rownames(data2) <- data2[,1]
data2 = data2[,-1]
data2[is.na(data2)] <- 0
write.csv(data2, file = "data2_ab.csv", row.names = F)

#### rarefaction curve ####

# rarefaction curve
data2 = read.csv("data2_ab.csv")
data2 = colSums(data2) # with data2_ab
rarefaction_result <- iNEXT(data2, q = 0, datatype = "abundance", endpoint = 6000)

ggiNEXT(rarefaction_result, type=1, se=TRUE, grey=F)+
  ggtitle('Rarefaction curve')+
  theme_bw()




#### test the site effect ####
data$site = ifelse(data$sector == "Bisous_dans_le_cou" | data$sector == "Rognon", "Site1", "Site2")
ggplot(data, aes(x = site, y = rel_abu, fill = site)) +
  geom_boxplot() +
  labs(x = "Site", y = "Relative abundance") +
  theme_bw()+
  theme(legend.position = "none")

ggplot(data, aes(x = site, y = SR, fill = site)) +
  geom_boxplot()+
  labs(x = "Site", y = "Species richness") +
  theme_bw() +
  theme(legend.position = "none")



data$abu_afc = ifelse(data$rel_abu > 0, 1, 0)
data_afc = data[,c("route", "sp_code", "abu_afc")]
data_afc <- unique(data_afc)
data_afc = data_afc %>% pivot_wider(names_from = "sp_code", values_from = "abu_afc")
data_afc[is.na(data_afc)] <- 0
data_afc <- as.data.frame(data_afc)
rownames(data_afc) <- data_afc[,1]
data_afc <-  data_afc[,-1]

data_freq <- data[,c("route", "site")]
data_freq <- unique(data_freq)
data_afc$site <- data_freq$site
data_afc$site = as.factor(data_afc$site)


afc = dudi.coa(data_afc, scannf = F, nf = 4)
s.label(afc$co)
s.class(afc$li, data_afc$site,  csta = 0, cgrid = 0, 
        clab = 0.5, col = c("#5ae66a", "#6a8fe6", "#e66369"), axesell = F,
        grid = F, xlim = c(0,0.3), ylim = c(-1.5,1.5), cellipse = 2)

#### Morphology test ####

#SR
datag = data %>% group_by(qual_freq, morpho) %>% mutate(SR = n_distinct(sp_code)) %>%
  distinct(qual_freq, morpho, SR)

ggplot(datag, aes(x = reorder(morpho, -SR), y = SR, fill = qual_freq)) +
  geom_bar(position = "dodge", stat = "identity", width = 0.75) +
  theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.25),
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        panel.background = element_rect(fill = 'white'),
        panel.grid.major = element_line(color = 'lightgray')) +
  labs(x = "", y = "Species richness") +
  scale_fill_manual(values = c("#5ae66a", "#6a8fe6", "#e66369"))+
  facet_wrap(vars(qual_freq), scales = "free_x") 

# abu
datag = data %>% group_by(qual_freq, morpho) %>% mutate(ab_per_mor = sum(rel_abu)) %>%
  distinct(qual_freq, morpho, ab_per_mor)

ggplot(datag, aes(x = reorder(morpho, -ab_per_mor), y = ab_per_mor, fill = qual_freq)) +
  geom_bar(position = "dodge", stat = "identity", width = 0.75) +
  theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.25),
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        panel.background = element_rect(fill = 'white'),
        panel.grid.major = element_line(color = 'lightgray')) +
  labs(x = "", y = "Relative abundance") +
  scale_fill_manual(values = c("#5ae66a", "#e66369"))+
  facet_wrap(vars(qual_freq), scales = "free_x") 

#### Analyses Raffinées ####
# Ressortir le graphique avec toutes les sp lorsque le dataset est nettoyé 

#graph alpha

ggplot(data3, aes(x = reorder(sp_code, -rel_abu), y = rel_abu, fill = qual_freq)) +
  geom_bar(position = "dodge", stat = "identity", width = 0.75) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) +
  labs(x = "", y = "Relative abundance") +
  scale_fill_manual(values = c("#5ae66a", "#e66369"))+
  guides(fill=guide_legend(title="")) +
  theme(axis.text.x = element_text(face = "bold", size = 10),panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

datag = data3 %>% group_by(qual_freq, sp_code) %>% mutate(ab_per_qual = sum(rel_abu)) %>%
  select(qual_freq, sp_code, ab_per_qual) %>% distinct(qual_freq, sp_code, ab_per_qual)

ggplot(datag, aes(x = reorder(sp_code, -ab_per_qual), y = ab_per_qual, fill = qual_freq)) +
  geom_bar(position = "dodge", stat = "identity", width = 0.75) +
  theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.25),
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        panel.background = element_rect(fill = 'white'),
        panel.grid.major = element_line(color = 'lightgray')) +
  labs(x = "", y = "Relative abundance") +
  scale_fill_manual(values = c("#5ae66a", "#e66369"))+
  facet_wrap(vars(qual_freq), scales = "free_x") 

#beta
data3$abu_afc = ifelse(data3$rel_abu > 0, 1, 0)
data_afc = data3[,c("route", "sp_code", "abu_afc")]
data_afc <- unique(data_afc)
data_afc = data_afc %>% pivot_wider(names_from = "sp_code", values_from = "abu_afc")
data_afc[is.na(data_afc)] <- 0
data_afc <- as.data.frame(data_afc)
rownames(data_afc) <- data_afc[,1]
data_afc <-  data_afc[,-1]

data_freq <- data3[,c("route", "qual_freq")]
data_freq <- unique(data_freq)
data_afc$qual_freq <- data_freq$qual_freq

afc = dudi.coa(data_afc, scannf = F, nf = 4)
s.label(afc$co)
s.class(afc$li, data_afc$qual_freq,  csta = 0, cgrid = 0, 
        clab = 0, col = c("#5ae66a", "#e66369"), axesell = F,
        grid = F, xlim = c(0,0.3), ylim = c(-1.5,1.5), cellipse = 2)


#specialist richness ratio
data3$ri_spe = ifelse(data3$function. == "specialiste", 1, 0)
data3 = data3 %>% group_by(qual_freq, route) %>% 
  mutate(ratio_ri_spe = sum(ri_spe) / SR)

ggplot(data3, aes(x = qual_freq, y = ratio_ri_spe, fill = qual_freq)) + 
  geom_half_boxplot(center=TRUE, errorbar.draw=FALSE, width=0.5) +
  geom_half_violin(side="r", nudge  = 0.05) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#5ae66a", "#6a8fe6", "#e66369")) +
  labs(title = "Specialist ratio by frequentation modalities", y = "Specialist richness ratio",
       x = "") 

# specialist abundance ratio
data3$ab_spe = ifelse(data3$function. == "specialiste", data$rel_abu, 0)
data3 = data3 %>% group_by(qual_freq, route) %>% 
  mutate(ratio_ab_spe = sum(ab_spe) / sum(rel_abu))

ggplot(data3, aes(x = qual_freq, y = ratio_ab_spe, fill = qual_freq)) + 
  geom_half_boxplot(center=TRUE, errorbar.draw=FALSE, width=0.5) +
  geom_half_violin(side="r", nudge  = 0.05) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#5ae66a", "#6a8fe6", "#e66369")) +
  labs(title = "Specialist ratio by frequentation modalities", y = "Specialist abundance ratio",
       x = "") 



# Shannon and equitability for each frequency
data3 = data3 %>%
  group_by(qual_freq) %>%
  mutate(SR_freq = n_distinct(sp_code),
         shannon = diversity(rel_abu, index = "shannon"), 
         simpson = diversity(rel_abu, index = "simpson"),
         hill = (simpson^(-1))/exp(shannon),
         equitability = shannon / SR_freq)

data3 %>% group_by(qual_freq) %>% distinct( SR = SR_freq, shannon, simpson, hill, equitability) #some difference but not that much

# Shannon and equitability and SR for each route

data3 = data3 %>%
  group_by(route) %>%
  mutate(shannon2 = diversity(rel_abu, index = "shannon"), 
         simpson2 = diversity(rel_abu, index = "simpson"),
         hill2 = (simpson2^(-1))/exp(shannon2),
         equitability2 = shannon2 / SR)

data3 %>% group_by(qual_freq) %>% summarise( mean(rel_abu), mean(SR),mean(shannon2), mean(simpson2), mean(hill2), mean(equitability2)) #some difference but not that much



# result graph
data2 = data %>% group_by(sp_code) %>% mutate(appear = n()) %>%
  filter(appear > 1) %>% group_by(qual_freq, sp_code) %>% mutate(rel_abu = sum(rel_abu)) %>%
  distinct(qual_freq, sp_code, rel_abu)


ggplot(data2, aes(x = reorder(sp_code, -rel_abu), y = rel_abu, fill = qual_freq)) +
  geom_bar(position = "dodge", stat = "identity", width = 0.75) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0, hjust= 0.7)) +
  labs(x = "", y = "Relative abundance") +
  scale_fill_manual(values = c("#FFBF00", "#2274A5", "#E83F6F"))+
  guides(fill=guide_legend(title="")) +
  theme(axis.text.x = element_text(face = "bold", size = 10, vjust=0.35),panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.position = "none")



# pdf with abundance diff
data = data3
filtre = data %>% ungroup() %>% distinct(qual_freq, sp_code) %>% group_by(sp_code) %>%
  mutate(n = n()) %>% filter(n > 1) %>% distinct(sp_code)
filtre = filtre$sp_code
data = data %>% filter(sp_code %in% filtre)


pdf("abundance test.pdf")
for(s in unique(data$sp_code)) {
  data2 = data %>% filter(sp_code == s)
  
  p_value <- round(summary(aov(data2$rel_abu ~ data2$qual_freq))[[1]]$`Pr(>F)`[1], digits = 10)

  titre = paste(s, "  p_value =", p_value)
  
  g = ggplot(data2, aes(x = qual_freq, y = rel_abu)) +
    geom_boxplot() +
    labs(x = "", y = "Relative abundance", title = titre)+
    theme_bw()
    
  print(g)
}
dev.off()

