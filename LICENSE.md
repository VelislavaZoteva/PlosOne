# plot extra-axial tumor formation, logistic regression and a Cox proportional hazards model.
# --------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------
library(tidyverse)
library(effects)
library(survival)
library(readxl)
library(graphics)
library(ggplot2)
library(lubridate)
library(ggpubr)
library(rstatix)
library(ggstatsplot)

# logistic regression
rat <- readxl::read_excel("Statistiek.xlsx", sheet = 1)

rat <- rat %>% 
  rename(verloop_ratio = `Verloop ratio`,
         Group = `Group Ratio`) %>% 
  mutate(Group = str_replace_all(Group, " ", "")) %>% 
  mutate(cells = as.numeric(str_extract(Group, "\\d+")))
  
rat <- rat %>% 
  mutate(verloop_bin = ifelse(verloop_ratio == "gunstig", 0, 1))

mod <- glm(verloop_bin ~ Group, data = rat, family = binomial(link = "logit") )
summary(mod)

mod0 <- glm(verloop_bin ~ 1, data = rat, family = "binomial")
drop1(mod, test = "LRT")
anova(mod0, mod, test = "LRT")

rat %>% 
  count(Group, verloop_bin)
table(rat$Group, rat$verloop_bin) 
table(rat$Group, rat$verloop_bin) %>% prop.table(1) %>% round(2)

mod_num <- glm(verloop_bin ~ cells, data = rat, family = binomial(link = "logit") )
summary(mod_num)
plot(effects::predictorEffects(mod_num), type = "response", xlab = 'F98 cells', ylab='Probability', main='Prevalence of EA tumor growth')

# Cox proportional hazards model
rat <- readxl::read_excel("Statistiek.xlsx", sheet = 1)

rat <- rat %>% 
  rename(verloop_ratio = `Verloop ratio`,
         Group = `Group Ratio`) %>% 
  mutate(Group = str_replace_all(Group, " ", "")) %>% 
  mutate(cells = as.numeric(str_extract(Group, "\\d+")))
#create long data
# ----------------
cells <- rat %>% select(ID, Group, cells)
rat <- readxl::read_excel("Statistiek.xlsx", sheet = 2)
rat <- left_join(rat, cells)
# Data set creation for interval based data, one event type, multiple events per subject possible
# -----------------------------------------------------------------------------------------------
# add event
rat <- rat %>% 
  mutate(Event = ifelse(Verloop == "ongunstig", 1, 0)) %>% 
  arrange(ID, Tijd)

# add indicator
rat <- rat %>% 
  group_by(ID) %>% 
  mutate(indicator = ifelse((sum(Event) == 0 & Tijd == max(Tijd)) | # Alleen gunstig --> laatste waarde van Tijd
                              (sum(Event) > 0  & Event - lag(Event, default = 0) == 1 ), # Eerste verandering naar ongunstig
                            1, 0)) 

# add tstart and tstop interval
rat <- rat %>% 
  mutate(tstart = lag(Tijd,default=-1),
         tstop = Tijd) %>% 
  arrange(ID, Tijd)

rat %>% print(n = 50)

# Data set creation for interval based data, one event type, first events per subject
# -----------------------------------------------------------------------------------

# filter by indicator
survrat <- rat %>% 
  filter(indicator == 1) 

# data file for 
survrat
# Differences in survival curves between groups?
# ----------------------------------------------

# Log rank test:
# --------------
survdiff(Surv(Tijd, Event) ~ cells, survrat)

# Cox PH model (cells categorical)
# --------------------------------

# Create a factor 
survrat$cellsc <- factor(survrat$cells, levels = dput(unique(survrat$cells)))
is.factor(survrat$cellsc)

# fit Cox PH model
options(show.signif.stars=T) # display statistical intelligence
cfit1 <- coxph(Surv(Tijd, Event) ~ cellsc, data=survrat)
print(cfit1, digits=3)
summary(cfit1, digits=3)
anova(cfit1)

# plot fitted model
dummy <- expand.grid(cellsc=c(500, 1000,5000,10000,20000))
dummy$cellsc <- factor(dummy$cells, levels = dput(unique(survrat$cells)))
dummy
csurv1 <- survfit(cfit1, newdata=dummy)
dim(csurv1)

# visualisation fitted model
par(mar=c(5, 4, 4, 8), xpd=TRUE)  # crëert meer ruimte rechts van plot
plot(csurv1, col=1:5, xlab="Days from GB confirmation (day 0)", ylab="Survival", xaxs="S")
legend(40, .99, legend =  c("500 cells", "1000 cells","5000 cells", "10000 cells", "20000 cells"),
      col=1:5, lwd=0.3, bty='n')
#plot(csurv1, col=1:5, xlab="Days since extra tumor appearance", ylab="Survival", xaxs="S",conf.times = c(3, 6, 9, 12, 15, 18, 21, 24, 27, 30))

# check PH assumption

zp1 <- cox.zph(cfit1)
zp1

# No indication of violation of PH assumption.

# plot for illustration of extracranial tumor formation probabilities
# -------------------------------------------------------------------
# -------------------------------------------------------------------

data <- read_excel("ECT.xlsx",sheet = 2)
# binairy scoring for present/absent extracranial tumor
yes <- c(sum(data$`500 cells`[1:5] == 'yes'), sum(data$`1000 cells` == 'yes'), sum(data$`5000 cells` == 'yes'),sum(data$`10 000 cells` == 'yes'), sum(data$`20 000 cells`[1:5] == 'yes') )
no <- c(sum(data$`500 cells`[1:5] == 'no'), sum(data$`1000 cells` == 'no'), sum(data$`5000 cells` == 'no'),sum(data$`10 000 cells` == 'no'), sum(data$`20 000 cells`[1:5] == 'no') )

percentage <- (yes/(yes+no))*100
groups <- c('500 cells', '1000 cells', '5000 cells', '10 000 cells', '20 000 cells')

ECT <- data.frame(groups,yes,no,percentage)
View(ECT)

stat <- read_excel('ECT.xlsx', sheet = 3)
stat <- as.data.frame(stat)
View(stat)

test <- fisher.test(table(stat$group,stat$ECT))
p <- round(test$p.value, digit = 4)
p
ggplot(ECT, aes(x=groups,y=percentage))+geom_bar(stat='identity')
plot <- ggplot(ECT, aes(x=groups,y=percentage))+
  geom_bar(stat='identity', fill="gold")+
  theme_minimal()+
  geom_text(ECT, mapping = aes(label = paste(round(percentage, digits = 1), '%'),hjust = -0.4, fontface = "bold")) +
  scale_x_discrete(limits=c('500 cells','1000 cells', '5000 cells', '10 000 cells', '20 000 cells'))+
  annotate("text",x=1.5,y=70,label= paste0("Fisher's exact test, p-value = ", p))
plot <- plot + ylim(0,100) #y-as limiet naar 100% krijgen
plot <- plot+coord_flip()
plot <- print(plot + ggtitle("EC prevalence")+theme(plot.title = element_text(hjust = 0.5))+labs(x= "", y = "%"))
plot
