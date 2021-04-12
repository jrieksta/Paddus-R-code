# -------------------------------------------------------------
# LMM PADDUS DATA example to Tao 12042021
# Jolanta Rieksta

# -------------------------------------------------------------

#############Load libraries################
.libPaths()
.libPaths("C:/R Packages")
Sys.setenv(LANGUAGE='en')
suppressPackageStartupMessages({
  library(dplyr)
  library(HistData)
  library(vegan)
  library(devtools)
  library(tidyverse)
  library(scales)
  library(ggplot2)
  library(LabApplStat)
  library(emmeans)
  library(multcompView)
  library(Rmisc)
  library(sciplot)
  library(ggpubr)
  library(lme4)
  library(nlme)
  library(wesanderson)
  library(multcomp)
  library(readxl)
  library(plyr)
  library(RColorBrewer)
  library(dplyr)
  library(car)
  library(e1071)
  library(lme4)
  library(bbmle)
  library(generics)
  library(caret)
  library(leaps)
  library(emmeans)
  library(devtools)
  library(LabApplStat)
  library(sjmisc)
  library(sjPlot)
  library(ggpubr)
  library(viridisLite)
  library(viridis)
  library(treemap)
  library(shiny)
  library(esquisse)
  library(remotes)
  library(ggpattern)
})
#################Load data##########

setwd("C:/Users/lcm767/OneDrive - K?benhavns Universitet/Desktop/Projects/Paddus 2018")
d1<- read_excel("Paddus to Tao_TAO_1412202020.xlsx", sheet = "Totals")
d1$log(sum+1) <- log(d1$sum+1) 
str(d1)
d1$Date[d1$Date=="30.06.2018"] <- "June 30"
d1$Date[d1$Date=="05.07.2018"] <- "July 5"
d1$Date[d1$Date=="11.07.2018"] <- "July 11"

names(d1)[names(d1) == "warming"] <- "treatment"



#######################Figure parameters - level order, theme#####################
level_order <- c('C','S', "W", 'C+H','S+H' ,"W+H" )
GLV1$Date1 = factor(d1$Date, levels=c('June 30','July 5','July 11'))


my_theme=theme(text = element_text(size=14),
               axis.text.x = element_text(size = 12),
               axis.text.y = element_text(size = 12),
               strip.text.x = element_text(size = 14),
               strip.text.y  = element_text(size = 14),
               panel.border = element_rect(colour = "black", fill=NA, size=0.5),
               
               strip.background = element_blank(),
               axis.title.x=element_text(size = 14),
               axis.title.y=element_text(size = 14),
               
) 


######################## LMM #######################
# SUBSET DATA  non std

GLV<-subset(d1,(Compound_group=="green leaf volatile")) #Y
HT<-subset(d1,(Compound_group=="HT")) #Y
HC<-subset(d1,(Compound_group=="hydrocarbons"))
ISO<-subset(d1,(Compound_group=="isoprenoids"))
MT<-subset(d1,(Compound_group=="MT")) #Y
other<-subset(d1,(Compound_group=="other"))
OVOC<-subset(d1,(Compound_group=="oxygenated VOCs"))
SQT<-subset(d1,(Compound_group=="SQT")) #Y
total<-subset(d1,(Compound_group=="Total")) #Y


test<-ddply(MT,.(herbivory,Date, treatment),summarize, N=length(sum), mean=mean(sum), sd=sd(sum),se=sd/sqrt(N) )



############# LOG TRANSF vs original values #################
total$Block=as.factor(total$Block)
total$Date=as.factor(total$Date)
total$Plot=as.factor(total$Plot)
total$treatment=as.factor(total$treatment)
total$herbivory=as.factor(total$herbivory)
total$Date=as.factor(total$Date)

names(total)
a<-lmer(sum~treatment*herbivory*Date+(1|Block/Plot),data=total)
plot(a)
plot(a, factor(Block)~resid(.), abline=0)
plot(a, factor(Plot)~resid(.), abline=0)
c<-qqnorm(resid(a))
d<-qqline(resid(a))
str(total)
b<-lmer(log(sum+1)~treatment*herbivory*Date+(1|Block/Plot),data=total)

plot(b)
e<-qqnorm(resid(b))
f<-qqline(resid(b))
plot(b, factor(Block)~resid(.), abline=0)
plot(b, factor(Plot)~resid(.), abline=0)




#install.packages("insight")
#library(insight)
#vars <- insight::get_variance(m4)
#r2_marginal <- vars$var.fixed / (vars$var.fixed + vars$var.random + vars$var.residual)
#r2_conditional <- (vars$var.fixed + vars$var.random) / (vars$var.fixed + vars$var.random + vars$var.residual)

#So, the marginal R2 is the fixed effects variance, divided by the total variance (i.e. fixed + random + residual). This value indicates how much of the "model variance" is explained by the fixed effects part only.
#The conditional R2 is the fixed+random effects variance divided by the total variance, and indicates how much of the "model variance" is explained by your "complete" model.
#If you would like to know how much of the proportion of variance can be explained by the random effects only, then you have the ICC. This is simply:
#icc_adjusted <- vars$var.random / (vars$var.random + vars$var.residual)




##################################    MT      ################################
plot(DD(sum~treatment*herbivory*Date,random=~Block/Plot,data=MT),circle="MSS")

MT$herbivory <- factor(MT$herbivory)
MT$treatment <- factor(MT$treatment)
MT$Block <- factor(MT$Block)
MT$Plot <- factor(MT$Plot)
MT$Date <- factor(MT$Date)
# 
# Validate initial model
# 

# Fit LMM (linear mixed effects model)
m0<-lmer(log(sum+1)~treatment*herbivory*Date+(1|Block/Plot),data=MT)
anova(m0)

# Validate model: remember that we have two sets of random effects in our model
plot(m0)
qqnorm(residuals(m0))
qqnorm(ranef(m0)$'Block'[,1])
qqnorm(ranef(m0)$'Plot'[,1])

qqline(resid(m0))
plot(m0, factor(Block)~resid(.), abline=0)
plot(m0, factor(Plot)~resid(.), abline=0)

# 
# Is there an effect?
# 

# Do model selection starting from m0
drop1(m0,test="Chisq")
drop1(m1 <- update(m0,.~.-treatment:herbivory:Date),test="Chisq")
drop1(m2 <- update(m1,.~.-treatment:herbivory),test="Chisq")
#validate 
anova(m2)
plot(m2)
qqnorm(residuals(m2))
r.rmse <- sqrt(mean(residuals(m2)^2))
plot(fitted(m2),total$log(sum+1))

tab_model(m2)
summary(m2)

# 
# Where is the effect? : emmeans
# 
# emmeans for 

#treatment
sigma(m2)
emmeans(m2, specs=trt.vs.ctrl~treatment,type="response", bias.adj=TRUE, sigma=0.6643695)


#Date
Date_total <- emmeans(m2,~Date, cont="pairwise", int.adjust="tukey", type="response", alpha=0.05, bias.adj=TRUE, sigma=0.6643695)
Date_total

#herbivory
emmeans(m2,~herbivory, cont="pairwise", int.adjust="tukey", type="response", alpha=0.05, bias.adj=TRUE, sigma=0.6643695)
#herbivory x date
emmeans(m2,~herbivory|Date, cont="pairwise", int.adjust="tukey", type="response", alpha=0.05, bias.adj=TRUE, sigma=0.6643695)

#treatment x date
emmeans(m2, specs=trt.vs.ctrl~treatment|Date,type="response", bias.adj=TRUE, sigma=0.6643695)

