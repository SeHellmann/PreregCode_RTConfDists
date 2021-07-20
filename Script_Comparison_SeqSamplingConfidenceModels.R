###########################################################################
#####      Comparison Sequential Sampling Confidence Models        #######
###########################################################################

# Sebastian Hellmann, 20.07.2021

# 1) Read in data and aggregate data for later visualization
# 2) Fit the models (dynWEV, 2DSD, IRM(t) and PCRM(t))  and predict rating and rt distribution
# 3) Analysis of Criteria and Computation of BayesFactors
# 4) Visualisation of observations and predictions
# 5) Analyse covariances between parameters



###   Preamble and imports    ####
rm(list = ls())
# # insert path by hand:
# script_path <- paste0("C:/Users/PPA859/Documents", # insert right path, here
#                       "/AccumulatorModels")
# or use RStudio to find script file path
script_path <- rstudioapi::getSourceEditorContext()$path
setwd(script_path)

library(plyr)
library(snow)
library(doSNOW)
library(BayesFactor)
library(tidyverse)
library(dynWEV)
library(RColorBrewer)
library(gridExtra)
source("functions/ROCs.R")

#=============================================================
#=============================================================
###   Load results from previous analysis   ####
### (Includes results from sections 1) & 2) ###
load("collected_fitsNpredicts.RData")
#=============================================================
#=============================================================
########### 1) Read, Preprocess, Aggregate Data ##############
#### Load Data   #####
# 1) Read in data
data_folder <- "dataRausch2018Exp2/"
source("functions/gather_data.R")
load("dataRausch2018.RData")
nRatings <- length(unique(Data$rating))
nConds <- length(unique(Data$condition))
cond_levels <- sort(unique(Data$condition))
Data <- Data %>% 
  mutate(condition = factor(condition,levels = cond_levels, 
                                           labels = paste(c("8.3", "16.7", "33.3", "66.7", "133.3"), "")))

#### Compute confidence rating distribution of the Data    #####
Data_RatingDist_part <- Data %>% 
  group_by(stimulus, condition, rating, response, correct, participant) %>% 
  summarise(p = n()/(mean(nrows))) %>% 
  full_join(y = expand.grid(stimulus = unique(Data$stimulus), 
                          condition = unique(Data$condition),
                          rating = 1:nRatings, 
                          correct = c(0,1), 
                          participant = unique(Data$participant))) %>%
  mutate(p = ifelse(is.na(p), 0, p)) %>%
  group_by(correct, condition, rating, participant) %>% 
  summarise(p = mean(p))
### Sanity Checks:
# sum(Data_RatingDist_part$p)
# table((Data_RatingDist_part %>% group_by(condition, participant) %>% summarise(p = sum(p)))$p)

# For the plots we won't differentiate between 
# stimulus directions and participants
Data_RatingDist_corr_cond <- Data_RatingDist_part %>% 
  group_by(correct, condition, rating) %>% 
  summarise(p = mean(p))
#sum(Data_RatingDist_corr_cond$p)

#### Compute Mean Rating of the Data    #####
Data_MRating_corr_cond <- Data_RatingDist_part %>% 
  group_by(condition, participant, correct) %>%
  summarise(MRating = sum(rating*p)/sum(p)) %>%
  group_by(condition, correct) %>%
  summarise(SER = sd(MRating, na.rm = TRUE)/sqrt(n()),
            MRating = mean(MRating, na.rm = TRUE))

#### Compute Reaction Time Quantiles of the Data with different grouping factors  #####
Data_RTQuants_cond <- Data %>% 
  group_by(condition) %>%
  summarise(p=c(.1,.3,.5,.7,.9), q = quantile(rt, probs = c(.1,.3,.5,.7,.9))) %>%
  left_join(summarise(group_by(Data_RatingDist_corr_cond, condition), p_correct=sum(p)), by=c("condition")) 
Data_RTQuants_corr_cond <- Data  %>%
  group_by(condition, correct) %>%
  summarise(p=c(.1,.3,.5,.7,.9), q = quantile(rt, probs = c(.1,.3,.5,.7,.9))) %>%
  left_join(summarise(Data_RatingDist_corr_cond, p_correct=sum(p)), by=c("correct","condition"))
Data_RTQuants_rating <- Data %>%
  group_by(rating) %>%
  summarise(p=c(.1,.3,.5,.7,.9), q = quantile(rt, probs = c(.1,.3,.5,.7,.9))) %>%
  left_join(summarise(group_by(Data_RatingDist_corr_cond, rating), p_correct=sum(p)/nConds), by=c("rating"))%>%
  mutate(rating = factor(rating,
                         levels=1:nRatings, 
                         labels = paste0(c("Guessing", "Very Unsure", "Unsure", "Rather Sure", "Sure"), "")))
Data_RTQuants_corr_rating <- Data %>%
  group_by(rating, correct) %>%
  summarise(p=c(.1,.3,.5,.7,.9), q = quantile(rt, probs = c(.1,.3,.5,.7,.9))) %>%
  left_join(summarise(group_by(Data_RatingDist_corr_cond, rating, correct), p_correct=mean(p)), by=c("correct","rating"))
### Compute ROC curves in Data
Data_ROC_cond <- Data %>% group_by(condition) %>%
  mutate(stimulus = as.factor(stimulus)) %>%
  summarise(getROCcoordFromData(cur_data()$rating, cur_data()$stimulus, cur_data()$correct)) %>%
  mutate(HitRate = qnorm(HitRate), FalseAlarmRate = qnorm(FalseAlarmRate))


####################### 2) Fit Models and Predict Distributions   #################################
n.cores <- parallel::detectCores() - 1
cl <- makeCluster(n.cores, "SOCK", outfile = "")
registerDoSNOW(cl)
clusterEvalQ(cl, library(dynWEV))

######### Start fitting with models including postdecisional accumulation of evidence #############
fitData <- cbind(rbind(Data, Data), model = rep(c("2DSD", "WEVmu"), each = nrow(Data)))
clusterExport(cl, "fitData")

t00 <- Sys.time()
fits_WEVmodels <- ddply(fitData,.(participant, model),
                              function(df) fitWEV(df, df$model[1], logging = TRUE, nRatings = 5),
                              .parallel = F) # .parallel refers to an internal parallelization 
# within the fitting (we parallelize over participants, here, 
# and not within the fitting process)

dir.create("saved_fits")
save(file = "saved_fits/fits_2DSD_WEV.RData", fits_WEVmodels)
print(paste("Fitting 2DSD and WEV took...",
            as.character(round(as.double(difftime(Sys.time(),t00,units = "mins")), 2)),
            " mins"))

##### Fit Race Models/ bounded accumulation models  #######################
fitData <- cbind(rbind(Data, Data, Data, Data), 
                 model = rep(c("IRM", "PCRM", "IRMt", "PCRMt"), each = nrow(Data)))
clusterExport(cl, "fitData")

t00 <- Sys.time()
fits_RMmodels <- ddply(fitData,.(participant, model),
                       function(df) fitRM(df, df$model[1], time_scaled = FALSE,
                                          logging = TRUE, nRatings = 5),
                       .parallel = F)
dir.create("saved_fits")
save(file = "saved_fits/fits_RacingModels.RData", fits_RMmodels)
print(paste("Fitting both Racing Models took...",
            as.character(round(as.double(difftime(Sys.time(),t00,units = "mins")), 2)),
            " mins"))

#=============================================================
#####   Compute model predictions with fitted parameters 
#####   Prediction of response and rating distributions  #####
clusterExport(cl, c("fits_WEVmodels"))
clusterExport(cl, c("fits_RMmodels"))
preds_WEV <- ddply(fits_WEVmodels,.(participant, model),
                    function(df) predictWEV_Conf(df, model = df$model[1], subdivisions = 1000),
                    .parallel = T)
table(preds_WEV$info)
preds_RM <- ddply(fits_RMmodels,.(participant, model),
               function(df) predictRM_Conf(df, model = df$model[1], subdivisions = 1000, time_scaled = FALSE),
               .parallel = TRUE)
table(preds_RM$info)


#####          Prediction of RT densities              #######
maxrt <- 6  ## For visualization a maximum RT of 6 is enough
compute_RTdens_from_fits <- function(Confpred_data) {
  if (Confpred_data$model[1] %in% c("IRM","PCRM","IRMt", "PCRMt")) {
    paramDf <- subset(fits_RMmodels, 
                      model == Confpred_data$model[1] & participant == Confpred_data$participant[1])
    res <- predictRM_RT(paramDf = paramDf, 
                        model = Confpred_data$model[1], 
                        maxrt = maxrt, subdivisions = 300, minrt = min(fits_WEVmodels$t0, fits_RMmodels$t0),
                        scaled = TRUE, DistConf = Confpred_data,
                        .progress = FALSE)
  } else {
    paramDf <- subset(fits_WEVmodels, 
                      model == Confpred_data$model[1] & participant == Confpred_data$participant[1])
    res <- predictWEV_RT(paramDf = paramDf, 
                         model = Confpred_data$model[1], 
                         maxrt = maxrt, subdivisions = 300, minrt = min(fits_WEVmodels$t0, fits_RMmodels$t0), 
                         scaled = TRUE, DistConf = Confpred_data,
                         .progress = FALSE)
  }
  return(res)
}

clusterExport(cl, c("preds_WEV", "preds_RM", "maxrt", "compute_RTdens_from_fits"))
RT_dist_RM <- ddply(preds_RM,.(participant, model),
                    .fun = compute_RTdens_from_fits, 
                    .parallel = T)
RT_dist_RM <- RT_dist_RM %>% group_by(model, participant, correct, rating, condition, rt) %>%
  summarise(dens = mean(dens), densscaled = mean(densscaled))

RT_dist_WEV <- ddply(preds_WEV,.(participant, model),
                     .fun = compute_RTdens_from_fits, 
                     .parallel = T)
RT_dist_WEV <- RT_dist_WEV %>% group_by(model, participant, correct, rating, condition, rt) %>%
  summarise(dens = mean(dens), densscaled = mean(densscaled))

stopCluster(cl)
rm(list = c("cl"))


#=============================================================
#####   Aggregating Predictions for Visualization     ########
# For the visualization we aggregate always over stimulus and
# response category, because we are only interested in correct
# vs. incorrect responses;
# Aggregation is always within each model (keeping them distinct)

##   Combine and Aggregate Confidence Rating Distribution   ##
Preds_RatingDist_corr_cond <- rbind(preds_WEV, preds_RM) %>%
  group_by(model, rating, correct, condition) %>%
  summarise(p = mean(p)) %>%
  mutate(model = factor(model, levels = c("WEVmu", "2DSD", "IRMt", "IRM", "PCRMt", "PCRM"),
                        labels = c("WEVmu", "2DSD", "IRMt", "IRM", "PCRMt", "PCRM")), 
         condition = factor(condition,levels = 1:nConds, 
                            labels = paste(c("8.3", "16.7", "33.3", "66.7", "133.3"), "")))


# # # Sanity checks:
# Preds_RatingDist_corr_cond %>% filter(model %in% c("IRM", "IRMt", "PCRM", "PCRMt")) %>%
#   group_by(model, condition) %>% summarise(p = sum(p))
# Preds_RatingDist_corr_cond %>% group_by(model, condition) %>%
#   summarise(p=sum(p))


##    Compute Mean Rating Accross Conditions 
# This is good to visualize a folded-X- and double-increase-pattern 
Preds_MRating_corr_cond <- Preds_RatingDist_corr_cond %>% group_by(model, condition, correct) %>%
  summarise(MRating = sum(p*rating)/sum(p))

##   Compute ROC curves for model predictions 
Preds_ROC_cond <- rbind(preds_RM, preds_WEV) %>%
  mutate(model = factor(model, levels = c("WEVmu", "2DSD", "IRMt", "IRM", "PCRMt", "PCRM"),
                        labels = c("WEVmu", "2DSD", "IRMt", "IRM", "PCRMt", "PCRM"))) %>%
  group_by(model, rating, correct, response, condition, stimulus) %>%
  summarise(p = mean(p)) %>% group_by(model) %>%
  summarise(getROCcoordFromPred(cur_data())) %>%
  mutate(HitRate = qnorm(HitRate), FalseAlarmRate = qnorm(FalseAlarmRate),
         condition = factor(condition, levels = 1:5, labels = paste(c("8.3", "16.7", "33.3", "66.7", "133.3"), "")))


##   Combine and Aggregate RT densities 
Preds_RTdens_corr_cond_rating <- rbind(RT_dist_RM, RT_dist_WEV) %>% ungroup() %>%
  select(-participant) %>%
  group_by(rating, condition, model, correct, rt) %>%
  summarise(dens = mean(dens), densscaled = mean(densscaled)) %>%
  mutate(model = factor(model, levels = c("WEVmu", "2DSD", "IRMt", "IRM", "PCRMt", "PCRM"),
                      labels = c("WEVmu", "2DSD", "IRMt", "IRM", "PCRMt", "PCRM")), 
         condition = factor(condition,levels = 1:nConds, 
                            labels = paste(c("8.3", "16.7", "33.3", "66.7", "133.3"), "")))  


####    Computation of RT quantiles    
Preds_RTQuants_corr_cond_rating <- Preds_RTdens_corr_cond_rating %>% select(-"densscaled") %>%
  RTDensityToQuantiles(c(.1,.3,.5,.7,.9))
Preds_RTQuants_cond <- Preds_RTdens_corr_cond_rating %>% group_by(model, rt, condition) %>%
  summarise(dens = sum(dens))%>%
  mutate(model = factor(model, levels=c("WEVmu", "2DSD", "IRMt", "IRM", "PCRMt", "PCRM"),
                        labels=c("WEVmu", "2DSD", "IRMt", "IRM", "PCRMt", "PCRM"))) %>%
  RTDensityToQuantiles(.) %>% 
  left_join(summarise(group_by(Preds_RatingDist_corr_cond, 
                               model, condition), 
                      p_correct=sum(p)))
Preds_RTQuants_corr_cond <- Preds_RTdens_corr_cond_rating %>% group_by(model, rt, correct, condition) %>%
  summarise(dens = sum(dens))%>%
  mutate(model = factor(model, levels=c("WEVmu", "2DSD", "IRMt", "IRM", "PCRMt", "PCRM"),
                        labels=c("WEVmu", "2DSD", "IRMt", "IRM", "PCRMt", "PCRM"))) %>% 
  group_by(model) %>% 
  RTDensityToQuantiles(.) %>% 
  left_join(summarise(group_by(Preds_RatingDist_corr_cond, 
                               model, condition, correct), 
                      p_correct=sum(p)))
Preds_RTQuants_rating <- Preds_RTdens_corr_cond_rating %>% group_by(model, rt, rating) %>%
  summarise(dens = sum(dens)/nConds)%>%
  mutate(model = factor(model, levels=c("WEVmu", "2DSD", "IRMt", "IRM", "PCRMt", "PCRM"),
                        labels=c("WEVmu", "2DSD", "IRMt", "IRM", "PCRMt", "PCRM"))) %>%
  RTDensityToQuantiles(.) %>% 
  left_join(summarise(group_by(Preds_RatingDist_corr_cond, 
                               model, rating), 
                      p_correct=sum(p)/nConds)) %>%
  mutate(rating = factor(rating,
                         levels=1:nRatings, 
                         labels = paste0(c("Guessing", "Very Unsure", "Unsure", "Rather Sure", "Sure"), "")))
Preds_RTQuants_corr_rating <- Preds_RTdens_corr_cond_rating %>% group_by(model, rt, correct, rating) %>%
  summarise(dens = mean(dens))%>%
  mutate(model = factor(model, levels=c("WEVmu", "2DSD", "IRMt", "IRM", "PCRMt", "PCRM"),
                        labels=c("WEVmu", "2DSD", "IRMt", "IRM", "PCRMt", "PCRM"))) %>%
  group_by(model) %>% 
  RTDensityToQuantiles(.) %>% 
  left_join(summarise(group_by(Preds_RatingDist_corr_cond, 
                               model, rating, correct), 
                      p_correct=mean(p)))

#=============================================================
#####          Clean Up and Save Results           ###########
gc()

save(file="collected_fitsNpredicts.RData", 
     Data , nConds, nRatings,  maxrt, cond_levels,
     Data_RatingDist_part ,  Data_RatingDist_corr_cond ,  Data_MRating_corr_cond ,  Data_ROC_cond , 
     fits_WEVmodels ,  fits_RMmodels ,
     preds_RM ,  preds_WEV ,
     Preds_RatingDist_corr_cond ,  Preds_MRating_corr_cond ,  Preds_ROC_cond ,
     RT_dist_RM ,  RT_dist_WEV ,  RT_dist ,   Preds_RTdens_corr_cond_rating  , 
     Preds_RTQuants_corr_cond_rating, Preds_RTQuants_corr_cond, Preds_RTQuants_cond,
     Preds_RTQuants_corr_rating, Preds_RTQuants_rating,
     Data_RTQuants_corr_cond, Data_RTQuants_cond,
     Data_RTQuants_corr_rating, Data_RTQuants_rating)

###############################################################

#=============================================================
#=============================================================
##### 3) Analysis of Criteria and Computation of BayesFactors  #####
##### Figure XYZ: Frequency of differences w.r.t. BIC     ##########
plotBIC_data <- rbind(fits_RMmodels[,c("participant", "model", "BIC")],
                      fits_WEVmodels[,c("participant", "model", "BIC")]) %>%
  filter(model !="WEVmu") %>%
  left_join(rename(fits_WEVmodels[fits_WEVmodels$model=="WEVmu" ,c("participant", "BIC")], BICWEV = BIC)) %>%
  mutate(model= factor(model, levels=c("PCRMt", "PCRM", "IRMt", "IRM",  "2DSD")),
         deltaBIC = BIC-BICWEV,
         binnedDeltaBIC = as.integer(cut(deltaBIC, breaks=c(-Inf,-100, -10,-2,2, 10, 100, Inf))))
table(plotBIC_data$binnedDeltaBIC, plotBIC_data$model)

Plot_deltaBIC_Confidence_2 <- 
  ggplot(plotBIC_data, 
         aes(x=factor(binnedDeltaBIC, levels=1:6, labels = 1:6), 
             fill= factor(binnedDeltaBIC, levels=1:7, labels = 1:7))) + 
  facet_wrap( ~ model, ncol=3, dir = "v",  scales = "fixed", as.table = F) + 
  geom_bar(colour = "black") + 
  scale_y_continuous(name = "Number of participants") + 
  scale_x_discrete(drop=FALSE, 
    name = expression(Delta~BIC~vs.~dynWEV~model))  + 
  scale_fill_brewer(
    type = "div", palette=5, direction = -1, 
    labels = c(expression(paste(Delta~BIC, " < ", "\u2212",100)),
               expression(paste("\u2212",100," < ", Delta~BIC, " < ", "\u2212",10)), 
               expression(paste("\u2212",10," < ", Delta~BIC, " < ", "\u2212",2)), 
               expression(paste("\u2212",2," < ", Delta~BIC, " < ", 2)),
               expression(paste(2," < ", Delta~BIC, " < ", 10)), 
               expression(paste(10, " < ", Delta~BIC, " < ", 100)), 
               expression(paste(Delta~BIC, " > ", 100))), 
    name = expression(Delta~BIC), drop=FALSE) + 
  ggtitle("Bayes information criterion") + 
  theme_bw() + 
  theme(text = element_text(size=10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.title.y = element_text( angle=90),# switch off minor gridlines
        legend.text.align = 0, 
        panel.grid.minor = element_blank(),
        #  plot.margin = unit(c(.5,.5,.5,.5), "lines"),
        panel.spacing = unit(.1, "lines"), 
        legend.position = c(.125, .25), 
        legend.title = element_blank())
Plot_deltaBIC_Confidence_2
dir.create("figures")
ggsave(file = "figures/BICs.png",
       dpi=600,
       height = 10, width = 14, 
       units = "cm")
##### Bayes t-test between WEVmu and all other models ########
compute_bf <- function(x,y) {
  drop <- is.na(x) | is.na(y)
  x <- x[!drop]
  y <- y[!drop]
  bf <- ttestBF(x,y,  paired = TRUE, rscale=1)
  res <-  data.frame(bf = as.numeric(extractBF(bf)["bf"]), Lower = NA, Upper = NA)
  if (! is.na(res$bf)) {
    posterior <- quantile(posterior(bf, iterations=100000)[,"delta"], probs = c(.025, .975))
    res$Lower <- posterior[1]
    res$Upper <- posterior[2]
  }
  res
}

criteria_comparison <- rbind(select(fits_WEVmodels, c("participant","model","BIC", "AICc", "AIC")),
                             select(fits_RMmodels, c("participant","model","BIC", "AICc", "AIC"))) %>% 
  pivot_longer(cols=c("BIC", "AICc", "AIC"), names_to="criteria", values_to="value") %>%
  #pivot_wider(id_cols = c("participant", "criteria"), names_from = model, values_from = c("value")) %>%
  filter(model != "WEVmu") %>%
  group_by(criteria, model) %>%
  summarise(compute_bf(x=value, 
                       y=subset(fits_WEVmodels, model =="WEVmu")[,cur_data_all()$criteria[1]]),
            M_diff=mean(value-subset(fits_WEVmodels, model =="WEVmu")[,cur_data_all()$criteria[1]], na.rm=TRUE)) %>%
  rename(BF_model_vs_WEVmu = bf)
criteria_comparison <- arrange(criteria_comparison, model, criteria)
criteria_comparison
knitr::kable(criteria_comparison, digits=3)




#=============================================================
#=============================================================
##### 4) Visualisation of observations and predictions  ######
#=============================================================
######   Plots of Accuracy and Rating Distributions    #######
##### Figure XYZ: Plot of Mean Ratings     ###################
pd <- position_dodge(0.05)
p_MRating <- ggplot(mutate(Data_MRating_corr_cond, condition = as.integer(condition)),
                    aes(x=condition, y=MRating, group = as.factor(correct), shape=as.factor(correct))) +
  geom_line(data=Preds_MRating_corr_cond, aes(color=as.factor(correct)), size=1)+
  geom_point(position = pd)+
  geom_line(linetype="dashed", alpha=0.5,position = pd)+
  geom_errorbar(aes(ymin=MRating-SER, ymax=MRating+SER), colour="black", width=.1,position =pd) +
  facet_wrap(.~model, nrow = 2, dir="v")+
  ylab("Mean Confidence Rating")+
  scale_x_discrete(name="Stimulus-onset-asynchrony [ms]")+  
  scale_color_manual(values= c("red", "green4"),
                     name = "Correctness",
                     labels=c("Wrong", "Correct")) +
  scale_shape_manual(values=c(24,21),
                     name = "Correctness",
                     labels=c("Wrong", "Correct")) +
  theme_bw() + 
  theme(text = element_text(size=10),
        axis.title.y = element_text( angle=90),
        panel.grid.minor = element_blank(),  # switch off minor gridlines
        panel.grid.major = element_blank(),
        strip.text = element_text(size=8),
        legend.position = c(.005, .995), legend.justification = c(0,1),
        legend.title = element_blank(), 
        legend.key = element_blank(),
        legend.key.width=unit(2.5,"line"))
p_MRating
ggsave("figures/meanRatingmodels.png", width = 19, height=13, units="cm",dpi=900)



##### Figure XYZ: Plot of Ratings Dist.    ###################
p_ratingdist <- ggplot(data=Data_RatingDist_corr_cond, aes(x=rating, y=p))+
  geom_bar( aes(fill=as.factor(correct)), stat = "identity", show.legend = FALSE)+
  geom_point(data=Preds_RatingDist_corr_cond,aes(shape=model), position=position_dodge(1), col="black")+
  facet_grid(cols=vars(correct), rows = vars(condition), 
             labeller = labeller(correct=c("0"="Wrong", "1"="Correct"), 
                                 condition=function(x) paste(x, "ms")))+
  scale_x_continuous(
    name = "Identification confidence [% scale width]", 
    breaks = 1:5, 
    labels = c("0-20","20-40", "40-60", "60-80", "80-100")) +
  scale_y_continuous(name = "Probability", breaks = c(0, 0.4, 0.8))+ 
  #ggtitle("Confidence Distribution: Data vs. Predictions; Fitted per participants (aggregated for plotting)") +
  # scale_fill_brewer(type = "qual", palette=6, 
  #                   name = "Model prediction" ) +
  scale_shape_discrete(name = "Model prediction") +
  theme_bw() + 
  theme(text = element_text(size=10),
        axis.text = element_text(color = "black"), 
        panel.grid.minor = element_blank(),  # switch off minor gridlines
        legend.key = element_blank(),
        strip.text.y = element_text(angle = 0),
        panel.spacing = unit(.1, "lines"))
p_ratingdist
ggsave("figures/distributionRatingmodels.png", width = 21, height=13, units="cm",dpi=600)


#=============================================================
############## Visualization of ROC curves   #################
ggplot(ungroup(Preds_ROC_cond), aes(x = FalseAlarmRate, y = HitRate)) +
  geom_line(aes(col = "Data")) +
  #geom_point(col="red", shape=2)+
  geom_point(data = ungroup(Data_ROC_cond), aes(col = "Model \nFits")) +
  facet_grid(cols = vars(condition), rows = vars(model),
             labeller = labeller(condition=function(x) paste(x, "ms"))) +
  scale_color_manual(values = c("red", "black")) +
  theme_bw() + 
  theme(text = element_text(size = 10),
        axis.title.y = element_text( angle = 90),
        panel.grid.minor = element_blank(),  # switch off minor gridlines
        #panel.grid.major = element_blank(),
        strip.text = element_text(size = 8),
        #legend.position = c(.995, .005), legend.justification = c(1,0),
        legend.title = element_blank(), 
        legend.key = element_blank(),
        legend.key.width = unit(2.5,"line"))
ggsave("figures/ROCfits.png", width = 27, height = 19, units = "cm",dpi = 900)






#=============================================================
########  Plots of RT Densities    ####################
### Comparison of overall RTs
# left_join(Preds_RTdens_corr_cond_rating, Preds_RatingDist_corr_cond) %>% 
#   group_by(model, rt, correct, condition) %>%

plot_RT_pred <- Preds_RTdens_corr_cond_rating %>% group_by(model, rt, correct) %>%
  summarise(dens = sum(dens)/nConds)%>%
  mutate(model = factor(model, levels=c("WEVmu", "2DSD", "IRMt", "IRM", "PCRMt", "PCRM"),
                        labels=c("WEVmu", "2DSD", "IRMt", "IRM", "PCRMt", "PCRM")))
plot_RT_pred <- left_join(plot_RT_pred, plot_RT_pred %>% group_by(model, correct) %>% 
                            summarise(p = sum(dens)*(rt[5]-rt[4])) %>% 
                            group_by(model) %>%
                            summarise(p=sum(p))) %>%
  mutate(dens=dens/p)

RTplot_overall <- ggplot()+
  geom_density(data=Data, 
                 aes(x=rt, y=..count../nrow(Data), 
                     col=as.factor(correct), group=correct, 
                     linetype="Data"), col="black")+  
  geom_line(data = plot_RT_pred, 
            aes(x=rt,y=dens, col=as.factor(correct), linetype="Model fit"),
            size=1.2)+
  xlim(c(0, maxrt))+ xlab("Reaction Time [s]") +
  scale_linetype_discrete(name="")+
  scale_color_manual(name="", breaks=c(0,1), 
                     values=RColorBrewer::brewer.pal(3, "RdBu")[c(1,3)], labels=c("Wrong", "Right"))+
  facet_wrap(.~model, nrow = 2, dir="v")
RTplot_overall
ggsave("figures/distributionRToverall.png", width = 22, height=13, units="cm",dpi=600)

  

Data <- Data %>% mutate(col_dens = rating * (-1)^(1-correct))
plot_maxt <- 4.3
plotRTpred <- Preds_RTdens_corr_cond_rating  %>% 
  mutate(col_dens = rating * (-1)^(1 - correct)) %>%
  group_by(model, rt, correct, rating, col_dens) %>%
  summarise(dens=mean(dens), densscaled=mean(densscaled))

pal <- RColorBrewer::brewer.pal(10, "RdBu")
pal[1:5] <- rev(pal[1:5])
col_breaks <- c(paste(c(1:5, 1:5), rep(0:1, each=5), sep="."))
type = "seq"
{
p1 <- ggplot()+
  geom_density(data=cbind(Data, strip="Observed data"), aes(x=rt, fill=interaction(as.factor(rating), as.factor(correct)), group=col_dens, y=..count../(nrow(Data))),
               position ="stack", bw=0.07)+
  xlim(c(0,plot_maxt))+ scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1, 1.25))+
  xlab("Reaction time [s]")+ylab("Density")+
  theme_bw()+
  theme(legend.position = "none",        
        panel.grid.minor = element_blank(),  # switch off minor gridlines
        panel.grid.major = element_blank())+
  scale_fill_manual(breaks=col_breaks,
                    values = pal)+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  facet_wrap(.~strip)
p2 <- ggplot()+
  geom_density(data=Data, aes(x=rt,  fill=interaction(as.factor(rating), as.factor(correct)),group=col_dens, y=..count../(nrow(Data))),
               position ="fill", bw=0.1)+
  theme_bw()+
  xlim(c(0,plot_maxt))+
  xlab("Reaction time [s]")+ylab("Relative proportion")+
  scale_fill_manual(breaks=col_breaks,
                    values = pal,labels=c("Unsure", "","", "", "Sure", rep("", 5)),
                    guide=guide_legend(title="Wrong\n\nCorrect", byrow=TRUE, 
                                       label=TRUE,label.position = "top",
                                       label.theme = element_text()))+
  theme(legend.position = "bottom",        
        panel.grid.minor = element_blank(),  # switch off minor gridlines
        panel.grid.major = element_blank(),
        legend.spacing.x = unit(0, 'cm'),
        legend.spacing.y = unit(0, 'cm'))
p3 <- ggplot(plotRTpred)+
  geom_area(aes(x=rt, y=dens, fill=interaction( as.factor(rating), as.factor(correct)), group=col_dens),
            col="black",
            position = "fill", stat = "identity")+
  theme_bw()+
  xlim(c(0,plot_maxt))+
  facet_wrap(.~model, nrow=2,dir="v")+
  scale_fill_manual(breaks=col_breaks,
                    values = pal)+
  theme(legend.position = "none",        
        panel.grid.minor = element_blank(),  # switch off minor gridlines
        panel.grid.major = element_blank())+
  xlab("Reaction time [s]")+ylab("Relative Proportion")

}
grid.arrange(p1,p2,p3, widths=c(0.3, 0.7), heights=c(0.35, 0.65),
             layout_matrix=rbind(c(1,3),
                                 c(2,3)))
ggsave(plot = grid.arrange(p1,p2,p3, widths=c(0.3, 0.7), heights=c(0.35, 0.65),
                                      layout_matrix=rbind(c(1,3),
                                                          c(2,3))),
       file = "figures/relativeRatingRT1.png",
       width = 21, height=13, units="cm",dpi=600)


plotRTpred <- Preds_RTdens_corr_cond_rating %>%
  group_by(model, rt, correct, condition) %>%
  summarise(dens=mean(dens), densscaled=mean(densscaled)) %>%
  mutate(condition = as.numeric(condition),
         col_dens = condition * (-1)^(1-correct))
Data_RTplot <- Data %>% mutate(condition = as.numeric(condition),
                               col_dens = condition * (-1)^(1-correct))
pal <- RColorBrewer::brewer.pal(10, "RdYlBu")
pal[1:5] <- rev(pal[1:5])
col_breaks <- c(paste(c(1:5, 1:5), rep(0:1, each=5), sep="."))
type = "seq"
{
  p1 <- ggplot()+
    geom_density(data=cbind(Data_RTplot, strip="Observed data"), aes(x=rt, fill=interaction(as.factor(condition), as.factor(correct)), 
                                                              group=col_dens, y=..count../(nrow(Data))),
                 position ="stack", bw=0.07)+
    xlim(c(0,plot_maxt))+scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1, 1.25))+
    xlab("Reaction time [s]")+ylab("Density")+
    theme_bw()+
    theme(legend.position = "none",        
          panel.grid.minor = element_blank(),  # switch off minor gridlines
          panel.grid.major = element_blank())+
    scale_fill_manual(breaks=col_breaks,
                      values = pal)+
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank())+
    facet_wrap(.~strip)
  p2 <- ggplot()+
    geom_density(data=Data_RTplot, aes(x=rt,  fill=interaction(as.factor(condition), as.factor(correct)),
                                group=col_dens, y=..count../(nrow(Data))),
                 position ="fill", bw=0.13)+
    theme_bw()+
    xlim(c(0,plot_maxt))+
    xlab("Reaction time [s]")+ylab("Relative proportion")+
    scale_fill_manual(breaks=col_breaks,name="SOA [ms]", 
                      values = pal,labels=c("8.3", "16.7", "33.3", "66.7", "133.3", rep("", 5)),
                      guide=guide_legend(title="Wrong\n\nCorrect", byrow=TRUE, 
                                         label=TRUE,label.position = "top",
                                         label.theme = element_text()))+
    theme(legend.position = "bottom",        
          panel.grid.minor = element_blank(),  # switch off minor gridlines
          panel.grid.major = element_blank(),
          legend.spacing.x = unit(0, 'cm'),
          legend.spacing.y = unit(0, 'cm'))
  p3 <- ggplot(plotRTpred)+
    geom_area(aes(x=rt, y=dens, fill=interaction( as.factor(condition), as.factor(correct)), 
                  group=col_dens),
              col="black",
              position = "fill", stat = "identity")+
    theme_bw()+
    xlim(c(0,plot_maxt))+
    facet_wrap(.~model, nrow=2,dir="v")+
    scale_fill_manual(breaks=col_breaks,
                      values = pal)+
    theme(legend.position = "none",        
          panel.grid.minor = element_blank(),  # switch off minor gridlines
          panel.grid.major = element_blank())+
    xlab("Reaction time [s]")+ylab("Relative proportion")
}
grid.arrange(p1,p2,p3, widths=c(0.3, 0.7),heights=c(0.35, 0.65),
             layout_matrix=rbind(c(1,3),
                                 c(2,3)))
ggsave(plot = grid.arrange(p1,p2,p3, widths=c(0.3, 0.7), heights=c(0.35, 0.65),
                                      layout_matrix=rbind(c(1,3),
                                                          c(2,3))),
       file = "figures/relativeCorrCondRT.png",
       width = 21, height=13, units="cm",dpi=600)




plotRTpred <- Preds_RTdens_corr_cond_rating %>%
  group_by(model, rt, condition) %>%
  summarise(dens=mean(dens), densscaled=mean(densscaled)) %>%
  mutate(condition=as.numeric(condition))
pal <- RColorBrewer::brewer.pal(5, "YlGn")
col_breaks <- 1:5
type = "seq"
{
  p1 <- ggplot()+
    geom_density(data=cbind(Data_RTplot, strip="Observed data"), aes(x=rt, fill=as.factor(condition), 
                                                              group=as.factor(condition), y=..count../(nrow(Data))),
                 position ="stack", bw=0.07)+
    xlim(c(0,plot_maxt))+scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1, 1.25))+
    xlab("Reaction time [s]")+ylab("Density")+
    theme_bw()+
    theme(legend.position = "none",        
          panel.grid.minor = element_blank(),  # switch off minor gridlines
          panel.grid.major = element_blank())+
    scale_fill_manual(breaks=col_breaks,
                      values = pal)+
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank())+
    facet_wrap(.~strip)
  p2 <- ggplot()+
    geom_density(data=Data_RTplot, aes(x=rt,  fill=as.factor(condition),
                                group=as.factor(condition), y=..count../(nrow(Data))),
                 position ="fill", bw=0.13)+
    theme_bw()+
    xlim(c(0,plot_maxt))+
    xlab("Reaction time [s]")+ylab("Relative proportion")+
    scale_fill_manual(breaks=col_breaks,
                      values = pal,labels=c("8.3", "16.7", "33.3", "66.7", "133.3"),
                      guide=guide_legend(title="SOA [ms]", byrow=TRUE, 
                                         label=TRUE,label.position = "top",
                                         label.theme = element_text(),
                                         title.position = "bottom"))+
    theme(legend.position = "bottom",        
          panel.grid.minor = element_blank(),  # switch off minor gridlines
          panel.grid.major = element_blank(),
          legend.spacing.x = unit(0.2, 'cm'),
          legend.key.width = unit(1, "cm"),
          legend.spacing.y = unit(0.1, 'cm'))
  p3 <- ggplot(plotRTpred)+
    geom_area(aes(x=rt, y=dens, fill=as.factor(condition), 
                  group=as.factor(condition)),
              col="black",
              position = "fill", stat = "identity")+
    theme_bw()+
    xlim(c(0,plot_maxt))+
    facet_wrap(.~model, nrow=2,dir="v")+
    scale_fill_manual(breaks=col_breaks,
                      values = pal)+
    theme(legend.position = "none",        
          panel.grid.minor = element_blank(),  # switch off minor gridlines
          panel.grid.major = element_blank())+
    xlab("Reaction time [s]")+ylab("Density")
}
grid.arrange(p1,p2,p3, widths=c(0.3, 0.7),heights=c(0.35, 0.65),
             layout_matrix=rbind(c(1,3),
                                 c(2,3)))
ggsave(plot = grid.arrange(p1,p2,p3, widths=c(0.3, 0.7), heights=c(0.35, 0.65),
                           layout_matrix=rbind(c(1,3),
                                               c(2,3))),
       file = "figures/relativeCondRT.png",
       width = 21, height=13, units="cm",dpi=600)




#=============================================================
##############   Plots of RT-Quantiles      ##################
##### Figure XYZ: RTQuantiles accross conditions    ##########
### Similar to Fig.5 in Ratcliff & Starns (2009)
ggplot()+
  geom_line(data=Data_RTQuants_cond, aes(x=condition, y=q, group=as.factor(p), linetype="Data"))+
  geom_point(data=Data_RTQuants_cond, aes(x=condition, y=q, shape="Data"),
             size=2, alpha=0.4)+
  geom_point( data=Preds_RTQuants_cond, aes(x=condition, y=q,shape="Model fit"),
              size=3, alpha=0.4) +
  geom_line(data=Preds_RTQuants_cond, aes(x=condition, y=q, group=as.factor(p), linetype="Model fit"))+
  scale_shape_manual(name="class", values=c(17, 2))+ scale_linetype_discrete(name="class")+
  scale_x_discrete(name="SOA [ms]")+
  scale_y_continuous(breaks = c(1,2,3), name="Reaction Time Quantiles [s]")+
  #scale_linetype_manual(values= c("solid", "dashed", "dotted"), breaks = c("Data", "WEVmu", "2DSD"))+
  facet_wrap(.~model, nrow=2,dir="v")+
  theme_bw() + 
  theme(text = element_text(size=10),
        axis.title.y = element_text( angle=90),
        panel.grid.minor = element_blank(),  # switch off minor gridlines
        #panel.grid.major = element_blank(),
        strip.text = element_text(size=8),
        legend.position = c(.005, .995), legend.justification = c(0,1),
        legend.title = element_blank(), 
        legend.key = element_blank(),
        legend.key.width=unit(2.5,"line"))
ggsave("figures/RTQuants_cond.png", 
       width=19, height=14, dpi=1000, units="cm")



##### Figure XYZ: RTQuantiles accross correct X conditions    ##########
ggplot()+
  geom_line(data=Data_RTQuants_corr_cond, aes(x=p_correct, y=q, group=as.factor(p), linetype="Data"))+
  geom_point(data=Data_RTQuants_corr_cond, aes(x=p_correct, y=q, shape="Data", col=as.factor(correct)),
             size=3, alpha=0.6)+
  geom_point( data=Preds_RTQuants_corr_cond, aes(x=p_correct, y=q,shape="Model fit", col=as.factor(correct)),
              size=3, alpha=0.6) +
  geom_line(data=Preds_RTQuants_corr_cond, aes(x=p_correct, y=q, group=as.factor(p), linetype="Model fit"))+
  scale_shape_manual(name="class", values=c(17, 2))+ scale_linetype_discrete(name="class")+
  scale_color_manual(name="Correctness", breaks=c(0,1), 
                     values=RColorBrewer::brewer.pal(3, "RdBu")[c(1,3)], labels=c("Wrong", "Right"))+
  scale_x_continuous(breaks = c(.25, .5, .75), name="Choice Probability")+
  scale_y_continuous(breaks = c(1,2,3), name="Reaction Time Quantiles [s]")+
  #scale_linetype_manual(values= c("solid", "dashed", "dotted"), breaks = c("Data", "WEVmu", "2DSD"))+
  facet_wrap(.~model, nrow=2,dir="v")+
  theme_bw() + 
  theme(text = element_text(size=10),
        axis.title.y = element_text( angle=90),
        panel.grid.minor = element_blank(),  # switch off minor gridlines
        #panel.grid.major = element_blank(),
        strip.text = element_text(size=8),
        #legend.box = "horizontal",
        #legend.position = c(.005, .995), legend.justification = c(0,1),
        legend.title = element_blank(), 
        legend.key = element_blank(),
        legend.key.width=unit(2.5,"line"))
ggsave("figures/RTQuants_corr_cond.png", 
       width=25, height=18, dpi=600, units="cm")




Data_RTQuants_corr_cond <- Data_RTQuants_corr_cond %>%
  mutate(condition = as.numeric(condition),
         X = factor(paste(condition, correct, sep="."),
                    levels = c(paste(c(5:1, 1:5), rep(0:1, each=5), sep=".")),
                    labels = c(paste(c(5:1, 1:5), rep(0:1, each=5), sep="."))))
Preds_RTQuants_corr_cond <- Preds_RTQuants_corr_cond %>%
  mutate(condition = as.numeric(condition),
         X = factor(paste(condition, correct, sep="."),
                    levels = c(paste(c(5:1, 1:5), rep(0:1, each=5), sep=".")),
                    labels = c(paste(c(5:1, 1:5), rep(0:1, each=5), sep="."))))
ggplot()+
  geom_line(data=Data_RTQuants_corr_cond, aes(x=X, y=q, group=as.factor(p), linetype="Data"))+
  geom_point(data=Data_RTQuants_corr_cond, aes(x=X, y=q, shape="Data", col=as.factor(correct)),
             size=2)+
  geom_point( data=Preds_RTQuants_corr_cond, aes(x=X, y=q,shape="Model fit", col=as.factor(correct)),
              size=3) +
  geom_line(data=Preds_RTQuants_corr_cond, aes(x=X, y=q, group=as.factor(p), linetype="Model fit"))+
  scale_shape_manual(name="class", values=c(17, 2))+ scale_linetype_discrete(name="class")+
  scale_color_manual(name="Correctness", breaks=c(0,1), 
                     values=RColorBrewer::brewer.pal(3, "RdBu")[c(1,3)], labels=c("Wrong", "Right"))+
  scale_x_discrete(name="Response x Condition")+
  scale_y_continuous(breaks = c(1,2,3), name="Reaction Time Quantiles [s]")+
  #scale_linetype_manual(values= c("solid", "dashed", "dotted"), breaks = c("Data", "WEVmu", "2DSD"))+
  facet_wrap(.~model, nrow=2,dir="v")+
  theme_bw() + 
  theme(text = element_text(size=10),
        axis.title.y = element_text( angle=90),
        panel.grid.minor = element_blank(),  # switch off minor gridlines
        #panel.grid.major = element_blank(),
        strip.text = element_text(size=8),
        #legend.box = "horizontal",
        #legend.position = c(.005, .995), legend.justification = c(0,1),
        legend.title = element_blank(), 
        legend.key = element_blank(),
        legend.key.width=unit(2.5,"line"))
ggsave("figures/RTQuants_corr_cond2.png", 
       width=25, height=18, dpi=600, units="cm")









##### Figure XYZ: RTQuantiles accross rating    ##########
ggplot()+
  geom_line(data=Data_RTQuants_rating, aes(x=rating, y=q, group=as.factor(p), linetype="Data"))+
  geom_point(data=Data_RTQuants_rating, aes(x=rating, y=q, shape="Data"),
             size=2, alpha=0.4)+
  geom_point( data=Preds_RTQuants_rating, aes(x=rating, y=q,shape="Model fit"),
              size=3, alpha=0.4) +
  geom_line(data=Preds_RTQuants_rating, aes(x=rating, y=q, group=as.factor(p), linetype="Model fit"))+
  scale_shape_manual(name="class", values=c(17, 2))+ scale_linetype_discrete(name="class")+
  scale_x_discrete(name="Confidence Rating")+
  scale_y_continuous(breaks = c(1,2,3), name="Reaction Time Quantiles [s]")+
  #scale_linetype_manual(values= c("solid", "dashed", "dotted"), breaks = c("Data", "WEVmu", "2DSD"))+
  facet_wrap(.~model, nrow=2,dir="v")+
  theme_bw() + 
  theme(text = element_text(size=10),
        axis.title.y = element_text( angle=90),
        panel.grid.minor = element_blank(),  # switch off minor gridlines
        #panel.grid.major = element_blank(),
        strip.text = element_text(size=8),
        legend.position = c(.005, .995), legend.justification = c(0,1),
        legend.title = element_blank(), 
        legend.key = element_blank(),
        legend.key.width=unit(2.5,"line"))
ggsave("figures/RTQuants_rating.png", 
       width=26, height=14, dpi=1000, units="cm")






##### Figure XYZ: RTQuantiles accross correct X rating    ##########
ggplot()+
  geom_line(data=Data_RTQuants_corr_rating, aes(x=p_correct, y=q, group=as.factor(p), linetype="Data"))+
  geom_point(data=Data_RTQuants_corr_rating, aes(x=p_correct, y=q, shape="Data", col=as.factor(correct)),
             size=3, alpha=0.6)+
  geom_point( data=Preds_RTQuants_corr_rating, aes(x=p_correct, y=q,shape="Model fit", col=as.factor(correct)),
              size=3, alpha=0.6) +
  geom_line(data=Preds_RTQuants_corr_rating, aes(x=p_correct, y=q, group=as.factor(p), linetype="Model fit"))+
  scale_shape_manual(name="class", values=c(17, 2))+ scale_linetype_discrete(name="class")+
  scale_color_manual(name="Correctness", breaks=c(0,1), 
                     values=RColorBrewer::brewer.pal(3, "RdBu")[c(1,3)], labels=c("Wrong", "Right"))+
  scale_x_continuous(breaks = c(0, .1,.2, .3), name="Choice Probability")+
  scale_y_continuous(breaks = c(1,2,3), name="Reaction Time Quantiles [s]")+
  #scale_linetype_manual(values= c("solid", "dashed", "dotted"), breaks = c("Data", "WEVmu", "2DSD"))+
  facet_wrap(.~model, nrow=2,dir="v")+
  theme_bw() + 
  theme(text = element_text(size=10),
        axis.title.y = element_text( angle=90),
        panel.grid.minor = element_blank(),  # switch off minor gridlines
        #panel.grid.major = element_blank(),
        strip.text = element_text(size=8),
        #legend.box = "horizontal",
        #legend.position = c(.005, .995), legend.justification = c(0,1),
        legend.title = element_blank(), 
        legend.key = element_blank(),
        legend.key.width=unit(2.5,"line"))
ggsave("figures/RTQuants_corr_rating.png", 
       width=25, height=18, dpi=600, units="cm")


Data_RTQuants_corr_rating <- Data_RTQuants_corr_rating %>%
  mutate(X = factor(paste(rating, correct, sep="."),
                    levels = c(paste(c(5:1, 1:5), rep(0:1, each=5), sep=".")),
                    labels = c(paste(c(5:1, 1:5), rep(0:1, each=5), sep="."))))
Preds_RTQuants_corr_rating <- Preds_RTQuants_corr_rating %>%
  mutate(X = factor(paste(rating, correct, sep="."),
                    levels = c(paste(c(5:1, 1:5), rep(0:1, each=5), sep=".")),
                    labels = c(paste(c(5:1, 1:5), rep(0:1, each=5), sep="."))))
ggplot()+
  geom_line(data=Data_RTQuants_corr_rating, aes(x=X, y=q, group=as.factor(p), linetype="Data"))+
  geom_point(data=Data_RTQuants_corr_rating, aes(x=X, y=q, shape="Data", col=as.factor(correct)),
             size=2)+
  geom_point( data=Preds_RTQuants_corr_rating, aes(x=X, y=q,shape="Model fit", col=as.factor(correct)),
              size=3) +
  geom_line(data=Preds_RTQuants_corr_rating, aes(x=X, y=q, group=as.factor(p), linetype="Model fit"))+
  scale_shape_manual(name="class", values=c(17, 2))+ scale_linetype_discrete(name="class")+
  scale_color_manual(name="Correctness", breaks=c(0,1), 
                     values=RColorBrewer::brewer.pal(3, "RdBu")[c(1,3)], labels=c("Wrong", "Right"))+
  scale_x_discrete(name="Rating x Response", labels=c("Sure", "", "", "", "         Unsure", "", "", "", "", "Sure"))+
  scale_y_continuous(breaks = c(1,2,3), name="Reaction Time Quantiles [s]")+
  #scale_linetype_manual(values= c("solid", "dashed", "dotted"), breaks = c("Data", "WEVmu", "2DSD"))+
  facet_wrap(.~model, nrow=2,dir="v")+
  theme_bw() + 
  theme(text = element_text(size=10),
        axis.title.y = element_text( angle=90),
        panel.grid.minor = element_blank(),  # switch off minor gridlines
        #panel.grid.major = element_blank(),
        strip.text = element_text(size=8),
        #legend.box = "horizontal",
        #legend.position = c(.005, .995), legend.justification = c(0,1),
        legend.title = element_blank(), 
        legend.key = element_blank(),
        legend.key.width=unit(2.5,"line"))
ggsave("figures/RTQuants_corr_rating2.png", 
       width=25, height=18, dpi=600, units="cm")

## Maybe interesting:
ggplot()+
  geom_line(data=Data_RTQuants_corr_rating, aes(x=X, y=q, group=as.factor(p), linetype="Data"))+
  geom_point(data=Data_RTQuants_corr_rating, aes(x=X, y=q, shape="Data", size=p_correct), alpha=0.3)+
  geom_point( data=Preds_RTQuants_corr_rating, aes(x=X, y=q,shape="Model fit", size=p_correct)) +
  geom_line(data=Preds_RTQuants_corr_rating, aes(x=X, y=q, group=as.factor(p), linetype="Model fit"))+
  scale_shape_manual(name="class", values=c(17, 2))+ scale_linetype_discrete(name="class")+
  scale_color_manual(name="Correctness", breaks=c(0,1), 
                     values=RColorBrewer::brewer.pal(3, "RdBu")[c(1,3)], labels=c("Wrong", "Right"))+
  scale_x_discrete(name="Rating x Response", labels=c("Wrong/Sure", "", "", "", "         Unsure", "", "", "", "", "Right/Sure"))+
  # scale_x_continuous(name="Correctness and Confidence", 
  #                    breaks=c(1, 3.5, 6), minor_breaks = 1:6, 
  #                    labels = c("Falsch/Hohe Konfidenz", "Niedrige Konfidenz", "Richtig/Hohe Konfidenz"))+
  scale_y_continuous(breaks = c(1,2,3), name="Reaction Time Quantiles [s]")+
  #scale_linetype_manual(values= c("solid", "dashed", "dotted"), breaks = c("Data", "WEVmu", "2DSD"))+
  facet_wrap(.~model, nrow=2,dir="v")+
  theme_bw() + 
  theme(text = element_text(size=10),
        axis.title.y = element_text( angle=90),
        panel.grid.minor = element_blank(),  # switch off minor gridlines
        #panel.grid.major = element_blank(),
        strip.text = element_text(size=8),
        #legend.box = "horizontal",
        #legend.position = c(.005, .995), legend.justification = c(0,1),
        legend.title = element_blank(), 
        legend.key = element_blank(),
        legend.key.width=unit(2.5,"line"))
ggsave("figures/RTQuants_corr_rating3.png", 
       width=25, height=18, dpi=600, units="cm")


###### REST ##################









#save.image("DynamicalConfidenceResults.RData")
# load("DynamicalConfidenceResults.RData")


#rm(list=ls())
