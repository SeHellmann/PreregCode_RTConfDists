### Gather data from data folder ####
library(tidyverse)
library(BayesFactor)

cat("\014") 

all_csv <- list.files(data_folder,
                      recursive = TRUE, 
                      full.names = TRUE, 
                      pattern = ".csv")
#file.info(all_csv)
Data <- data.frame()
pb <- txtProgressBar(style=3, min=0, max=length(all_csv))
for (i in 1:length(all_csv)) {
  oldw <- getOption("warn")
  options(warn=-1)
  temp <- read_csv(all_csv[i], col_types = cols())  
  options(warn=oldw)
  temp <- temp %>% filter(!is.na(blocks.thisRepN)) %>%
    select(participant, gender, age, SOA, Orientrierung, expectedAnswer, SenktOderWaag.keys,
           SenktOderWaag.corr, SenktOderWaag.rt, Rating1, Rating1RT) %>%
    rename(Orientierung=Orientrierung, 
           rt = SenktOderWaag.rt, 
           correct = SenktOderWaag.corr, 
           response = SenktOderWaag.keys,
           stimulus = expectedAnswer,
           condition = SOA) %>%
    mutate(rating = as.numeric(as.factor(cut(Rating1, breaks = c(-1, -.6, -.2, .2, .6, 1), include.lowest = TRUE))),
           rating3 = as.numeric(as.factor(cut(Rating1, breaks = c(-1, -.6, .6, 1), include.lowest = TRUE))))
  
  Data <- rbind(Data, temp)
  setTxtProgressBar(pb, i)
  
}
Data <- Data %>% group_by(participant)%>%
  mutate(participant = as.numeric(participant))
Ntotal <-  Data %>%
  summarise(Ntot = sum(n())) 


Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

CheckParticipants <- Data %>%
  summarise(
    nTrials = length(correct),
    Performance = mean(correct),
    AboveChance = as.numeric(try(extractBF(proportionBF(y=sum(correct), N=length(correct), p=1/2))$bf)),
    WhichModeRating1 = Mode(Rating1),
    NumModeRating1 = sum(Rating1  == Mode(Rating1))/length(Rating1))

BadSubjects <- CheckParticipants %>%
  filter(Performance < .5 | AboveChance < 3 | 
           NumModeRating1 > .9)

BadSubjects

Data <- subset(Data, !participant %in% BadSubjects$participant)
Data <- Data %>%
  filter(rt < mean(rt) + 4*sd(rt)  & rt > .3)


Data <- Data %>%
  filter(min(table(condition, stimulus)) >= 20 )

Nanalysis <-  left_join(summarise(Data, Nana = sum(n())), Ntotal)
cat(paste("Excluded Participants based on Accuracy:", paste(BadSubjects$participant, collapse = ", "), "\n"))
cat("Summary of Proportion of Excluded Trials per (non-excluded) Participant :\n")
print(summary(1- Nanalysis$Nana/Nanalysis$Ntot))


save(Data, BadSubjects, file="dataRausch2018.RData")