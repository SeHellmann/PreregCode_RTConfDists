getROCcoordFromData <- 
  function(rating, stimulus, correct){
    if (is.numeric(rating)) {
      if (all(rating == floor(rating))) {
        if (min(rating) == 0) {
          Rlevels = 0:max(rating)
        } else {
          Rlevels = 1:max(rating)
        }
        rating <- factor(rating, levels = Rlevels)
      }
    }
    
    if (!is.integer(rating)) {
      if (!is.factor(rating)) {
        rating <- as.factor(rating)
        warning("rating converted to factor!")
      } 
    }
    nRatings <- length(levels(rating))
    
    if (!is.factor(stimulus)) {
       warning("stimulus converted to factor!")
       stimulus <- as.factor(stimulus)
    }
    if (length(levels(stimulus)) != 2) stop("stimulus should have 2 levels")
    if (!all(correct %in% c(0,1))) stop("correct should be 1 or 0")
    
    A <- levels(stimulus)[1] 
    B <- levels(stimulus)[2]

    
    N_SA_RA <- table(rating[stimulus == A & correct == 1]) ##+ .001
    N_SA_RB <- table(rating[stimulus == A & correct == 0]) #+ .001
    N_SB_RA <- table(rating[stimulus == B & correct == 0]) #+ 001
    N_SB_RB <- table(rating[stimulus == B & correct == 1]) #+ 001
    
    HR <- cumsum(c(rev(N_SB_RB),N_SB_RA))/sum(c(N_SB_RB,N_SB_RA))
    FA <- cumsum(c(rev(N_SA_RB),N_SA_RA))/sum(c(N_SA_RB,N_SA_RA))
    HitRates <- HR[-length(HR)]
    FalseAlarms <- FA[-length(FA)]
    res <- data.frame(HitRate = as.vector(t(HitRates)), 
                      FalseAlarmRate = as.vector(t(FalseAlarms)),
                      criterion = 1:(2 * nRatings - 1))
    
    res
    
  }    

getROCcoordFromPred <- function(predDf){
  if (!"rating" %in% colnames(predDf)) {
    stop("there should be a column named rating!")
  }
  if (!"stimulus" %in% colnames(predDf)) {
    stop("there should be a column named stimulus!")
  }
  if (!"response" %in% colnames(predDf)) {
    stop("there should be a column named response!")
  }
  if (!"condition" %in% colnames(predDf)) {
    stop("there should be a column named condition!")
  }
  if (!"p" %in% colnames(predDf)) {
    stop("there should be a column named p")
  }
  # 
  # require(plyr)
  # 
  # res <- ddply(predDf, ~ condition, 
  #              function(df) {
  #                p_SA_RA <- df$p[df$stimulus == -1 & df$response == -1] ##+ .001
  #                p_SA_RB <- df$p[df$stimulus == -1 & df$response == 1] #+ .001
  #                p_SB_RA <- df$p[df$stimulus == 1 & df$response == -1] #+ 001
  #                p_SB_RB <- df$p[df$stimulus == 1 & df$response == 1]
  #                
  #                HR <- cumsum(c(rev(p_SB_RB),p_SB_RA))/sum(c(p_SB_RB,p_SB_RA))
  #                FA <- cumsum(c(rev(p_SA_RB),p_SA_RA))/sum(c(p_SA_RB,p_SA_RA))
  #                
  #                HitRates <- HR[-length(HR)]
  #                FalseAlarms <- FA[-length(FA)]
  #                res <- data.frame(HitRate = as.vector(t(HitRates)), 
  #                                  FalseAlarmRate = as.vector(t(FalseAlarms)))
  #              }, .progress = "text")
  # res
  # 
  if (2 %in% predDf$response) {
    predDf <- predDf %>% mutate(response = if_else(response == 2, -1, 1), 
                      stimulus = if_else(stimulus == 2, -1, 1))
  }
  # nRatings <- max(predDf$rating)
  # res <- predDf %>% arrange(condition, stimulus, response, rating) %>% group_by(condition, stimulus) %>%
  #   dplyr::summarise(positives = cumsum(c(rev(p[1:nRatings]), p[(nRatings+1):n()]))[-n()]) %>%
  #   pivot_wider(names_from=stimulus, values_from = positives) 
  
  res <- predDf %>% arrange(rating) %>% group_by(condition) %>%
    do({
        df <- .
        p_SA_RA <- df$p[df$stimulus == -1 & df$response == -1] ##+ .001
        p_SA_RB <- df$p[df$stimulus == -1 & df$response == 1] #+ .001
        p_SB_RA <- df$p[df$stimulus == 1 & df$response == -1] #+ 001
        p_SB_RB <- df$p[df$stimulus == 1 & df$response == 1]
        
        HR <- cumsum(c(rev(p_SB_RB),p_SB_RA))/sum(c(p_SB_RB,p_SB_RA))
        FA <- cumsum(c(rev(p_SA_RB),p_SA_RA))/sum(c(p_SA_RB,p_SA_RA))
        HitRates <- HR[-length(HR)]
        FalseAlarms <- FA[-length(FA)]
        res <- data.frame(HitRate = as.vector(t(HitRates)), 
                          FalseAlarmRate = as.vector(t(FalseAlarms)))
        res
      })
  res
  
}
