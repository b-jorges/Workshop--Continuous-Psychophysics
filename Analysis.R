#install.packages(c("dplyr", "ggplot2", "cowplot", "ggdist"))

require(dplyr)
require(ggplot2)
require(cowplot)
require(ggdist)
theme_set(theme_cowplot())

#load all datafiles (focussing on the .csv files)
flist <- list.files(paste0(dirname(rstudioapi::getSourceEditorContext()$path), "/Experiment/ContPsyTracking/data/"))[c(grep(".csv", list.files(paste0(dirname(rstudioapi::getSourceEditorContext()$path), "//Experiment/ContPsyTracking/data"), full.names = TRUE)))]

Responses_Wide = c()
Data2 = c()
j = 0
for (i in flist){
  j = j+1
  print(j)
  if (j == 1){
    Responses_Wide  = read.csv(paste0(dirname(rstudioapi::getSourceEditorContext()$path), "/Experiment/ContPsyTracking/data/",i), header = TRUE, row.names = NULL)
  } else {
    Data2  = read.csv(paste0(dirname(rstudioapi::getSourceEditorContext()$path), "/Experiment/ContPsyTracking/data/",i), header = TRUE, row.names = NULL)
  }
  Responses_Wide = rbind(Responses_Wide,Data2)
}

Responses_Wide = Responses_Wide %>%
  #transform recorded mouse data into same unit as displayed polygon (cm from center of screen)
                      mutate(mouse_2.x = as.numeric(substr(mouse_2.x[1], 2, 6)),
                             mouse_2.y = as.numeric(substr(mouse_2.y[1], 2, 6)),
                             mouse_x_cm = x_coord_mouse*(1/mouse_2.x),
                             mouse_y_cm = y_coord_mouse*(1/mouse_2.y),
                             
                             mouse_speed.x = (mouse_x_cm-lag(mouse_x_cm))/(Responses_Wide$time_in_run - 
                                                                             lag(Responses_Wide$time_in_run, 1)),
                             mouse_speed.y = (mouse_y_cm-lag(mouse_y_cm))/(Responses_Wide$time_in_run - 
                                                                             lag(Responses_Wide$time_in_run, 1)),
                             target_speed.x = (x_coord_target-lag(x_coord_target))/(Responses_Wide$time_in_run - 
                                                                             lag(Responses_Wide$time_in_run, 1)),
                             target_speed.y = (y_coord_target-lag(y_coord_target))/(Responses_Wide$time_in_run - 
                                                                             lag(Responses_Wide$time_in_run, 1))) %>%
  
  #filter out rows with no data (because PsychoPy .csv files be a little messy otherwise)
                      filter(!is.na(opacity)) %>% 
                      group_by(participant) %>%
  
  #get median duration of frame for each participant
                      mutate(FrameDuration = median(Responses_Wide$time_in_run - lag(Responses_Wide$time_in_run, 1), na.rm = TRUE))

##################
#some basic plots#
##################
#target versus mouse (x)
ggplot(Responses_Wide, aes(time_in_run, x_coord_target)) +
  geom_point() +
  geom_point(aes(time_in_run, mouse_x_cm), color = "red") +
  ylab("x position") +
  xlab("Time (s)") +
  facet_grid(participant~opacity) +
  theme(axis.title=element_text(size=18),
        axis.text = element_text(size=18),
        strip.text = element_text(size = 18))

#target versus mouse (y)
ggplot(Responses_Wide, aes(time_in_run, y_coord_target)) +
  geom_point() +
  geom_point(aes(time_in_run, mouse_y_cm), color = "red") +
  ylab("y position") +
  xlab("Time (s)") +
  facet_grid(participant~opacity) +
  theme(axis.title=element_text(size=18),
        axis.text = element_text(size=18),
        strip.text = element_text(size = 18))

###################
#Crosscorrelations#
###################

#how many frames in 1s
HowManyFrames = 1/median(Responses_Wide$FrameDuration)

CCG_Frame = data.frame()

#create copies of the dataframe where we lag target and response positions
#by 1 to HowManyFrames in steps of 4 (only 4 to keep size of this dataframe manageable)
for (i in seq(1,HowManyFrames,4)){
  CCG_Frame = rbind(CCG_Frame, Responses_Wide %>%
                      ungroup() %>%
                      select(participant,opacity,x_coord_target,mouse_x_cm,
                             y_coord_target,FrameDuration,mouse_y_cm, 
                             mouse_speed.x, mouse_speed.y, target_speed.x, target_speed.y) %>%
                      group_by(participant, opacity) %>%
                      mutate(x_coord_target_lagged = lag(x_coord_target, i),
                             y_coord_target_lagged = lag(y_coord_target, i),
                             
                             target_speed.x_lagged = lag(target_speed.x, i),
                             target_speed.y_lagged = lag(target_speed.y, i),
                             
                             lag = i))
}

CCG_Frame %>% filter(opacity == 0.05 & lag == 1)

CCG_Frame = CCG_Frame %>%
  group_by(participant, opacity, lag) %>%
  
  #Calculate the correlation between target and response positions for each lag
  #(separately per participant and condition and also separately for x/y directions)
  mutate(Correlation_x = cor.test(mouse_x_cm, x_coord_target_lagged)[4]$estimate[[1]],
         Correlation_y = cor.test(mouse_y_cm, y_coord_target_lagged)[4]$estimate[[1]],
         
         Correlation_x_Speed = cor.test(mouse_speed.x, target_speed.x_lagged)[4]$estimate[[1]],
         Correlation_y_Speed = cor.test(mouse_speed.y, target_speed.y_lagged)[4]$estimate[[1]],
         Correlation_Overall_Speed = (Correlation_x_Speed + Correlation_y_Speed)/2) %>% 
  group_by(participant, opacity) %>%
  
  #get the maximum correlation and the lag at which this maximum correlation is located
  mutate(MaxCorr_x = max(Correlation_x),
         Time_MaxCorr_x = lag[which.max(Correlation_x)],
         MaxCorr_y = max(Correlation_y),
         Time_MaxCorr_y = lag[which.max(Correlation_y)],
         
         MaxCorr_Speed_x = max(Correlation_x_Speed),
         Time_MaxCorr_Speed_x = lag[which.max(Correlation_x)],
         MaxCorr_Speed_y = max(Correlation_y_Speed),
         Time_MaxCorr_Speed_y = lag[which.max(Correlation_y)],
         
  #take the average maximum correlation and time lag across x and y directions
         MaxCorr_Overall = (MaxCorr_x + MaxCorr_y)/2,
         Time_MaxCorr_Overall = (Time_MaxCorr_x+Time_MaxCorr_y)/2,
  
         MaxCorr_Overall_Speed = max(Correlation_Overall_Speed),
         Time_MaxCorr_Overall_Speed = lag[which.max(Correlation_Overall_Speed)]
  )

# save(CCG_Frame, file = paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/SavedVariables/CCG_Frame.RData"))
# load(file = paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/SavedVariables/CCG_Frame.RData"))

#########################
#Cross Correlation Plots#
#########################
ggplot(CCG_Frame %>% 
         group_by(participant,opacity,lag) %>% 
         slice(1),
       aes(lag*FrameDuration,Correlation_x)) +
  geom_line(linewidth = 2) +
  xlab("Lag (s)") + 
  ylab("Correlation") +
  facet_grid(participant ~ opacity) +
  geom_vline(aes(xintercept = Time_MaxCorr_x*FrameDuration)) +
  geom_hline(aes(yintercept = MaxCorr_x), linetype = 2, linewidth = 1.5) +
  theme(axis.title=element_text(size=18),
        axis.text = element_text(size=18),
        strip.text = element_text(size = 18))

ggplot(CCG_Frame %>% 
         group_by(participant,opacity,lag) %>% 
         slice(1),
       aes(lag*FrameDuration,Correlation_y)) +
  geom_line(linewidth = 2) +
  xlab("Lag (s)") + 
  ylab("Correlation") +
  facet_grid(participant ~ opacity) +
  geom_vline(aes(xintercept = Time_MaxCorr_y*FrameDuration)) +
  geom_hline(aes(yintercept = MaxCorr_y), linetype = 2, linewidth = 1.5) +
  theme(axis.title=element_text(size=18),
        axis.text = element_text(size=18),
        strip.text = element_text(size = 18))

ggplot(CCG_Frame %>% 
         group_by(participant,opacity,lag) %>% 
         slice(1),
       aes(lag*FrameDuration,Correlation_Overall_Speed)) +
  geom_line(linewidth = 2) +
  xlab("Lag (s)") + 
  ylab("Correlation") +
  facet_grid(participant ~ opacity) +
  geom_vline(aes(xintercept = Time_MaxCorr_Overall_Speed*FrameDuration)) +
  geom_hline(aes(yintercept = MaxCorr_Overall_Speed), linetype = 2, linewidth = 1.5) +
  theme(axis.title=element_text(size=18),
        axis.text = element_text(size=18),
        strip.text = element_text(size = 18))

##################
#Outlier Analysis#
##################
#get those conditions (per participant) where the best lag is at the lower or the upper end of the interval of time lags
Outliers = unique((CCG_Frame %>% group_by(participant, opacity) %>%
                     slice(1) %>%
                     select(participant, opacity, Time_MaxCorr_Overall_Speed, FrameDuration) %>%
                     # filter(Time_MaxCorr_x*FrameDuration > 0.98 | Time_MaxCorr_x*FrameDuration > 0.98 | 
                     #        Time_MaxCorr_y*FrameDuration < 0.02 | Time_MaxCorr_y*FrameDuration < 0.02))$participant)
                     filter(Time_MaxCorr_Overall_Speed*FrameDuration > 0.98 | 
                            Time_MaxCorr_Overall_Speed*FrameDuration < 0.02))$participant)

###############
#Summary Plots#
###############
require(ggdist)

#Time
ggplot(CCG_Frame %>%
         filter(!(participant %in% Outliers)) %>% 
         group_by(participant,opacity) %>%
         slice(1) %>%
         group_by(opacity) %>%
         mutate(Mean = mean(Time_MaxCorr_Overall_Speed),
                SD = sd(Time_MaxCorr_Overall_Speed)), 
       aes(as.factor(opacity),Time_MaxCorr_Overall_Speed*FrameDuration)) +
  geom_point(aes(as.factor(opacity),Mean*FrameDuration), position = position_dodge(width = 0.2), size = 5) +
  geom_errorbar(aes(ymin = Mean*FrameDuration-SD*FrameDuration, 
                    ymax = Mean*FrameDuration+SD*FrameDuration), position = position_dodge(width = 0.2), width = 0.2, linewidth = 1.5) +
  stat_dots(side = "right", justification = -0.2, size = 2, alpha = 0.33) +
  theme(axis.text=element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.title = element_text(size=18),
        legend.text = element_text(size=18)) +
  ylab("Time Lag of Maximum Correlation (s)") +
  scale_x_discrete(name = "Opacity")

#maximum correlation
ggplot(CCG_Frame %>%
         filter(!(participant %in% Outliers)) %>% 
         group_by(participant,opacity) %>%
         slice(1) %>%
         group_by(opacity) %>%
         mutate(Mean = mean(MaxCorr_Overall_Speed),
                SD = sd(MaxCorr_Overall_Speed)), 
       aes(as.factor(opacity),MaxCorr_Overall_Speed)) +
  geom_point(aes(as.factor(opacity),Mean), position = position_dodge(width = 0.2), size = 5) +
  geom_errorbar(aes(ymin = Mean-SD, 
                    ymax = Mean+SD), position = position_dodge(width = 0.2), width = 0.2, linewidth = 1.5) +
  stat_dots(side = "right", justification = -0.2, size = 2, alpha = 0.33) +
  theme(axis.text=element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.title = element_text(size=18),
        legend.text = element_text(size=18)) +
  ylab("Maximum Correlation") +
  scale_x_discrete(name = "Opacity")

#Time
ggplot(CCG_Frame %>%
         filter(!(participant %in% Outliers)) %>% 
         group_by(participant,opacity) %>%
         slice(1) %>%
         group_by(opacity) %>%
         mutate(Mean = mean(Time_MaxCorr_Overall),
                SD = sd(Time_MaxCorr_Overall)), 
       aes(as.factor(opacity),Time_MaxCorr_Overall*FrameDuration)) +
  geom_point(aes(as.factor(opacity),Mean*FrameDuration), position = position_dodge(width = 0.2), size = 5) +
  geom_errorbar(aes(ymin = Mean*FrameDuration-SD*FrameDuration, 
                    ymax = Mean*FrameDuration+SD*FrameDuration), position = position_dodge(width = 0.2), width = 0.2, linewidth = 1.5) +
  stat_dots(side = "right", justification = -0.2, size = 2, alpha = 0.33) +
  theme(axis.text=element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.title = element_text(size=18),
        legend.text = element_text(size=18)) +
  ylab("Time Lag of Maximum Correlation (s)") +
  scale_x_discrete(name = "Opacity")

#maximum correlation
ggplot(CCG_Frame %>%
         filter(!(participant %in% Outliers)) %>% 
         group_by(participant,opacity) %>%
         slice(1) %>%
         group_by(opacity) %>%
         mutate(Mean = mean(MaxCorr_Overall),
                SD = sd(MaxCorr_Overall)), 
       aes(as.factor(opacity),MaxCorr_Overall)) +
  geom_point(aes(as.factor(opacity),Mean), position = position_dodge(width = 0.2), size = 5) +
  geom_errorbar(aes(ymin = Mean-SD, 
                    ymax = Mean+SD), position = position_dodge(width = 0.2), width = 0.2, linewidth = 1.5) +
  stat_dots(side = "right", justification = -0.2, size = 2, alpha = 0.33) +
  theme(axis.text=element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.title = element_text(size=18),
        legend.text = element_text(size=18)) +
  ylab("Maximum Correlation") +
  scale_x_discrete(name = "Opacity")

###############
#Kalman Filter#
###############

#####Kalman filter (as per Bonnen 2015)
KalmanFilter = function(R){
  
  Responses_Wide_Frame_Temp = Responses_Wide %>% filter(participant == j & opacity == k)
  BestLag = round(0.3/Responses_Wide_Frame_Temp$FrameDuration[1]) #pick a reasonable-ish lag
  
  #x
  x_t.x = lag(Responses_Wide_Frame_Temp$x_coord_target,BestLag)
  x_hat_t.x = Responses_Wide_Frame_Temp$mouse_x_cm
  #y
  x_t.y = lag(Responses_Wide_Frame_Temp$y_coord_target,BestLag)
  x_hat_t.y = Responses_Wide_Frame_Temp$mouse_y_cm
  
  Q = 0.3^2 #variance (in stimulus change), not standard deviation
  P <- Q / 2 * (sqrt(1 + 4 * R / Q) - 1)  #posterior variance
  K <- (P + Q)/(P + Q + R)  #Kalman Gain
  
  # compute residual; in Bonnen 2015's script, this is d, 
  #but whenever you print d to console, it is a vector of N times K. So - simplified. 
  #I am assuming this is because matlab just really loves matrices
  temp.x <- x_hat_t.x - K * x_t.x  
  temp.x = temp.x[!is.na(temp.x)] #remove those residuals that come out as NA
  temp.y <- x_hat_t.y - K * x_t.y
  temp.y = temp.y[!is.na(temp.y)] #remove those residuals that come out as NA
  
  temp = c(temp.x, temp.y)
  
  N = length(temp)
  nLL = - (-1/(2*(K^2)*R)*sum(temp^2) - N/2*log(R) - N*log(K))

  # #get the negative log likelihood that the residuals (temp) are from a normal distribution with a mean of 0 and a standard deviation of K^2*R
  # nLL = -sum(LogProb)
  nLL
  # 
  # nLL = -sum(log(dnorm(temp, 0, K^2*R))) #get the negative log likelihood that the residuals (temp) are from a normal distribution with a mean of 0 and a standard deviation of K^2*R
  # nLL
}

#Optimize for each condition per participant
OptimResults = data.frame()
for (j in unique((Responses_Wide$participant))){
  for (k in unique((Responses_Wide %>% filter(participant == j))$opacity)){
    
    print(j)
    print(k)
    
    #keep in mind that R is the VARIANCE
    MLE_Result = stats4::mle(KalmanFilter, start = list(R=0.3^2), lower = 0.005, upper = 10000, method = "Brent")
    
    OptimResults = rbind(OptimResults,
                         data.frame(participant = j,
                                    opacity = k,
                                    Estimate = MLE_Result@coef[1],
                                    LogLike = MLE_Result@details$value))
    print("Done")
  }
}

ggplot(OptimResults %>%
         group_by(participant,opacity) %>% 
         dplyr::slice(1) %>%
         group_by(opacity) %>%
         mutate(Mean_Per_Condition = mean(Estimate^0.5),
                SD_Per_Condition = sd(Estimate^0.5)), 
       aes(as.factor(opacity), Estimate^0.5)) + #take the square root because we want to report standard deviation not variance
  geom_point(aes(as.factor(opacity),Mean_Per_Condition), position = position_dodge(width = 0.2), size = 5) +
  geom_errorbar(aes(ymin = Mean_Per_Condition-SD_Per_Condition, 
                    ymax = Mean_Per_Condition+SD_Per_Condition), 
                          position = position_dodge(width = 0.2), width = 0.2, linewidth = 1.5) +
  stat_dots(side = "right", justification = -0.2, size = 2, alpha = 0.33) +
  theme(axis.text=element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.title = element_text(size=18),
        legend.text = element_text(size=18)) +
  ylab("Sensory Noise Parameter") +
  scale_x_discrete(name = "Opacity")
ggsave(paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/Figures/Results of Kalman Filter.jpg"), w = 14, h = 5)