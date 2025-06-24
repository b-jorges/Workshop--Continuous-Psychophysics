require(dplyr)
require(ggplot2)
require(cowplot)
require(ggdist)
theme_set(theme_cowplot())

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
                      mutate(mouse_2.x = as.numeric(substr(mouse_2.x[1], 2, 7)),
                             mouse_2.y = as.numeric(substr(mouse_2.y[1], 2, 7)),
                             mouse_x_cm = x_coord_mouse*(1/mouse_2.x),
                             mouse_y_cm = y_coord_mouse*(1/mouse_2.y)) %>% 
                      filter(!is.na(opacity)) %>% 
                      group_by(participant, opacity) %>% 
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

HowManyFrames = 1/median(Responses_Wide$FrameDuration)
CCG_Frame = data.frame()
for (i in seq(1,HowManyFrames,4)){
  CCG_Frame = rbind(CCG_Frame, Responses_Wide %>%
                      ungroup() %>%
                      select(participant,opacity,x_coord_target,mouse_x_cm,y_coord_target,FrameDuration,mouse_y_cm) %>%
                      group_by(participant, opacity) %>%
                      mutate(x_coord_target_lagged = lag(x_coord_target, i),
                             y_coord_target_lagged = lag(y_coord_target, i),
                             lag = i))
}

#Calculate the mean and median difference for each lag (separately per participant and condition)
CCG_Frame = CCG_Frame %>%
  group_by(participant, opacity, lag) %>%
  mutate(Correlation_x = cor.test(mouse_x_cm, x_coord_target_lagged)[4]$estimate[[1]],
         Correlation_y = cor.test(mouse_y_cm, y_coord_target_lagged)[4]$estimate[[1]]) %>% 
  group_by(participant, opacity) %>%
  mutate(MaxCorr_x = max(Correlation_x),
         Time_MaxCorr_x = lag[which.max(Correlation_x)],
         MaxCorr_y = max(Correlation_y),
         Time_MaxCorr_y = lag[which.max(Correlation_y)])

save(CCG_Frame, file = paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/SavedVariables/CCG_Frame.RData"))
load(file = paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/SavedVariables/CCG_Frame.RData"))

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

##################
#Outlier Analysis#
##################

#get those conditions (per participant) where the best lag is at the lower or the upper end of the interval of time lags
Outliers = unique((CCG_Frame %>% group_by(participant, opacity) %>%
                     slice(1) %>%
                     select(participant, opacity, Time_MaxCorr_x, Time_MaxCorr_y, FrameDuration) %>%
                     filter(Time_MaxCorr_x*FrameDuration > 0.98 | Time_MaxCorr_x*FrameDuration > 0.98 | 
                            Time_MaxCorr_y*FrameDuration < 0.02 | Time_MaxCorr_y*FrameDuration < 0.02))$participant)

###############
#Summary Plots#
###############
require(ggdist)

#Time // x
ggplot(CCG_Frame %>%
         filter(!(participant %in% Outliers)) %>% 
         group_by(participant,opacity) %>%
         slice(1) %>%
         group_by(opacity) %>%
         mutate(Mean = mean(Time_MaxCorr_x),
                SD = sd(Time_MaxCorr_x)), 
       aes(as.factor(opacity),Time_MaxCorr_x*FrameDuration)) +
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

#Time // y
ggplot(CCG_Frame %>%
         filter(!(participant %in% Outliers)) %>% 
         group_by(participant,opacity) %>%
         slice(1) %>%
         group_by(opacity) %>%
         mutate(Mean = mean(Time_MaxCorr_y),
                SD = sd(Time_MaxCorr_y)), 
       aes(as.factor(opacity),Time_MaxCorr_y*FrameDuration)) +
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

#maximum correlation // x
ggplot(CCG_Frame %>%
         filter(!(participant %in% Outliers)) %>% 
         group_by(participant,opacity) %>%
         slice(1) %>%
         group_by(opacity) %>%
         mutate(Mean = mean(MaxCorr_x),
                SD = sd(MaxCorr_x)), 
       aes(as.factor(opacity),MaxCorr_x)) +
  geom_point(aes(as.factor(opacity),Mean), position = position_dodge(width = 0.2), size = 5) +
  geom_errorbar(aes(ymin = Mean-SD, 
                    ymax = Mean+SD), position = position_dodge(width = 0.2), width = 0.2, linewidth = 1.5) +
  stat_dots(side = "right", justification = -0.2, size = 2, alpha = 0.33) +
  theme(axis.text=element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.title = element_text(size=18),
        legend.text = element_text(size=18)) +
  ylab("Time Lag of Maximum Correlation (s)") +
  scale_x_discrete(name = "Opacity")

#maximum correlation // y
ggplot(CCG_Frame %>%
         filter(!(participant %in% Outliers)) %>% 
         group_by(participant,opacity) %>%
         slice(1) %>%
         group_by(opacity) %>%
         mutate(Mean = mean(MaxCorr_y),
                SD = sd(MaxCorr_y)), 
       aes(as.factor(opacity),MaxCorr_y)) +
  geom_point(aes(as.factor(opacity),Mean), position = position_dodge(width = 0.2), size = 5) +
  geom_errorbar(aes(ymin = Mean-SD, 
                    ymax = Mean+SD), position = position_dodge(width = 0.2), width = 0.2, linewidth = 1.5) +
  stat_dots(side = "right", justification = -0.2, size = 2, alpha = 0.33) +
  theme(axis.text=element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.title = element_text(size=18),
        legend.text = element_text(size=18)) +
  ylab("Time Lag of Maximum Correlation (s)") +
  scale_x_discrete(name = "Opacity")

###############
#Kalman Filter#
###############

#####Kalman filter (as per Bonnen 2015)
KalmanFilter = function(R){
  
  BestLag = round(0.3/FrameDuration) #pick a reasonable-ish lag
  
  Responses_Wide_Frame_Temp = Responses_Wide %>% filter(participant == j & opacity == k)

  #x
  x_t.x = lag(Responses_Wide_Frame_Temp$x_coord_target,BestLag)
  x_hat_t.x = Responses_Wide_Frame_Temp$x_coord_mouse
  x_t_plus_1.x = x_t.x    #angle on next step
  v_t.x = rnorm(length(x_t.x),0,R)    #v_t is noise is perceived angle
  y_t.x = lag(x_t.x,1) + v_t.x    #y_t is perceived angle (actual angle + noise)
  
  #y
  x_t.y = lag(Responses_Wide_Frame_Temp$y_coord_target,BestLag)
  x_hat_t.y = Responses_Wide_Frame_Temp$y_coord_mouse
  x_t_plus_1.x = x_t.x    #angle on next step
  v_t.x = rnorm(length(x_t.x),0,R)    #v_t is noise is perceived angle
  y_t.x = lag(x_t.x,1) + v_t.x    #y_t is perceived angle (actual angle + noise)
  
  #Updating fraction
  Q = 0.3^2 #Q is actual variability in angle (as variance, that's why we square it)
  P = Q/2*((1+4*R/Q)^0.5 - 1) #posterior variance
  K_t = (Q + P)/(Q + P + R) #Kalman gain, used to update the old angle by a percentage of the new perceived angle
  
  #as per Bonnen 2015, the error between predicted and observed values should come from a normal distribution within mean = 0 and sd = K^2 * R
  #x
  Errors.x = x_hat_t.x - K_t*x_t.x
  #y
  Errors.y = x_hat_t.y - K_t*x_t.y
  
  #Get the log likelihood that the observed errors come from this distribution
  #x
  Likelihood.x = log(dnorm(Errors.x, 0, K_t^2 * R))
  #y
  Likelihood.y = log(dnorm(Errors.y, 0, K_t^2 * R))
  
  #combine errors from x and y
  Likelihood = c(Likelihood.x, Likelihood.y)
  
  #Sometimes the values are so likely that R gives us "-Inf" as output.
  #We just crudely substitute those with -500, which is smaller than the smallest log likelihood we can get from R
  Likelihood[Likelihood == -Inf] = -500
  
  #sum up the negative log likelihoods
  -sum(Likelihood, na.rm = TRUE)
}

#Optimize for each condition per participant
OptimResults = data.frame()
for (j in unique((Responses_Wide$participant))){
  for (k in unique((Responses_Wide %>% filter(participant == j))$opacity)){
    
    print(j)
    print(k)
    
    #keep in mind that R is the VARIANCE
    MLE_Result = stats4::mle(KalmanFilter, start = list(R=1), lower = 0, upper = 3, method = "Brent")
    
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