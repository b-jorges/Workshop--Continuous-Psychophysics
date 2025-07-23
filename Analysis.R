# install.packages(c("dplyr", "ggplot2", "cowplot", "ggdist"))

require(dplyr)
require(ggplot2)
require(cowplot)
require(ggdist)
theme_set(theme_cowplot())

#load all datafiles (we use the .csv files)
flist <- list.files(paste0(dirname(rstudioapi::getSourceEditorContext()$path), "/Experiment/ContPsyTracking/data/"))[c(grep(".csv", list.files(paste0(dirname(rstudioapi::getSourceEditorContext()$path), "//Experiment/ContPsyTracking/data"), full.names = TRUE)))]

Responses_Wide = c()
Data2 = c()
j = 0
for (i in flist){
  j = j+1
  print(i)
  if (j == 1){
    
    #make a new dataframe based on the first .csv we load
    Responses_Wide  = read.csv(paste0(dirname(rstudioapi::getSourceEditorContext()$path), "/Experiment/ContPsyTracking/data/",i), header = TRUE, row.names = NULL)

    #only get those columns that we need because PsychoPy loves to just save e v e r y t h i n g
    Responses_Wide = Responses_Wide %>% select(mouse_2.x, mouse_2.y, time_in_run, x_coord_mouse, y_coord_mouse, x_coord_target, 
                                               y_coord_target, opacity, participant, frameRate, session)
  } else {
    
    #load subsequent datafiles
    Data2  = read.csv(paste0(dirname(rstudioapi::getSourceEditorContext()$path), "/Experiment/ContPsyTracking/data/",i), header = TRUE, row.names = NULL)
    Data2 = Data2 %>% select(mouse_2.x, mouse_2.y, time_in_run, x_coord_mouse, y_coord_mouse, x_coord_target, 
                                               y_coord_target, opacity, participant, frameRate, session)
  }
  
  #add every new .csv to the first dataframe
  Responses_Wide = rbind(Responses_Wide,Data2)
}

Responses_Wide$
Responses_Wide = Responses_Wide %>%
  group_by(participant, session) %>%
  #transform recorded mouse data into same unit as displayed polygon (cm from center of screen)
                      mutate(mouse_2.x = as.numeric(substr(mouse_2.x[1], 2, 5)),
                             mouse_2.y = as.numeric(substr(mouse_2.y[1], 2, 5)),
                             mouse_x_cm = x_coord_mouse*(5/mouse_2.x),
                             mouse_y_cm = y_coord_mouse*(5/mouse_2.y),
                             
                             #to take away the correlation between target position at t and target position at t + 1,
                             #we calculate the x and y speeds of target and mouse cursor
                             mouse_speed.x = (mouse_x_cm-lag(mouse_x_cm))/(time_in_run - lag(time_in_run, 1)),
                             mouse_speed.y = (mouse_y_cm-lag(mouse_y_cm))/(time_in_run - lag(time_in_run, 1)),
                             target_speed.x = (x_coord_target-lag(x_coord_target))/(time_in_run - lag(time_in_run, 1)),
                             target_speed.y = (y_coord_target-lag(y_coord_target))/(time_in_run - lag(time_in_run, 1))) %>%
  
  #filter out rows with no data (because PsychoPy .csv files can be a little messy otherwise)
                      filter(!is.na(opacity)) %>% 
                      group_by(participant) %>%
  
  #get median duration of frame for each participant
                      mutate(FrameDuration = median(Responses_Wide$time_in_run - lag(Responses_Wide$time_in_run, 1), na.rm = TRUE))

##################
#some basic plots#
##################
#target versus mouse (x)
ggplot(Responses_Wide %>% filter(participant == unique(Responses_Wide$participant)[2] & session == 1), 
       aes(time_in_run, x_coord_target)) +
  geom_point() +
  geom_point(aes(time_in_run, mouse_x_cm), color = "red") +
  ylab("x position") +
  xlab("Time (s)") +
  facet_grid(.~opacity) +
  theme(axis.title=element_text(size=18),
        axis.text = element_text(size=18),
        strip.text = element_text(size = 18))

#target versus mouse (y)
ggplot(Responses_Wide %>% filter(participant == unique(Responses_Wide$participant)[1] & session == 1), 
       aes(time_in_run, y_coord_target)) +
  geom_point() +
  geom_point(aes(time_in_run, mouse_y_cm), color = "red") +
  ylab("y position") +
  xlab("Time (s)") +
  facet_grid(participant~opacity) +
  theme(axis.title=element_text(size=18),
        axis.text = element_text(size=18),
        strip.text = element_text(size = 18))

#x versus y
ggplot(Responses_Wide %>% filter(participant == unique(Responses_Wide$participant)[1] & 
                                   session == 1 &
                                   time_in_run > 10 &  time_in_run < 20), 
       aes(x_coord_target, y_coord_target)) +
  geom_point() +
  geom_point(aes(mouse_x_cm, mouse_y_cm), color = "red") +
  ylab("y position (cm)") +
  xlab("x position (cm)") +
  facet_grid(participant~opacity) +
  theme(axis.title=element_text(size=18),
        axis.text = element_text(size=18),
        strip.text = element_text(size = 18))

###################
#Crosscorrelations#
###################

#how many frames in 1s
CCG_Frame_Temp = data.frame()

for (j in unique(Responses_Wide$participant)){
  
  Responses_Wide_Temp = Responses_Wide %>% filter(participant == j)
  HowManyFrames = 1/median(Responses_Wide_Temp$FrameDuration)
  
  #create copies of the dataframe where we lag target and response positions
  #by 1 to HowManyFrames in steps of 4 (only 4 to keep size of this dataframe manageable)
  for (i in seq(1,HowManyFrames,4)){
    CCG_Frame_Temp = rbind(CCG_Frame_Temp, Responses_Wide_Temp %>%
                        ungroup() %>%
                        
                        #select only those columns that we need for further analysis
                        select(participant, opacity, FrameDuration, session,
                               mouse_speed.x, mouse_speed.y, 
                               target_speed.x, target_speed.y) %>%
                        group_by(participant, opacity, session) %>%
                        
                        #lag the target speed by i frames (such as to cover a range of no lag to 1s lag)
                        mutate(target_speed.x_lagged = lag(target_speed.x, i),
                               target_speed.y_lagged = lag(target_speed.y, i),
                               
                               lag = i))
  }
}

CCG_Frame = CCG_Frame_Temp %>%
  group_by(participant, opacity, session, lag) %>%
  
  #Calculate the correlation between target and response speeds for each lag
  #(separately per participant, condition and session, for x/y directions)
  mutate(Correlation_x_Speed = cor.test(mouse_speed.x, target_speed.x_lagged)[4]$estimate[[1]],
         Correlation_y_Speed = cor.test(mouse_speed.y, target_speed.y_lagged)[4]$estimate[[1]],
         
         #calculate the mean correlation of x and y correlations because we don't think this should differ between
         #x and y directions
         Correlation_Overall_Speed = (Correlation_x_Speed + Correlation_y_Speed)/2) %>%
  
  #take the average maximum correlation and time lag of the maximum correlation (Correlation_Overall_Speed)
  #for each participant, opacity level and session. 
  #These are our main dependent variables
  group_by(participant, opacity, session) %>%
  mutate(MaxCorr_Overall_Speed = max(Correlation_Overall_Speed),
        Time_MaxCorr_Overall_Speed = lag[which.max(Correlation_Overall_Speed)]) %>%
  
  #pair the size of this dataframe down to free up some space
  group_by(participant, opacity, session, lag) %>% 
  slice(1)

#save this data frame so we don't have to run the cross-correlation again
save(CCG_Frame, file = paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/SavedVariables/CCG_Frame.RData"))
# load(file = paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/SavedVariables/CCG_Frame.RData"))
#########################

#########################
#Cross Correlation Plots#
#########################
ggplot(CCG_Frame %>% 
         group_by(participant,opacity,lag) %>% 
         slice(1) %>% filter(participant == unique(CCG_Frame$participant)[1]),
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

##################
#Outlier Analysis#
##################
#get those conditions (per participant) where the best lag is at the lower or the upper end of the interval of time lags
Outliers = unique((CCG_Frame %>% group_by(participant, opacity) %>%
                     slice(1) %>%
                     select(participant, opacity, Time_MaxCorr_Overall_Speed, FrameDuration) %>%
                     filter(Time_MaxCorr_Overall_Speed*FrameDuration > 0.98 | 
                            Time_MaxCorr_Overall_Speed*FrameDuration < 0.02))$participant)
##################

###############
#Summary Plots#
###############

#Time lag of maximum correlation
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
###############


###############
#Kalman Filter#
###############
KalmanFilter = function(R){
  
  #take the condition that we are interested in from the data frame
  Responses_Wide_Frame_Temp = Responses_Wide %>% filter(participant == j & session == k & opacity == m)
  
  #pick a reasonable-ish lag: here I chose 0.3s, based on the crosscorrelations
  BestLag = round(0.3/Responses_Wide_Frame_Temp$FrameDuration[1]) 

  #get the actual target positions over time in x (x_t.x) and y (x_t.y) directions, 
  #offset by the "BestLag" set above
  #also get the participant responses over time in x (x_hat_t.x) and y (x_hat_t.y) directions
  #x
  x_t.x = lag(Responses_Wide_Frame_Temp$x_coord_target,BestLag)
  x_hat_t.x = Responses_Wide_Frame_Temp$mouse_x_cm
  #y
  x_t.y = lag(Responses_Wide_Frame_Temp$y_coord_target,BestLag)
  x_hat_t.y = Responses_Wide_Frame_Temp$mouse_y_cm
  
  #Q is the variance (in stimulus change), not the standard deviation. 
  #We assume that participants learn the distribution used for the random walk very quickly 
  #(as per Bonnen 2015) and use it as a prior
  Q = 0.3^2 
  
  #via... math, the posterior variance is calculated from the variance of the prior (Q) 
  #and the variance of the likelihood (R). R represents the sensory uncertainty that we are interested in
  #and which will be fitted.
  P = (Q / 2) * (sqrt(1 + 4 * R / Q) - 1)  #posterior variance
  
  #the Kalman gain (K) can then be computed based on the variance of the prior (Q), the variance of
  #the posterior (P) and the variance of the likelihood (R)
  K = (P + Q)/(P + Q + R)
  
  #here we beginn fitting. For each sensory uncertainty (R) we try out, we generate predictions 50 times
  #and compare them to what we actually observe from the participants.
  #the goal is to minimize these residuals.
  #Here we deviate from Bonnen et al. 2015 to implement an RMSE minimization rather than the arguably
  #harder-to-follow maximum likelihood estimation employed in Bonnen et al.
  Residuals = c()
  
  for (i in 1:50){
  #y_t.x and y_t.x are the perceived location of the stimulus in x and y directions. 
  #We assume that there is no bias, so that would be the actual position (x_t.x) 
  #with noise added (via the rnorm() function)
  y_t.x = x_t.x + rnorm(length(x_t.x),0,R)
  
  #The Kalman filter weighs the internal representation of the object location (x_hat_t.x) at t-1
  #versus the uncertainty with which the stimulus is perceived at t+0 (we use the lag() function)
  #to achieve the t-1 // t+0 offset.
  Predictions.x = (1 - K)*x_hat_t.x + K*lag(y_t.x,1)
  
  #same thing but in y direction
  y_t.y = x_t.y + rnorm(length(x_t.y),0,R)
  Predictions.y = (1 - K)*x_hat_t.y + K*lag(y_t.y,1)
  
  #get the residuals by comparing predictions (Predictions.x, Predictions.y) to the observed
  #data (x_hat_t.x, x_hat_t.y). We use lag() to line them up properly
  Residuals = c(Residuals, lag(Predictions.x, 1) - x_hat_t.x, lag(Predictions.y, 1) - x_hat_t.y)}

  #get the Root Mean Squared Error from the residuals
  (mean(Residuals^2, na.rm = TRUE))^0.5
}

#We loop through each participant and condition separately and get the best fit for R
OptimResults = data.frame()
for (j in unique((Responses_Wide$participant))){
  for (k in unique((Responses_Wide %>% filter(participant == j))$session)){
    for(m in unique((Responses_Wide %>% filter(participant == j & session == k))$opacity)){
      
      #keep track of where we are:
      print(paste0("Participant: ", j))
      print(paste0("Session: ", k))
      print(paste0("Opacity: ", m))
      
      #we use the optimize() function to minimize the RMSE in the KalmanFilter() function set up above.
      #lower and upper are the limits for the R values optimize() will try out
      Fit = optimize(KalmanFilter, lower = 0.005, upper = 10)
      #keep in mind that R is the VARIANCE
      
      #save the participant/condition, the fitted R, as well as the RMSE corresponding to the fitted R
      OptimResults = rbind(OptimResults,
                           data.frame(participant = j,
                                      opacity = m,
                                      session = k,
                                      Estimate = Fit$minimum,
                                      RMSE = Fit$objective))
      
      print("Done")
    }
  }
}


#########################
###plot the fitted Rs####
#########################
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
  ylab("Sensory Noise Parameter (cm)") +
  scale_x_discrete(name = "Opacity")
#########################