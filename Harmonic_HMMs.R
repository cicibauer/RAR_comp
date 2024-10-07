#' Given the output of depmix (Harmonic HMMs), this function summarises some useful results for further analysis and plotting
#' @param HMM: 24-hour circadian harmonic HMM fitting results (the return of depmix function). Three activity states are assumed:
# 1: inactive state; 2: moderately active state; 3: highly active state
#' @return ML_states: maximal likelihood states at each time point (results of local decoding)
#' @return prob_ML_states: probability of ML_states at each time point (results of local decoding)
#' @return transition_prob (oscillates with 24-h period): transition probability
#' @return circadian_states_prob (oscillates with 24-h period): state probability #' @return AIC and BIC
#' @return obs_density_sq: mean and standard deviation of observation densities (square root), conditioned on 3 states
#' @return obs_density: 5%, 50%, 95% quantile of observation densities (original scale), conditioned on 3 states
#' 
Harmonic_HMMs <-function(HMM,sin_part,cos_part){
  ##in HMM, the order of 3 states are random. we need to re-arrange them according to their mean values
  #note: all the following results are based on 3 states!
  #for other number of states, the parameters order in HMM will change and one should carefully extract them from getpars(HMM)
  obs_params<-data.frame(mean=summary(HMM)[1:3],sd=summary(HMM)[4:6])
  state1_lable<-which.min(obs_params$mean) 
  state3_lable<-which.max(obs_params$mean) 
  state2_lable<-which(obs_params$mean==median(obs_params$mean))
  if(length(state2_lable) >1) {state2_lable <- which(1:3 %ni% c(state1_lable, state3_lable))}
  
  ########################################################## 
  ###local decoding
  ########################################################## 
  e1<-forwardbackward(HMM)$gamma
  L<-nrow(e1)
  state1<-e1[,state1_lable]
  state2<-e1[,state2_lable] 
  state3<-e1[,state3_lable]
  states<-rep(NA,L)
  for(i in 1:L){ states[i]<-ifelse((state1[i]>state2[i] &
                                      state1[i]>state3[i]),0,ifelse(state2[i]>state3[i],1,2)) }
  #ML_states: maximal likelihood of current states 
  ML_states<-states+1
  #ML_states has 3 numbers where
  # 1: inactive; 2: moderately active; 3: highly active
  #probability of each states 
  prob_ML_states<-cbind(state1,state2,state3)
  ########################################################## 
  ###computue the time-varing transition probabilities: transition_prob
  ########################################################## 
  trans1_1<-rep(NA,L) # transition from state 1 to 1
  trans1_2<-rep(NA,L) # transition from state 1 to 2
  trans1_3<-rep(NA,L)
  trans2_1<-rep(NA,L) 
  trans2_2<-rep(NA,L) 
  trans2_3<-rep(NA,L) 
  trans3_1<-rep(NA,L) 
  trans3_2<-rep(NA,L) 
  trans3_3<-rep(NA,L)
  ## in HMM, the order of 3 states are random. we need to re-arrange them 
  A<-getpars(HMM)
  # parameters of the time-varing transition matrix. see eq(3) of the paper. 
  tmat_all<-array(NA,c(3,3,3))
  tmat_all[1,1,]<-A[4:6] 
  tmat_all[1,2,]<-A[7:9] 
  tmat_all[1,3,]<-A[10:12]
  tmat_all[2,1,]<-A[13:15] 
  tmat_all[2,2,]<-A[16:18] 
  tmat_all[2,3,]<-A[19:21]
  tmat_all[3,1,]<-A[22:24] 
  tmat_all[3,2,]<-A[25:27] 
  tmat_all[3,3,]<-A[28:30]
  tmat1<-tmat_all[state1_lable,,] 
  tmat2<-tmat_all[state2_lable,,] 
  tmat3<-tmat_all[state3_lable,,]
  transition_from_one_state<-function(a11,a22,a33){ a123<-c(a11,a22,a33)
  a1<-a123[state1_lable]
  a2<-a123[state2_lable]
  a3<-a123[state3_lable]
  trans<-c(a1,a2,a3)/sum(c(a1,a2,a3))
  trans1<-replace(trans, is.na(trans), (1-sum(trans[!is.na(trans)]))/
                    sum(is.na(trans)) ) 
  return(trans1) }
  for (i in 1:L){ 
    sin_i<-sin_part[i] 
    cos_i<-cos_part[i]
  # transition from state 1 
  a11<-exp(c(1,sin_i,cos_i)%*% tmat1[,1]) 
  a22<-exp(c(1,sin_i,cos_i)%*% tmat1[,2]) 
  a33<-exp(c(1,sin_i,cos_i)%*% tmat1[,3]) 
  cons<-a11+a22+a33 
  y1<-transition_from_one_state(a11,a22,a33) 
  trans1_1[i]<-y1[1]
  trans1_2[i]<-y1[2]
  trans1_3[i]<-y1[3]
  # transition from state 2 
  a11<-exp(c(1,sin_i,cos_i)%*% tmat2[,1]) 
  a22<-exp(c(1,sin_i,cos_i)%*% tmat2[,2]) 
  a33<-exp(c(1,sin_i,cos_i)%*% tmat2[,3]) 
  cons<-a11+a22+a33 
  y2<-transition_from_one_state(a11,a22,a33) 
  trans2_1[i]<-y2[1]
  trans2_2[i]<-y2[2]
  trans2_3[i]<-y2[3]
  # transition from state 3 
  a11<-exp(c(1,sin_i,cos_i)%*% tmat3[,1]) 
  a22<-exp(c(1,sin_i,cos_i)%*% tmat3[,2]) 
  a33<-exp(c(1,sin_i,cos_i)%*% tmat3[,3])
  cons<-a11+a22+a33 
  y3<-transition_from_one_state(a11,a22,a33) 
  trans3_1[i]<-y3[1]
  trans3_2[i]<-y3[2]
  trans3_3[i]<-y3[3]
  }
  transition_prob<- data.frame(trans1_1=trans1_1,trans1_2=trans1_2,trans1_3=trans1_3,
                               trans2_1=trans2_1,trans2_2=trans2_2,trans2_3=trans2_3, trans3_1=trans3_1,trans3_2=trans3_2,trans3_3=trans3_3)
  ########################################################## 
  ###compute the probability of each states (time-varying): prob_states 
  ##########################################################
  P_1<-rep(NA,L) 
  P_2<-rep(NA,L) 
  P_3<-rep(NA,L) 
  P_sum<-rep(NA,L)
  P_1[1]<-0;P_2[1]<-0;P_3[1]<-0;P_sum[1]<-1 
  initial_state<-A[1:3] 
  states_0_guess<-which.max(initial_state)
  #states_0<-which(states_0_guess==c(1,2,3))
  if(state1_lable==states_0_guess){states_0<-1;P_1[1]<-1} 
  if(state2_lable==states_0_guess){states_0<-2;P_2[1]<-1} 
  if(state3_lable==states_0_guess){states_0<-3;P_3[1]<-1}
  for (i in 2:L){
    trans_matrix<-matrix(NA,nrow=3,ncol=3) 
    trans_matrix[1,]<-c(trans1_1[i],trans1_2[i],trans1_3[i]) 
    trans_matrix[2,]<-c(trans2_1[i],trans2_2[i],trans2_3[i]) 
    trans_matrix[3,]<-c(trans3_1[i],trans3_2[i],trans3_3[i])
    P_pervious<-c(P_1[i-1],P_2[i-1],P_3[i-1])
    P_1[i]<-P_pervious %*% trans_matrix[,1] 
    P_2[i]<-P_pervious %*% trans_matrix[,2] 
    P_3[i]<-P_pervious %*% trans_matrix[,3] 
    P_sum[i]<-P_1[i]+P_2[i]+P_3[i]
  } 
  circadian_states_prob<-data.frame(state_1=P_1,state_2=P_2,state_3=P_3)
  ###############observation densities################# 
  mean_obs<-obs_params$mean
  mean_sq<- c(mean_obs[state1_lable],mean_obs[state2_lable],mean_obs[state3_lable])
  sd_sq<- c(obs_params$sd[state1_lable],obs_params$sd[state2_lable],obs_params$sd[state3_lable])
  nc_prams<-(mean_sq/sd_sq)^2
  var<-sd_sq^2 
  obs_density<-matrix(NA,3,3) 
  for (i in 1:3){
    obs_density[i,]<- qchisq(p=c(0.05,0.5,0.95), df=1, ncp = nc_prams[i])*var[i] 
    }
  obs_density_sq<-data.frame(mean=mean_sq,sd=sd_sq)
  return(list(ML_states=ML_states, prob_ML_states=prob_ML_states, transition_prob=transition_prob, circadian_states_prob=circadian_states_prob, AIC=AIC(HMM),BIC=BIC(HMM), obs_density=obs_density,obs_density_sq=obs_density_sq
              ,obs_params=obs_params
              ))
}