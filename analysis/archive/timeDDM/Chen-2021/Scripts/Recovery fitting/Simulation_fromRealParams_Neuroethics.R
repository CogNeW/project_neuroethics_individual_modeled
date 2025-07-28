library(dplyr)

## Fit the tSSM in parallel for multiple subjects
rm(list=ls())

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# CHANGE THIS WITH THE PATH WHERE YOU SAVE ALL THE SCRIPTS AND THE DATA
base_dir = "C:/Users/Brian Kim/Documents/OneDrive - Drexel University/Documents/classes/cognew/Neuroethics/IndividualDifferenceRiskBenefit/DATA_ANALYSIS"
codeFolder = paste0(base_dir, "/analysis/timeDDM/Chen-2021/Scripts/Recovery fitting")
dataFolder = paste0(base_dir, "/analysis/timeDDM/Chen-2021/Scripts/ModelResults/Neuroethics/")
main_folder <- file.path(codeFolder)

# Load Data
behavior_path = paste(base_dir, "/data/Cleaned_Pilot/exclusions.2020.10.21.Rdata", sep="")
load(behavior_path)
risk_data = Neuroethics_Judgement
risk_data <- risk_data %>% 
  # filter(rt <= 5000) %>% 
  mutate(rt = rt / 1000)
risk_data <- risk_data %>%
  mutate(
    response = ifelse(experimental_treatment_selected == "Yes to treatment", 1, 0),
    RT = if_else(response == 0, -rt, rt), # Multiply by -1 if response is 0
    risk = risk * 0.01,                   # Scale risk by 0.01
    gain = gain * 0.01,                    # Scale gain by 0.01
    no_effect = no_effect * 0.01
  )
risk_data$idxP <- as.numeric(as.factor(risk_data$user_id))

# Load stDDM Results
results = read.csv(paste0(dataFolder, "20240313STDDMResults.csv"), header = TRUE, sep = ",")

# Install required packages if necessary:
want = c("DEoptim", "Rcpp", "plyr", "parallel", "RcppParallel","stats4","pracma")
have = want %in% rownames(installed.packages())
if ( any(!have) ) { install.packages( want[!have] ) }
# Now load them all
suppressPackageStartupMessages(lapply(want, require, character.only = TRUE))

RcppParallel::setThreadOptions(numThreads = 1) #this is critical for mclapply

### -----
#load the c++ file containing the functions to simulate the delayed DDM
# will give some errors that can be ignored
sourceCpp(file.path(main_folder, "2ddm_r_cpp_M1.cpp"))

# Need to load the correct data

DDataSim<- data.frame(hd = risk_data$gain, td = risk_data$risk, subject = risk_data$idxP )
### Simulations ###

# parameters to recover with the fitting

DataSim<-DDataSim
subj <- unique(DataSim$subject)
countSub=1
for (s in subj){
  
  wt = subset(results, X == paste0("b2.p[", toString(s), "]"))$Mean
  wh = subset(results, X == paste0("b1.p[", toString(s), "]"))$Mean
  boundary = subset(results, X == paste0("alpha.p[", toString(s), "]"))$Mean
  nDT = subset(results, X == paste0("theta.p[", toString(s), "]"))$Mean
  tHin = subset(results, X == paste0("time[", toString(s), "]"))$Mean
  bias = subset(results, X == paste0("bias[", toString(s), "]"))$Mean
  
  DataSim$wt[DataSim$subject==s] = wt
  DataSim$wh[DataSim$subject==s] = wh
  DataSim$boundary[DataSim$subject==s] = boundary
  DataSim$nDT[DataSim$subject==s] = nDT
  DataSim$tHin[DataSim$subject==s] = tHin
  DataSim$bias[DataSim$subject==s] = bias
  
  vd <- (DataSim$td[DataSim$subject==s])
  hd <- (DataSim$hd[DataSim$subject==s])
  for (i in 1:length(vd)) {
    #ddm2_parallel(d_v, d_h, thres, nDT, tIn, bias, vd, hd, 1, N)
    rts = ddm2_parallel(wt,wh,boundary,nDT,tHin,bias,vd[i],hd[i],1,1)
    if (rts<0){
      DataSim$choice[DataSim$subject==s & DataSim$td ==vd[i] & DataSim$hd ==hd[i]] <- 0
      DataSim$rt[DataSim$subject==s & DataSim$td ==vd[i] & DataSim$hd ==hd[i]] <- -rts
    }else{
      DataSim$choice[DataSim$subject==s & DataSim$td ==vd[i] & DataSim$hd ==hd[i]] <- 1
      DataSim$rt[DataSim$subject==s & DataSim$td ==vd[i] & DataSim$hd ==hd[i]] <- rts
    }
  }
  countSub=countSub+1
}

# DataSim$rt[DataSim$rt>4]<-NaN
# DataSim$rt[DataSim$rt<0.2]<-NaN
# DataSim<-DataSim[!(is.nan(DataSim$rt)),]


save(DataSim,file=file.path(codeFolder, paste("NeuroethicsSimulation.RData", sep = "")))
write.csv(DataSim, file=paste(c(file.path(codeFolder,"NeuroethicsSimulation"), ".csv"), collapse = "") ) 


