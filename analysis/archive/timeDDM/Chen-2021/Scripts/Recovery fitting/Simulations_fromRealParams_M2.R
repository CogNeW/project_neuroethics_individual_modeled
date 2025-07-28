## Fit the tSSM in parallel for multiple subjects
rm(list=ls())

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# CHANGE THIS WITH THE PATH WHERE YOU SAVE ALL THE SCRIPTS AND THE DATA
codeFolder <- 'folder path'  #code folder
simFolder <- 'folder path'  #simulated data folder

# DATA FILE NAME

# load example data of a bynary food choice task
main_folder <- file.path(codeFolder)


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
sourceCpp(file.path(main_folder, "2ddm_r_cpp_M2.cpp"))

### ----


###  format choice & RT data ----

#Need a dataframe concatenating all subjects' data with the following variables. ----
#The data should be cleaned to remove missed trials or implausibly fast responses < 200 ms

# subject = indicator for each individual subject
# condition = indicator for different conditions if applicable
# choice = [1,0] indicator for yes/no or left/right choices
# rt = reaction times for every choice 
# td = difference in taste ratings between the two options 
# generally, these should be z-scored (compute mean + sd from all ratings in experiment)
# if using the 426-point rating scales, consider binning ratings into 10 or 11 bins to speed up fitting <- zscore after
# hd = difference in healthiness ratings between the two options (same details as above)
# ----

source(file.path(main_folder, "getBehData.R"))
Data=getBehData()

Data$Group[Data$Subject>200000]="OA"
Data$Group[Data$Subject<200000]="YA"

for (gg in c("YA","OA")){
    for (cc in c("C","A")){
        
    #### NB! change this to only use one age group of subjects
    Data_partial=Data[Data$Condition==cc & Data$Group==gg,]
    
    ### NB!! replace results_HtSSM_FIT_17052019_condA.RData with the name of the file .RData of the fitting you already run for condition A and specific age group
    load(file.path(pathToFolder,paste("code/results/addAccumulation/results_HtSSM_FIT_addAccumReward_210617_10000_H2_cond",cc,"_",gg,".RData", sep = "")))
    chain=rbind(results$mcmc[[1]], results$mcmc[[2]], results$mcmc[[3]])
    
    DDataSim<- data.frame(hd = Data_partial$hd, td = Data_partial$td, ad = Data_partial$ad, subject = Data_partial$Subject )
    ### Simulations ###
    
    # parameters to recover with the fitting
    
    DataSim<-DDataSim
    subj <- unique(DataSim$subject)
    countSub=1
    for (s in subj){
        
        wt = mean(chain[,c( paste( c("b2.p[",toString(countSub),"]"), collapse = ""))])
        wh = mean(chain[,c( paste( c("b1.p[",toString(countSub),"]"), collapse = ""))])
        wa = mean(chain[,c( paste( c("b3.p[",toString(countSub),"]"), collapse = ""))])
        boundary = mean(chain[,c( paste( c("alpha.p[",toString(countSub),"]"), collapse = ""))])
        nDT = mean(chain[,c( paste( c("theta.p[",toString(countSub),"]"), collapse = ""))])
        tHin = mean(chain[,c( paste( c("time[",toString(countSub),"]"), collapse = ""))])
        bias = mean(chain[,c( paste( c("bias[",toString(countSub),"]"), collapse = ""))])

        DataSim$wt[DataSim$subject==s] = wt
        DataSim$wh[DataSim$subject==s] = wh
        DataSim$wa[DataSim$subject==s] = wa
        DataSim$boundary[DataSim$subject==s] = boundary
        DataSim$nDT[DataSim$subject==s] = nDT
        DataSim$tHin[DataSim$subject==s] = tHin
        DataSim$bias[DataSim$subject==s] = bias
        
        vd <- (DataSim$td[DataSim$subject==s])
        hd <- (DataSim$hd[DataSim$subject==s])
        ad <- (DataSim$ad[DataSim$subject==s])
        for (i in 1:length(vd)) {
            #ddm2_parallel(d_v, d_h, thres, nDT, tIn, bias, vd, hd, 1, N)
            rts = ddm2_parallel(wt,wh,wa,boundary,nDT,tHin,bias,vd[i],hd[i],ad[i],1,1)
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
    
    DataSim$rt[DataSim$rt>4]<-NaN
    DataSim$rt[DataSim$rt<0.2]<-NaN
    DataSim<-DataSim[!(is.nan(DataSim$rt)),]
    
    
    save(DataSim,file=file.path(simFolder, paste("DataSim_beh_M2_cond",cc,"_",gg,".RData", sep = "")))
    }
}


