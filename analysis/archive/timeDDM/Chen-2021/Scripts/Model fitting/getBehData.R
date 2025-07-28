getBehData<-function () {
    # CHANGE THIS WITH THE PATH WHERE YOU SAVE ALL THE SCRIPTS AND THE DATA
  baseFolder = getwd()
  pathToFolder <- baseFolder
  
    ##################################

  main_folder <- file.path(pathToFolder,"Data")
  filelist = list.files(path = main_folder, pattern = ".*.txt")
  filelist = file.path(main_folder, filelist)
  datalist = lapply(filelist, FUN=read.table, header=TRUE) 
  datafr = do.call("rbind", datalist) 
  datafr = datafr[datafr$Timeout ==0 & datafr$RT > 200,]
  datafr$rt = datafr$RT/1000
  datafr$choice = datafr$Acceptance
  #d = w1*(GainProb - LossProb) + w2*Magnitude
  datafr$td = as.numeric(scale(I(datafr$WinProb - datafr$LossProb), center = F))
  datafr$hd = as.numeric(scale(datafr$Magnitude, center = F))
  datafr$adlog = c(rep(NA,length(datafr$AccumReward)))
  
  for (ss in unique(datafr$Subject)){
    ind_accum = datafr$AccumReward[datafr$Subject==ss]
    accumre = c(rep(NA,length(ind_accum)))
    accumre[1]= 0
    accumre[2:length(accumre)] = ind_accum[1:length(ind_accum)-1]
    log_indaccum = c(rep(NA,length(ind_accum)))
    log_indaccum[which(accumre>0)] = log(accumre[which(accumre>0)])
    log_indaccum[which(accumre<0)] = log(abs(accumre[which(accumre<0)]))*-1
    log_indaccum[is.na(log_indaccum)] =0
    
    datafr$adlog[datafr$Subject==ss] = as.numeric(log_indaccum)
  }
  datafr$ad = as.numeric(scale(datafr$adlog, center = F))
  return(datafr)
}
