rm(list=ls())

pack = "runjags"
have = pack %in% rownames(installed.packages())
if ( any(!have) ) { install.packages( want[!have] ) }
# Now load it
suppressPackageStartupMessages(lapply(pack, require, character.only = TRUE))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# CHANGE THIS WITH THE PATH WHERE YOU SAVE ALL THE SCRIPTS AND THE DATA
pathToFolder <- 'folder path'  #getwd()

# DATA FILE NAME

# load example data of a bynary food choice task
main_folder <- file.path(pathToFolder,"code")
source( file.path(main_folder, "getBehData.R" ))

Data=getBehData()
Data$Group[Data$Subject<200000] = 'YA'
Data$Group[Data$Subject>200000] = 'OA'

# aa <- data.frame('ID'=ordered(unique(Data$Subject)))
# write.csv(aa, file=paste(c(file.path(main_folder,"results","SubjectID"), ".csv"), collapse = "") )

for (gg in unique(Data$Group)){
    for (cc in unique(Data$Condition)){
        
        data=Data[Data$Group==gg & Data$Condition==cc,]
    
    
        idx = which(data$choice==0)
        data$RT <- data$rt
        data$RT[idx] = data$rt[idx] * -1
        
        idxP = as.numeric(ordered(data$Subject)) #makes a sequentially numbered subj index
        
        hd = (data$hd) # magnitude
        td = (data$td) # probability
    
        rtpos = data$rt
        
        y= data$RT
        
        
        N = length(y)
        
        ns = length(unique(idxP))
        
        #--------------------------------------------
        # fit the model
        
        dat <- dump.format(list(N=N, y=y, idxP=idxP, hd=hd, td=td, ns=ns)) #no RT enter
        
        
    
        inits3 <- dump.format(list( alpha.mu=0.5, alpha.pr=0.05,  theta.mu=0.1,
                                    theta.pr=0.05,  b1.mu=0.3, b1.pr=0.05, b2.mu=0.01, b2.pr=0.05, bias.mu=0.4,
                                    bias.kappa=1, y_pred=y,  .RNG.name="base::Super-Duper", .RNG.seed=99999))
        
        
        inits2 <- dump.format(list( alpha.mu=0.1, alpha.pr=0.05, theta.mu=0.01,
                                    theta.pr=0.05, b1.mu=0.3, b1.pr=0.05, b2.mu=0.1, b2.pr=0.05, bias.mu=0.5,
                                    bias.kappa=1, y_pred=y,  .RNG.name="base::Wichmann-Hill", .RNG.seed=1234))
        
        inits1 <- dump.format(list( alpha.mu=0.2, alpha.pr=0.05, theta.mu=0.15,
                                    theta.pr=0.05, b1.mu=0.1, b1.pr=0.05, b2.mu=0.05, b2.pr=0.05, bias.mu=0.3,
                                    bias.kappa=1, y_pred=y, .RNG.name="base::Mersenne-Twister", .RNG.seed=6666 ))
        
        monitor = c(
            "b1.mu", "b2.mu", "theta.mu", "alpha.mu","bias.mu", "b1.p", "b2.p",  "theta.p", "bias","alpha.p", "deviance","pd","dic" )
        
        
        results <- run.jags(model=file.path(main_folder,"tssmHT_modelWIENER_standardDDM.txt"), 
                            monitor=monitor, data=dat, n.chains=3, inits=c(inits1,inits2, inits3), 
                            monitor.pd = TRUE, monitor.pd.i = TRUE, monitor.deviance = TRUE,
                            plots = TRUE, module="wiener", burnin=5000, sample=1000, thin=10)
        suuum<-summary(results)
        save(results,data, file=paste(c(file.path(main_folder,"results","results_HtSSM_FIT_noTime_addDIC_210214_cond"),cc,'_',gg, ".RData"), collapse = ""))
        write.csv(suuum, file=paste(c(file.path(main_folder,"results","results_HtSSM_FIT_noTime_addDIC_210214_cond"),cc,'_',gg, ".csv"), collapse = "") ) 
    }
}

