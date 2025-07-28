rm(list=ls())

pack = "runjags"
have = pack %in% rownames(installed.packages())
if ( any(!have) ) { install.packages( want[!have] ) }
# Now load it
suppressPackageStartupMessages(lapply(pack, require, character.only = TRUE))


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# CHANGE THIS WITH THE PATH WHERE YOU SAVE ALL THE SCRIPTS AND THE DATA
codeFolder <- 'folder path'  #code folder
simFolder <- 'folder path'  #simulated data folder
resultFolder <- 'folder path'  #result folder

# DATA FILE NAME

# load example data of a bynary food choice task
main_folder <- file.path(codeFolder)
source( file.path(codeFolder, "getBehData.R" ))

aa = 1
age = c('OA','YA')
cond= c('C','A')

for (gg in age){
    for (cc in cond){
 

        load(file.path(simFolder, paste("DataSim_beh_M1_cond",cc,"_",gg,".RData",sep="")))
        Data=DataSim
    
    
        idx = which(Data$choice==0)
        Data$RT <- Data$rt
        Data$RT[idx] = Data$rt[idx] * -1
        
        idxP = as.numeric(ordered(Data$subject)) #makes a sequentially numbered subj index
        
        hd = (Data$hd) # magnitude
        td = (Data$td) # probability
    
        rtpos = Data$rt
        
        y= Data$RT
        
        
        N = length(y)
        
        ns = length(unique(idxP))
    
        #--------------------------------------------
        # fit the model
        
        dat <- dump.format(list(N=N, y=y, idxP=idxP, hd=hd, td=td, rt=rtpos, ns=ns))
    
        
        inits3 <- dump.format(list( alpha.mu=2, alpha.pr=0.05,  time.mu=0.5, time.pr= 0.05, theta.mu=0.1,
                                    theta.pr=0.05,  b1.mu=0.3, b1.pr=0.05, b2.mu=0.01, b2.pr=0.05, bias.mu=0.4,
                                    bias.kappa=1, y_pred=y,  .RNG.name="base::Super-Duper", .RNG.seed=99999))
        
        
        inits2 <- dump.format(list( alpha.mu=2.2, alpha.pr=0.05, time.mu=-0.1, time.pr= 0.05, theta.mu=0.01,
                                    theta.pr=0.05, b1.mu=0.3, b1.pr=0.05, b2.mu=0.1, b2.pr=0.05, bias.mu=0.5,
                                    bias.kappa=1, y_pred=y,  .RNG.name="base::Wichmann-Hill", .RNG.seed=1234))
        
        inits1 <- dump.format(list( alpha.mu=2.4, alpha.pr=0.05, time.mu=0, time.pr= 0.05, theta.mu=0.15,
                                    theta.pr=0.05, b1.mu=0.1, b1.pr=0.05, b2.mu=0.05, b2.pr=0.05, bias.mu=0.3,
                                    bias.kappa=1, y_pred=y, .RNG.name="base::Mersenne-Twister", .RNG.seed=6666 ))
        
        monitor = c(
            "b1.mu", "b2.mu", "theta.mu", "alpha.mu","bias.mu","time.mu", "b1.p", "time.pr","time", "b2.p",  "theta.p", "bias","alpha.p", "deviance" )
        
        
        results <- run.jags(model=file.path(pathToFolder,"tssmHT_model_priors_M1.txt"), 
                            monitor=monitor, data=dat, n.chains=3, inits=c(inits1,inits2, inits3), 
                            plots = TRUE, method="parallel", module="wiener", burnin=50000, sample=10000, thin=10)
        
    
        suuum<-summary(results)
        save(results,Data, file= file.path(resultFolder,paste("results_HtSSM_RecovFitRealParams_M1_",gg,"_Cond",cc,".RData",sep="")) )
        write.csv(suuum, file= file.path(resultFolder,paste("results_HtSSM_RecovFitRealParams_M1_",gg,"_Cond",cc,".csv",sep="")))
    }
}

