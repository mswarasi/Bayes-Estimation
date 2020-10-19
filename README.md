This repository contains R and FORTRAN programs associated with the article, “Estimating the prevalence of multiple diseases from two-stage hierarchical pooling.”

Dual.GT.Dorfman.Bayes64.R – main R program that fits the proposed Bayesian model to estimate the prevalence of multiple infections.

multigibbs64.dll – compiled FORTRAN subroutine that works under the main R program.

multigibbs64.f90 – source code of the compiled FORTRAN subroutine.

Documentation and user instructions are provided within the R function: Dual.GT.Dorfman.Bayes64.R.

Reference:

Warasi, M., Tebbs, J., McMahan, C., and Bilder, C. (2016). Estimating the prevalence of multiple diseases from two-stage hierarchical pooling. Statistics in Medicine, 35, 3851-3864.



############### Simulation Examples ###############


# Specify the working directory
setwd(dir = "C:/Users/msarker/Desktop/TwoStageMultDis_GibHub_101920")

# Load the "MCMCpack" package, used in the function: Dual.GT.Dorfman.Bayes( )
library(MCMCpack)	


# Simuating IPP group testing data using the function "Dual.IPP.data" and then fit our Bayesian model using "Dual.GT.Dorfman.Bayes". 
 
N <- 1000

p <- c(.80,.10,.09,.01)

Se <- c(.95,.95)

Sp <- c(.99,.99)

c <- 3


IppData <- Dual.IPP.data(N,p,Se,Sp,c)$data      # data simulation step

indiv_diag <- IppData[ ,c(5,6)]					# individual diagnoses

group_resp <- IppData[ ,c(3,4)][seq(1,N,c), ]	# group testing responses


# Note: To start iterating, the diagnosed statuses of the individuals are used as their true statuses. One can, however, specify only 0's.
indiv_true <- indiv_diag
		

# The model is fit with the Dirichlet prior with the power parameter a0 = 0.5 and the blat beta priors for Se and Sp.
param <- Dual.GT.Dorfman.Bayes(indiv.diag=indiv_diag,group.resp=group_resp,indiv.true=indiv_true,gSize=c,p0=p,N0=N,a0=0.5,TP1=0,FN1=0,TN1=0,FP1=0,TP2=0,FN2=0,TN2=0,FP2=0,GI=6000,burn=1000,thin=5)


# Some estimation results
apply(param$prevalence, 2, median)# Median estimates for prevalence parameters

        p00         p10         p01         p11 
        
0.809034290 0.096036529 0.087255741 0.007241692 

apply(param$prevalence, 2, sd)# Standard deviations for prevalence parameters

        p00         p10         p01         p11 
        
0.011143304 0.008817562 0.007876397 0.002248644 


apply(param$accuracy, 2, median)# Median estimates for accuracy parameters

      Se1       Se2       Sp1       Sp2 
      
0.9545391 0.9582054 0.9863735 0.9955384 

apply(param$accuracy, 2, sd)   # Standard deviations for accuracy parameters

        Se1         Se2         Sp1         Sp2 
        
0.032557391 0.026275380 0.007875868 0.004515715 



