
#****************************************************************************************#
# Gibbs sampler to make Bayesian inference using group testing data 
# under the Iowa IPP pooling algorithm for two diseases (chlamydia and gonorrhea)
# 
# Note: This R function makes use of an external FORTRAN DLL - "multigibbs64.dll"
#       The .dll only works with a 64-bit R
#
# Note: This R function requires the "MCMCpack" package to sample from a Dirichlet distribution
#		> Install the MCMCpack package before running the code
#
#
# Input:
# indiv.diag = N by 2 matrix of individual diagnoses (i.e., Stage 2 test results)
#			 = Matrix(Y_{i1k},Y_{i2k})
# group.resp = K by 2 matrix of pool testing responses, where K is the number of pools
#			 = Matrix(Z_{1k},Z_{2k})
# indiv.true = N by 2 matrix of individual true statuses
#			 = Matrix(\widetilde{Y}_{i1k},\widetilde{Y}_{i2k})
# gSize 	 = Pool size
# p0   		 = An estimate of p from historical data
# N0   		 = Historical data sample size
# a0   		 = Precision parameter in [0,1] 
# TP1 		 = Number of truly positive test results for disease 1 (chlamydia)
# FN1		 = Number of falsely negative test results for disease 1 (chlamydia)
# TN1		 = Number of truly negative test results for disease 1 (chlamydia)
# FP1		 = Number of falsely positive test results for disease 1 (chlamydia)
# TP2 	 	 = Number of truly positive test results for disease 2 (gonorrhea)
# FN2		 = Number of falsely negative test results for disease 2 (gonorrhea)
# TN2		 = Number of truly negative test results for disease 2 (gonorrhea)
# FP2		 = Number of falsely positive test results for disease 2 (gonorrhea)
# GI		 = Number of Gibbs iterates
# burn	     = Burn-in period
# thin		 = Thinning parameter which picks every "thin"th sample point from posterior samples after burn-in period.
#			
#		 	
# Output: 
# A list-object of posterior samples for prevalence and accuracy parameters
# prevalence = Matrix(p00,p10,p01,p11)
# accuracy	 = Matrix(Se1,Se2,Sp1,Sp2)
#
Dual.GT.Dorfman.Bayes <- function(indiv.diag,group.resp,indiv.true=indiv.diag,gSize,p0=c(.25,.25,.25,.25),N0=100,a0=0,TP1=0,FN1=0,TN1=0,FP1=0,TP2=0,FN2=0,TN2=0,FP2=0,GI=2500,burn=500,thin=1){
	#
	# The dll "multigibbs64.dll" samples \widetilde{Y}_{i1k} and 
	# \widetilde{Y}_{i2k} from a posterior multinomial distribution.
	#
	dyn.load("multigibbs64.dll")		# Specify the source directory 

	data <- IPP.structured.data(indiv.diag,group.resp,indiv.true,gSize)	
	N <- dim(data)[1]		# Sample size
	Y1til <- data[ ,1]		# Individual true statuses for disease 1
	Y2til <- data[ ,2]		# Individual true statuses for disease 2
	Y1 <- data[ ,5]			# Individual retesting results for disease 1
	Y2 <- data[ ,6]			# Individual retesting results for disease 2
	
	# This block is necessary if N/gSize > 0
	if((N%%gSize) > 0){
		Y1til <- c(Y1til,array(0,(gSize-N%%gSize)))
		Y2til <- c(Y2til,array(0,(gSize-N%%gSize)))
		Y1 <- c(Y1,array(0,(gSize-N%%gSize)))
		Y2 <- c(Y2,array(0,(gSize-N%%gSize)))
	}
	
	Z1til <- apply(matrix(Y1til,nrow=gSize),2,max) 	# True pool statuses for disease 1
	Z2til <- apply(matrix(Y2til,nrow=gSize),2,max) 	# True pool statuses for disease 2
	seq.vec <- seq(1,N,gSize)
	Z1 <- data[ ,3][seq.vec]			# Pool testing results for disease 1
	Z2 <- data[ ,4][seq.vec]			# Pool testing results for disease 2
	Zj.sum <- apply(cbind(Z1,Z2),1,max) # Zj.sum is the indicator fn I(Z_1k+Z_2k>0)

	# Dirichlet prior for p
	p.prior <- p0*N0*a0
	alpha1 <- 1 + p.prior[1] 
	alpha2 <- 1 + p.prior[2]
	alpha3 <- 1 + p.prior[3]
	alpha4 <- 1 + p.prior[4]

	# Beta priors for Se1 and Se2
	beta11 <- 1 + TP1
	beta12 <- 1 + FN1 
	beta21 <- 1 + TP2
	beta22 <- 1 + FN2

	# Beta priors for Sp1 and Sp2
	gamma11 <- 1 + TN1 
	gamma12 <- 1 + FP1
	gamma21 <- 1 + TN2
	gamma22 <- 1 + FP2

	# Define matrices to save estimates before thinning
	p.save <- matrix(-9,nrow=GI,ncol=4)
	SeSp.save <- matrix(-9,nrow=GI,ncol=4)
	
	# Start iterating from the "full conditional" distributions
	#
	for(g in 1:GI){
	
		# Sample p from Dirichlet
		a1 <- alpha1 + sum((1-Y1til)*(1-Y2til))
		a2 <- alpha2 + sum(Y1til*(1-Y2til)) 
		a3 <- alpha3 + sum((1-Y1til)*Y2til) 
		a4 <- alpha4 + sum(Y1til*Y2til) 
		p <- rdirichlet(1,c(a1,a2,a3,a4))
		p.save[g, ] <- p

		# Sample Se1 from beta
		temp1 <- colSums(matrix(Y1*Y1til,nrow=gSize)) 
		tau11 <- beta11+sum(Z1*Z1til)+sum(Zj.sum*temp1)
		temp1comp <- colSums(matrix((1-Y1)*Y1til,nrow=gSize))
		tau12 <- beta12+sum((1-Z1)*Z1til)+sum(Zj.sum*temp1comp)
		Se1 <- rbeta(1,tau11,tau12)

		# Sample Se2 from beta
		temp2 <- colSums(matrix(Y2*Y2til,nrow=gSize)) 
		tau21 <- beta21+sum(Z2*Z2til)+sum(Zj.sum*temp2)
		temp2comp <- colSums(matrix((1-Y2)*Y2til,nrow=gSize))
		tau22 <- beta22+sum((1-Z2)*Z2til)+sum(Zj.sum*temp2comp)
		Se2 <- rbeta(1,tau21,tau22)

		# Sample Sp1 from beta
		temp3 <- colSums(matrix((1-Y1)*(1-Y1til),nrow=gSize))
		lambda11 <- gamma11+sum((1-Z1)*(1-Z1til))+sum(Zj.sum*temp3)
		temp3comp <- colSums(matrix(Y1*(1-Y1til),nrow=gSize))
		lambda12 <- gamma12+sum(Z1*(1-Z1til))+sum(Zj.sum*temp3comp)
		Sp1 <- rbeta(1,lambda11,lambda12)

		# Sample Sp2 from beta
		temp4 <- colSums(matrix((1-Y2)*(1-Y2til),nrow=gSize))
		lambda21 <- gamma21+sum((1-Z2)*(1-Z2til))+sum(Zj.sum*temp4)
		temp4comp <- colSums(matrix(Y2*(1-Y2til),nrow=gSize))
		lambda22 <- gamma22+sum(Z2*(1-Z2til))+sum(Zj.sum*temp4comp)
		Sp2 <- rbeta(1,lambda21,lambda22)
		
		SeSp.save[g, ] <- as.matrix(c(Se1,Se2,Sp1,Sp2))	
		Se <- c(Se1,Se2)
		Sp <- c(Sp1,Sp2)
	
		# Now sample Y1til and Y2til from multinomial.
		# Note that W is a N by 4 matrix; N is the sample size. 
		# W is not contributing anything to the estimation;
		# however, we need to specify W to run the code.
		# W was used in the program for a different purpose.
		#
		W <- matrix(0,nrow=N,ncol=4)
		u <- matrix(runif(N),nrow=N,ncol=1)
		res <- .C("multigibbs64",as.integer(N),as.integer(data),as.integer(dim(data)[2]),as.integer(1),as.integer(0),Se,Sp,p,as.integer(W),u)
		data <- matrix(res[[2]], nrow=N, ncol=dim(data)[2])

		Y1til <- data[ ,1]
		Y2til <- data[ ,2]
		if((N%%gSize) > 0){
			Y1til <- c(Y1til,array(0,(gSize-N%%gSize)))
			Y2til <- c(Y2til,array(0,(gSize-N%%gSize)))
		}

		Z1til <- apply(matrix(Y1til,nrow=gSize),2,max) 
		Z2til <- apply(matrix(Y2til,nrow=gSize),2,max)
	}
	
	par.sample <- cbind(p.save, SeSp.save)
	estm <- par.sample[(burn + 1):dim(par.sample)[1], ]	# discards the burn-in period
	seqn <- seq(1, dim(estm)[1], thin)		
	theta <- estm[seqn, ] 	# picks every "thin"th sample to avoid auto-correlation	
	colnames(theta) <- c("p00","p10","p01","p11","Se1","Se2","Sp1","Sp2")
	return(list("prevalence"=theta[ ,1:4], "accuracy"=theta[ ,5:8]))		
}

#****************************************************************************************#
# It's a support function for Dual.GT.Dorfman.Bayes()
# This function creates a special data structure required to 
# run the Dual.GT.Dorfman.Bayes() function
#
#
# Input:
# indiv.diag = N by 2 matrix of individual diagnoses (i.e., Stage 2 test results)
#			 = Matrix(Y_{i1k},Y_{i2k})
# group.resp = K by 2 matrix of pool testing responses, where K is the number of pools
#			 = Matrix(Z_{1k},Z_{2k})
# indiv.true = N by 2 matrix of individual true statuses
#			 = Matrix(\widetilde{Y}_{i1k},\widetilde{Y}_{i2k})
# c = pool size
#
#
# Output: 
# data = N by (7+c) data matrix 
#
IPP.structured.data <- function(indiv.diag,group.resp,indiv.true,c){
	N <- dim(indiv.diag)[1]
	data <- matrix(-99,nrow=N,ncol=(7+c))
	data[ ,c(1,2)] <- indiv.true
	data[ ,c(5,6)] <- indiv.diag	
	for(j in 1:ceiling(N/c)){
		upper <- min((j*c),N)
		lower <- ((j-1)*c+1)
		ind <- lower:upper
		psz <- length(ind)		
		data[ind,3] <- group.resp[j,1]	
		data[ind,4] <- group.resp[j,2]			
		data[ind,7] <- psz
		for(i in ind){
			data[i,(8:(7+psz))] <- ind
		}		
	}
	return(data)
}
#
##########################################################################################


#****************************************************************************************#
# This function simulates the IPP group testing data for two diseases.
# We provide this function ONLY to illustrate how our program works 
#
#
# Input: 
# N		= Sample size
# p 	= c(p00, p10, p01, p11) = Cell probabilities
# Se 	= c(Se1, Se2) = Sensitivities
# Sp 	= c(Sp1, Sp2) = Specificities
# c		= Pool size.
#
# Output: 
# data	= Matrix(\widetilde{Y}_{i1k},\widetilde{Y}_{i2k},Z_{1k},Z_{2k},Y_{i1k},Y_{i2k},c_k, c(Vector that identifies the other individuals in the kth pool))
# cost 	= Number of tests expended.
#
Dual.IPP.data <- function(N,p,Se,Sp,c){

	# Generate true individual statuses
	Y <- matrix(-99,nrow=N,ncol=2)
	for(i in 1:N){
		stat <- rmultinom(1,1,p)
		if(stat[1]==1){Y[i, ] <- c(0,0)}
		if(stat[2]==1){Y[i, ] <- c(1,0)}
		if(stat[3]==1){Y[i, ] <- c(0,1)}
		if(stat[4]==1){Y[i, ] <- c(1,1)}
	}
	
	# Create IPP group testing data
	data <- matrix(-99,nrow=N,ncol=(7+c))
	data[ ,c(1,2)] <- Y
	num.test <- 0
	for(j in 1:ceiling(N/c)){
		upper <- min((j*c),N)
		lower <- ((j-1)*c+1)
		ind <- lower:upper
		psz <- length(ind)			# Pool size
		prob1 <- ifelse(sum(Y[ind,1])>0,Se[1],(1-Sp[1]))
		prob2 <- ifelse(sum(Y[ind,2])>0,Se[2],(1-Sp[2]))
		Z1 <- rbinom(1,1,prob1)		# Pool testing result for disease 1
		Z2 <- rbinom(1,1,prob2)		# Pool testing result for disease 2
		num.test <- num.test+1
		data[ind,3] <- Z1
		data[ind,4] <- Z2
		data[ind,7] <- psz
		for(i in ind){
			data[i,(8:(7+psz))] <- ind
		}
		
		# Test results for all individuals of a pool are recorded as negative 
		# for both diseases if the pool tests negative for both diseases.
		#
		if(Z1==0 & Z2==0){
			data[ind,c(5,6)] <- 0
		}
		
		if(length(ind)==1){
			data[ind,c(5,6)] <- c(Z1,Z2)
		}
		
		# Retest individuals of a pool that tests positive for either disease.
		#
		if(length(ind)>1){
			if(Z1==1 | Z2==1){
				for(i in ind){
					prob1 <- ifelse(Y[i,1]>0,Se[1],(1-Sp[1]))
					prob2 <- ifelse(Y[i,2]>0,Se[2],(1-Sp[2]))
					Y1 <- rbinom(1,1,prob1) # Individual testing result for disease 1
					Y2 <- rbinom(1,1,prob2) # Individual testing result for disease 2
					num.test <- num.test + 1 
					data[i,c(5,6)] <- c(Y1,Y2)
				}
			}
		}
	}
	
	return(list("data"=data, "cost"=num.test))
}

#*********************************************************************
# For illustration, we simulate IPP group testing data 
# using the function "Dual.IPP.data"
# and then fit our Bayesian model using "Dual.GT.Dorfman.Bayes".
#
#
N <- 1000
p <- c(.80,.10,.09,.01)
Se <- c(.95,.95)
Sp <- c(.99,.99)
c <- 3

# Simulate data
#
IppData <- Dual.IPP.data(N,p,Se,Sp,c)$data
indiv_diag <- IppData[ ,c(5,6)]					# individual diagnoses
group_resp <- IppData[ ,c(3,4)][seq(1,N,c), ]	# group testing responses

# Note: to start iterating, we use individual diagnosed
# statuses as individual true statuses
#
indiv_true <- indiv_diag

# Set R's working directory
#
setwd(dir = "C:\\programs")

# Load the "MCMCpack" package which is required for the Dual.GT.Dorfman.Bayes() function
#
library(MCMCpack)			

# Fit the model with the following priors:
# Informative Dirichlet prior with a0 = 0.5 
# Flat priors for Se/Sp; that is, Se1, Se2, Sp1, Sp2 ~ beta (1, 1)
#
param <- Dual.GT.Dorfman.Bayes(indiv.diag=indiv_diag,group.resp=group_resp,indiv.true=indiv_true,gSize=c,p0=p,N0=N,a0=0.5,TP1=0,FN1=0,TN1=0,FP1=0,TP2=0,FN2=0,TN2=0,FP2=0,GI=6000,burn=1000,thin=5)

# Some estimation results
#
apply(param$prevalence, 2, median)		# Median estimates for prevalence parameters
apply(param$prevalence, 2, sd)			# Standard deviations for prevalence parameters

apply(param$accuracy, 2, median)		# Median estimates for accuracy parameters
apply(param$accuracy, 2, sd)		   	# Standard deviations for accuracy parameters

