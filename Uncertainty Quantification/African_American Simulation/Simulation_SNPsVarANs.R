#!/usr/bin/env Rscript

# Load required libraries 

require(Summix)
require(tidyverse)
require(data.table)
require(R.utils)
require(nloptr)
require(latex2exp)
require(gridExtra)

# Set working directory on server 

setwd("/newhome/agboolol/Gnome_Data")
# import files
files <- list.files(path="/newhome/agboolol/Gnome_Data", pattern="*.txt.gz")
data_gnom <- bind_rows(lapply(files, fread)) # row bind files
rm(files) # remove files



#######################################################################
#####  Ancestry estimation function
#######################################################################

#

# Loaded in manually so each core on cluster would have access to it

ancestr = function(refmatrix, obsvector){
  
  testmatrix = cbind(refmatrix, obsvector)
  starting = numeric(ncol(refmatrix))
  
  for (i in 1:(ncol(refmatrix))){
    starting[i] = 1/ncol(refmatrix)
  }
  
  fn.ancmix = function(x){
    minfunc = 0
    for (i in 1:ncol(refmatrix)){
      minfunc = minfunc + x[i]*testmatrix[i]
    }
    minfunc = minfunc - testmatrix[ncol(refmatrix) + 1]
    minfunc = sum((minfunc)**2)
    return(minfunc)
  }
  
  
  
  gr.ancmix <- function(x){
    gradvec = matrix(0,ncol(refmatrix),1)
    gradfunc = 0
    for (i in 1:ncol(refmatrix)){
      gradfunc = gradfunc + x[i]*testmatrix[i]
    }
    
    gradfunc = gradfunc - testmatrix[ncol(refmatrix) + 1]
    for (i in 1:ncol(refmatrix)){
      gradvec[i] = sum(2 * testmatrix[i] * gradfunc)
    }
    return(gradvec)
  }
  
  
  
  heq.ancmix = function(x){
    equality = 0
    for (i in 1:ncol(refmatrix)){
      equality = equality + x[i]
    }
    return(equality - 1)
  }
  
  hin.ancmix <- function(x){
    h = numeric(ncol(refmatrix))
    for (i in 1:ncol(refmatrix)){
      h[i] = x[i]
    }
    return(h)
  }
  
  start_time = Sys.time()
  
  suppressMessages(
    {S = slsqp(starting,
               fn = fn.ancmix,
               gr = gr.ancmix,
               hin = hin.ancmix,
               heq = heq.ancmix)
    }
  )
  end_time = Sys.time()
  ttime = end_time - start_time
  val = c( S$par,
           S$value,
           S$iter,
           ttime )
  
  return(val)
}


# transform to data frame
data_gnom <- as.data.frame(unclass(data_gnom))
data_red <- data_gnom[,c(1:8,92,94)] # extract required columns 

## SIMULATION ##

set.seed(23456) # Set seed
samp = 1000 # number of SNPs to be extracted - change to 10,000 and 100,000 to get desired plots

# Select 1K SNPS from reference data 
# Extract required columns 
refdat = data_red %>% 
  sample_n(samp) %>% 
  select(ref_AF_afr_1000G, ref_AF_eur_1000G, gnomad_AF_nfe, gnomad_AF_afr) %>% 
  rename(ref_AFR = ref_AF_afr_1000G, ref_EUR = ref_AF_eur_1000G, obs_EUR = gnomad_AF_nfe, obs_AFR = gnomad_AF_afr)


########################################################################
##### For N = 10, 20, 50, 100, 500. (25 possible pairwise combinations)
########################################################################

# Generate allele counts for population from binomial distribution for N1 = N2 = 10, P = Allele Frequency
# Repeat process for all pairwise combinations of N = 10, 20, 50, 100, 500 - comments not repeated for other combinations
# AF_ref_afr and AF_ref_eur are reference allele frequencies simulated or AF*
# refdat[,1] and refdat[,2] are reference allele frequencies from 1000G or AFobs
# N varies for each ancestry group

N1 = 10 # number of people (--> allele number = 2 * N1)
N2 = 10 # number of people (--> allele number = 2 * N2)
rep = 1000 # number of replicates 
AF_ref_afr <- t(sapply(refdat[,1], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})) 
AF_ref_eur <- t(sapply(refdat[,2], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)}))

## Creating storage for outputs
## Array stores multiple matrices - each matrix represents a dataset

Res_10a <- matrix(data = NA, ncol = 5, nrow = rep)
DAT_arr_AFR <- array(data=NA,dim = c(samp,3,rep))

for (i in 1:rep){
  
  DAT_arr_AFR[ , ,i] = cbind(AF_ref_afr[,i],AF_ref_eur[,i],refdat[,4]) # create a matrix of 3 columns - two ref AFs and one obs AF
  # Each matrix is a dataset to be run in summix. So we have 1000 datasets
  
  Res_10a[i,] = ancestr(as.data.frame(DAT_arr_AFR[,1:2,i]), as.data.frame(DAT_arr_AFR[,3,i])) # run summix on each dataset in the array
  # first entry is a matrix of ref AFs, second entry is a vector of observed AF
  
}

setwd("/newhome/agboolol")

Res10_afr <- as.data.frame(cbind(Anc_Prop = Res_10a[,1], N1 = rep("N=10",1000), N2 = rep("N=10",1000))) # extract ancestry proportion estimates for the two ancestry groups: AFR and EUR

LSEs_10a_1000 <- as.data.frame(cbind(LS = Res_10a[,3], N1 = rep("N=10",1000), N2 = rep("N=10",1000))) # extract least square values for each replicate

Diff10 <- as.data.frame(cbind(AFR=AF_ref_afr[,1]-refdat[,1], N1 = rep("N=10",1000), N2 = rep("N=10",1000))) # create dataframe for difference between simulated and observed AF_AFR


##########################################################################
# Repeat process above for N1=10, N2=20 
# Ditto N1=10, N2=50, ... , N1=10, N2=500 

N1 = 10 # number of people (--> allele number = 2 * N1)
N2 = 20 # number of people (--> allele number = 2 * N2)
rep = 1000 # number of replicates 
AF_ref_afr <- t(sapply(refdat[,1], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})) 
AF_ref_eur <- t(sapply(refdat[,2], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)}))

## Creating storage for outputs
## Array stores multiple matrices - each matrix represents a dataset

Res_1020a <- matrix(data = NA, ncol = 5, nrow = rep)
DAT_arr_AFR <- array(data=NA,dim = c(samp,3,rep))

for (i in 1:rep){
  
  DAT_arr_AFR[ , ,i] = cbind(AF_ref_afr[,i],AF_ref_eur[,i],refdat[,4]) # create a matrix of 3 columns - two ref AFs and one obs AF
  # Each matrix is a dataset to be run in summix. So we have 1000 datasets
  
  Res_1020a[i,] = ancestr(as.data.frame(DAT_arr_AFR[,1:2,i]), as.data.frame(DAT_arr_AFR[,3,i])) # run summix on each dataset in the array
  # first entry is a matrix of ref AFs, second entry is a vector of observed AF
  
}

Res1020_afr <- as.data.frame(cbind(Anc_Prop = Res_1020a[,1], N1 = rep("N=10",1000), N2 = rep("N=20",1000)))

LSEs_1020a_1000 <- as.data.frame(cbind(LS = Res_1020a[,3], N1 = rep("N=10",1000), N2 = rep("N=20",1000)))

Diff1020 <- as.data.frame(cbind(AFR=AF_ref_afr[,1]-refdat[,1], N1 = rep("N=10",1000), N2 = rep("N=20",1000))) 


##########################################################################

N1 = 10 # number of people (--> allele number = 2 * N1)
N2 = 50 # number of people (--> allele number = 2 * N2)
rep = 1000 # number of replicates 
AF_ref_afr <- t(sapply(refdat[,1], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})) 
AF_ref_eur <- t(sapply(refdat[,2], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)}))

## Creating storage for outputs
## Array stores multiple matrices - each matrix represents a dataset

Res_1050a <- matrix(data = NA, ncol = 5, nrow = rep)
DAT_arr_AFR <- array(data=NA,dim = c(samp,3,rep))

for (i in 1:rep){
  
  DAT_arr_AFR[ , ,i] = cbind(AF_ref_afr[,i],AF_ref_eur[,i],refdat[,4]) # create a matrix of 3 columns - two ref AFs and one obs AF
  # Each matrix is a dataset to be run in summix. So we have 1000 datasets
  
  Res_1050a[i,] = ancestr(as.data.frame(DAT_arr_AFR[,1:2,i]), as.data.frame(DAT_arr_AFR[,3,i])) # run summix on each dataset in the array
  # first entry is a matrix of ref AFs, second entry is a vector of observed AF
  
}


Res1050_afr <- as.data.frame(cbind(Anc_Prop = Res_1050a[,1], N1 = rep("N=10",1000), N2 = rep("N=50",1000)))

LSEs_1050a_1000 <- as.data.frame(cbind(LS = Res_1050a[,3], N1 = rep("N=10",1000), N2 = rep("N=50",1000)))

Diff1050 <- as.data.frame(cbind(AFR=AF_ref_afr[,1]-refdat[,1], N1 = rep("N=10",1000), N2 = rep("N=50",1000)))



##########################################################################

N1 = 10 # number of people (--> allele number = 2 * N1)
N2 = 100 # number of people (--> allele number = 2 * N2)
rep = 1000 # number of replicates 
AF_ref_afr <- t(sapply(refdat[,1], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})) 
AF_ref_eur <- t(sapply(refdat[,2], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)}))

## Creating storage for outputs
## Array stores multiple matrices - each matrix represents a dataset

Res_10100a <- matrix(data = NA, ncol = 5, nrow = rep)
DAT_arr_AFR <- array(data=NA,dim = c(samp,3,rep))

for (i in 1:rep){
  
  DAT_arr_AFR[ , ,i] = cbind(AF_ref_afr[,i],AF_ref_eur[,i],refdat[,4]) # create a matrix of 3 columns - two ref AFs and one obs AF
  # Each matrix is a dataset to be run in summix. So we have 1000 datasets
  
  Res_10100a[i,] = ancestr(as.data.frame(DAT_arr_AFR[,1:2,i]), as.data.frame(DAT_arr_AFR[,3,i])) # run summix on each dataset in the array
  # first entry is a matrix of ref AFs, second entry is a vector of observed AF
  
}


Res10100_afr <- as.data.frame(cbind(Anc_Prop = Res_10100a[,1], N1 = rep("N=10",1000), N2 = rep("N=100",1000)))

LSEs_10100a_1000 <- as.data.frame(cbind(LS = Res_10100a[,3], N1 = rep("N=10",1000), N2 = rep("N=100",1000)))

Diff10100 <- as.data.frame(cbind(AFR=AF_ref_afr[,1]-refdat[,1], N1 = rep("N=10",1000), N2 = rep("N=100",1000)))


##########################################################################

N1 = 10 # number of people (--> allele number = 2 * N1)
N2 = 500 # number of people (--> allele number = 2 * N2)
rep = 1000 # number of replicates 
AF_ref_afr <- t(sapply(refdat[,1], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})) 
AF_ref_eur <- t(sapply(refdat[,2], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)}))

## Creating storage for outputs
## Array stores multiple matrices - each matrix represents a dataset

Res_10500a <- matrix(data = NA, ncol = 5, nrow = rep)
DAT_arr_AFR <- array(data=NA,dim = c(samp,3,rep))

for (i in 1:rep){
  
  DAT_arr_AFR[ , ,i] = cbind(AF_ref_afr[,i],AF_ref_eur[,i],refdat[,4]) # create a matrix of 3 columns - two ref AFs and one obs AF
  # Each matrix is a dataset to be run in summix. So we have 1000 datasets
  
  Res_10500a[i,] = ancestr(as.data.frame(DAT_arr_AFR[,1:2,i]), as.data.frame(DAT_arr_AFR[,3,i])) # run summix on each dataset in the array
  # first entry is a matrix of ref AFs, second entry is a vector of observed AF
  
}


Res10500_afr <- as.data.frame(cbind(Anc_Prop = Res_10500a[,1], N1 = rep("N=10",1000), N2 = rep("N=500",1000)))

LSEs_10500a_1000 <- as.data.frame(cbind(LS = Res_10500a[,3], N1 = rep("N=10",1000), N2 = rep("N=500",1000)))

Diff10500 <- as.data.frame(cbind(AFR=AF_ref_afr[,1]-refdat[,1], N1 = rep("N=10",1000), N2 = rep("N=500",1000)))


##########################################################################


##########################################################################
# Repeat process above for N1=20, N2=10 
# Ditto N1=20, N2=20, ... , N1=20, N2=500 

N1 = 20 # number of people (--> allele number = 2 * N1)
N2 = 10 # number of people (--> allele number = 2 * N2)
rep = 1000 # number of replicates 
AF_ref_afr <- t(sapply(refdat[,1], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})) 
AF_ref_eur <- t(sapply(refdat[,2], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)}))

## Creating storage for outputs
## Array stores multiple matrices - each matrix represents a dataset

Res_2010a <- matrix(data = NA, ncol = 5, nrow = rep)
DAT_arr_AFR <- array(data=NA,dim = c(samp,3,rep))

for (i in 1:rep){
  
  DAT_arr_AFR[ , ,i] = cbind(AF_ref_afr[,i],AF_ref_eur[,i],refdat[,4]) # create a matrix of 3 columns - two ref AFs and one obs AF
  # Each matrix is a dataset to be run in summix. So we have 1000 datasets
  
  Res_2010a[i,] = ancestr(as.data.frame(DAT_arr_AFR[,1:2,i]), as.data.frame(DAT_arr_AFR[,3,i])) # run summix on each dataset in the array
  # first entry is a matrix of ref AFs, second entry is a vector of observed AF
  
}

Res2010_afr <- as.data.frame(cbind(Anc_Prop = Res_2010a[,1], N1 = rep("N=20",1000), N2 = rep("N=10",1000)))

LSEs_2010a_1000 <- as.data.frame(cbind(LS = Res_2010a[,3], N1 = rep("N=20",1000), N2 = rep("N=10",1000)))

Diff2010 <- as.data.frame(cbind(AFR=AF_ref_afr[,1]-refdat[,1], N1 = rep("N=20",1000), N2 = rep("N=10",1000))) 


##########################################################################

N1 = 20 # number of people (--> allele number = 2 * N1)
N2 = 20 # number of people (--> allele number = 2 * N2)
rep = 1000 # number of replicates 
AF_ref_afr <- t(sapply(refdat[,1], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})) 
AF_ref_eur <- t(sapply(refdat[,2], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)}))

## Creating storage for outputs
## Array stores multiple matrices - each matrix represents a dataset

Res_2020a <- matrix(data = NA, ncol = 5, nrow = rep)
DAT_arr_AFR <- array(data=NA,dim = c(samp,3,rep))

for (i in 1:rep){
  
  DAT_arr_AFR[ , ,i] = cbind(AF_ref_afr[,i],AF_ref_eur[,i],refdat[,4]) # create a matrix of 3 columns - two ref AFs and one obs AF
  # Each matrix is a dataset to be run in summix. So we have 1000 datasets
  
  Res_2020a[i,] = ancestr(as.data.frame(DAT_arr_AFR[,1:2,i]), as.data.frame(DAT_arr_AFR[,3,i])) # run summix on each dataset in the array
  # first entry is a matrix of ref AFs, second entry is a vector of observed AF
  
}

Res2020_afr <- as.data.frame(cbind(Anc_Prop = Res_2020a[,1], N1 = rep("N=20",1000), N2 = rep("N=20",1000)))

LSEs_2020a_1000 <- as.data.frame(cbind(LS = Res_2020a[,3], N1 = rep("N=20",1000), N2 = rep("N=20",1000)))

Diff2020 <- as.data.frame(cbind(AFR=AF_ref_afr[,1]-refdat[,1], N1 = rep("N=20",1000), N2 = rep("N=20",1000))) 


##########################################################################

N1 = 20 # number of people (--> allele number = 2 * N1)
N2 = 50 # number of people (--> allele number = 2 * N2)
rep = 1000 # number of replicates 
AF_ref_afr <- t(sapply(refdat[,1], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})) 
AF_ref_eur <- t(sapply(refdat[,2], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)}))

## Creating storage for outputs
## Array stores multiple matrices - each matrix represents a dataset

Res_2050a <- matrix(data = NA, ncol = 5, nrow = rep)
DAT_arr_AFR <- array(data=NA,dim = c(samp,3,rep))

for (i in 1:rep){
  
  DAT_arr_AFR[ , ,i] = cbind(AF_ref_afr[,i],AF_ref_eur[,i],refdat[,4]) # create a matrix of 3 columns - two ref AFs and one obs AF
  # Each matrix is a dataset to be run in summix. So we have 1000 datasets
  
  Res_2050a[i,] = ancestr(as.data.frame(DAT_arr_AFR[,1:2,i]), as.data.frame(DAT_arr_AFR[,3,i])) # run summix on each dataset in the array
  # first entry is a matrix of ref AFs, second entry is a vector of observed AF
  
}

Res2050_afr <- as.data.frame(cbind(Anc_Prop = Res_2050a[,1], N1 = rep("N=20",1000), N2 = rep("N=50",1000)))

LSEs_2050a_1000 <- as.data.frame(cbind(LS = Res_2050a[,3], N1 = rep("N=20",1000), N2 = rep("N=50",1000)))

Diff2050 <- as.data.frame(cbind(AFR=AF_ref_afr[,1]-refdat[,1], N1 = rep("N=20",1000), N2 = rep("N=50",1000))) 


##########################################################################

N1 = 20 # number of people (--> allele number = 2 * N1)
N2 = 100 # number of people (--> allele number = 2 * N2)
rep = 1000 # number of replicates 
AF_ref_afr <- t(sapply(refdat[,1], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})) 
AF_ref_eur <- t(sapply(refdat[,2], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)}))

## Creating storage for outputs
## Array stores multiple matrices - each matrix represents a dataset

Res_20100a <- matrix(data = NA, ncol = 5, nrow = rep)
DAT_arr_AFR <- array(data=NA,dim = c(samp,3,rep))

for (i in 1:rep){
  
  DAT_arr_AFR[ , ,i] = cbind(AF_ref_afr[,i],AF_ref_eur[,i],refdat[,4]) # create a matrix of 3 columns - two ref AFs and one obs AF
  # Each matrix is a dataset to be run in summix. So we have 1000 datasets
  
  Res_20100a[i,] = ancestr(as.data.frame(DAT_arr_AFR[,1:2,i]), as.data.frame(DAT_arr_AFR[,3,i])) # run summix on each dataset in the array
  # first entry is a matrix of ref AFs, second entry is a vector of observed AF
  
}

Res20100_afr <- as.data.frame(cbind(Anc_Prop = Res_20100a[,1], N1 = rep("N=20",1000), N2 = rep("N=100",1000)))

LSEs_20100a_1000 <- as.data.frame(cbind(LS = Res_20100a[,3], N1 = rep("N=20",1000), N2 = rep("N=100",1000)))

Diff20100 <- as.data.frame(cbind(AFR=AF_ref_afr[,1]-refdat[,1], N1 = rep("N=20",1000), N2 = rep("N=100",1000))) 


##########################################################################


N1 = 20 # number of people (--> allele number = 2 * N1)
N2 = 500 # number of people (--> allele number = 2 * N2)
rep = 1000 # number of replicates 
AF_ref_afr <- t(sapply(refdat[,1], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})) 
AF_ref_eur <- t(sapply(refdat[,2], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)}))

## Creating storage for outputs
## Array stores multiple matrices - each matrix represents a dataset

Res_20500a <- matrix(data = NA, ncol = 5, nrow = rep)
DAT_arr_AFR <- array(data=NA,dim = c(samp,3,rep))

for (i in 1:rep){
  
  DAT_arr_AFR[ , ,i] = cbind(AF_ref_afr[,i],AF_ref_eur[,i],refdat[,4]) # create a matrix of 3 columns - two ref AFs and one obs AF
  # Each matrix is a dataset to be run in summix. So we have 1000 datasets
  
  Res_20500a[i,] = ancestr(as.data.frame(DAT_arr_AFR[,1:2,i]), as.data.frame(DAT_arr_AFR[,3,i])) # run summix on each dataset in the array
  # first entry is a matrix of ref AFs, second entry is a vector of observed AF
  
}

Res20500_afr <- as.data.frame(cbind(Anc_Prop = Res_20500a[,1], N1 = rep("N=20",1000), N2 = rep("N=500",1000)))

LSEs_20500a_1000 <- as.data.frame(cbind(LS = Res_20500a[,3], N1 = rep("N=20",1000), N2 = rep("N=500",1000)))

Diff20500 <- as.data.frame(cbind(AFR=AF_ref_afr[,1]-refdat[,1], N1 = rep("N=20",1000), N2 = rep("N=500",1000))) 


##########################################################################


##########################################################################
# Repeat process above for N1=50, N2=10 
# Ditto N1=50, N2=20, ... , N1=50, N2=500 

N1 = 50 # number of people (--> allele number = 2 * N1)
N2 = 10 # number of people (--> allele number = 2 * N2)
rep = 1000 # number of replicates 
AF_ref_afr <- t(sapply(refdat[,1], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})) 
AF_ref_eur <- t(sapply(refdat[,2], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)}))

## Creating storage for outputs
## Array stores multiple matrices - each matrix represents a dataset

Res_5010a <- matrix(data = NA, ncol = 5, nrow = rep)
DAT_arr_AFR <- array(data=NA,dim = c(samp,3,rep))

for (i in 1:rep){
  
  DAT_arr_AFR[ , ,i] = cbind(AF_ref_afr[,i],AF_ref_eur[,i],refdat[,4]) # create a matrix of 3 columns - two ref AFs and one obs AF
  # Each matrix is a dataset to be run in summix. So we have 1000 datasets
  
  Res_5010a[i,] = ancestr(as.data.frame(DAT_arr_AFR[,1:2,i]), as.data.frame(DAT_arr_AFR[,3,i])) # run summix on each dataset in the array
  # first entry is a matrix of ref AFs, second entry is a vector of observed AF
  
}

Res5010_afr <- as.data.frame(cbind(Anc_Prop = Res_5010a[,1], N1 = rep("N=50",1000), N2 = rep("N=10",1000)))

LSEs_5010a_1000 <- as.data.frame(cbind(LS = Res_5010a[,3], N1 = rep("N=50",1000), N2 = rep("N=10",1000)))

Diff5010 <- as.data.frame(cbind(AFR=AF_ref_afr[,1]-refdat[,1], N1 = rep("N=50",1000), N2 = rep("N=10",1000))) 


##########################################################################

N1 = 50 # number of people (--> allele number = 2 * N1)
N2 = 20 # number of people (--> allele number = 2 * N2)
rep = 1000 # number of replicates 
AF_ref_afr <- t(sapply(refdat[,1], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})) 
AF_ref_eur <- t(sapply(refdat[,2], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)}))

## Creating storage for outputs
## Array stores multiple matrices - each matrix represents a dataset

Res_5020a <- matrix(data = NA, ncol = 5, nrow = rep)
DAT_arr_AFR <- array(data=NA,dim = c(samp,3,rep))

for (i in 1:rep){
  
  DAT_arr_AFR[ , ,i] = cbind(AF_ref_afr[,i],AF_ref_eur[,i],refdat[,4]) # create a matrix of 3 columns - two ref AFs and one obs AF
  # Each matrix is a dataset to be run in summix. So we have 1000 datasets
  
  Res_5020a[i,] = ancestr(as.data.frame(DAT_arr_AFR[,1:2,i]), as.data.frame(DAT_arr_AFR[,3,i])) # run summix on each dataset in the array
  # first entry is a matrix of ref AFs, second entry is a vector of observed AF
  
}

Res5020_afr <- as.data.frame(cbind(Anc_Prop = Res_5020a[,1], N1 = rep("N=50",1000), N2 = rep("N=20",1000)))

LSEs_5020a_1000 <- as.data.frame(cbind(LS = Res_5020a[,3], N1 = rep("N=50",1000), N2 = rep("N=20",1000)))

Diff5020 <- as.data.frame(cbind(AFR=AF_ref_afr[,1]-refdat[,1], N1 = rep("N=50",1000), N2 = rep("N=20",1000))) 


##########################################################################

N1 = 50 # number of people (--> allele number = 2 * N1)
N2 = 50 # number of people (--> allele number = 2 * N2)
rep = 1000 # number of replicates 
AF_ref_afr <- t(sapply(refdat[,1], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})) 
AF_ref_eur <- t(sapply(refdat[,2], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)}))

## Creating storage for outputs
## Array stores multiple matrices - each matrix represents a dataset

Res_5050a <- matrix(data = NA, ncol = 5, nrow = rep)
DAT_arr_AFR <- array(data=NA,dim = c(samp,3,rep))

for (i in 1:rep){
  
  DAT_arr_AFR[ , ,i] = cbind(AF_ref_afr[,i],AF_ref_eur[,i],refdat[,4]) # create a matrix of 3 columns - two ref AFs and one obs AF
  # Each matrix is a dataset to be run in summix. So we have 1000 datasets
  
  Res_5050a[i,] = ancestr(as.data.frame(DAT_arr_AFR[,1:2,i]), as.data.frame(DAT_arr_AFR[,3,i])) # run summix on each dataset in the array
  # first entry is a matrix of ref AFs, second entry is a vector of observed AF
  
}

Res5050_afr <- as.data.frame(cbind(Anc_Prop = Res_5050a[,1], N1 = rep("N=50",1000), N2 = rep("N=50",1000)))

LSEs_5050a_1000 <- as.data.frame(cbind(LS = Res_5050a[,3], N1 = rep("N=50",1000), N2 = rep("N=50",1000)))

Diff5050 <- as.data.frame(cbind(AFR=AF_ref_afr[,1]-refdat[,1], N1 = rep("N=50",1000), N2 = rep("N=50",1000))) 


##########################################################################

N1 = 50 # number of people (--> allele number = 2 * N1)
N2 = 100 # number of people (--> allele number = 2 * N2)
rep = 1000 # number of replicates 
AF_ref_afr <- t(sapply(refdat[,1], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})) 
AF_ref_eur <- t(sapply(refdat[,2], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)}))

## Creating storage for outputs
## Array stores multiple matrices - each matrix represents a dataset

Res_50100a <- matrix(data = NA, ncol = 5, nrow = rep)
DAT_arr_AFR <- array(data=NA,dim = c(samp,3,rep))

for (i in 1:rep){
  
  DAT_arr_AFR[ , ,i] = cbind(AF_ref_afr[,i],AF_ref_eur[,i],refdat[,4]) # create a matrix of 3 columns - two ref AFs and one obs AF
  # Each matrix is a dataset to be run in summix. So we have 1000 datasets
  
  Res_50100a[i,] = ancestr(as.data.frame(DAT_arr_AFR[,1:2,i]), as.data.frame(DAT_arr_AFR[,3,i])) # run summix on each dataset in the array
  # first entry is a matrix of ref AFs, second entry is a vector of observed AF
  
}

Res50100_afr <- as.data.frame(cbind(Anc_Prop = Res_50100a[,1], N1 = rep("N=50",1000), N2 = rep("N=100",1000)))

LSEs_50100a_1000 <- as.data.frame(cbind(LS = Res_50100a[,3], N1 = rep("N=50",1000), N2 = rep("N=100",1000)))

Diff50100 <- as.data.frame(cbind(AFR=AF_ref_afr[,1]-refdat[,1], N1 = rep("N=50",1000), N2 = rep("N=100",1000))) 


##########################################################################

N1 = 50 # number of people (--> allele number = 2 * N1)
N2 = 500 # number of people (--> allele number = 2 * N2)
rep = 1000 # number of replicates 
AF_ref_afr <- t(sapply(refdat[,1], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})) 
AF_ref_eur <- t(sapply(refdat[,2], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)}))

## Creating storage for outputs
## Array stores multiple matrices - each matrix represents a dataset

Res_50500a <- matrix(data = NA, ncol = 5, nrow = rep)
DAT_arr_AFR <- array(data=NA,dim = c(samp,3,rep))

for (i in 1:rep){
  
  DAT_arr_AFR[ , ,i] = cbind(AF_ref_afr[,i],AF_ref_eur[,i],refdat[,4]) # create a matrix of 3 columns - two ref AFs and one obs AF
  # Each matrix is a dataset to be run in summix. So we have 1000 datasets
  
  Res_50500a[i,] = ancestr(as.data.frame(DAT_arr_AFR[,1:2,i]), as.data.frame(DAT_arr_AFR[,3,i])) # run summix on each dataset in the array
  # first entry is a matrix of ref AFs, second entry is a vector of observed AF
  
}

Res50500_afr <- as.data.frame(cbind(Anc_Prop = Res_50500a[,1], N1 = rep("N=50",1000), N2 = rep("N=500",1000)))

LSEs_50500a_1000 <- as.data.frame(cbind(LS = Res_50500a[,3], N1 = rep("N=50",1000), N2 = rep("N=500",1000)))

Diff50500 <- as.data.frame(cbind(AFR=AF_ref_afr[,1]-refdat[,1], N1 = rep("N=50",1000), N2 = rep("N=500",1000))) 


##########################################################################



##########################################################################
# Repeat process above for N1=100, N2=10 
# Ditto N1=100, N2=20, ... , N1=100, N2=500 

N1 = 100 # number of people (--> allele number = 2 * N1)
N2 = 10 # number of people (--> allele number = 2 * N2)
rep = 1000 # number of replicates 
AF_ref_afr <- t(sapply(refdat[,1], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})) 
AF_ref_eur <- t(sapply(refdat[,2], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)}))

## Creating storage for outputs
## Array stores multiple matrices - each matrix represents a dataset

Res_10010a <- matrix(data = NA, ncol = 5, nrow = rep)
DAT_arr_AFR <- array(data=NA,dim = c(samp,3,rep))

for (i in 1:rep){
  
  DAT_arr_AFR[ , ,i] = cbind(AF_ref_afr[,i],AF_ref_eur[,i],refdat[,4]) # create a matrix of 3 columns - two ref AFs and one obs AF
  # Each matrix is a dataset to be run in summix. So we have 1000 datasets
  
  Res_10010a[i,] = ancestr(as.data.frame(DAT_arr_AFR[,1:2,i]), as.data.frame(DAT_arr_AFR[,3,i])) # run summix on each dataset in the array
  # first entry is a matrix of ref AFs, second entry is a vector of observed AF
  
}

Res10010_afr <- as.data.frame(cbind(Anc_Prop = Res_10010a[,1], N1 = rep("N=100",1000), N2 = rep("N=10",1000)))

LSEs_10010a_1000 <- as.data.frame(cbind(LS = Res_10010a[,3], N1 = rep("N=100",1000), N2 = rep("N=10",1000)))

Diff10010 <- as.data.frame(cbind(AFR=AF_ref_afr[,1]-refdat[,1], N1 = rep("N=100",1000), N2 = rep("N=10",1000))) 


##########################################################################


N1 = 100 # number of people (--> allele number = 2 * N1)
N2 = 20 # number of people (--> allele number = 2 * N2)
rep = 1000 # number of replicates 
AF_ref_afr <- t(sapply(refdat[,1], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})) 
AF_ref_eur <- t(sapply(refdat[,2], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)}))

## Creating storage for outputs
## Array stores multiple matrices - each matrix represents a dataset

Res_10020a <- matrix(data = NA, ncol = 5, nrow = rep)
DAT_arr_AFR <- array(data=NA,dim = c(samp,3,rep))

for (i in 1:rep){
  
  DAT_arr_AFR[ , ,i] = cbind(AF_ref_afr[,i],AF_ref_eur[,i],refdat[,4]) # create a matrix of 3 columns - two ref AFs and one obs AF
  # Each matrix is a dataset to be run in summix. So we have 1000 datasets
  
  Res_10020a[i,] = ancestr(as.data.frame(DAT_arr_AFR[,1:2,i]), as.data.frame(DAT_arr_AFR[,3,i])) # run summix on each dataset in the array
  # first entry is a matrix of ref AFs, second entry is a vector of observed AF
  
}

Res10020_afr <- as.data.frame(cbind(Anc_Prop = Res_10020a[,1], N1 = rep("N=100",1000), N2 = rep("N=20",1000)))

LSEs_10020a_1000 <- as.data.frame(cbind(LS = Res_10020a[,3], N1 = rep("N=100",1000), N2 = rep("N=20",1000)))

Diff10020 <- as.data.frame(cbind(AFR=AF_ref_afr[,1]-refdat[,1], N1 = rep("N=100",1000), N2 = rep("N=20",1000))) 


##########################################################################

N1 = 100 # number of people (--> allele number = 2 * N1)
N2 = 50 # number of people (--> allele number = 2 * N2)
rep = 1000 # number of replicates 
AF_ref_afr <- t(sapply(refdat[,1], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})) 
AF_ref_eur <- t(sapply(refdat[,2], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)}))

## Creating storage for outputs
## Array stores multiple matrices - each matrix represents a dataset

Res_10050a <- matrix(data = NA, ncol = 5, nrow = rep)
DAT_arr_AFR <- array(data=NA,dim = c(samp,3,rep))

for (i in 1:rep){
  
  DAT_arr_AFR[ , ,i] = cbind(AF_ref_afr[,i],AF_ref_eur[,i],refdat[,4]) # create a matrix of 3 columns - two ref AFs and one obs AF
  # Each matrix is a dataset to be run in summix. So we have 1000 datasets
  
  Res_10050a[i,] = ancestr(as.data.frame(DAT_arr_AFR[,1:2,i]), as.data.frame(DAT_arr_AFR[,3,i])) # run summix on each dataset in the array
  # first entry is a matrix of ref AFs, second entry is a vector of observed AF
  
}

Res10050_afr <- as.data.frame(cbind(Anc_Prop = Res_10050a[,1], N1 = rep("N=100",1000), N2 = rep("N=50",1000)))

LSEs_10050a_1000 <- as.data.frame(cbind(LS = Res_10050a[,3], N1 = rep("N=100",1000), N2 = rep("N=50",1000)))

Diff10050 <- as.data.frame(cbind(AFR=AF_ref_afr[,1]-refdat[,1], N1 = rep("N=100",1000), N2 = rep("N=50",1000))) 


##########################################################################

N1 = 100 # number of people (--> allele number = 2 * N1)
N2 = 100 # number of people (--> allele number = 2 * N2)
rep = 1000 # number of replicates 
AF_ref_afr <- t(sapply(refdat[,1], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})) 
AF_ref_eur <- t(sapply(refdat[,2], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)}))

## Creating storage for outputs
## Array stores multiple matrices - each matrix represents a dataset

Res_100100a <- matrix(data = NA, ncol = 5, nrow = rep)
DAT_arr_AFR <- array(data=NA,dim = c(samp,3,rep))

for (i in 1:rep){
  
  DAT_arr_AFR[ , ,i] = cbind(AF_ref_afr[,i],AF_ref_eur[,i],refdat[,4]) # create a matrix of 3 columns - two ref AFs and one obs AF
  # Each matrix is a dataset to be run in summix. So we have 1000 datasets
  
  Res_100100a[i,] = ancestr(as.data.frame(DAT_arr_AFR[,1:2,i]), as.data.frame(DAT_arr_AFR[,3,i])) # run summix on each dataset in the array
  # first entry is a matrix of ref AFs, second entry is a vector of observed AF
  
}

Res100100_afr <- as.data.frame(cbind(Anc_Prop = Res_100100a[,1], N1 = rep("N=100",1000), N2 = rep("N=100",1000)))

LSEs_100100a_1000 <- as.data.frame(cbind(LS = Res_100100a[,3], N1 = rep("N=100",1000), N2 = rep("N=100",1000)))

Diff100100 <- as.data.frame(cbind(AFR=AF_ref_afr[,1]-refdat[,1], N1 = rep("N=100",1000), N2 = rep("N=100",1000))) 


##########################################################################

N1 = 100 # number of people (--> allele number = 2 * N1)
N2 = 500 # number of people (--> allele number = 2 * N2)
rep = 1000 # number of replicates 
AF_ref_afr <- t(sapply(refdat[,1], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})) 
AF_ref_eur <- t(sapply(refdat[,2], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)}))

## Creating storage for outputs
## Array stores multiple matrices - each matrix represents a dataset

Res_100500a <- matrix(data = NA, ncol = 5, nrow = rep)
DAT_arr_AFR <- array(data=NA,dim = c(samp,3,rep))

for (i in 1:rep){
  
  DAT_arr_AFR[ , ,i] = cbind(AF_ref_afr[,i],AF_ref_eur[,i],refdat[,4]) # create a matrix of 3 columns - two ref AFs and one obs AF
  # Each matrix is a dataset to be run in summix. So we have 1000 datasets
  
  Res_100500a[i,] = ancestr(as.data.frame(DAT_arr_AFR[,1:2,i]), as.data.frame(DAT_arr_AFR[,3,i])) # run summix on each dataset in the array
  # first entry is a matrix of ref AFs, second entry is a vector of observed AF
  
}

Res100500_afr <- as.data.frame(cbind(Anc_Prop = Res_100500a[,1], N1 = rep("N=100",1000), N2 = rep("N=500",1000)))

LSEs_100500a_1000 <- as.data.frame(cbind(LS = Res_100500a[,3], N1 = rep("N=100",1000), N2 = rep("N=500",1000)))

Diff100500 <- as.data.frame(cbind(AFR=AF_ref_afr[,1]-refdat[,1], N1 = rep("N=100",1000), N2 = rep("N=500",1000))) 


##########################################################################


##########################################################################
# Repeat process above for N1=100, N2=10 
# Ditto N1=100, N2=20, ... , N1=100, N2=500 

N1 = 500 # number of people (--> allele number = 2 * N1)
N2 = 10 # number of people (--> allele number = 2 * N2)
rep = 1000 # number of replicates 
AF_ref_afr <- t(sapply(refdat[,1], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})) 
AF_ref_eur <- t(sapply(refdat[,2], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)}))

## Creating storage for outputs
## Array stores multiple matrices - each matrix represents a dataset

Res_50010a <- matrix(data = NA, ncol = 5, nrow = rep)
DAT_arr_AFR <- array(data=NA,dim = c(samp,3,rep))

for (i in 1:rep){
  
  DAT_arr_AFR[ , ,i] = cbind(AF_ref_afr[,i],AF_ref_eur[,i],refdat[,4]) # create a matrix of 3 columns - two ref AFs and one obs AF
  # Each matrix is a dataset to be run in summix. So we have 1000 datasets
  
  Res_50010a[i,] = ancestr(as.data.frame(DAT_arr_AFR[,1:2,i]), as.data.frame(DAT_arr_AFR[,3,i])) # run summix on each dataset in the array
  # first entry is a matrix of ref AFs, second entry is a vector of observed AF
  
}

Res50010_afr <- as.data.frame(cbind(Anc_Prop = Res_50010a[,1], N1 = rep("N=500",1000), N2 = rep("N=10",1000)))

LSEs_50010a_1000 <- as.data.frame(cbind(LS = Res_50010a[,3], N1 = rep("N=500",1000), N2 = rep("N=10",1000)))

Diff50010 <- as.data.frame(cbind(AFR=AF_ref_afr[,1]-refdat[,1], N1 = rep("N=500",1000), N2 = rep("N=10",1000))) 


##########################################################################

N1 = 500 # number of people (--> allele number = 2 * N1)
N2 = 20 # number of people (--> allele number = 2 * N2)
rep = 1000 # number of replicates 
AF_ref_afr <- t(sapply(refdat[,1], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})) 
AF_ref_eur <- t(sapply(refdat[,2], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)}))

## Creating storage for outputs
## Array stores multiple matrices - each matrix represents a dataset

Res_50020a <- matrix(data = NA, ncol = 5, nrow = rep)
DAT_arr_AFR <- array(data=NA,dim = c(samp,3,rep))

for (i in 1:rep){
  
  DAT_arr_AFR[ , ,i] = cbind(AF_ref_afr[,i],AF_ref_eur[,i],refdat[,4]) # create a matrix of 3 columns - two ref AFs and one obs AF
  # Each matrix is a dataset to be run in summix. So we have 1000 datasets
  
  Res_50020a[i,] = ancestr(as.data.frame(DAT_arr_AFR[,1:2,i]), as.data.frame(DAT_arr_AFR[,3,i])) # run summix on each dataset in the array
  # first entry is a matrix of ref AFs, second entry is a vector of observed AF
  
}

Res50020_afr <- as.data.frame(cbind(Anc_Prop = Res_50020a[,1], N1 = rep("N=500",1000), N2 = rep("N=20",1000)))

LSEs_50020a_1000 <- as.data.frame(cbind(LS = Res_50020a[,3], N1 = rep("N=500",1000), N2 = rep("N=20",1000)))

Diff50020 <- as.data.frame(cbind(AFR=AF_ref_afr[,1]-refdat[,1], N1 = rep("N=500",1000), N2 = rep("N=20",1000))) 


##########################################################################

N1 = 500 # number of people (--> allele number = 2 * N1)
N2 = 50 # number of people (--> allele number = 2 * N2)
rep = 1000 # number of replicates 
AF_ref_afr <- t(sapply(refdat[,1], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})) 
AF_ref_eur <- t(sapply(refdat[,2], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)}))

## Creating storage for outputs
## Array stores multiple matrices - each matrix represents a dataset

Res_50050a <- matrix(data = NA, ncol = 5, nrow = rep)
DAT_arr_AFR <- array(data=NA,dim = c(samp,3,rep))

for (i in 1:rep){
  
  DAT_arr_AFR[ , ,i] = cbind(AF_ref_afr[,i],AF_ref_eur[,i],refdat[,4]) # create a matrix of 3 columns - two ref AFs and one obs AF
  # Each matrix is a dataset to be run in summix. So we have 1000 datasets
  
  Res_50050a[i,] = ancestr(as.data.frame(DAT_arr_AFR[,1:2,i]), as.data.frame(DAT_arr_AFR[,3,i])) # run summix on each dataset in the array
  # first entry is a matrix of ref AFs, second entry is a vector of observed AF
  
}

Res50050_afr <- as.data.frame(cbind(Anc_Prop = Res_50050a[,1], N1 = rep("N=500",1000), N2 = rep("N=50",1000)))

LSEs_50050a_1000 <- as.data.frame(cbind(LS = Res_50050a[,3], N1 = rep("N=500",1000), N2 = rep("N=50",1000)))

Diff50050 <- as.data.frame(cbind(AFR=AF_ref_afr[,1]-refdat[,1], N1 = rep("N=500",1000), N2 = rep("N=50",1000))) 


##########################################################################


N1 = 500 # number of people (--> allele number = 2 * N1)
N2 = 100 # number of people (--> allele number = 2 * N2)
rep = 1000 # number of replicates 
AF_ref_afr <- t(sapply(refdat[,1], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})) 
AF_ref_eur <- t(sapply(refdat[,2], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)}))

## Creating storage for outputs
## Array stores multiple matrices - each matrix represents a dataset

Res_500100a <- matrix(data = NA, ncol = 5, nrow = rep)
DAT_arr_AFR <- array(data=NA,dim = c(samp,3,rep))

for (i in 1:rep){
  
  DAT_arr_AFR[ , ,i] = cbind(AF_ref_afr[,i],AF_ref_eur[,i],refdat[,4]) # create a matrix of 3 columns - two ref AFs and one obs AF
  # Each matrix is a dataset to be run in summix. So we have 1000 datasets
  
  Res_500100a[i,] = ancestr(as.data.frame(DAT_arr_AFR[,1:2,i]), as.data.frame(DAT_arr_AFR[,3,i])) # run summix on each dataset in the array
  # first entry is a matrix of ref AFs, second entry is a vector of observed AF
  
}

Res500100_afr <- as.data.frame(cbind(Anc_Prop = Res_500100a[,1], N1 = rep("N=500",1000), N2 = rep("N=100",1000)))

LSEs_500100a_1000 <- as.data.frame(cbind(LS = Res_500100a[,3], N1 = rep("N=500",1000), N2 = rep("N=100",1000)))

Diff500100 <- as.data.frame(cbind(AFR=AF_ref_afr[,1]-refdat[,1], N1 = rep("N=500",1000), N2 = rep("N=100",1000))) 


##########################################################################


N1 = 500 # number of people (--> allele number = 2 * N1)
N2 = 500 # number of people (--> allele number = 2 * N2)
rep = 1000 # number of replicates 
AF_ref_afr <- t(sapply(refdat[,1], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})) 
AF_ref_eur <- t(sapply(refdat[,2], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)}))

## Creating storage for outputs
## Array stores multiple matrices - each matrix represents a dataset

Res_500500a <- matrix(data = NA, ncol = 5, nrow = rep)
DAT_arr_AFR <- array(data=NA,dim = c(samp,3,rep))

for (i in 1:rep){
  
  DAT_arr_AFR[ , ,i] = cbind(AF_ref_afr[,i],AF_ref_eur[,i],refdat[,4]) # create a matrix of 3 columns - two ref AFs and one obs AF
  # Each matrix is a dataset to be run in summix. So we have 1000 datasets
  
  Res_500500a[i,] = ancestr(as.data.frame(DAT_arr_AFR[,1:2,i]), as.data.frame(DAT_arr_AFR[,3,i])) # run summix on each dataset in the array
  # first entry is a matrix of ref AFs, second entry is a vector of observed AF
  
}

Res500500_afr <- as.data.frame(cbind(Anc_Prop = Res_500500a[,1], N1 = rep("N=500",1000), N2 = rep("N=500",1000)))

LSEs_500500a_1000 <- as.data.frame(cbind(LS = Res_500500a[,3], N1 = rep("N=500",1000), N2 = rep("N=500",1000)))

Diff500500 <- as.data.frame(cbind(AFR=AF_ref_afr[,1]-refdat[,1], N1 = rep("N=500",1000), N2 = rep("N=500",1000))) 


##########################################################################

#########################################################################
## Combine all data frames (outputs) for every combination of N  
## for ancestry proportions, least squares values, and difference in AFs
#########################################################################

df.combine <- rbind(Res10_afr, Res1020_afr, Res1050_afr, Res10100_afr, Res10500_afr,
                    Res2010_afr, Res2020_afr, Res2050_afr, Res20100_afr, Res20500_afr,
                    Res5010_afr, Res5020_afr, Res5050_afr, Res50100_afr, Res50500_afr,
                    Res10010_afr, Res10020_afr, Res10050_afr, Res100100_afr, Res100500_afr,
                    Res50010_afr, Res50020_afr, Res50050_afr, Res500100_afr, Res500500_afr) # Combining all ancestry proportions estimates for varying N
df.combine$N1 <- factor(df.combine$N1, levels = c("N=10", "N=20", "N=50", "N=100", "N=500"))
df.combine$N2 <- factor(df.combine$N2, levels = c("N=10", "N=20", "N=50", "N=100", "N=500"))

LS.combine <- rbind(LSEs_10a_1000, LSEs_1020a_1000, LSEs_1050a_1000, LSEs_10100a_1000, LSEs_10500a_1000,
                    LSEs_2010a_1000, LSEs_2020a_1000, LSEs_2050a_1000, LSEs_20100a_1000, LSEs_20500a_1000,
                    LSEs_5010a_1000, LSEs_5020a_1000, LSEs_5050a_1000, LSEs_50100a_1000, LSEs_50500a_1000,
                    LSEs_10010a_1000, LSEs_10020a_1000, LSEs_10050a_1000, LSEs_100100a_1000, LSEs_100500a_1000,
                    LSEs_50010a_1000, LSEs_50020a_1000, LSEs_50050a_1000, LSEs_500100a_1000, LSEs_500500a_1000) # combining all least square values for varying N
LS.combine$N1 <- factor(LS.combine$N1, levels = c("N=10", "N=20", "N=50", "N=100", "N=500"))
LS.combine$N2 <- factor(LS.combine$N2, levels = c("N=10", "N=20", "N=50", "N=100", "N=500"))

DF.combine <- rbind(Diff10, Diff1020, Diff1050, Diff10100, Diff10500,
                    Diff2010, Diff2020, Diff2050, Diff20100, Diff20500,
                    Diff5010, Diff5020, Diff5050, Diff50100, Diff50500,
                    Diff10010, Diff10020, Diff10050, Diff100100, Diff100500,
                    Diff50010, Diff50020, Diff50050, Diff500100, Diff500500) # Combining difference in AFs for varying N
DF.combine$N1 <- factor(DF.combine$N1, levels = c("N=10", "N=20", "N=50", "N=100", "N=500"))
DF.combine$N2 <- factor(DF.combine$N2, levels = c("N=10", "N=20", "N=50", "N=100", "N=500"))

## Create plot for ancestry proportion estimates for varying N

pdf('Rplot_SNP1000_Anc_VarN.pdf')
ggplot(df.combine, aes(x = as.factor(N1), y = as.numeric(Anc_Prop), color = as.factor(N2))) + 
  geom_boxplot() + labs(x = "AFR: Allele Number (2 * N)", y = "Ancestry Proportion")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Variation of Summix Ancestry Proportion Estimate (SNP = 1,000) \n African/African American GnomAD")+
  scale_x_discrete(limits=c("N=10","N=20","N=50","N=100","N=500"))+
  labs(color = "EUR: \nAllele Number (2 * N)")
dev.off()


## Create plot for least square values for varying N

pdf('Rplot_SNP1000_LS_VarN.pdf')
ggplot(LS.combine, aes(x = as.factor(N1), y = as.numeric(LS), color = as.factor(N2))) + 
  geom_boxplot() + labs(x = "AFR: Allele Number (2 * N)", y = "Least Square Values")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Distribution of Least Square Values (SNP = 1,000) \n African/African American GnomAD")+
  scale_x_discrete(limits=c("N=10","N=20","N=50","N=100","N=500"))+
  labs(color = "EUR: \nAllele Number (2 * N)")
dev.off()


## Create plot for difference in AFs for varying N

pdf('Rplot_SNP1000_Dif_VarN.pdf')
ggplot(DF.combine, aes(x = as.factor(N1), y = as.numeric(AFR), color = as.factor(N2))) + 
  geom_boxplot() + labs(x = "AFR: Allele Number (2 * N)", y = "Difference in Allele Frequencies")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Distribution of Difference in Allele Frequencies (SNP = 1,000) \n African/African American GnomAD")+
  scale_x_discrete(limits=c("N=10","N=20","N=50","N=100","N=500"))+
  labs(color = "EUR: \nAllele Number (2 * N)")
dev.off()
