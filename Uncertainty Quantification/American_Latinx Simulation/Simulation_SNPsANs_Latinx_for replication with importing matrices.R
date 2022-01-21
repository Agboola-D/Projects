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


# select relevant columns (observed and reference data) and preparing a dataset (refdat) for Latinx analyses
data_red <- data_gnom[,c(1:11, 93)]

## SIMULATION ##

# Set Seed
set.seed(23456)

samp = 1000 #number of SNPs to be extracted

refdat = data_red %>% 
          sample_n(samp) %>% 
          select(ref_AF_eur_1000G,ref_AF_afr_1000G, ref_AF_sas_1000G,ref_AF_eas_1000G,ref_AF_iam_1000G, gnomad_AF_amr) %>% 
          
          rename(ref_AFR = ref_AF_afr_1000G, ref_EUR = ref_AF_eur_1000G, ref_SAS = ref_AF_sas_1000G, ref_EAS = ref_AF_eas_1000G,
                                                                         ref_IAM = ref_AF_iam_1000G, obs_AMR = gnomad_AF_amr)


##########################################################################
##### Codes for N = 5, 10, 15, 20, 25, 10, 50, 100, 500, 1000
##########################################################################

##########################################################################
# Generate allele counts for population from binomial distribution 
# for N = 5, P = Allele Frequency
# 
# Repeat process for each N - comments not repeated for other Ns
# AF_ref_*** are reference allele frequency simulated for the reference 
# ancestry groups or AF*
# refdat[,i] are reference allele frequencies from 1000G or AFobs
# N is fixed for each ancestry group
#
# There was a little issue with the function used below that we could not resolve with the Math server,
# so we worked around it.
# We needed 1000 x 1000 matrices, but the server produced 1 x 1000 vectors.
# Our R produced matrices but the server produced vectors.
# So we collected the matrices from ours, saved them as datasets,
# and used them on the server because the server runs way faster than our machines
############################################################################################################

N = 5 # number of people (--> allele number = 2 * N)
rep = 1000 # number of replicates 

# Created the AFs and saved them as datasets on our machine
# Comments not repeated for each N

# AF_ref_eur_5_1000 <- t(sapply(refdat[,1], function(x){rbinom(rep, 2 * N, x)/ (2 * N)})); write.table(AF_ref_eur_5_1000,file = "AF_ref_eur_5_1000.txt", sep = ' ')
# AF_ref_afr_5_1000 <- t(sapply(refdat[,2], function(x){rbinom(rep, 2 * N, x)/ (2 * N)})); write.table(AF_ref_afr_5_1000,file = "AF_ref_afr_5_1000.txt", sep = ' ') 
# AF_ref_sas_5_1000 <- t(sapply(refdat[,3], function(x){rbinom(rep, 2 * N, x)/ (2 * N)})); write.table(AF_ref_sas_5_1000,file = "AF_ref_sas_5_1000.txt", sep = ' ')
# AF_ref_eas_5_1000 <- t(sapply(refdat[,4], function(x){rbinom(rep, 2 * N, x)/ (2 * N)})); write.table(AF_ref_eas_5_1000,file = "AF_ref_eas_5_1000.txt", sep = ' ')
# AF_ref_iam_5_1000 <- t(sapply(refdat[,5], function(x){rbinom(rep, 2 * N, x)/ (2 * N)})); write.table(AF_ref_iam_5_1000,file = "AF_ref_iam_5_1000.txt", sep = ' ')

# importing the matrices as datasets - see OneDrive for "AF_matrices_Latinx" folder. It contains the matrices
# comment not repeated for each N

AF_ref_eur_5_1000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_eur_5_1000.txt", sep=" ")) 
AF_ref_afr_5_1000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_afr_5_1000.txt", sep=" "))
AF_ref_sas_5_1000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_sas_5_1000.txt", sep=" "))
AF_ref_eas_5_1000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_eas_5_1000.txt", sep=" "))
AF_ref_iam_5_1000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_iam_5_1000.txt", sep=" "))


## Creating storage for outputs
Res_5a <- matrix(data = NA, ncol = 8, nrow = rep)

## Array stores multiple matrices - each matrix represents a dataset
DAT_arr_AMR <- array(data=NA,dim = c(samp,6,rep))

##
for (i in 1:rep){
  # create a matrix of 6 columns - five ref AFs and one obs AF
  # Each matrix is a dataset to be run in summix. So we have 1000 datasets
  DAT_arr_AMR[ , ,i] = cbind(AF_ref_eur_5_1000[,i], AF_ref_afr_5_1000[,i], AF_ref_sas_5_1000[,i], AF_ref_eas_5_1000[,i],  AF_ref_iam_5_1000[,i], refdat[,6])
  
  # run summix on each dataset in the array
  # first entry is a matrix of ref AFs, second entry is a vector of observed AF
  Res_5a[i,] = ancestr(as.data.frame(DAT_arr_AMR[,1:5,i]), as.data.frame(DAT_arr_AMR[,6,i])) 
}

# extract ancestry proportion estimates for the five ancestry groups and the rename columns
Res5_amr <- as.data.frame(cbind(Res_5a[,1:5])) 
colnames(Res5_amr) <- c("Anc_Prop_EUR","Anc_Prop_AFR", "Anc_Prop_SAS", "Anc_Prop_EAS", "Anc_Prop_IAM")

## extract least square values for each replicate
LSEs_5amr_1000 <- as.data.frame(cbind(LS = Res_5a[,6], N = rep("N=5", 5000))) # extract least square values for each replicate


## Creating dataset for plots (convert data from wide to long format and create column to specify N) 
df5 <- reshape(Res5_amr, times = c("EUR", "AFR","SAS", "EAS","IAM"), timevar = "Anc_Grp",
               varying = list(names(Res5_amr)), direction = "long") 

df5$size <- rep("N=5",5000)


## Dataframe for AFstar Minus AFobs (difference between simulated and observed AF for each ancestry group)
Diff <- as.data.frame(cbind(EUR=AF_ref_eur_5_1000[,1]-refdat[,1],   AFR=AF_ref_afr_5_1000[,1]-refdat[,2],  SAS=AF_ref_sas_5_1000[,1]-refdat[,3], 
                            EAS=AF_ref_eas_5_1000[,1] - refdat[,4],  IAM=AF_ref_iam_5_1000[,1]- refdat[,5]))

## convert data from wide to long format and create column to specify N
DF1 <- reshape(Diff, times = c("EUR","AFR","SAS","EAS","IAM"), timevar = "Anc_Grp",
               varying = list(names(Diff)), direction = "long") # convert wide to long. 
DF1$size <- rep("N=5", 5000)


##########################################################################
# Generate allele counts for population using binom for N = 10
# All other manipulations are similar to the above
##########################################################################
N = 10 # number of people
rep = 1000 # number of replicates 

# AF_ref_eur_10_1000 <- t(sapply(refdat[,1], function(x){rbinom(rep, 2 * N, x)/ (2 * N)})); write.table(AF_ref_eur_10_1000,file = "AF_ref_eur_10_1000.txt", sep = ' ')
# AF_ref_afr_10_1000 <- t(sapply(refdat[,2], function(x){rbinom(rep, 2 * N, x)/ (2 * N)})); write.table(AF_ref_afr_10_1000,file = "AF_ref_afr_10_1000.txt", sep = ' ')
# AF_ref_sas_10_1000 <- t(sapply(refdat[,3], function(x){rbinom(rep, 2 * N, x)/ (2 * N)})); write.table(AF_ref_sas_10_1000,file = "AF_ref_sas_10_1000.txt", sep = ' ')
# AF_ref_eas_10_1000 <- t(sapply(refdat[,4], function(x){rbinom(rep, 2 * N, x)/ (2 * N)})); write.table(AF_ref_eas_10_1000,file = "AF_ref_eas_10_1000.txt", sep = ' ')
# AF_ref_iam_10_1000 <- t(sapply(refdat[,5], function(x){rbinom(rep, 2 * N, x)/ (2 * N)})); write.table(AF_ref_iam_10_1000,file = "AF_ref_iam_10_1000.txt", sep = ' ')

AF_ref_eur_10_1000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_eur_10_1000.txt", sep=" ")) 
AF_ref_afr_10_1000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_afr_10_1000.txt", sep=" "))
AF_ref_sas_10_1000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_sas_10_1000.txt", sep=" "))
AF_ref_eas_10_1000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_eas_10_1000.txt", sep=" "))
AF_ref_iam_10_1000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_iam_10_1000.txt", sep=" "))


Res_10a <- matrix(data = NA, ncol = 8, nrow = rep)

DAT_arr_AMR <- array(data=NA,dim = c(samp,6,rep))

##
for (i in 1:rep){
  
  DAT_arr_AMR[ , ,i] = cbind(AF_ref_eur_10_1000[,i], AF_ref_afr_10_1000[,i], AF_ref_sas_10_1000[,i], AF_ref_eas_10_1000[,i],  AF_ref_iam_10_1000[,i], refdat[,6])
  
  Res_10a[i,] = ancestr(as.data.frame(DAT_arr_AMR[,1:5,i]), as.data.frame(DAT_arr_AMR[,6,i])) 
}


Res10_amr <- as.data.frame(cbind(Res_10a[,1:5]))
colnames(Res10_amr) <- c("Anc_Prop_EUR","Anc_Prop_AFR", "Anc_Prop_SAS", "Anc_Prop_EAS", "Anc_Prop_IAM")

### Output of least square errors
LSEs_10amr_1000 <- as.data.frame(cbind(LS = Res_10a[,6], N = rep("N=10", 5000))) # extract least square values for each replicate

## Creating datasets for plots
df10 <- reshape(Res10_amr, times = c("EUR", "AFR","SAS", "EAS","IAM"), timevar = "Anc_Grp",
               varying = list(names(Res10_amr)), direction = "long") # convert wide to long. 

df10$size <- rep("N=10",5000)


##### Data Frame for AFstar Minus AFobs
Diff <- as.data.frame(cbind(EUR=AF_ref_eur_10_1000[,1]-refdat[,1],   AFR=AF_ref_afr_10_1000[,1]-refdat[,2],  SAS=AF_ref_sas_10_1000[,1]-refdat[,3], 
                            EAS=AF_ref_eas_10_1000[,1] - refdat[,4],  IAM=AF_ref_iam_10_1000[,1]- refdat[,5]))

DF2 <- reshape(Diff, times = c("EUR","AFR","SAS","EAS","IAM"), timevar = "Anc_Grp",
               varying = list(names(Diff)), direction = "long") # convert wide to long. 
DF2$size <- rep("N=10", 5000)





######################################################################
# Generate allele counts for population using binom for N = 15
######################################################################
N = 15 # number of people
rep = 1000 # number of replicates 

# AF_ref_eur_15_1000 <- t(sapply(refdat[,1], function(x){rbinom(rep, 2 * N, x)/ (2 * N)})); write.table(AF_ref_eur_15_1000,file = "AF_ref_eur_15_1000.txt", sep = ' ')
# AF_ref_afr_15_1000 <- t(sapply(refdat[,2], function(x){rbinom(rep, 2 * N, x)/ (2 * N)})); write.table(AF_ref_afr_15_1000,file = "AF_ref_afr_15_1000.txt", sep = ' ') 
# AF_ref_sas_15_1000 <- t(sapply(refdat[,3], function(x){rbinom(rep, 2 * N, x)/ (2 * N)})); write.table(AF_ref_sas_15_1000,file = "AF_ref_sas_15_1000.txt", sep = ' ')
# AF_ref_eas_15_1000 <- t(sapply(refdat[,4], function(x){rbinom(rep, 2 * N, x)/ (2 * N)})); write.table(AF_ref_eas_15_1000,file = "AF_ref_eas_15_1000.txt", sep = ' ')
# AF_ref_iam_15_1000 <- t(sapply(refdat[,5], function(x){rbinom(rep, 2 * N, x)/ (2 * N)})); write.table(AF_ref_iam_15_1000,file = "AF_ref_iam_15_1000.txt", sep = ' ')

AF_ref_eur_15_1000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_eur_15_1000.txt", sep=" ")) 
AF_ref_afr_15_1000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_afr_15_1000.txt", sep=" "))
AF_ref_sas_15_1000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_sas_15_1000.txt", sep=" "))
AF_ref_eas_15_1000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_eas_15_1000.txt", sep=" "))
AF_ref_iam_15_1000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_iam_15_1000.txt", sep=" "))


## Creating matrices to pass to the Summix function
Res_15a <- matrix(data = NA, ncol = 8, nrow = rep)

DAT_arr_AMR <- array(data=NA,dim = c(samp,6,rep))

for (i in 1:rep){
  
  DAT_arr_AMR[ , ,i] = cbind(AF_ref_eur_15_1000[,i], AF_ref_afr_15_1000[,i], AF_ref_sas_15_1000[,i], AF_ref_eas_15_1000[,i],  AF_ref_iam_15_1000[,i], refdat[,6])
  
  Res_15a[i,] = ancestr(as.data.frame(DAT_arr_AMR[,1:5,i]), as.data.frame(DAT_arr_AMR[,6,i])) 
}


Res15_amr <- as.data.frame(cbind(Res_15a[,1:5]))
colnames(Res15_amr) <- c("Anc_Prop_EUR","Anc_Prop_AFR", "Anc_Prop_SAS", "Anc_Prop_EAS", "Anc_Prop_IAM")

### Output of least square errors
LSEs_15amr_1000 <- as.data.frame(cbind(LS = Res_15a[,6], N = rep("N=15", 5000))) # extract least square values for each replicate


## Creating datasets for plots
df15 <- reshape(Res15_amr, times = c("EUR", "AFR","SAS", "EAS","IAM"), timevar = "Anc_Grp",
                varying = list(names(Res15_amr)), direction = "long") # convert wide to long. 

df15$size <- rep("N=15",5000)


##### Data Frame for AFstar Minus AFobs
Diff <- as.data.frame(cbind(EUR=AF_ref_eur_15_1000[,1]-refdat[,1],   AFR=AF_ref_afr_15_1000[,1]- refdat[,2],  SAS=AF_ref_sas_15_1000[,1]-refdat[,3], 
                            EAS=AF_ref_eas_15_1000[,1] - refdat[,4],  IAM=AF_ref_iam_15_1000[,1]- refdat[,5]))

DF3 <- reshape(Diff, times = c("EUR","AFR","SAS","EAS","IAM"), timevar = "Anc_Grp",
               varying = list(names(Diff)), direction = "long") # convert wide to long. 
DF3$size <- rep("N=15", 5000)



##########################################################################
# Generate allele counts for population using binom for N = 20
##########################################################################
N = 20 # number of people
rep = 1000 # number of replicates 

# AF_ref_eur_20_1000 <- t(sapply(refdat[,1], function(x){rbinom(rep, 2 * N, x)/ (2 * N)})); write.table(AF_ref_eur_20_1000,file = "AF_ref_eur_20_1000.txt", sep = ' ')
# AF_ref_afr_20_1000 <- t(sapply(refdat[,2], function(x){rbinom(rep, 2 * N, x)/ (2 * N)})); write.table(AF_ref_afr_20_1000,file = "AF_ref_afr_20_1000.txt", sep = ' ') 
# AF_ref_sas_20_1000 <- t(sapply(refdat[,3], function(x){rbinom(rep, 2 * N, x)/ (2 * N)})); write.table(AF_ref_sas_20_1000,file = "AF_ref_sas_20_1000.txt", sep = ' ')
# AF_ref_eas_20_1000 <- t(sapply(refdat[,4], function(x){rbinom(rep, 2 * N, x)/ (2 * N)})); write.table(AF_ref_eas_20_1000,file = "AF_ref_eas_20_1000.txt", sep = ' ')
# AF_ref_iam_20_1000 <- t(sapply(refdat[,5], function(x){rbinom(rep, 2 * N, x)/ (2 * N)})); write.table(AF_ref_iam_20_1000,file = "AF_ref_iam_20_1000.txt", sep = ' ')


AF_ref_eur_20_1000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_eur_20_1000.txt", sep=" ")) 
AF_ref_afr_20_1000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_afr_20_1000.txt", sep=" "))
AF_ref_sas_20_1000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_sas_20_1000.txt", sep=" "))
AF_ref_eas_20_1000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_eas_20_1000.txt", sep=" "))
AF_ref_iam_20_1000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_iam_20_1000.txt", sep=" "))


## Creating matrices to pass to the Summix function
Res_20a <- matrix(data = NA, ncol = 8, nrow = rep)

DAT_arr_AMR <- array(data=NA,dim = c(samp,6,rep))

##
for (i in 1:rep){
  
  DAT_arr_AMR[ , ,i] = cbind(AF_ref_eur_20_1000[,i], AF_ref_afr_20_1000[,i], AF_ref_sas_20_1000[,i], AF_ref_eas_20_1000[,i],  AF_ref_iam_20_1000[,i], refdat[,6])
  
  Res_20a[i,] = ancestr(as.data.frame(DAT_arr_AMR[,1:5,i]), as.data.frame(DAT_arr_AMR[,6,i])) 
}


Res20_amr <- as.data.frame(cbind(Res_20a[,1:5]))
colnames(Res20_amr) <- c("Anc_Prop_EUR","Anc_Prop_AFR", "Anc_Prop_SAS", "Anc_Prop_EAS", "Anc_Prop_IAM")

### Output of least square errors
LSEs_20amr_1000 <- as.data.frame(cbind(LS = Res_20a[,6], N = rep("N=20", 5000))) # extract least square values for each replicate

## Creating datasets for plots
df20 <- reshape(Res20_amr, times = c("EUR", "AFR","SAS", "EAS","IAM"), timevar = "Anc_Grp",
                varying = list(names(Res20_amr)), direction = "long") # convert wide to long. 

df20$size <- rep("N=20",5000)


##### Data Frame for AFstar Minus AFobs
Diff <- as.data.frame(cbind(EUR=AF_ref_eur_20_1000[,1]-refdat[,1],   AFR=AF_ref_afr_20_1000[,1]-refdat[,2],  SAS=AF_ref_sas_20_1000[,1]-refdat[,3], 
                            EAS=AF_ref_eas_20_1000[,1] - refdat[,4],  IAM=AF_ref_iam_20_1000[,1]- refdat[,5]))

DF4 <- reshape(Diff, times = c("EUR","AFR","SAS","EAS","IAM"), timevar = "Anc_Grp",
               varying = list(names(Diff)), direction = "long") # convert wide to long. 
DF4$size <- rep("N=20", 5000)




##########################################################################
# Generate allele counts for population using binom for N = 25
##########################################################################
N = 25 # number of people
rep = 1000 # number of replicates 

# AF_ref_eur_25_1000 <- t(sapply(refdat[,1], function(x){rbinom(rep, 2 * N, x)/ (2 * N)})); write.table(AF_ref_eur_25_1000,file = "AF_ref_eur_25_1000.txt", sep = ' ')
# AF_ref_afr_25_1000 <- t(sapply(refdat[,2], function(x){rbinom(rep, 2 * N, x)/ (2 * N)})); write.table(AF_ref_afr_25_1000,file = "AF_ref_afr_25_1000.txt", sep = ' ') 
# AF_ref_sas_25_1000 <- t(sapply(refdat[,3], function(x){rbinom(rep, 2 * N, x)/ (2 * N)})); write.table(AF_ref_sas_25_1000,file = "AF_ref_sas_25_1000.txt", sep = ' ')
# AF_ref_eas_25_1000 <- t(sapply(refdat[,4], function(x){rbinom(rep, 2 * N, x)/ (2 * N)})); write.table(AF_ref_eas_25_1000,file = "AF_ref_eas_25_1000.txt", sep = ' ')
# AF_ref_iam_25_1000 <- t(sapply(refdat[,5], function(x){rbinom(rep, 2 * N, x)/ (2 * N)})); write.table(AF_ref_iam_25_1000,file = "AF_ref_iam_25_1000.txt", sep = ' ')

AF_ref_eur_25_1000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_eur_25_1000.txt", sep=" ")) 
AF_ref_afr_25_1000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_afr_25_1000.txt", sep=" "))
AF_ref_sas_25_1000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_sas_25_1000.txt", sep=" "))
AF_ref_eas_25_1000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_eas_25_1000.txt", sep=" "))
AF_ref_iam_25_1000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_iam_25_1000.txt", sep=" "))


## Creating matrices to pass to the Summix function
Res_25a <- matrix(data = NA, ncol = 8, nrow = rep)

DAT_arr_AMR <- array(data=NA,dim = c(samp,6,rep))

for (i in 1:rep){
  
  DAT_arr_AMR[ , ,i] = cbind(AF_ref_eur_25_1000[,i], AF_ref_afr_25_1000[,i], AF_ref_sas_25_1000[,i], AF_ref_eas_25_1000[,i],  AF_ref_iam_25_1000[,i], refdat[,6])
  
  Res_25a[i,] = ancestr(as.data.frame(DAT_arr_AMR[,1:5,i]), as.data.frame(DAT_arr_AMR[,6,i])) 
}


Res25_amr <- as.data.frame(cbind(Res_25a[,1:5]))
colnames(Res25_amr) <- c("Anc_Prop_EUR","Anc_Prop_AFR", "Anc_Prop_SAS", "Anc_Prop_EAS", "Anc_Prop_IAM")

### Output of least square errors

LSEs_25amr_1000 <- as.data.frame(cbind(LS = Res_25a[,6], N = rep("N=25", 5000))) # extract least square values for each replicate

## Creating datasets for plots
df25 <- reshape(Res25_amr, times = c("EUR", "AFR","SAS", "EAS","IAM"), timevar = "Anc_Grp",
                varying = list(names(Res25_amr)), direction = "long") # convert wide to long. 

df25$size <- rep("N=25",5000)


##### Data Frame for AFstar Minus AFobs
Diff <- as.data.frame(cbind(EUR=AF_ref_eur_25_1000[,1]-refdat[,1],   AFR=AF_ref_afr_25_1000[,1]-refdat[,2],  SAS=AF_ref_sas_25_1000[,1]-refdat[,3], 
                            EAS=AF_ref_eas_25_1000[,1] - refdat[,4],  IAM=AF_ref_iam_25_1000[,1]- refdat[,5]))

DF5 <- reshape(Diff, times = c("EUR","AFR","SAS","EAS","IAM"), timevar = "Anc_Grp",
               varying = list(names(Diff)), direction = "long") # convert wide to long. 
DF5$size <- rep("N=25", 5000)



##########################################################################
# Generate allele counts for population using binom for N = 50
##########################################################################
N = 50 # number of people
rep = 1000 # number of replicates 

# AF_ref_eur_50_1000 <- t(sapply(refdat[,1], function(x){rbinom(rep, 2 * N, x)/ (2 * N)})); write.table(AF_ref_eur_50_1000,file = "AF_ref_eur_50_1000.txt", sep = ' ')
# AF_ref_afr_50_1000 <- t(sapply(refdat[,2], function(x){rbinom(rep, 2 * N, x)/ (2 * N)})); write.table(AF_ref_afr_50_1000,file = "AF_ref_afr_50_1000.txt", sep = ' ') 
# AF_ref_sas_50_1000 <- t(sapply(refdat[,3], function(x){rbinom(rep, 2 * N, x)/ (2 * N)})); write.table(AF_ref_sas_50_1000,file = "AF_ref_sas_50_1000.txt", sep = ' ')
# AF_ref_eas_50_1000 <- t(sapply(refdat[,4], function(x){rbinom(rep, 2 * N, x)/ (2 * N)})); write.table(AF_ref_eas_50_1000,file = "AF_ref_eas_50_1000.txt", sep = ' ')
# AF_ref_iam_50_1000 <- t(sapply(refdat[,5], function(x){rbinom(rep, 2 * N, x)/ (2 * N)})); write.table(AF_ref_iam_50_1000,file = "AF_ref_iam_50_1000.txt", sep = ' ')


AF_ref_eur_50_1000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_eur_50_1000.txt", sep=" ")) 
AF_ref_afr_50_1000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_afr_50_1000.txt", sep=" "))
AF_ref_sas_50_1000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_sas_50_1000.txt", sep=" "))
AF_ref_eas_50_1000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_eas_50_1000.txt", sep=" "))
AF_ref_iam_50_1000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_iam_50_1000.txt", sep=" "))


## Creating matrices to pass to the Summix function
Res_50a <- matrix(data = NA, ncol = 8, nrow = rep)

DAT_arr_AMR <- array(data=NA,dim = c(samp,6,rep))

for (i in 1:rep){
  
  DAT_arr_AMR[ , ,i] = cbind(AF_ref_eur_50_1000[,i], AF_ref_afr_50_1000[,i], AF_ref_sas_50_1000[,i], AF_ref_eas_50_1000[,i],  AF_ref_iam_50_1000[,i], refdat[,6])
  
  Res_50a[i,] = ancestr(as.data.frame(DAT_arr_AMR[,1:5,i]), as.data.frame(DAT_arr_AMR[,6,i])) 
}


Res50_amr <- as.data.frame(cbind(Res_50a[,1:5]))
colnames(Res50_amr) <- c("Anc_Prop_EUR","Anc_Prop_AFR", "Anc_Prop_SAS", "Anc_Prop_EAS", "Anc_Prop_IAM")

### Output of least square errors
LSEs_50amr_1000 <- as.data.frame(cbind(LS = Res_50a[,6], N = rep("N=50", 5000))) # extract least square values for each replicate

## Creating datasets for plots
df50 <- reshape(Res50_amr, times = c("EUR", "AFR","SAS", "EAS","IAM"), timevar = "Anc_Grp",
                varying = list(names(Res50_amr)), direction = "long") # convert wide to long. 

df50$size <- rep("N=50",5000)


##### Data Frame for AFstar Minus AFobs

Diff <- as.data.frame(cbind(EUR=AF_ref_eur_50_1000[,1]-refdat[,1],   AFR=AF_ref_afr_50_1000[,1]-refdat[,2],  SAS=AF_ref_sas_50_1000[,1]-refdat[,3], 
                            EAS=AF_ref_eas_50_1000[,1] - refdat[,4],  IAM=AF_ref_iam_50_1000[,1]- refdat[,5]))

DF6 <- reshape(Diff, times = c("EUR","AFR","SAS","EAS","IAM"), timevar = "Anc_Grp",
               varying = list(names(Diff)), direction = "long") # convert wide to long. 
DF6$size <- rep("N=50", 5000)



##########################################################################
# Generate allele counts for population using binom for N = 100
##########################################################################
N = 100 # number of people
rep = 1000 # number of replicates 

# AF_ref_eur_100_1000 <- t(sapply(refdat[,1], function(x){rbinom(rep, 2 * N, x)/ (2 * N)})); write.table(AF_ref_eur_100_1000,file = "AF_ref_eur_100_1000.txt", sep = ' ')
# AF_ref_afr_100_1000 <- t(sapply(refdat[,2], function(x){rbinom(rep, 2 * N, x)/ (2 * N)})); write.table(AF_ref_afr_100_1000,file = "AF_ref_afr_100_1000.txt", sep = ' ') 
# AF_ref_sas_100_1000 <- t(sapply(refdat[,3], function(x){rbinom(rep, 2 * N, x)/ (2 * N)})); write.table(AF_ref_sas_100_1000,file = "AF_ref_sas_100_1000.txt", sep = ' ')
# AF_ref_eas_100_1000 <- t(sapply(refdat[,4], function(x){rbinom(rep, 2 * N, x)/ (2 * N)})); write.table(AF_ref_eas_100_1000,file = "AF_ref_eas_100_1000.txt", sep = ' ')
# AF_ref_iam_100_1000 <- t(sapply(refdat[,5], function(x){rbinom(rep, 2 * N, x)/ (2 * N)})); write.table(AF_ref_iam_100_1000,file = "AF_ref_iam_100_1000.txt", sep = ' ')


AF_ref_eur_100_1000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_eur_100_1000.txt", sep=" ")) 
AF_ref_afr_100_1000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_afr_100_1000.txt", sep=" "))
AF_ref_sas_100_1000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_sas_100_1000.txt", sep=" "))
AF_ref_eas_100_1000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_eas_100_1000.txt", sep=" "))
AF_ref_iam_100_1000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_iam_100_1000.txt", sep=" "))


## Creating matrices to pass to the Summix function
Res_100a <- matrix(data = NA, ncol = 8, nrow = rep)

DAT_arr_AMR <- array(data=NA,dim = c(samp,6,rep))

##
for (i in 1:rep){
  
  DAT_arr_AMR[ , ,i] = cbind(AF_ref_eur_100_1000[,i], AF_ref_afr_100_1000[,i], AF_ref_sas_100_1000[,i], AF_ref_eas_100_1000[,i],  AF_ref_iam_100_1000[,i], refdat[,6])
  
  Res_100a[i,] = ancestr(as.data.frame(DAT_arr_AMR[,1:5,i]), as.data.frame(DAT_arr_AMR[,6,i])) 
}


Res100_amr <- as.data.frame(cbind(Res_100a[,1:5]))
colnames(Res100_amr) <- c("Anc_Prop_EUR","Anc_Prop_AFR", "Anc_Prop_SAS", "Anc_Prop_EAS", "Anc_Prop_IAM")

### Output of least square errors
LSEs_100amr_1000 <- as.data.frame(cbind(LS = Res_100a[,6], N = rep("N=100", 5000))) # extract least square values for each replicate

## Creating datasets for plots
df100 <- reshape(Res100_amr, times = c("EUR", "AFR","SAS", "EAS","IAM"), timevar = "Anc_Grp",
                 varying = list(names(Res100_amr)), direction = "long") # convert wide to long. 

df100$size <- rep("N=100",5000)


##### Data Frame for AFstar Minus AFobs
Diff <- as.data.frame(cbind(EUR=AF_ref_eur_100_1000[,1]-refdat[,1],   AFR=AF_ref_afr_100_1000[,1]-refdat[,2],  SAS=AF_ref_sas_100_1000[,1]-refdat[,3], 
                            EAS=AF_ref_eas_100_1000[,1]-refdat[,4],  IAM=AF_ref_iam_100_1000[,1]-refdat[,5]))

DF7 <- reshape(Diff, times = c("EUR","AFR","SAS","EAS","IAM"), timevar = "Anc_Grp",
               varying = list(names(Diff)), direction = "long") # convert wide to long. 
DF7$size <- rep("N=100", 5000)




##########################################################################
# Generate allele counts for population using binom for N = 500
##########################################################################
N = 500 # number of people
rep = 1000 # number of replicates 

# AF_ref_eur_500_1000 <- t(sapply(refdat[,1], function(x){rbinom(rep, 2 * N, x)/ (2 * N)})); write.table(AF_ref_eur_500_1000,file = "AF_ref_eur_500_1000.txt", sep = ' ')
# AF_ref_afr_500_1000 <- t(sapply(refdat[,2], function(x){rbinom(rep, 2 * N, x)/ (2 * N)})); write.table(AF_ref_afr_500_1000,file = "AF_ref_afr_500_1000.txt", sep = ' ') 
# AF_ref_sas_500_1000 <- t(sapply(refdat[,3], function(x){rbinom(rep, 2 * N, x)/ (2 * N)})); write.table(AF_ref_sas_500_1000,file = "AF_ref_sas_500_1000.txt", sep = ' ')
# AF_ref_eas_500_1000 <- t(sapply(refdat[,4], function(x){rbinom(rep, 2 * N, x)/ (2 * N)})); write.table(AF_ref_eas_500_1000,file = "AF_ref_eas_500_1000.txt", sep = ' ')
# AF_ref_iam_500_1000 <- t(sapply(refdat[,5], function(x){rbinom(rep, 2 * N, x)/ (2 * N)})); write.table(AF_ref_iam_500_1000,file = "AF_ref_iam_500_1000.txt", sep = ' ')


AF_ref_eur_500_1000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_eur_500_1000.txt", sep=" ")) 
AF_ref_afr_500_1000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_afr_500_1000.txt", sep=" "))
AF_ref_sas_500_1000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_sas_500_1000.txt", sep=" "))
AF_ref_eas_500_1000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_eas_500_1000.txt", sep=" "))
AF_ref_iam_500_1000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_iam_500_1000.txt", sep=" "))


## Creating matrices to pass to the Summix function
Res_500a <- matrix(data = NA, ncol = 8, nrow = rep)

DAT_arr_AMR <- array(data=NA,dim = c(samp,6,rep))

##
for (i in 1:rep){
  
  DAT_arr_AMR[ , ,i] = cbind(AF_ref_eur_500_1000[,i], AF_ref_afr_500_1000[,i], AF_ref_sas_500_1000[,i], AF_ref_eas_500_1000[,i],  AF_ref_iam_500_1000[,i], refdat[,6])
  
  Res_500a[i,] = ancestr(as.data.frame(DAT_arr_AMR[,1:5,i]), as.data.frame(DAT_arr_AMR[,6,i])) 
}


Res500_amr <- as.data.frame(cbind(Res_500a[,1:5]))
colnames(Res500_amr) <- c("Anc_Prop_EUR","Anc_Prop_AFR", "Anc_Prop_SAS", "Anc_Prop_EAS", "Anc_Prop_IAM")

### Output of least square errors
LSEs_500amr_1000 <- as.data.frame(cbind(LS = Res_500a[,6], N = rep("N=500", 5000))) # extract least square values for each replicate

## Creating datasets for plots
df500 <- reshape(Res500_amr, times = c("EUR", "AFR","SAS", "EAS","IAM"), timevar = "Anc_Grp",
                 varying = list(names(Res500_amr)), direction = "long") # convert wide to long. 

df500$size <- rep("N=500",5000)


##### Data Frame for AFstar Minus AFobs

Diff <- as.data.frame(cbind(EUR=AF_ref_eur_500_1000[,1]-refdat[,1],   AFR=AF_ref_afr_500_1000[,1]-refdat[,2],  SAS=AF_ref_sas_500_1000[,1]-refdat[,3], 
                            EAS=AF_ref_eas_500_1000[,1] - refdat[,4],  IAM=AF_ref_iam_500_1000[,1]- refdat[,5]))

DF8 <- reshape(Diff, times = c("EUR","AFR","SAS","EAS","IAM"), timevar = "Anc_Grp",
               varying = list(names(Diff)), direction = "long") # convert wide to long. 
DF8$size <- rep("N=500", 5000)



##########################################################################
# Generate allele counts for population using binom for N = 1000
##########################################################################
N = 1000 # number of people
rep = 1000 # number of replicates 

# AF_ref_eur_1000_1000 <- t(sapply(refdat[,1], function(x){rbinom(rep, 2 * N, x)/ (2 * N)})); write.table(AF_ref_eur_1000_1000,file = "AF_ref_eur_1000_1000.txt", sep = ' ')
# AF_ref_afr_1000_1000 <- t(sapply(refdat[,2], function(x){rbinom(rep, 2 * N, x)/ (2 * N)})); write.table(AF_ref_afr_1000_1000,file = "AF_ref_afr_1000_1000.txt", sep = ' ') 
# AF_ref_sas_1000_1000 <- t(sapply(refdat[,3], function(x){rbinom(rep, 2 * N, x)/ (2 * N)})); write.table(AF_ref_sas_1000_1000,file = "AF_ref_sas_1000_1000.txt", sep = ' ')
# AF_ref_eas_1000_1000 <- t(sapply(refdat[,4], function(x){rbinom(rep, 2 * N, x)/ (2 * N)})); write.table(AF_ref_eas_1000_1000,file = "AF_ref_eas_1000_1000.txt", sep = ' ')
# AF_ref_iam_1000_1000 <- t(sapply(refdat[,5], function(x){rbinom(rep, 2 * N, x)/ (2 * N)})); write.table(AF_ref_iam_1000_1000,file = "AF_ref_iam_1000_1000.txt", sep = ' ')

AF_ref_eur_1000_1000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_eur_1000_1000.txt", sep=" ")) 
AF_ref_afr_1000_1000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_afr_1000_1000.txt", sep=" "))
AF_ref_sas_1000_1000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_sas_1000_1000.txt", sep=" "))
AF_ref_eas_1000_1000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_eas_1000_1000.txt", sep=" "))
AF_ref_iam_1000_1000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_iam_1000_1000.txt", sep=" "))


## Creating matrices to pass to the Summix function
Res_1000a <- matrix(data = NA, ncol = 8, nrow = rep)

DAT_arr_AMR <- array(data=NA,dim = c(samp,6,rep))

##
for (i in 1:rep){
  
  DAT_arr_AMR[ , ,i] = cbind(AF_ref_eur_1000_1000[,i], AF_ref_afr_1000_1000[,i], AF_ref_sas_1000_1000[,i], AF_ref_eas_1000_1000[,i],  AF_ref_iam_1000_1000[,i], refdat[,6])
  
  Res_1000a[i,] = ancestr(as.data.frame(DAT_arr_AMR[,1:5,i]), as.data.frame(DAT_arr_AMR[,6,i])) 
}


Res1000_amr <- as.data.frame(cbind(Res_1000a[,1:5]))
colnames(Res1000_amr) <- c("Anc_Prop_EUR","Anc_Prop_AFR", "Anc_Prop_SAS", "Anc_Prop_EAS", "Anc_Prop_IAM")

### Output of least square errors
LSEs_1000amr_1000 <- as.data.frame(cbind(LS = Res_1000a[,6], N = rep("N=1000", 5000))) # extract least square values for each replicate

## Creating datasets for plots
df1000 <- reshape(Res1000_amr, times = c("EUR", "AFR","SAS", "EAS","IAM"), timevar = "Anc_Grp",
                  varying = list(names(Res1000_amr)), direction = "long") # convert wide to long. 

df1000$size <- rep("N=1000",5000)


##### Data Frame for AFstar Minus AFobs

Diff <- as.data.frame(cbind(EUR=AF_ref_eur_1000_1000[,1]-refdat[,1],   AFR=AF_ref_afr_1000_1000[,1]-refdat[,2],  SAS=AF_ref_sas_1000_1000[,1]-refdat[,3], 
                            EAS=AF_ref_eas_1000_1000[,1] - refdat[,4],  IAM=AF_ref_iam_1000_1000[,1]- refdat[,5]))

DF9 <- reshape(Diff, times = c("EUR","AFR","SAS","EAS","IAM"), timevar = "Anc_Grp",
               varying = list(names(Diff)), direction = "long") # convert wide to long. 
DF9$size <- rep("N=1000", 5000)



##################################################
####Combined data
##################################################
# Combining all ancestry proportions estimates for each N
df.combined <- rbind(df5,df10, df15, df20, df25, df50, df100, df500, df1000)

# Combining difference in AFs for each N
DF.combine <- rbind(DF1, DF2, DF3, DF4, DF5, DF6, DF7, DF8, DF9)

# combining all least square values for each N
LS.combine <- rbind(LSEs_5amr_1000, LSEs_10amr_1000, LSEs_15amr_1000, LSEs_20amr_1000, LSEs_25amr_1000, LSEs_50amr_1000,
                    LSEs_100amr_1000, LSEs_500amr_1000, LSEs_1000amr_1000) 


pdf('Rplot_SNP1000_facA1.pdf')
ggplot(df.combined, aes(x=size, y=Anc_Prop_EUR, color=Anc_Grp)) + #x = N, y = ancestry prop, fill = ancestry groups
  geom_boxplot() + labs(x = "Allele Number (2 * N)", y = "Ancestry Proportion")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Variation of Summix Ancestry Proportion Estimate (SNP = 1,000) \n American/ Latinx GnomAD")+
  scale_x_discrete(limits=c("N=5","N=10","N=15","N=20","N=25","N=50","N=100","N=500","N=1000"))+
  facet_wrap(~Anc_Grp, scales = "free")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


pdf('Rplot_SNP1000_diff1D1.pdf')
ggplot(DF.combine, aes(x=size, y=EUR, color=Anc_Grp)) + #x = N, y = ancestry prop, fill = ancestry groups
  geom_boxplot() + labs(x = "Allele Number (2 * N)", y = "Difference in Allele Frequencies")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Distribution of Difference in Allele Frequencies (SNP = 1,000) \n American/ Latinx GnomAD")+
  scale_x_discrete(limits=c("N=5","N=10","N=15","N=20","N=25","N=50","N=100","N=500","N=1000"))+
  facet_wrap(~Anc_Grp, scales = "free")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


## Create plot for least square values for each N
pdf('Rplot_SNP1000Latinx_LS.pdf')
ggplot(LS.combine, aes(x=N, y=as.numeric(LS))) +
  geom_boxplot() + labs(x = "Allele Number (2 * N)", y = "Least Square Values")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Distribution of Least Square Values (SNP = 100,000) \n American/ Latinx GnomAD")+
  scale_x_discrete(limits=c("N=5","N=10","N=15","N=20","N=25","N=50","N=100","N=500","N=1000"))+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


df.combined <- df.combined %>% 
  mutate(N = ifelse(size=="N=5",5,
                    ifelse(size == "N=10", 10,
                           ifelse(size == "N=15", 15,
                                  ifelse(size == "N=20", 20,
                                         ifelse(size == "N=25", 25,
                                                ifelse(size == "N=50", 50, 
                                                       ifelse(size == "N=100", 100, 
                                                              ifelse(size == "N=500", 500, 1000)))))))))

dfA <- df.combined %>%
  group_by(Anc_Grp,N) %>%
  summarize(mean = mean(Anc_Prop_EUR),
            std = sd(Anc_Prop_EUR),
            minimum =min(Anc_Prop_EUR), 
            Q1=quantile(Anc_Prop_EUR, .25), 
            median = median(Anc_Prop_EUR),
            Q3=quantile(Anc_Prop_EUR, .75),
            maximum = max(Anc_Prop_EUR), .groups = 'keep')

write.table(dfA, file = "SNPs_1000A1.txt", sep = ' ')


DF.combine <- DF.combine %>% 
  mutate(N = ifelse(size=="N=5",5,
                    ifelse(size == "N=10", 10,
                           ifelse(size == "N=15", 15,
                                  ifelse(size == "N=20", 20,
                                         ifelse(size == "N=25", 25,
                                                ifelse(size == "N=50", 50, 
                                                       ifelse(size == "N=100", 100, 
                                                              ifelse(size == "N=500", 500, 1000)))))))))


dfA1 <- DF.combine %>%
  group_by(Anc_Grp,N) %>%
  summarize(mean = mean(EUR),
            std = sd(EUR),
            minimum =min(EUR), 
            Q1=quantile(EUR, .25), 
            median = median(EUR),
            Q3=quantile(EUR, .75),
            maximum = max(EUR), .groups = 'keep')

write.table(dfA1, file = "SNPs_1000D1_1.txt", sep = ' ')



# Create plot for variability in simulated and observed AFs for each 1/N and 1/N^2
pdf('Rplot_SNP10000Latinx_Vardiff_Ninv.pdf')
ggplot(dfA1, aes(x=as.numeric(1/N), y=as.numeric(std^2), color=Anc_Grp)) +
  geom_point() + labs(x = TeX("$\\frac{1}{N}$"), y = "Variance")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Variability of Difference in Allele Frequencies (SNP = 1000) \n American/ Latinx GnomAD")+
  facet_wrap(~Anc_Grp, scales = "free")
dev.off()

pdf('Rplot_SNP10000Latinx_Vardiff_Nsqinv.pdf')
ggplot(dfA1, aes(x=as.numeric(1/(N^2)), y=as.numeric(std^2), color=Anc_Grp)) +
  geom_point() + labs(x = TeX("$\\frac{1}{N^2}$"), y = "Variance")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Variability of Difference in Allele Frequencies (SNP = 1000) \n American/ Latinx GnomAD")+
  facet_wrap(~Anc_Grp, scales = "free")
dev.off()


##########################################################################
##### Selecting 10K SNPS from reference data(SNP = 10,000)
##### For N = 5, 10, 15, 20, 25, 10, 50, 100, 500, 1000
##########################################################################

# Set Seed
set.seed(122356)
samp1 = 10000

refdat1 = data_red %>% 
  sample_n(samp1) %>% 
  select(ref_AF_eur_1000G,ref_AF_afr_1000G, ref_AF_sas_1000G,ref_AF_eas_1000G,ref_AF_iam_1000G, gnomad_AF_amr) %>% 
  
  rename(ref_AFR = ref_AF_afr_1000G, ref_EUR = ref_AF_eur_1000G, ref_SAS = ref_AF_sas_1000G, ref_EAS = ref_AF_eas_1000G,
         ref_IAM = ref_AF_iam_1000G, obs_AMR = gnomad_AF_amr)


##########################################################################
# Generate allele counts for population using binom for N = 5
##########################################################################

N1 = 5 # number of people (--> allele number = 2 * N)
rep = 1000 # number of replicates 

# AF_ref_eur1_5_10000 <- t(sapply(refdat1[,1], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})); write.table(AF_ref_eur1_5_10000,file = "AF_ref_eur1_5_10000.txt", sep = ' ')
# AF_ref_afr1_5_10000 <- t(sapply(refdat1[,2], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})); write.table(AF_ref_afr1_5_10000,file = "AF_ref_afr1_5_10000.txt", sep = ' ') 
# AF_ref_sas1_5_10000 <- t(sapply(refdat1[,3], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})); write.table(AF_ref_sas1_5_10000,file = "AF_ref_sas1_5_10000.txt", sep = ' ')
# AF_ref_eas1_5_10000 <- t(sapply(refdat1[,4], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})); write.table(AF_ref_eas1_5_10000,file = "AF_ref_eas1_5_10000.txt", sep = ' ')
# AF_ref_iam1_5_10000 <- t(sapply(refdat1[,5], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})); write.table(AF_ref_iam1_5_10000,file = "AF_ref_iam1_5_10000.txt", sep = ' ')


AF_ref_eur1_5_10000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_eur1_5_10000.txt", sep=" ")) 
AF_ref_afr1_5_10000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_afr1_5_10000.txt", sep=" "))
AF_ref_sas1_5_10000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_sas1_5_10000.txt", sep=" "))
AF_ref_eas1_5_10000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_eas1_5_10000.txt", sep=" "))
AF_ref_iam1_5_10000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_iam1_5_10000.txt", sep=" "))


## Creating matrices to pass to the Summix function
Res_5a1 <- matrix(data = NA, ncol = 8, nrow = rep)

DAT_arr_AMR1 <- array(data=NA,dim = c(samp1,6,rep))


for (i in 1:rep){
  
  DAT_arr_AMR1[ , ,i] = cbind(AF_ref_eur1_5_10000[,i], AF_ref_afr1_5_10000[,i], AF_ref_sas1_5_10000[,i], AF_ref_eas1_5_10000[,i],
                              AF_ref_iam1_5_10000[,i], refdat1[,6])
  
  Res_5a1[i,] = ancestr(as.data.frame(DAT_arr_AMR1[,1:5,i]), as.data.frame(DAT_arr_AMR1[,6,i])) 
}

setwd("/newhome/agboolol/")

Res5_amr1 <- as.data.frame(cbind(Res_5a1[,1:5]))
colnames(Res5_amr1) <- c("Anc_Prop_EUR","Anc_Prop_AFR", "Anc_Prop_SAS", "Anc_Prop_EAS", "Anc_Prop_IAM")

### Output of least square errors
LSEs_5amr1_10000 <- as.data.frame(cbind(LS = Res_5a1[,6], N = rep("N=5", 50000))) # extract least square values for each replicate

## Creating datasets for plots
df5_1 <- reshape(Res5_amr1, times = c("EUR", "AFR","SAS", "EAS","IAM"), timevar = "Anc_Grp",
                 varying = list(names(Res5_amr1)), direction = "long") # convert wide to long. 


df5_1$size <- rep("N=5",5000)

##### Data Frame for AFstar Minus AFobs
Diff <- as.data.frame(cbind(EUR=AF_ref_eur1_5_10000[,1]-refdat1[,1],   AFR=AF_ref_afr1_5_10000[,1]-refdat1[,2],  SAS=AF_ref_sas1_5_10000[,1]-refdat1[,3], 
                            EAS=AF_ref_eas1_5_10000[,1] - refdat1[,4],  IAM=AF_ref_iam1_5_10000[,1]- refdat1[,5]))

DF1_1 <- reshape(Diff, times = c("EUR","AFR","SAS","EAS","IAM"), timevar = "Anc_Grp",
                 varying = list(names(Diff)), direction = "long") # convert wide to long. 

DF1_1$size <- rep("N=5", 50000)



##########################################################################
# Generate allele counts for population using binom for N = 10
##########################################################################
N1 = 10 # number of people (--> allele number = 2 * N)
rep = 1000 # number of replicates 

# AF_ref_eur1_10_10000 <- t(sapply(refdat1[,1], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})); write.table(AF_ref_eur1_10_10000,file = "AF_ref_eur1_10_10000.txt", sep = ' ')
# AF_ref_afr1_10_10000 <- t(sapply(refdat1[,2], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})); write.table(AF_ref_afr1_10_10000,file = "AF_ref_afr1_10_10000.txt", sep = ' ')
# AF_ref_sas1_10_10000 <- t(sapply(refdat1[,3], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})); write.table(AF_ref_sas1_10_10000,file = "AF_ref_sas1_10_10000.txt", sep = ' ')
# AF_ref_eas1_10_10000 <- t(sapply(refdat1[,4], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})); write.table(AF_ref_eas1_10_10000,file = "AF_ref_eas1_10_10000.txt", sep = ' ')
# AF_ref_iam1_10_10000 <- t(sapply(refdat1[,5], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})); write.table(AF_ref_iam1_10_10000,file = "AF_ref_iam1_10_10000.txt", sep = ' ')


AF_ref_eur1_10_10000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_eur1_10_10000.txt", sep=" ")) 
AF_ref_afr1_10_10000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_afr1_10_10000.txt", sep=" "))
AF_ref_sas1_10_10000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_sas1_10_10000.txt", sep=" "))
AF_ref_eas1_10_10000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_eas1_10_10000.txt", sep=" "))
AF_ref_iam1_10_10000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_iam1_10_10000.txt", sep=" "))


## Creating matrices to pass to the Summix function
Res_10a1 <- matrix(data = NA, ncol = 8, nrow = rep)

DAT_arr_AMR1 <- array(data=NA,dim = c(samp1,6,rep))

##
for (i in 1:rep){
  
  DAT_arr_AMR1[ , ,i] = cbind(AF_ref_eur1_10_10000[,i], AF_ref_afr1_10_10000[,i], AF_ref_sas1_10_10000[,i], AF_ref_eas1_10_10000[,i],  AF_ref_iam1_10_10000[,i], refdat1[,6])
  
  Res_10a1[i,] = ancestr(as.data.frame(DAT_arr_AMR1[,1:5,i]), as.data.frame(DAT_arr_AMR1[,6,i])) 
}


Res10_amr1 <- as.data.frame(cbind(Res_10a1[,1:5]))
colnames(Res10_amr1) <- c("Anc_Prop_EUR","Anc_Prop_AFR", "Anc_Prop_SAS", "Anc_Prop_EAS", "Anc_Prop_IAM")

### Output of least square errors
LSEs_10amr1_10000 <- as.data.frame(cbind(LS = Res_10a1[,6], N = rep("N=10", 50000))) # extract least square values for each replicate

## Creating datasets for plots
df10_1 <- reshape(Res10_amr1, times = c("EUR", "AFR","SAS", "EAS","IAM"), timevar = "Anc_Grp",
                  varying = list(names(Res10_amr1)), direction = "long") # convert wide to long. 


df10_1$size <- rep("N=10",50000)


##### Data Frame for AFstar Minus AFobs
Diff <- as.data.frame(cbind(EUR=AF_ref_eur1_10_10000[,1]-refdat1[,1],   AFR=AF_ref_afr1_10_10000[,1]-refdat1[,2],  SAS=AF_ref_sas1_10_10000[,1]-refdat1[,3], 
                            EAS=AF_ref_eas1_10_10000[,1] - refdat1[,4],  IAM=AF_ref_iam1_10_10000[,1]- refdat1[,5]))

DF2_1 <- reshape(Diff, times = c("EUR","AFR","SAS","EAS","IAM"), timevar = "Anc_Grp",
                 varying = list(names(Diff)), direction = "long") # convert wide to long. 

DF2_1$size <- rep("N=10", 50000)




##########################################################################
# Generate allele counts for population using binom for N = 15
##########################################################################
N1 = 15 # number of people (--> allele number = 2 * N)
rep = 1000 # number of replicates 

# AF_ref_eur1_15_10000 <- t(sapply(refdat1[,1], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})); write.table(AF_ref_eur1_15_10000,file = "AF_ref_eur1_15_10000.txt", sep = ' ')
# AF_ref_afr1_15_10000 <- t(sapply(refdat1[,2], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})); write.table(AF_ref_afr1_15_10000,file = "AF_ref_afr1_15_10000.txt", sep = ' ') 
# AF_ref_sas1_15_10000 <- t(sapply(refdat1[,3], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})); write.table(AF_ref_sas1_15_10000,file = "AF_ref_sas1_15_10000.txt", sep = ' ')
# AF_ref_eas1_15_10000 <- t(sapply(refdat1[,4], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})); write.table(AF_ref_eas1_15_10000,file = "AF_ref_eas1_15_10000.txt", sep = ' ')
# AF_ref_iam1_15_10000 <- t(sapply(refdat1[,5], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})); write.table(AF_ref_iam1_15_10000,file = "AF_ref_iam1_15_10000.txt", sep = ' ')

AF_ref_eur1_15_10000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_eur1_15_10000.txt", sep=" ")) 
AF_ref_afr1_15_10000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_afr1_15_10000.txt", sep=" "))
AF_ref_sas1_15_10000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_sas1_15_10000.txt", sep=" "))
AF_ref_eas1_15_10000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_eas1_15_10000.txt", sep=" "))
AF_ref_iam1_15_10000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_iam1_15_10000.txt", sep=" "))


## Creating matrices to pass to the Summix function
Res_15a1 <- matrix(data = NA, ncol = 8, nrow = rep)

DAT_arr_AMR1 <- array(data=NA,dim = c(samp1,6,rep))

##
for (i in 1:rep){
  
  DAT_arr_AMR1[ , ,i] = cbind(AF_ref_eur1_15_10000[,i], AF_ref_afr1_15_10000[,i], AF_ref_sas1_15_10000[,i], AF_ref_eas1_15_10000[,i],  AF_ref_iam1_15_10000[,i], refdat1[,6])
  
  Res_15a1[i,] = ancestr(as.data.frame(DAT_arr_AMR1[,1:5,i]), as.data.frame(DAT_arr_AMR1[,6,i])) 
}


Res15_amr1 <- as.data.frame(cbind(Res_15a1[,1:5]))
colnames(Res15_amr1) <- c("Anc_Prop_EUR","Anc_Prop_AFR", "Anc_Prop_SAS", "Anc_Prop_EAS", "Anc_Prop_IAM")

### Output of least square errors
LSEs_15amr1_10000 <- as.data.frame(cbind(LS = Res_15a1[,6], N = rep("N=15", 50000))) # extract least square values for each replicate

## Creating datasets for plots
df15_1 <- reshape(Res15_amr1, times = c("EUR", "AFR","SAS", "EAS","IAM"), timevar = "Anc_Grp",
                  varying = list(names(Res15_amr1)), direction = "long") # convert wide to long. 


df15_1$size <- rep("N=15",5000)


##### Data Frame for AFstar Minus AFobs
Diff <- as.data.frame(cbind(EUR=AF_ref_eur1_15_10000[,1]-refdat1[,1],   AFR=AF_ref_afr1_15_10000[,1]-refdat1[,2],  SAS=AF_ref_sas1_15_10000[,1]-refdat1[,3], 
                            EAS=AF_ref_eas1_15_10000[,1] - refdat1[,4],  IAM=AF_ref_iam1_15_10000[,1]- refdat1[,5]))



DF3_1 <- reshape(Diff, times = c("EUR","AFR","SAS","EAS","IAM"), timevar = "Anc_Grp",
                 varying = list(names(Diff)), direction = "long") # convert wide to long. 

DF3_1$size <- rep("N=15", 50000)



##########################################################################
# Generate allele counts for population using binom for N = 20
##########################################################################
N1 = 20 # number of people (--> allele number = 2 * N)
rep = 1000 # number of replicates 

# AF_ref_eur1_20_10000 <- t(sapply(refdat1[,1], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})); write.table(AF_ref_eur1_20_10000,file = "AF_ref_eur1_20_10000.txt", sep = ' ')
# AF_ref_afr1_20_10000 <- t(sapply(refdat1[,2], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})); write.table(AF_ref_afr1_20_10000,file = "AF_ref_afr1_20_10000.txt", sep = ' ') 
# AF_ref_sas1_20_10000 <- t(sapply(refdat1[,3], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})); write.table(AF_ref_sas1_20_10000,file = "AF_ref_sas1_20_10000.txt", sep = ' ')
# AF_ref_eas1_20_10000 <- t(sapply(refdat1[,4], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})); write.table(AF_ref_eas1_20_10000,file = "AF_ref_eas1_20_10000.txt", sep = ' ')
# AF_ref_iam1_20_10000 <- t(sapply(refdat1[,5], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})); write.table(AF_ref_iam1_20_10000,file = "AF_ref_iam1_20_10000.txt", sep = ' ')

AF_ref_eur1_20_10000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_eur1_20_10000.txt", sep=" ")) 
AF_ref_afr1_20_10000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_afr1_20_10000.txt", sep=" "))
AF_ref_sas1_20_10000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_sas1_20_10000.txt", sep=" "))
AF_ref_eas1_20_10000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_eas1_20_10000.txt", sep=" "))
AF_ref_iam1_20_10000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_iam1_20_10000.txt", sep=" "))

## Creating matrices to pass to the Summix function
Res_20a1 <- matrix(data = NA, ncol = 8, nrow = rep)

DAT_arr_AMR1 <- array(data=NA,dim = c(samp1,6,rep))

##
for (i in 1:rep){
  
  DAT_arr_AMR1[ , ,i] = cbind(AF_ref_eur1_20_10000[,i], AF_ref_afr1_20_10000[,i], AF_ref_sas1_20_10000[,i], AF_ref_eas1_20_10000[,i],  AF_ref_iam1_20_10000[,i], refdat1[,6])
  
  Res_20a1[i,] = ancestr(as.data.frame(DAT_arr_AMR1[,1:5,i]), as.data.frame(DAT_arr_AMR1[,6,i])) 
}


Res20_amr1 <- as.data.frame(cbind(Res_20a1[,1:5]))
colnames(Res20_amr1) <- c("Anc_Prop_EUR","Anc_Prop_AFR", "Anc_Prop_SAS", "Anc_Prop_EAS", "Anc_Prop_IAM")

### Output of least square errors
LSEs_20amr1_10000 <- as.data.frame(cbind(LS = Res_20a1[,6], N = rep("N=20", 50000))) # extract least square values for each replicate

## Creating datasets for plots
df20_1 <- reshape(Res20_amr1, times = c("EUR", "AFR","SAS", "EAS","IAM"), timevar = "Anc_Grp",
                  varying = list(names(Res20_amr1)), direction = "long") # convert wide to long. 


df20_1$size <- rep("N=20",5000)


##### Data Frame for AFstar Minus AFobs
Diff <- as.data.frame(cbind(EUR=AF_ref_eur1_20_10000[,1]-refdat1[,1],   AFR=AF_ref_afr1_20_10000[,1]-refdat1[,2],  SAS=AF_ref_sas1_20_10000[,1]-refdat1[,3], 
                            EAS=AF_ref_eas1_20_10000[,1] - refdat1[,4],  IAM=AF_ref_iam1_20_10000[,1]- refdat1[,5]))

DF4_1 <- reshape(Diff, times = c("EUR","AFR","SAS","EAS","IAM"), timevar = "Anc_Grp",
                 varying = list(names(Diff)), direction = "long") # convert wide to long. 

DF4_1$size <- rep("N=20", 50000)




##########################################################################
# Generate allele counts for population using binom for N = 25
##########################################################################
N1 = 25 # number of people (--> allele number = 2 * N)
rep = 1000 # number of replicates 

# AF_ref_eur1_25_10000 <- t(sapply(refdat1[,1], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})); write.table(AF_ref_eur1_25_10000,file = "AF_ref_eur1_25_10000.txt", sep = ' ')
# AF_ref_afr1_25_10000 <- t(sapply(refdat1[,2], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})); write.table(AF_ref_afr1_25_10000,file = "AF_ref_afr1_25_10000.txt", sep = ' ') 
# AF_ref_sas1_25_10000 <- t(sapply(refdat1[,3], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})); write.table(AF_ref_sas1_25_10000,file = "AF_ref_sas1_25_10000.txt", sep = ' ')
# AF_ref_eas1_25_10000 <- t(sapply(refdat1[,4], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})); write.table(AF_ref_eas1_25_10000,file = "AF_ref_eas1_25_10000.txt", sep = ' ')
# AF_ref_iam1_25_10000 <- t(sapply(refdat1[,5], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})); write.table(AF_ref_iam1_25_10000,file = "AF_ref_iam1_25_10000.txt", sep = ' ')


AF_ref_eur1_25_10000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_eur1_25_10000.txt", sep=" ")) 
AF_ref_afr1_25_10000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_afr1_25_10000.txt", sep=" "))
AF_ref_sas1_25_10000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_sas1_25_10000.txt", sep=" "))
AF_ref_eas1_25_10000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_eas1_25_10000.txt", sep=" "))
AF_ref_iam1_25_10000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_iam1_25_10000.txt", sep=" "))

## Creating matrices to pass to the Summix function
Res_25a1 <- matrix(data = NA, ncol = 8, nrow = rep)

DAT_arr_AMR1 <- array(data=NA,dim = c(samp1,6,rep))

##
for (i in 1:rep){
  
  DAT_arr_AMR1[ , ,i] = cbind(AF_ref_eur1_25_10000[,i], AF_ref_afr1_25_10000[,i], AF_ref_sas1_25_10000[,i], AF_ref_eas1_25_10000[,i],  AF_ref_iam1_25_10000[,i], refdat1[,6])
  
  Res_25a1[i,] = ancestr(as.data.frame(DAT_arr_AMR1[,1:5,i]), as.data.frame(DAT_arr_AMR1[,6,i])) 
}


Res25_amr1 <- as.data.frame(cbind(Res_25a1[,1:5]))
colnames(Res25_amr1) <- c("Anc_Prop_EUR","Anc_Prop_AFR", "Anc_Prop_SAS", "Anc_Prop_EAS", "Anc_Prop_IAM")

### Output of least square errors
LSEs_25amr1_10000 <- as.data.frame(cbind(LS = Res_25a1[,6], N = rep("N=25", 50000))) # extract least square values for each replicate

## Creating datasets for plots
df25_1 <- reshape(Res25_amr1, times = c("EUR", "AFR","SAS", "EAS","IAM"), timevar = "Anc_Grp",
                  varying = list(names(Res25_amr1)), direction = "long") # convert wide to long. 

df25_1$size <- rep("N=25",5000)


##### Data Frame for AFstar Minus AFobs
Diff <- as.data.frame(cbind(EUR=AF_ref_eur1_25_10000[,1]-refdat1[,1],   AFR=AF_ref_afr1_25_10000[,1]-refdat1[,2],  SAS=AF_ref_sas1_25_10000[,1]-refdat1[,3], 
                            EAS=AF_ref_eas1_25_10000[,1] - refdat1[,4],  IAM=AF_ref_iam1_25_10000[,1]- refdat1[,5]))

DF5_1 <- reshape(Diff, times = c("EUR","AFR","SAS","EAS","IAM"), timevar = "Anc_Grp",
                 varying = list(names(Diff)), direction = "long") # convert wide to long. 

DF5_1$size <- rep("N=25", 50000)




##########################################################################
# Generate allele counts for population using binom for N = 50
##########################################################################
N1 = 50 # number of people (--> allele number = 2 * N)
rep = 1000 # number of replicates 

# AF_ref_eur1_50_10000 <- t(sapply(refdat1[,1], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})); write.table(AF_ref_eur1_50_10000,file = "AF_ref_eur1_50_10000.txt", sep = ' ')
# AF_ref_afr1_50_10000 <- t(sapply(refdat1[,2], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})); write.table(AF_ref_afr1_50_10000,file = "AF_ref_afr1_50_10000.txt", sep = ' ') 
# AF_ref_sas1_50_10000 <- t(sapply(refdat1[,3], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})); write.table(AF_ref_sas1_50_10000,file = "AF_ref_sas1_50_10000.txt", sep = ' ')
# AF_ref_eas1_50_10000 <- t(sapply(refdat1[,4], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})); write.table(AF_ref_eas1_50_10000,file = "AF_ref_eas1_50_10000.txt", sep = ' ')
# AF_ref_iam1_50_10000 <- t(sapply(refdat1[,5], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})); write.table(AF_ref_iam1_50_10000,file = "AF_ref_iam1_50_10000.txt", sep = ' ')


AF_ref_eur1_500_10000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_eur1_50_10000.txt", sep=" ")) 
AF_ref_afr1_500_10000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_afr1_50_10000.txt", sep=" "))
AF_ref_sas1_500_10000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_sas1_50_10000.txt", sep=" "))
AF_ref_eas1_500_10000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_eas1_50_10000.txt", sep=" "))
AF_ref_iam1_500_10000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_iam1_50_10000.txt", sep=" "))


## Creating matrices to pass to the Summix function
Res_50a1 <- matrix(data = NA, ncol = 8, nrow = rep)

DAT_arr_AMR1 <- array(data=NA,dim = c(samp1,6,rep))

##
for (i in 1:rep){
  
  DAT_arr_AMR1[ , ,i] = cbind(AF_ref_eur1_50_10000[,i], AF_ref_afr1_50_10000[,i], AF_ref_sas1_50_10000[,i], AF_ref_eas1_50_10000[,i],  AF_ref_iam1_50_10000[,i], refdat1[,6])
  
  Res_50a1[i,] = ancestr(as.data.frame(DAT_arr_AMR1[,1:5,i]), as.data.frame(DAT_arr_AMR1[,6,i])) 
}


Res50_amr1 <- as.data.frame(cbind(Res_50a1[,1:5]))
colnames(Res50_amr1) <- c("Anc_Prop_EUR","Anc_Prop_AFR", "Anc_Prop_SAS", "Anc_Prop_EAS", "Anc_Prop_IAM")

### Output of least square errors
LSEs_50amr1_10000 <- as.data.frame(cbind(LS = Res_50a1[,6], N = rep("N=50", 50000))) # extract least square values for each replicate

## Creating datasets for plots
df50_1 <- reshape(Res50_amr1, times = c("EUR", "AFR","SAS", "EAS","IAM"), timevar = "Anc_Grp",
                  varying = list(names(Res50_amr1)), direction = "long") # convert wide to long. 

df50_1$size <- rep("N=50",5000)


##### Data Frame for AFstar Minus AFobs
Diff <- as.data.frame(cbind(EUR=AF_ref_eur1_50_10000[,1]-refdat1[,1],   AFR=AF_ref_afr1_50_10000[,1]-refdat1[,2],  SAS=AF_ref_sas1_50_10000[,1]-refdat1[,3], 
                            EAS=AF_ref_eas1_50_10000[,1] - refdat1[,4],  IAM=AF_ref_iam1_50_10000[,1]- refdat1[,5]))

DF6_1 <- reshape(Diff, times = c("EUR","AFR","SAS","EAS","IAM"), timevar = "Anc_Grp",
                 varying = list(names(Diff)), direction = "long") # convert wide to long. 

DF6_1$size <- rep("N=50", 50000)



##########################################################################
# Generate allele counts for population using binom for N = 100
##########################################################################

N1 = 100 # number of people (--> allele number = 2 * N)
rep = 1000 # number of replicates 

# AF_ref_eur1_100_10000 <- t(sapply(refdat1[,1], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})); write.table(AF_ref_eur1_100_10000,file = "AF_ref_eur1_100_10000.txt", sep = ' ')
# AF_ref_afr1_100_10000 <- t(sapply(refdat1[,2], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})); write.table(AF_ref_afr1_100_10000,file = "AF_ref_afr1_100_10000.txt", sep = ' ') 
# AF_ref_sas1_100_10000 <- t(sapply(refdat1[,3], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})); write.table(AF_ref_sas1_100_10000,file = "AF_ref_sas1_100_10000.txt", sep = ' ')
# AF_ref_eas1_100_10000 <- t(sapply(refdat1[,4], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})); write.table(AF_ref_eas1_100_10000,file = "AF_ref_eas1_100_10000.txt", sep = ' ')
# AF_ref_iam1_100_10000 <- t(sapply(refdat1[,5], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})); write.table(AF_ref_iam1_100_10000,file = "AF_ref_iam1_100_10000.txt", sep = ' ')


AF_ref_eur1_500_10000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_eur1_100_10000.txt", sep=" ")) 
AF_ref_afr1_500_10000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_afr1_100_10000.txt", sep=" "))
AF_ref_sas1_500_10000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_sas1_100_10000.txt", sep=" "))
AF_ref_eas1_500_10000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_eas1_100_10000.txt", sep=" "))
AF_ref_iam1_500_10000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_iam1_100_10000.txt", sep=" "))


## Creating matrices to pass to the Summix function
Res_100a1 <- matrix(data = NA, ncol = 8, nrow = rep)

DAT_arr_AMR1 <- array(data=NA,dim = c(samp1,6,rep))

##
for (i in 1:rep){
  
  DAT_arr_AMR1[ , ,i] = cbind(AF_ref_eur1_100_10000[,i], AF_ref_afr1_100_10000[,i], AF_ref_sas1_100_10000[,i], AF_ref_eas1_100_10000[,i],  AF_ref_iam1_100_10000[,i], refdat1[,6])
  
  Res_100a1[i,] = ancestr(as.data.frame(DAT_arr_AMR1[,1:5,i]), as.data.frame(DAT_arr_AMR1[,6,i])) 
}


Res100_amr1 <- as.data.frame(cbind(Res_100a1[,1:5]))
colnames(Res100_amr1) <- c("Anc_Prop_EUR","Anc_Prop_AFR", "Anc_Prop_SAS", "Anc_Prop_EAS", "Anc_Prop_IAM")

### Output of least square errors
#write.table(Res_100a1[,6],file = "LSEs_100a1_10000.txt", sep = ' ')
LSEs_100amr1_10000 <- as.data.frame(cbind(LS = Res_100a1[,6], N = rep("N=100", 50000))) # extract least square values for each replicate
## Creating datasets for plots
df100_1 <- reshape(Res100_amr1, times = c("EUR", "AFR","SAS", "EAS","IAM"), timevar = "Anc_Grp",
                   varying = list(names(Res100_amr1)), direction = "long") # convert wide to long. 

df100_1$size <- rep("N=100",5000)


##### Data Frame for AFstar Minus AFobs
Diff <- as.data.frame(cbind(EUR=AF_ref_eur1_100_10000[,1]-refdat1[,1],   AFR=AF_ref_afr1_100_10000[,1]-refdat1[,2],  SAS=AF_ref_sas1_100_10000[,1]-refdat1[,3], 
                            EAS=AF_ref_eas1_100_10000[,1] - refdat1[,4],  IAM=AF_ref_iam1_100_10000[,1]- refdat1[,5]))

DF7_1 <- reshape(Diff, times = c("EUR","AFR","SAS","EAS","IAM"), timevar = "Anc_Grp",
                 varying = list(names(Diff)), direction = "long") # convert wide to long. 

DF7_1$size <- rep("N=100", 50000)



##########################################################################
# Generate allele counts for population using binom for N = 500
##########################################################################
N1 = 500 # number of people (--> allele number = 2 * N)
rep = 1000 # number of replicates 

# AF_ref_eur1_500_10000 <- t(sapply(refdat1[,1], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})); write.table(AF_ref_eur1_500_10000,file = "AF_ref_eur1_500_10000.txt", sep = ' ')
# AF_ref_afr1_500_10000 <- t(sapply(refdat1[,2], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})); write.table(AF_ref_afr1_500_10000,file = "AF_ref_afr1_500_10000.txt", sep = ' ') 
# AF_ref_sas1_500_10000 <- t(sapply(refdat1[,3], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})); write.table(AF_ref_sas1_500_10000,file = "AF_ref_sas1_500_10000.txt", sep = ' ')
# AF_ref_eas1_500_10000 <- t(sapply(refdat1[,4], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})); write.table(AF_ref_eas1_500_10000,file = "AF_ref_eas1_500_10000.txt", sep = ' ')
# AF_ref_iam1_500_10000 <- t(sapply(refdat1[,5], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})); write.table(AF_ref_iam1_500_10000,file = "AF_ref_iam1_500_10000.txt", sep = ' ')


AF_ref_eur1_500_10000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_eur1_500_10000.txt", sep=" ")) 
AF_ref_afr1_500_10000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_afr1_500_10000.txt", sep=" "))
AF_ref_sas1_500_10000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_sas1_500_10000.txt", sep=" "))
AF_ref_eas1_500_10000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_eas1_500_10000.txt", sep=" "))
AF_ref_iam1_500_10000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_iam1_500_10000.txt", sep=" "))


## Creating matrices to pass to the Summix function
Res_500a1 <- matrix(data = NA, ncol = 8, nrow = rep)

DAT_arr_AMR1 <- array(data=NA,dim = c(samp1,6,rep))

##
for (i in 1:rep){
  
  DAT_arr_AMR1[ , ,i] = cbind(AF_ref_eur1_500_10000[,i], AF_ref_afr1_500_10000[,i], AF_ref_sas1_500_10000[,i], AF_ref_eas1_500_10000[,i],  AF_ref_iam1_500_10000[,i], refdat1[,6])
  
  Res_500a1[i,] = ancestr(as.data.frame(DAT_arr_AMR1[,1:5,i]), as.data.frame(DAT_arr_AMR1[,6,i])) 
}


Res500_amr1 <- as.data.frame(cbind(Res_500a1[,1:5]))
colnames(Res500_amr1) <- c("Anc_Prop_EUR","Anc_Prop_AFR", "Anc_Prop_SAS", "Anc_Prop_EAS", "Anc_Prop_IAM")

### Output of least square errors
LSEs_500amr1_10000 <- as.data.frame(cbind(LS = Res_500a1[,6], N = rep("N=500", 50000))) # extract least square values for each replicate

## Creating datasets for plots
df500_1 <- reshape(Res500_amr1, times = c("EUR", "AFR","SAS", "EAS","IAM"), timevar = "Anc_Grp",
                   varying = list(names(Res500_amr1)), direction = "long") # convert wide to long. 

df500_1$size <- rep("N=500",5000)


##### Data Frame for AFstar Minus AFobs
Diff <- as.data.frame(cbind(EUR=AF_ref_eur1_500_10000[,1]-refdat1[,1],   AFR=AF_ref_afr1_500_10000[,1]-refdat1[,2],  SAS=AF_ref_sas1_500_10000[,1]-refdat1[,3], 
                            EAS=AF_ref_eas1_500_10000[,1] - refdat1[,4],  IAM=AF_ref_iam1_500_10000[,1]- refdat1[,5]))

DF8_1 <- reshape(Diff, times = c("EUR","AFR","SAS","EAS","IAM"), timevar = "Anc_Grp",
                 varying = list(names(Diff)), direction = "long") # convert wide to long. 

DF8_1$size <- rep("N=500", 50000)



##########################################################################
# Generate allele counts for population using binom for N = 1000
##########################################################################
N1 = 1000 # number of people (--> allele number = 2 * N)
rep = 1000 # number of replicates 

# AF_ref_eur1_1000_10000 <- t(sapply(refdat1[,1], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})); write.table(AF_ref_eur1_1000_10000,file = "AF_ref_eur1_1000_10000.txt", sep = ' ')
# AF_ref_afr1_1000_10000 <- t(sapply(refdat1[,2], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})); write.table(AF_ref_afr1_1000_10000,file = "AF_ref_afr1_1000_10000.txt", sep = ' ') 
# AF_ref_sas1_1000_10000 <- t(sapply(refdat1[,3], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})); write.table(AF_ref_sas1_1000_10000,file = "AF_ref_sas1_1000_10000.txt", sep = ' ')
# AF_ref_eas1_1000_10000 <- t(sapply(refdat1[,4], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})); write.table(AF_ref_eas1_1000_10000,file = "AF_ref_eas1_1000_10000.txt", sep = ' ')
# AF_ref_iam1_1000_10000 <- t(sapply(refdat1[,5], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})); write.table(AF_ref_iam1_1000_10000,file = "AF_ref_iam1_1000_10000.txt", sep = ' ')

AF_ref_eur1_1000_10000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_eur1_1000_10000.txt", sep=" ")) 
AF_ref_afr1_1000_10000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_afr1_1000_10000.txt", sep=" "))
AF_ref_sas1_1000_10000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_sas1_1000_10000.txt", sep=" "))
AF_ref_eas1_1000_10000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_eas1_1000_10000.txt", sep=" "))
AF_ref_iam1_1000_10000 <- as.matrix(read.table("/newhome/agboolol/AF_matrices_Latinx/AF_ref_iam1_1000_10000.txt", sep=" "))


## Creating matrices to pass to the Summix function
Res_1000a1 <- matrix(data = NA, ncol = 8, nrow = rep)

DAT_arr_AMR1 <- array(data=NA,dim = c(samp1,6,rep))

##
for (i in 1:rep){
  
  DAT_arr_AMR1[ , ,i] = cbind(AF_ref_eur1_1000_10000[,i], AF_ref_afr1_1000_10000[,i], AF_ref_sas1_1000_10000[,i], AF_ref_eas1_1000_10000[,i],  AF_ref_iam1_1000_10000[,i], refdat1[,6])
  
  Res_1000a1[i,] = ancestr(as.data.frame(DAT_arr_AMR1[,1:5,i]), as.data.frame(DAT_arr_AMR1[,6,i])) 
}


Res1000_amr1 <- as.data.frame(cbind(Res_1000a1[,1:5]))
colnames(Res1000_amr1) <- c("Anc_Prop_EUR","Anc_Prop_AFR", "Anc_Prop_SAS", "Anc_Prop_EAS", "Anc_Prop_IAM")

### Output of least square errors
LSEs_1000amr1_10000 <- as.data.frame(cbind(LS = Res_1000a1[,6], N = rep("N=1000", 50000))) # extract least square values for each replicate

## Creating datasets for plots
df1000_1 <- reshape(Res1000_amr1, times = c("EUR", "AFR","SAS", "EAS","IAM"), timevar = "Anc_Grp",
                    varying = list(names(Res1000_amr1)), direction = "long") # convert wide to long. 

df1000_1$size <- rep("N=1000",5000)


##### Data Frame for AFstar Minus AFobs
Diff <- as.data.frame(cbind(EUR=AF_ref_eur1_1000_10000[,1]-refdat1[,1],   AFR=AF_ref_afr1_1000_10000[,1]-refdat1[,2],  SAS=AF_ref_sas1_1000_10000[,1]-refdat1[,3], 
                            EAS=AF_ref_eas1_1000_10000[,1] - refdat1[,4],  IAM=AF_ref_iam1_1000_10000[,1]- refdat1[,5]))

DF9_1 <- reshape(Diff, times = c("EUR","AFR","SAS","EAS","IAM"), timevar = "Anc_Grp",
                 varying = list(names(Diff)), direction = "long") # convert wide to long. 

DF9_1$size  <-  rep("N=1000", 50000)



##################################################
####Combined data
##################################################
# Combining all ancestry proportions estimates for each N
df.combined1 <- rbind(df5_1,df10_1, df15_1, df20_1, df25_1, df50_1, df100_1, df500_1, df1000_1)

# Combining difference in AFs for each N
DF.combine1 <- rbind(DF1_1, DF2_1, DF3_1, DF4_1, DF5_1, DF6_1, DF7_1, DF8_1, DF9_1)

# combining all least square values for each N
LS.combine1 <- rbind(LSEs_5amr1_10000, LSEs_10amr1_10000, LSEs_15amr1_10000, LSEs_20amr1_10000, LSEs_25amr1_10000, LSEs_50amr1_10000,
                    LSEs_100amr1_10000, LSEs_500amr1_10000, LSEs_1000amr1_10000)


# dev.off()
pdf('Rplot_SNP10000_facA2.pdf')
# #M <- rep("N = 10",2000) # x variable
ggplot(df.combined1, aes(x=size, y=Anc_Prop_EUR, color=Anc_Grp)) + #x = N, y = ancestry prop, fill = ancestry groups
  geom_boxplot() + labs(x = "Allele Number (2 * N)", y = "Ancestry Proportion")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Variation of Summix Ancestry Proportion Estimate (SNP = 10,000) \n American/ Latinx GnomAD")+
  scale_x_discrete(limits=c("N=5","N=10","N=15","N=20","N=25","N=50","N=100","N=500","N=1000"))+
  facet_wrap(~Anc_Grp, scales = "free")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


pdf('Rplot_SNP10000_diff1D2.pdf')
# #M <- rep("N = 10",2000) # x variable
ggplot(DF.combine1, aes(x=size, y=EUR, color=Anc_Grp)) + #x = N, y = ancestry prop, fill = ancestry groups
  geom_boxplot() + labs(x = "Allele Number (2 * N)", y = "Difference in Allele Frequencies")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Distribution of Difference in Allele Frequencies (SNP = 10,000) \n American/ Latinx GnomAD")+
  scale_x_discrete(limits=c("N=5","N=10","N=15","N=20","N=25","N=50","N=100","N=500","N=1000"))+
  facet_wrap(~Anc_Grp, scales = "free")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf('Rplot_SNP10000Latinx_LS.pdf')
ggplot(LS.combine1, aes(x=N, y=as.numeric(LS))) +
  geom_boxplot() + labs(x = "Allele Number (2 * N)", y = "Least Square Values")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Distribution of Least Square Values (SNP = 10,000) \n American/ Latinx GnomAD")+
  scale_x_discrete(limits=c("N=5","N=10","N=15","N=20","N=25","N=50","N=100","N=500","N=1000"))+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()



df.combined1 <- df.combined1 %>% 
  mutate(N = ifelse(size=="N=5",5,
                    ifelse(size == "N=10", 10,
                           ifelse(size == "N=15", 15,
                                  ifelse(size == "N=20", 20,
                                         ifelse(size == "N=25", 25,
                                                ifelse(size == "N=50", 50,
                                                       ifelse(size == "N=100", 100,
                                                              ifelse(size == "N=500", 500, 1000)))))))))

dfB <- df.combined1 %>%
  group_by(Anc_Grp,N) %>%
  summarize(mean = mean(Anc_Prop_EUR),
            std = sd(Anc_Prop_EUR),
            minimum =min(Anc_Prop_EUR), 
            Q1=quantile(Anc_Prop_EUR, .25), 
            median = median(Anc_Prop_EUR),
            Q3=quantile(Anc_Prop_EUR, .75),
            maximum = max(Anc_Prop_EUR), .groups = 'keep')

write.table(dfB, file = "SNPs_10000A2.txt", sep = ' ')


DF.combine1 <- DF.combine1 %>% 
  mutate(N = ifelse(size=="N=5",5,
                    ifelse(size == "N=10", 10,
                           ifelse(size == "N=15", 15,
                                  ifelse(size == "N=20", 20,
                                         ifelse(size == "N=25", 25,
                                                ifelse(size == "N=50", 50, 
                                                       ifelse(size == "N=100", 100, 
                                                              ifelse(size == "N=500", 500, 1000)))))))))

dfB1 <- DF.combine1 %>%
  group_by(Anc_Grp,N) %>%
  summarize(mean = mean(EUR),
            std = sd(EUR),
            minimum =min(EUR), 
            Q1=quantile(EUR, .25), 
            median = median(EUR),
            Q3=quantile(EUR, .75),
            maximum = max(EUR), .groups = 'keep')

write.table(dfB1, file = "SNPs_10000D2_1.txt", sep = ' ')


# Create plot for variability in simulated and observed AFs for each 1/N and 1/N^2
 pdf('Rplot_SNP10000Latinx_Vardiff_Ninv.pdf')
 ggplot(dfB1, aes(x=as.numeric(1/N), y=as.numeric(std^2), color=Anc_Grp)) +
   geom_point() + labs(x = TeX("$\\frac{1}{N}$"), y = "Variance")+
   theme(plot.title = element_text(hjust = 0.5))+
   ggtitle("Variability of Difference in Allele Frequencies (SNP = 10,000) \n American/ Latinx GnomAD")+
   facet_wrap(~Anc_Grp, scales = "free")
 dev.off()

 pdf('Rplot_SNP10000Latinx_Vardiff_Nsqinv.pdf')
 ggplot(dfB1, aes(x=as.numeric(1/(N^2)), y=as.numeric(std^2), color=Anc_Grp)) +
   geom_point() + labs(x = TeX("$\\frac{1}{N^2}$"), y = "Variance")+
   theme(plot.title = element_text(hjust = 0.5))+
   ggtitle("Variability of Difference in Allele Frequencies (SNP = 10,000) \n American/ Latinx GnomAD")+
   facet_wrap(~Anc_Grp, scales = "free")
 dev.off()


 ############################################################################
 ## We were not able to replicate the above analyses for (SNP = 100,000).
 ## We run into an error Error:vector memory exhausted (limit reached?)
 ## trying to run this locally on our computer.
 ## 
 ## As a next step we may run the above analyses for SNP = 100,000. However,
 ## note that when SNP was increased from 1000 to 10000 we saw a decrease in 
 ## the variability of ancestry proportion estimates for the five ancestry 
 ## groups.
############################################################################
 
 # ##########################################################################
 # ##### Selecting 100K SNPS from reference data(SNP = 100,000)
 # ##### For N = 5, 10, 15, 20, 25, 10, 50, 100, 500, 1000
 # ##########################################################################
 # 
 # #Set Seed
 # set.seed(2656256)
 # samp2 = 100000
 # 
 # refdat2 = data_red %>% 
 #   sample_n(samp2) %>% 
 #   select(ref_AF_eur_1000G,ref_AF_afr_1000G, ref_AF_sas_1000G,ref_AF_eas_1000G,ref_AF_iam_1000G, gnomad_AF_amr) %>% 
 #   
 #   rename(ref_AFR = ref_AF_afr_1000G, ref_EUR = ref_AF_eur_1000G, ref_SAS = ref_AF_sas_1000G, ref_EAS = ref_AF_eas_1000G,
 #          ref_IAM = ref_AF_iam_1000G, obs_AMR = gnomad_AF_amr)
 # 
 # 
 # 
 # ##########################################################################
 # # Generate allele counts for population using binom for N = 5
 # ##########################################################################
 # 
 # N2 = 5 # number of people (--> allele number = 2 * N)
 # rep = 1000 # number of replicates 
 # 
 # AF_ref_eur2_5_1000 <- t(sapply(refdat2[,1], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})); write.table(AF_ref_eur2_5_1000,file = "AF_ref_eur2_5_1000.txt", sep = ' ')
 # AF_ref_afr2_5_1000 <- t(sapply(refdat2[,2], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})); write.table(AF_ref_afr2_5_1000,file = "AF_ref_afr2_5_1000.txt", sep = ' ') 
 # AF_ref_sas2_5_1000 <- t(sapply(refdat2[,3], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})); write.table(AF_ref_sas2_5_1000,file = "AF_ref_sas2_5_1000.txt", sep = ' ')
 # AF_ref_eas2_5_1000 <- t(sapply(refdat2[,4], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})); write.table(AF_ref_eas2_5_1000,file = "AF_ref_eas2_5_1000.txt", sep = ' ')
 # AF_ref_iam2_5_1000 <- t(sapply(refdat2[,5], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})); write.table(AF_ref_iam2_5_1000,file = "AF_ref_iam2_5_1000.txt", sep = ' ')
 # 
 # 
 # ## Creating matrices to pass to the Summix function
 # Res_5a2 <- matrix(data = NA, ncol = 8, nrow = rep)
 # 
 # DAT_arr_AMR2 <- array(data=NA,dim = c(samp2,6,rep))
 # 
 # ##
 # for (i in 1:rep){
 #   
 #   DAT_arr_AMR2[ , ,i] = cbind(AF_ref_eur2_5_1000[,i], AF_ref_afr2_5_1000[,i], AF_ref_sas2_5_1000[,i], AF_ref_eas2_5_1000[,i],  AF_ref_iam2_5_1000[,i], refdat2[,6])
 #   
 #   Res_5a2[i,] = ancestr(as.data.frame(DAT_arr_AMR2[,1:5,i]), as.data.frame(DAT_arr_AMR2[,6,i])) 
 # }
 # 
 # 
 # Res5_amr2 <- as.data.frame(cbind(Res_5a2[,1:5]))
 # colnames(Res5_amr2) <- c("Anc_Prop_EUR","Anc_Prop_AFR", "Anc_Prop_SAS", "Anc_Prop_EAS", "Anc_Prop_IAM")
 # 
 # ### Output of least square errors
 # write.table(Res_5a2[,6],file = "LSEs_5a2_100000.txt", sep = ' ')
 # 
 # ## Creating datasets for plots
 # df5_2 <- reshape(Res5_amr2, times = c("EUR", "AFR","SAS", "EAS","IAM"), timevar = "Anc_Grp",
 #                  varying = list(names(Res5_amr2)), direction = "long") # convert wide to long. 
 # 
 # df5_2$size <- rep("N=5",5000)
 # 
 # 
 # ##### Data Frame for AFstar Minus AFobs
 # Diff <- as.data.frame(cbind(EUR=AF_ref_eur2_5_1000[,1]-refdat2[,1],   AFR=AF_ref_afr2_5_1000[,1]-refdat2[,2],  SAS=AF_ref_sas2_5_1000[,1]-refdat2[,3], 
 #                             EAS=AF_ref_eas2_5_1000[,1] - refdat2[,4],  IAM=AF_ref_iam2_5_1000[,1]- refdat2[,5]))
 # 
 # DF1_2 <- reshape(Diff, times = c("EUR","AFR","SAS","EAS","IAM"), timevar = "Anc_Grp",
 #                  varying = list(names(Diff)), direction = "long") # convert wide to long. 
 # 
 # DF1_2$size <- rep("N=5", 500000)
 # 
 # 
 # 
 # ##########################################################################
 # # Generate allele counts for population using binom for N = 10
 # ##########################################################################
 # N2 = 10 # number of people (--> allele number = 2 * N)
 # rep = 1000 # number of replicates 
 # 
 # AF_ref_eur2_10_1000 <- t(sapply(refdat2[,1], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})); write.table(AF_ref_eur2_10_1000,file = "AF_ref_eur2_10_1000.txt", sep = ' ')
 # AF_ref_afr2_10_1000 <- t(sapply(refdat2[,2], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})); write.table(AF_ref_afr2_10_1000,file = "AF_ref_afr2_10_1000.txt", sep = ' ')
 # AF_ref_sas2_10_1000 <- t(sapply(refdat2[,3], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})); write.table(AF_ref_sas2_10_1000,file = "AF_ref_sas2_10_1000.txt", sep = ' ')
 # AF_ref_eas2_10_1000 <- t(sapply(refdat2[,4], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})); write.table(AF_ref_eas2_10_1000,file = "AF_ref_eas2_10_1000.txt", sep = ' ')
 # AF_ref_iam2_10_1000 <- t(sapply(refdat2[,5], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})); write.table(AF_ref_iam2_10_1000,file = "AF_ref_iam2_10_1000.txt", sep = ' ')
 # 
 # 
 # ## Creating matrices to pass to the Summix function
 # Res_10a2 <- matrix(data = NA, ncol = 8, nrow = rep)
 # 
 # DAT_arr_AMR2 <- array(data=NA,dim = c(samp2,6,rep))
 # 
 # ##
 # for (i in 1:rep){
 #   
 #   DAT_arr_AMR2[ , ,i] = cbind(AF_ref_eur2_10_1000[,i], AF_ref_afr2_10_1000[,i], AF_ref_sas2_10_1000[,i], AF_ref_eas2_10_1000[,i],  AF_ref_iam2_10_1000[,i], refdat2[,6])
 #   
 #   Res_10a2[i,] = ancestr(as.data.frame(DAT_arr_AMR2[,1:5,i]), as.data.frame(DAT_arr_AMR2[,6,i])) 
 # }
 # 
 # 
 # Res10_amr2 <- as.data.frame(cbind(Res_10a2[,1:5]))
 # colnames(Res10_amr2) <- c("Anc_Prop_EUR","Anc_Prop_AFR", "Anc_Prop_SAS", "Anc_Prop_EAS", "Anc_Prop_IAM")
 # 
 # ### Output of least square errors
 # write.table(Res_10a2[,6],file = "LSEs_10a2_100000.txt", sep = ' ')
 # 
 # ## Creating datasets for plots
 # df10_2 <- reshape(Res10_amr2, times = c("EUR", "AFR","SAS", "EAS","IAM"), timevar = "Anc_Grp",
 #                   varying = list(names(Res10_amr2)), direction = "long") # convert wide to long. 
 # 
 # df10_2$size <- rep("N=10",5000)
 # 
 # 
 # ##### Data Frame for AFstar Minus AFobs
 # Diff <- as.data.frame(cbind(EUR=AF_ref_eur2_10_1000[,1]-refdat2[,1],   AFR=AF_ref_afr2_10_1000[,1]-refdat2[,2],  SAS=AF_ref_sas2_10_1000[,1]-refdat2[,3], 
 #                             EAS=AF_ref_eas2_10_1000[,1] - refdat2[,4],  IAM=AF_ref_iam2_10_1000[,1]- refdat2[,5]))
 # 
 # DF2_2 <- reshape(Diff, times = c("EUR","AFR","SAS","EAS","IAM"), timevar = "Anc_Grp",
 #                  varying = list(names(Diff)), direction = "long") # convert wide to long. 
 # 
 # DF2_2$size <- rep("N=10", 500000)
 # 
 # 
 # 
 # 
 # 
 # ##########################################################################
 # # Generate allele counts for population using binom for N = 15
 # ##########################################################################
 # N2 = 15 # number of people (--> allele number = 2 * N)
 # rep = 1000 # number of replicates 
 # 
 # AF_ref_eur2_15_1000 <- t(sapply(refdat2[,1], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})); write.table(AF_ref_eur2_15_1000,file = "AF_ref_eur2_15_1000.txt", sep = ' ')
 # AF_ref_afr2_15_1000 <- t(sapply(refdat2[,2], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})); write.table(AF_ref_afr2_15_1000,file = "AF_ref_afr2_15_1000.txt", sep = ' ') 
 # AF_ref_sas2_15_1000 <- t(sapply(refdat2[,3], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})); write.table(AF_ref_sas2_15_1000,file = "AF_ref_sas2_15_1000.txt", sep = ' ')
 # AF_ref_eas2_15_1000 <- t(sapply(refdat2[,4], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})); write.table(AF_ref_eas2_15_1000,file = "AF_ref_eas2_15_1000.txt", sep = ' ')
 # AF_ref_iam2_15_1000 <- t(sapply(refdat2[,5], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})); write.table(AF_ref_iam2_15_1000,file = "AF_ref_iam2_15_1000.txt", sep = ' ')
 # 
 # ## Creating matrices to pass to the Summix function
 # Res_15a2 <- matrix(data = NA, ncol = 8, nrow = rep)
 # 
 # DAT_arr_AMR2 <- array(data=NA,dim = c(samp2,6,rep))
 # 
 # ##
 # for (i in 1:rep){
 #   
 #   DAT_arr_AMR2[ , ,i] = cbind(AF_ref_eur2_15_1000[,i], AF_ref_afr2_15_1000[,i], AF_ref_sas2_15_1000[,i], AF_ref_eas2_15_1000[,i],  AF_ref_iam2_15_1000[,i], refdat2[,6])
 #   
 #   Res_15a2[i,] = ancestr(as.data.frame(DAT_arr_AMR2[,1:5,i]), as.data.frame(DAT_arr_AMR2[,6,i])) 
 # }
 # 
 # 
 # Res15_amr2 <- as.data.frame(cbind(Res_15a2[,1:5]))
 # colnames(Res15_amr2) <- c("Anc_Prop_EUR","Anc_Prop_AFR", "Anc_Prop_SAS", "Anc_Prop_EAS", "Anc_Prop_IAM")
 # 
 # ### Output of least square errors
 # write.table(Res_15a2[,6],file = "LSEs_15a2_100000.txt", sep = ' ')
 # 
 # ## Creating datasets for plots
 # df15_2 <- reshape(Res15_amr2, times = c("EUR", "AFR","SAS", "EAS","IAM"), timevar = "Anc_Grp",
 #                   varying = list(names(Res15_amr2)), direction = "long") # convert wide to long. 
 # 
 # df15_2$size <- rep("N=15",5000)
 # 
 # 
 # ##### Data Frame for AFstar Minus AFobs
 # Diff <- as.data.frame(cbind(EUR=AF_ref_eur2_15_1000[,1]-refdat2[,1],   AFR=AF_ref_afr2_15_1000[,1]-refdat2[,2],  SAS=AF_ref_sas2_15_1000[,1]-refdat2[,3], 
 #                             EAS=AF_ref_eas2_15_1000[,1] - refdat2[,4],  IAM=AF_ref_iam2_15_1000[,1]- refdat2[,5]))
 # 
 # DF3_2 <- reshape(Diff, times = c("EUR","AFR","SAS","EAS","IAM"), timevar = "Anc_Grp",
 #                  varying = list(names(Diff)), direction = "long") # convert wide to long. 
 # 
 # DF3_2$size <- rep("N=15", 500000)
 # 
 # 
 # 
 # 
 # 
 # 
 # 
 # ##########################################################################
 # # Generate allele counts for population using binom for N = 20
 # ##########################################################################
 # N2 = 20 # number of people (--> allele number = 2 * N)
 # rep = 1000 # number of replicates 
 # 
 # AF_ref_eur2_20_1000 <- t(sapply(refdat2[,1], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})); write.table(AF_ref_eur2_20_1000,file = "AF_ref_eur2_20_1000.txt", sep = ' ')
 # AF_ref_afr2_20_1000 <- t(sapply(refdat2[,2], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})); write.table(AF_ref_afr2_20_1000,file = "AF_ref_afr2_20_1000.txt", sep = ' ') 
 # AF_ref_sas2_20_1000 <- t(sapply(refdat2[,3], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})); write.table(AF_ref_sas2_20_1000,file = "AF_ref_sas2_20_1000.txt", sep = ' ')
 # AF_ref_eas2_20_1000 <- t(sapply(refdat2[,4], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})); write.table(AF_ref_eas2_20_1000,file = "AF_ref_eas2_20_1000.txt", sep = ' ')
 # AF_ref_iam2_20_1000 <- t(sapply(refdat2[,5], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})); write.table(AF_ref_iam2_20_1000,file = "AF_ref_iam2_20_1000.txt", sep = ' ')
 # 
 # 
 # ## Creating matrices to pass to the Summix function
 # Res_15a2 <- matrix(data = NA, ncol = 8, nrow = rep)
 # 
 # DAT_arr_AMR2 <- array(data=NA,dim = c(samp2,6,rep))
 # 
 # ##
 # for (i in 1:rep){
 #   
 #   DAT_arr_AMR2[ , ,i] = cbind(AF_ref_eur2_20_1000[,i], AF_ref_afr2_20_1000[,i], AF_ref_sas2_20_1000[,i], AF_ref_eas2_20_1000[,i],  AF_ref_iam2_20_1000[,i], refdat2[,6])
 #   
 #   Res_20a2[i,] = ancestr(as.data.frame(DAT_arr_AMR2[,1:5,i]), as.data.frame(DAT_arr_AMR2[,6,i])) 
 # }
 # 
 # 
 # Res20_amr2 <- as.data.frame(cbind(Res_20a2[,1:5]))
 # colnames(Res20_amr2) <- c("Anc_Prop_EUR","Anc_Prop_AFR", "Anc_Prop_SAS", "Anc_Prop_EAS", "Anc_Prop_IAM")
 # 
 # ### Output of least square errors
 # write.table(Res_20a2[,6],file = "LSEs_20a2_100000.txt", sep = ' ')
 # 
 # ## Creating datasets for plots
 # df20_2 <- reshape(Res20_amr2, times = c("EUR", "AFR","SAS", "EAS","IAM"), timevar = "Anc_Grp",
 #                   varying = list(names(Res20_amr2)), direction = "long") # convert wide to long. 
 # 
 # df20_2$size <- rep("N=20",5000)
 # 
 # 
 # ##### Data Frame for AFstar Minus AFobs
 # Diff <- as.data.frame(cbind(EUR=AF_ref_eur2_20_1000[,1]-refdat2[,1],   AFR=AF_ref_afr2_20_1000[,1]-refdat2[,2],  SAS=AF_ref_sas2_20_1000[,1]-refdat2[,3], 
 #                             EAS=AF_ref_eas2_20_1000[,1] - refdat2[,4],  IAM=AF_ref_iam2_20_1000[,1]- refdat2[,5]))
 # 
 # DF4_2 <- reshape(Diff, times = c("EUR","AFR","SAS","EAS","IAM"), timevar = "Anc_Grp",
 #                  varying = list(names(Diff)), direction = "long") # convert wide to long. 
 # 
 # DF4_2$size <- rep("N=20", 500000)
 # 
 # 
 # 
 # 
 # 
 # ##########################################################################
 # # Generate allele counts for population using binom for N = 25
 # ##########################################################################
 # N2 = 25 # number of people (--> allele number = 2 * N)
 # rep = 1000 # number of replicates 
 # 
 # AF_ref_eur2_25_1000 <- t(sapply(refdat2[,1], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})); write.table(AF_ref_eur2_25_1000,file = "AF_ref_eur2_25_1000.txt", sep = ' ')
 # AF_ref_afr2_25_1000 <- t(sapply(refdat2[,2], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})); write.table(AF_ref_afr2_25_1000,file = "AF_ref_afr2_25_1000.txt", sep = ' ') 
 # AF_ref_sas2_25_1000 <- t(sapply(refdat2[,3], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})); write.table(AF_ref_sas2_25_1000,file = "AF_ref_sas2_25_1000.txt", sep = ' ')
 # AF_ref_eas2_25_1000 <- t(sapply(refdat2[,4], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})); write.table(AF_ref_eas2_25_1000,file = "AF_ref_eas2_25_1000.txt", sep = ' ')
 # AF_ref_iam2_25_1000 <- t(sapply(refdat2[,5], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})); write.table(AF_ref_iam2_25_1000,file = "AF_ref_iam2_25_1000.txt", sep = ' ')
 # 
 # 
 # ## Creating matrices to pass to the Summix function
 # Res_25a2 <- matrix(data = NA, ncol = 8, nrow = rep)
 # 
 # DAT_arr_AMR2 <- array(data=NA,dim = c(samp2,6,rep))
 # 
 # ##
 # for (i in 1:rep){
 #   
 #   DAT_arr_AMR2[ , ,i] = cbind(AF_ref_eur2_25_1000[,i], AF_ref_afr2_25_1000[,i], AF_ref_sas2_25_1000[,i], AF_ref_eas2_25_1000[,i],  AF_ref_iam2_25_1000[,i], refdat2[,6])
 #   
 #   Res_25a2[i,] = ancestr(as.data.frame(DAT_arr_AMR2[,1:5,i]), as.data.frame(DAT_arr_AMR2[,6,i])) 
 # }
 # 
 # 
 # Res25_amr2 <- as.data.frame(cbind(Res_25a2[,1:5]))
 # colnames(Res20_amr2) <- c("Anc_Prop_EUR","Anc_Prop_AFR", "Anc_Prop_SAS", "Anc_Prop_EAS", "Anc_Prop_IAM")
 # 
 # ### Output of least square errors
 # write.table(Res_25a2[,6],file = "LSEs_25a2_100000.txt", sep = ' ')
 # 
 # ## Creating datasets for plots
 # df25_2 <- reshape(Res25_amr2, times = c("EUR", "AFR","SAS", "EAS","IAM"), timevar = "Anc_Grp",
 #                   varying = list(names(Res25_amr2)), direction = "long") # convert wide to long. 
 # 
 # df25_2$size <- rep("N=25",5000)
 # 
 # 
 # ##### Data Frame for AFstar Minus AFobs
 # Diff <- as.data.frame(cbind(EUR=AF_ref_eur2_25_1000[,1]-refdat2[,1],   AFR=AF_ref_afr2_25_1000[,1]-refdat2[,2],  SAS=AF_ref_sas2_25_1000[,1]-refdat2[,3], 
 #                             EAS=AF_ref_eas2_25_1000[,1] - refdat2[,4],  IAM=AF_ref_iam2_25_1000[,1]- refdat2[,5]))
 # 
 # DF5_2 <- reshape(Diff, times = c("EUR","AFR","SAS","EAS","IAM"), timevar = "Anc_Grp",
 #                  varying = list(names(Diff)), direction = "long") # convert wide to long. 
 # 
 # DF5_2$size <- rep("N=25", 500000)
 # 
 # 
 # 
 # 
 # 
 # 
 # ##########################################################################
 # # Generate allele counts for population using binom for N = 50
 # ##########################################################################
 # N2 = 50 # number of people (--> allele number = 2 * N)
 # rep = 1000 # number of replicates 
 # 
 # AF_ref_eur2 <- t(sapply(refdat2[,1], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)}))
 # AF_ref_afr2 <- t(sapply(refdat2[,2], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})) 
 # AF_ref_sas2 <- t(sapply(refdat2[,3], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)}))
 # AF_ref_eas2 <- t(sapply(refdat2[,4], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)}))
 # AF_ref_iam2 <- t(sapply(refdat2[,5], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)}))
 # 
 # AF_ref_eur2_50_1000 <- t(sapply(refdat2[,1], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})); write.table(AF_ref_eur2_50_1000,file = "AF_ref_eur2_50_1000.txt", sep = ' ')
 # AF_ref_afr2_50_1000 <- t(sapply(refdat2[,2], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})); write.table(AF_ref_afr2_50_1000,file = "AF_ref_afr2_50_1000.txt", sep = ' ') 
 # AF_ref_sas2_50_1000 <- t(sapply(refdat2[,3], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})); write.table(AF_ref_sas2_50_1000,file = "AF_ref_sas2_50_1000.txt", sep = ' ')
 # AF_ref_eas2_50_1000 <- t(sapply(refdat2[,4], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})); write.table(AF_ref_eas2_50_1000,file = "AF_ref_eas2_50_1000.txt", sep = ' ')
 # AF_ref_iam2_50_1000 <- t(sapply(refdat2[,5], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})); write.table(AF_ref_iam2_50_1000,file = "AF_ref_iam2_50_1000.txt", sep = ' ')
 # 
 # 
 # ## Creating matrices to pass to the Summix function
 # Res_50a2 <- matrix(data = NA, ncol = 8, nrow = rep)
 # 
 # DAT_arr_AMR2 <- array(data=NA,dim = c(samp2,6,rep))
 # 
 # ##
 # for (i in 1:rep){
 #   
 #   DAT_arr_AMR2[ , ,i] = cbind(AF_ref_eur2_50_1000[,i], AF_ref_afr2_50_1000[,i], AF_ref_sas2_50_1000[,i], AF_ref_eas2_50_1000[,i],  AF_ref_iam2_50_1000[,i], refdat2[,6])
 #   
 #   Res_50a2[i,] = ancestr(as.data.frame(DAT_arr_AMR2[,1:5,i]), as.data.frame(DAT_arr_AMR2[,6,i])) 
 # }
 # 
 # 
 # Res50_amr2 <- as.data.frame(cbind(Res_50a2[,1:5]))
 # colnames(Res50_amr2) <- c("Anc_Prop_EUR","Anc_Prop_AFR", "Anc_Prop_SAS", "Anc_Prop_EAS", "Anc_Prop_IAM")
 # 
 # ### Output of least square errors
 # write.table(Res_50a2[,6],file = "LSEs_50a2_100000.txt", sep = ' ')
 # 
 # ## Creating datasets for plots
 # df50_2 <- reshape(Res50_amr2, times = c("EUR", "AFR","SAS", "EAS","IAM"), timevar = "Anc_Grp",
 #                   varying = list(names(Res50_amr2)), direction = "long") # convert wide to long. 
 # 
 # df50_2$size <- rep("N=50",5000)
 # 
 # 
 # ##### Data Frame for AFstar Minus AFobs
 # Diff <- as.data.frame(cbind(EUR=AF_ref_eur2_50_1000[,1]-refdat2[,1],   AFR=AF_ref_afr2_50_1000[,1]-refdat2[,2],  SAS=AF_ref_sas2_50_1000[,1]-refdat2[,3], 
 #                             EAS=AF_ref_eas2_50_1000[,1] - refdat2[,4],  IAM=AF_ref_iam2_50_1000[,1]- refdat2[,5]))
 # 
 # DF6_2 <- reshape(Diff, times = c("EUR","AFR","SAS","EAS","IAM"), timevar = "Anc_Grp",
 #                  varying = list(names(Diff)), direction = "long") # convert wide to long. 
 # 
 # DF6_2$size <- rep("N=50", 500000)
 # 
 # 
 # 
 # 
 # 
 # 
 # 
 # 
 # ##########################################################################
 # # Generate allele counts for population using binom for N = 100
 # ##########################################################################
 # N2 = 100 # number of people (--> allele number = 2 * N)
 # rep = 1000 # number of replicates 
 # 
 # AF_ref_eur2_100_1000 <- t(sapply(refdat2[,1], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})); write.table(AF_ref_eur2_100_1000,file = "AF_ref_eur2_100_1000.txt", sep = ' ')
 # AF_ref_afr2_100_1000 <- t(sapply(refdat2[,2], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})); write.table(AF_ref_afr2_100_1000,file = "AF_ref_afr2_100_1000.txt", sep = ' ') 
 # AF_ref_sas2_100_1000 <- t(sapply(refdat2[,3], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})); write.table(AF_ref_sas2_100_1000,file = "AF_ref_sas2_100_1000.txt", sep = ' ')
 # AF_ref_eas2_100_1000 <- t(sapply(refdat2[,4], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})); write.table(AF_ref_eas2_100_1000,file = "AF_ref_eas2_100_1000.txt", sep = ' ')
 # AF_ref_iam2_100_1000 <- t(sapply(refdat2[,5], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})); write.table(AF_ref_iam2_100_1000,file = "AF_ref_iam2_100_1000.txt", sep = ' ')
 # 
 # 
 # ## Creating matrices to pass to the Summix function
 # Res_100a2 <- matrix(data = NA, ncol = 8, nrow = rep)
 # 
 # DAT_arr_AMR2 <- array(data=NA,dim = c(samp2,6,rep))
 # 
 # ##
 # for (i in 1:rep){
 #   
 #   DAT_arr_AMR2[ , ,i] = cbind(AF_ref_eur2_100_1000[,i], AF_ref_afr2_100_1000[,i], AF_ref_sas2_100_1000[,i], AF_ref_eas2_100_1000[,i],  AF_ref_iam2_100_1000[,i], refdat2[,6])
 #   
 #   Res_100a2[i,] = ancestr(as.data.frame(DAT_arr_AMR2[,1:5,i]), as.data.frame(DAT_arr_AMR2[,6,i])) 
 # }
 # 
 # 
 # Res100_amr2 <- as.data.frame(cbind(Res_100a2[,1:5]))
 # colnames(Res100_amr2) <- c("Anc_Prop_EUR","Anc_Prop_AFR", "Anc_Prop_SAS", "Anc_Prop_EAS", "Anc_Prop_IAM")
 # 
 # ### Output of least square errors
 # write.table(Res_100a2[,6],file = "LSEs_100a2_100000.txt", sep = ' ')
 # 
 # ## Creating datasets for plots
 # df100_2 <- reshape(Res100_amr2, times = c("EUR", "AFR","SAS", "EAS","IAM"), timevar = "Anc_Grp",
 #                    varying = list(names(Res100_amr2)), direction = "long") # convert wide to long. 
 # 
 # df100_2$size <- rep("N=100",5000)
 # 
 # 
 # ##### Data Frame for AFstar Minus AFobs
 # Diff <- as.data.frame(cbind(EUR=AF_ref_eur2_100_1000[,1]-refdat2[,1],   AFR=AF_ref_afr2_100_1000[,1]-refdat2[,2],  SAS=AF_ref_sas2_100_1000[,1]-refdat2[,3], 
 #                             EAS=AF_ref_eas2_100_1000[,1] - refdat2[,4],  IAM=AF_ref_iam2_100_1000[,1]- refdat2[,5]))
 # 
 # DF7_2 <- reshape(Diff, times = c("EUR","AFR","SAS","EAS","IAM"), timevar = "Anc_Grp",
 #                  varying = list(names(Diff)), direction = "long") # convert wide to long. 
 # 
 # DF7_2$size <- rep("N=100", 500000)
 # 
 # 
 # 
 # 
 # 
 # 
 # 
 # ##########################################################################
 # # Generate allele counts for population using binom for N = 500
 # ##########################################################################
 # ###################################################################
 # 
 # N2 = 500 # number of people (--> allele number = 2 * N)
 # rep = 1000 # number of replicates 
 # 
 # AF_ref_eur2_500_1000 <- t(sapply(refdat2[,1], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})); write.table(AF_ref_eur2_500_1000,file = "AF_ref_eur2_500_1000.txt", sep = ' ')
 # AF_ref_afr2_500_1000 <- t(sapply(refdat2[,2], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})); write.table(AF_ref_afr2_500_1000,file = "AF_ref_afr2_500_1000.txt", sep = ' ') 
 # AF_ref_sas2_500_1000 <- t(sapply(refdat2[,3], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})); write.table(AF_ref_sas2_500_1000,file = "AF_ref_sas2_500_1000.txt", sep = ' ')
 # AF_ref_eas2_500_1000 <- t(sapply(refdat2[,4], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})); write.table(AF_ref_eas2_500_1000,file = "AF_ref_eas2_500_1000.txt", sep = ' ')
 # AF_ref_iam2_500_1000 <- t(sapply(refdat2[,5], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})); write.table(AF_ref_iam2_500_1000,file = "AF_ref_iam2_500_1000.txt", sep = ' ')
 # 
 # 
 # ## Creating matrices to pass to the Summix function
 # Res_500a2 <- matrix(data = NA, ncol = 8, nrow = rep)
 # 
 # DAT_arr_AMR2 <- array(data=NA,dim = c(samp2,6,rep))
 # 
 # ##
 # for (i in 1:rep){
 #   
 #   DAT_arr_AMR2[ , ,i] = cbind(AF_ref_eur2_500_1000[,i], AF_ref_afr2_500_1000[,i], AF_ref_sas2_500_1000[,i], AF_ref_eas2_500_1000[,i],  AF_ref_iam2_500_1000[,i], refdat2[,6])
 #   
 #   Res_500a2[i,] = ancestr(as.data.frame(DAT_arr_AMR2[,1:5,i]), as.data.frame(DAT_arr_AMR2[,6,i])) 
 # }
 # 
 # 
 # Res500_amr2 <- as.data.frame(cbind(Res_500a2[,1:5]))
 # colnames(Res500_amr2) <- c("Anc_Prop_EUR","Anc_Prop_AFR", "Anc_Prop_SAS", "Anc_Prop_EAS", "Anc_Prop_IAM")
 # 
 # ### Output of least square errors
 # write.table(Res_500a2[,6],file = "LSEs_500a2_100000.txt", sep = ' ')
 # 
 # ## Creating datasets for plots
 # df500_2 <- reshape(Res500_amr2, times = c("EUR", "AFR","SAS", "EAS","IAM"), timevar = "Anc_Grp",
 #                    varying = list(names(Res500_amr2)), direction = "long") # convert wide to long. 
 # 
 # df500_2$size <- rep("N=500",5000)
 # 
 # 
 # ##### Data Frame for AFstar Minus AFobs
 # Diff <- as.data.frame(cbind(EUR=AF_ref_eur2_500_1000[,1]-refdat2[,1],   AFR=AF_ref_afr2_500_1000[,1]-refdat2[,2],  SAS=AF_ref_sas2_500_1000[,1]-refdat2[,3], 
 #                             EAS=AF_ref_eas2_500_1000[,1] - refdat2[,4],  IAM=AF_ref_iam2_500_1000[,1]- refdat2[,5]))
 # 
 # DF8_2 <- reshape(Diff, times = c("EUR","AFR","SAS","EAS","IAM"), timevar = "Anc_Grp",
 #                  varying = list(names(Diff)), direction = "long") # convert wide to long. 
 # 
 # DF8_2$size <- rep("N=500", 500000)
 # 
 # 
 # 
 # 
 # 
 # ##########################################################################
 # # Generate allele counts for population using binom for N = 1000
 # ##########################################################################
 # N2 = 1000 # number of people (--> allele number = 2 * N)
 # rep = 1000 # number of replicates 
 # 
 # AF_ref_eur2 <- t(sapply(refdat2[,1], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)}))
 # AF_ref_afr2 <- t(sapply(refdat2[,2], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})) 
 # AF_ref_sas2 <- t(sapply(refdat2[,3], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)}))
 # AF_ref_eas2 <- t(sapply(refdat2[,4], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)}))
 # AF_ref_iam2 <- t(sapply(refdat2[,5], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)}))
 # 
 # AF_ref_eur2_1000_1000 <- t(sapply(refdat2[,1], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})); write.table(AF_ref_eur2_1000_1000,file = "AF_ref_eur2_1000_1000.txt", sep = ' ')
 # AF_ref_afr2_1000_1000 <- t(sapply(refdat2[,2], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})); write.table(AF_ref_afr2_1000_1000,file = "AF_ref_afr2_1000_1000.txt", sep = ' ') 
 # AF_ref_sas2_1000_1000 <- t(sapply(refdat2[,3], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})); write.table(AF_ref_sas2_1000_1000,file = "AF_ref_sas2_1000_1000.txt", sep = ' ')
 # AF_ref_eas2_1000_1000 <- t(sapply(refdat2[,4], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})); write.table(AF_ref_eas2_1000_1000,file = "AF_ref_eas2_1000_1000.txt", sep = ' ')
 # AF_ref_iam2_1000_1000 <- t(sapply(refdat2[,5], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})); write.table(AF_ref_iam2_1000_1000,file = "AF_ref_iam2_1000_1000.txt", sep = ' ')
 # 
 # 
 # ## Creating matrices to pass to the Summix function
 # Res_1000a2 <- matrix(data = NA, ncol = 8, nrow = rep)
 # 
 # DAT_arr_AMR2 <- array(data=NA,dim = c(samp2,6,rep))
 # 
 # ##
 # for (i in 1:rep){
 #   
 #   DAT_arr_AMR2[ , ,i] = cbind(AF_ref_eur2_1000_1000[,i], AF_ref_afr2_1000_1000[,i], AF_ref_sas2_1000_1000[,i], AF_ref_eas2_1000_1000[,i],  AF_ref_iam2_1000_1000[,i], refdat2[,6])
 #   
 #   Res_1000a2[i,] = ancestr(as.data.frame(DAT_arr_AMR2[,1:5,i]), as.data.frame(DAT_arr_AMR2[,6,i])) 
 # }
 # 
 # 
 # Res1000_amr2 <- as.data.frame(cbind(Res_1000a2[,1:5]))
 # colnames(Res1000_amr2) <- c("Anc_Prop_EUR","Anc_Prop_AFR", "Anc_Prop_SAS", "Anc_Prop_EAS", "Anc_Prop_IAM")
 # 
 # ### Output of least square errors
 # write.table(Res_1000a2[,6],file = "LSEs_1000a2_100000.txt", sep = ' ')
 # 
 # ## Creating datasets for plots
 # df1000_2 <- reshape(Res1000_amr2, times = c("EUR", "AFR","SAS", "EAS","IAM"), timevar = "Anc_Grp",
 #                     varying = list(names(Res1000_amr2)), direction = "long") # convert wide to long. 
 # 
 # df1000_2$size <- rep("N=1000",5000)
 # 
 # 
 # ##### Data Frame for AFstar Minus AFobs
 # Diff <- as.data.frame(cbind(EUR=AF_ref_eur2_1000_1000[,1]-refdat2[,1],   AFR=AF_ref_afr2_1000_1000[,1]-refdat2[,2],  SAS=AF_ref_sas2_1000_1000[,1]-refdat2[,3], 
 #                             EAS=AF_ref_eas2_1000_1000[,1] - refdat2[,4],  IAM=AF_ref_iam2_1000_1000[,1]- refdat2[,5]))
 # 
 # DF9_2 <- reshape(Diff, times = c("EUR","AFR","SAS","EAS","IAM"), timevar = "Anc_Grp",
 #                  varying = list(names(Diff)), direction = "long") # convert wide to long. 
 # 
 # DF9_2$size <- rep("N=500", 500000)
 # 
 # 
 # 
 # ##################################################
 # ####Combined data
 # ##################################################
 # 
 # df.combined2 <- rbind(df5_2,df10_2, df15_2, df20_2, df25_2, df50_2, df100_2, df500_2, df1000_2)
 # #View(df.combined)
 # 
 # DF.combine2 <- rbind(DF1_2, DF2_2, DF3_2, DF4_2, DF5_2, DF6_2, DF7_2, DF8_2, DF9_2)
 # 
 # 
 # 
 # pdf('Rplot_SNP100000_facA3.pdf')
 # # #M <- rep("N = 10",2000) # x variable
 # ggplot(df.combined2, aes(x=size, y=Anc_Prop_EUR, fill=Anc_Grp)) + #x = N, y = ancestry prop, fill = ancestry groups
 #   geom_boxplot() + labs(x = "Allele Number (2 * N)", y = "Ancestry Proportion")+
 #   theme(plot.title = element_text(hjust = 0.5))+
 #   ggtitle("Variation of Summix Ancestry Proportion Estimate (SNP = 100,000) \n American/ Latinx GnomAD")+
 #   scale_x_discrete(limits=c("N=5","N=10","N=15","N=20","N=25","N=50","N=100","N=500","N=1000"))+
 #   facet_wrap(~Anc_Grp, scales = "free") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
 # dev.off()
 # 
 # 
 # 
 # pdf('Rplot_SNP100000_diff1D3.pdf')
 # # #M <- rep("N = 10",2000) # x variable
 # ggplot(DF.combine2, aes(x=size, y=EUR, fill=Anc_Grp)) + #x = N, y = ancestry prop, fill = ancestry groups
 #   geom_boxplot() + labs(x = "Allele Number (2 * N)", y = "Difference in Allele Frequencies")+
 #   theme(plot.title = element_text(hjust = 0.5))+
 #   ggtitle("Distribution of Difference in Allele Frequencies (SNP = 100,000) \n American/ Latinx GnomAD")+
 #   scale_x_discrete(limits=c("N=5","N=10","N=15","N=20","N=25","N=50","N=100","N=500","N=1000"))+
 #   facet_wrap(~Anc_Grp, scales = "free")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
 # dev.off()
 # 
 # 
 # 
 # df.combined2 <- df.combined2 %>% 
 #   mutate(N = ifelse(size=="N=5",5,
 #                     ifelse(size == "N=10", 10,
 #                            ifelse(size == "N=15", 15,
 #                                   ifelse(size == "N=20", 20,
 #                                          ifelse(size == "N=25", 25, 
 #                                                 ifelse(size == "N=50", 50, 
 #                                                        ifelse(size == "N=100", 100, 
 #                                                               ifelse(size == "N=500", 500, 1000)))))))))
 # 
 # dfC <- df.combined2 %>%
 #   group_by(Anc_Grp,N) %>%
 #   summarize(mean = mean(Anc_Prop_EUR),
 #             std = sd(Anc_Prop_EUR),
 #             minimum =min(Anc_Prop_EUR), 
 #             Q1=quantile(Anc_Prop_EUR, .25), 
 #             median = median(Anc_Prop_EUR),
 #             Q3=quantile(Anc_Prop_EUR, .75),
 #             maximum = max(Anc_Prop_EUR), .groups = 'keep')
 # 
 # write.table(dfC, file = "SNPs_100000A3.txt", sep = ' ')
 # 
 # 
 # DF.combine2 <- DF.combine2 %>% 
 #   mutate(N = ifelse(size=="N=5",5,
 #                     ifelse(size == "N=10", 10,
 #                            ifelse(size == "N=15", 15,
 #                                   ifelse(size == "N=20", 20,
 #                                          ifelse(size == "N=25", 25,
 #                                                 ifelse(size == "N=50", 50, 
 #                                                        ifelse(size == "N=100", 100, 
 #                                                               ifelse(size == "N=500", 500, 1000)))))))))
 # 
 # dfC1 <- DF.combine2 %>%
 #   group_by(Anc_Grp,N) %>%
 #   summarize(mean = mean(EUR),
 #             std = sd(EUR),
 #             minimum =min(EUR), 
 #             Q1=quantile(EUR, .25), 
 #             median = median(EUR),
 #             Q3=quantile(EUR, .75),
 #             maximum = max(EUR), .groups = 'keep')
 # 
 # write.table(dfC1, file = "SNPs_100000D3_1.txt", sep = ' ')
 # 
 # 
 # 
