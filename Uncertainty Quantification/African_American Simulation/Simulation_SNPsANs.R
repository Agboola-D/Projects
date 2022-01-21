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
samp = 1000 # number of SNPs to be extracted

# Select 1K SNPS from reference data 
# Extract required columns 
refdat = data_red %>% 
  sample_n(samp) %>% 
  select(ref_AF_afr_1000G, ref_AF_eur_1000G, gnomad_AF_nfe, gnomad_AF_afr) %>% 
  rename(ref_AFR = ref_AF_afr_1000G, ref_EUR = ref_AF_eur_1000G, obs_EUR = gnomad_AF_nfe, obs_AFR = gnomad_AF_afr)


#######################################################################
##### For N = 5, 10, 15, 20, 25, 50, 100, 500, 1000
#######################################################################

# Generate allele counts for population from binomial distribution for N = 5, P = Allele Frequency
# Repeat process for each N - comments not repeated for other Ns
# AF_ref_afr and AF_ref_eur are reference allele frequency simulated or AF*
# refdat[,1] and refdat[,2] are reference allele frequencies from 1000G or AFobs
# N is fixed for each ancestry group

N = 5 # number of people (--> allele number = 2 * N)
rep = 1000 # number of replicates 
AF_ref_afr <- t(sapply(refdat[,1], function(x){rbinom(rep, 2 * N, x)/ (2 * N)})) 
AF_ref_eur <- t(sapply(refdat[,2], function(x){rbinom(rep, 2 * N, x)/ (2 * N)}))

## Creating storage for outputs
## Array stores multiple matrices - each matrix represents a dataset

Res_5a <- matrix(data = NA, ncol = 5, nrow = rep)
DAT_arr_AFR <- array(data=NA,dim = c(samp,3,rep))

for (i in 1:rep){
  
  DAT_arr_AFR[ , ,i] = cbind(AF_ref_afr[,i],AF_ref_eur[,i],refdat[,4]) # create a matrix of 3 columns - two ref AFs and one obs AF
  # Each matrix is a dataset to be run in summix. So we have 1000 datasets

  Res_5a[i,] = ancestr(as.data.frame(DAT_arr_AFR[,1:2,i]), as.data.frame(DAT_arr_AFR[,3,i])) # run summix on each dataset in the array
  # first entry is a matrix of ref AFs, second entry is a vector of observed AF
  
}

setwd("/newhome/agboolol")

Res5_afr <- as.data.frame(cbind(Res_5a[,1:2])) # extract ancestry proportion estimates for the two ancestry groups: AFR and EUR
colnames(Res5_afr) <- c("Anc_Prop_AFR","Anc_Prop_EUR") # rename columns 

LSEs_5a_1000 <- as.data.frame(cbind(LS = Res_5a[,3], N = rep("N=5", 1000))) # extract least square values for each replicate

df5 <- reshape(Res5_afr, times = c("AFR","EUR"), timevar = "Anc_Grp",
               varying = list(names(Res5_afr)), direction = "long") # convert from wide to long format for easier handling. 
df5$size <- rep("N=5",2000) # create column to specify sample size(s) 



#######################################################################
## Data Frame for AF* - AFobs (difference between simulated
##          and observed AF for each ancestry group)
#######################################################################

Diff <- as.data.frame(cbind(AFR=AF_ref_afr[,1]-refdat[,1],EUR=AF_ref_eur[,1]-refdat[,2])) # create dataframe for difference

DF1 <- reshape(Diff, times = c("AFR","EUR"), timevar = "Anc_Grp",
               varying = list(names(Diff)), direction = "long") # convert from wide to long format for easier handling. 
DF1$size <- rep("N=5", 2000) # create column to specify sample size(s)


##########################################################################
# Repeat process above for N=10 
# Ditto N=15,20, ... , 1000 

N = 10 # number of people (--> allele number = 2 * N)
rep = 1000 # number of replicates 
AF_ref_afr <- t(sapply(refdat[,1], function(x){rbinom(rep, 2 * N, x)/ (2 * N)})) 
AF_ref_eur <- t(sapply(refdat[,2], function(x){rbinom(rep, 2 * N, x)/ (2 * N)}))


Res_10a <- matrix(data = NA, ncol = 5, nrow = rep)
DAT_arr_AFR <- array(data=NA,dim = c(samp,3,rep))

for (i in 1:rep){
  
  DAT_arr_AFR[ , ,i] = cbind(AF_ref_afr[,i],AF_ref_eur[,i],refdat[,4])

  Res_10a[i,] = ancestr(as.data.frame(DAT_arr_AFR[,1:2,i]), as.data.frame(DAT_arr_AFR[,3,i])) # AFR
  
}

Res10_afr <- as.data.frame(cbind(Res_10a[,1:2]))
colnames(Res10_afr) <- c("Anc_Prop_AFR","Anc_Prop_EUR")

LSEs_10a_1000 <- as.data.frame(cbind(LS = Res_10a[,3], N = rep("N=10", 1000))) # extract least square values for each replicate

df10 <- reshape(Res10_afr, times = c("AFR","EUR"), timevar = "Anc_Grp",
                varying = list(names(Res10_afr)), direction = "long") # convert wide to long. 
df10$size <- rep("N=10",2000)



#######################################################################
## Data Frame for AF* - AFobs (difference between simulated
##          and observed AF for each ancestry group)
#######################################################################


Diff <- as.data.frame(cbind(AFR=AF_ref_afr[,1]-refdat[,1],EUR=AF_ref_eur[,1]-refdat[,2]))

DF2 <- reshape(Diff, times = c("AFR","EUR"), timevar = "Anc_Grp",
               varying = list(names(Diff)), direction = "long") # convert wide to long. 
DF2$size <- rep("N=10", 2000)


##########################################################################

N = 15 # number of people (--> allele number = 2 * N)
rep = 1000 # number of replicates 
AF_ref_afr <- t(sapply(refdat[,1], function(x){rbinom(rep, 2 * N, x)/ (2 * N)})) 
AF_ref_eur <- t(sapply(refdat[,2], function(x){rbinom(rep, 2 * N, x)/ (2 * N)}))

##  Creating Storage for Simulation Outputs ## 
Res_15a <- matrix(data = NA, ncol = 5, nrow = rep)
DAT_arr_AFR <- array(data=NA,dim = c(samp,3,rep))

for (i in 1:rep){
  
  DAT_arr_AFR[ , ,i] = cbind(AF_ref_afr[,i],AF_ref_eur[,i],refdat[,4])

  Res_15a[i,] = ancestr(as.data.frame(DAT_arr_AFR[,1:2,i]), as.data.frame(DAT_arr_AFR[,3,i])) # AFR
  
}

Res15_afr <- as.data.frame(cbind(Res_15a[,1:2]))
colnames(Res15_afr) <- c("Anc_Prop_AFR","Anc_Prop_EUR")

LSEs_15a_1000 <- as.data.frame(cbind(LS = Res_15a[,3], N = rep("N=15", 1000))) # extract least square values for each replicate

df15 <- reshape(Res15_afr, times = c("AFR","EUR"), timevar = "Anc_Grp",
                varying = list(names(Res15_afr)), direction = "long") # convert wide to long. 
df15$size <- rep("N=15",2000)



#######################################################################
## Data Frame for AF* - AFobs (difference between simulated
##          and observed AF for each ancestry group)
#######################################################################


Diff <- as.data.frame(cbind(AFR=AF_ref_afr[,1]-refdat[,1],EUR=AF_ref_eur[,1]-refdat[,2]))

DF3 <- reshape(Diff, times = c("AFR","EUR"), timevar = "Anc_Grp",
               varying = list(names(Diff)), direction = "long") # convert wide to long. 
DF3$size <- rep("N=15", 2000)


##########################################################################

N = 20 # number of people (--> allele number = 2 * N)
rep = 1000 # number of replicates 
AF_ref_afr <- t(sapply(refdat[,1], function(x){rbinom(rep, 2 * N, x)/ (2 * N)})) 
AF_ref_eur <- t(sapply(refdat[,2], function(x){rbinom(rep, 2 * N, x)/ (2 * N)}))

##  Creating Storage for Simulation Outputs ## 
Res_20a <- matrix(data = NA, ncol = 5, nrow = rep)
DAT_arr_AFR <- array(data=NA,dim = c(samp,3,rep))

for (i in 1:rep){
  
  DAT_arr_AFR[ , ,i] = cbind(AF_ref_afr[,i],AF_ref_eur[,i],refdat[,4])

  Res_20a[i,] = ancestr(as.data.frame(DAT_arr_AFR[,1:2,i]), as.data.frame(DAT_arr_AFR[,3,i])) # AFR

  
}

Res20_afr <- as.data.frame(cbind(Res_20a[,1:2]))
colnames(Res20_afr) <- c("Anc_Prop_AFR","Anc_Prop_EUR")

LSEs_20a_1000 <- as.data.frame(cbind(LS = Res_20a[,3], N = rep("N=20", 1000))) # extract least square values for each replicate

df20 <- reshape(Res20_afr, times = c("AFR","EUR"), timevar = "Anc_Grp",
                varying = list(names(Res20_afr)), direction = "long") # convert wide to long. 
df20$size <- rep("N=20",2000)


#######################################################################
## Data Frame for AF* - AFobs (difference between simulated
##          and observed AF for each ancestry group)
#######################################################################


Diff <- as.data.frame(cbind(AFR=AF_ref_afr[,1]-refdat[,1],EUR=AF_ref_eur[,1]-refdat[,2]))

DF4 <- reshape(Diff, times = c("AFR","EUR"), timevar = "Anc_Grp",
               varying = list(names(Diff)), direction = "long") # convert wide to long. 
DF4$size <- rep("N=20", 2000)


##########################################################################

N = 25 # number of people (--> allele number = 2 * N)
rep = 1000 # number of replicates 
AF_ref_afr <- t(sapply(refdat[,1], function(x){rbinom(rep, 2 * N, x)/ (2 * N)})) 
AF_ref_eur <- t(sapply(refdat[,2], function(x){rbinom(rep, 2 * N, x)/ (2 * N)}))

##  Creating Storage for Simulation Outputs ## 
Res_25a <- matrix(data = NA, ncol = 5, nrow = rep)
DAT_arr_AFR <- array(data=NA,dim = c(samp,3,rep))

for (i in 1:rep){
  
  DAT_arr_AFR[ , ,i] = cbind(AF_ref_afr[,i],AF_ref_eur[,i],refdat[,4])

  Res_25a[i,] = ancestr(as.data.frame(DAT_arr_AFR[,1:2,i]), as.data.frame(DAT_arr_AFR[,3,i])) # AFR

}

Res25_afr <- as.data.frame(cbind(Res_25a[,1:2]))
colnames(Res25_afr) <- c("Anc_Prop_AFR","Anc_Prop_EUR")

LSEs_25a_1000 <- as.data.frame(cbind(LS = Res_25a[,3], N = rep("N=25", 1000))) # extract least square values for each replicate

df25 <- reshape(Res25_afr, times = c("AFR","EUR"), timevar = "Anc_Grp",
                varying = list(names(Res25_afr)), direction = "long") # convert wide to long. 
df25$size <- rep("N=25",2000)


#######################################################################
## Data Frame for AF* - AFobs (difference between simulated
##          and observed AF for each ancestry group)
#######################################################################


Diff <- as.data.frame(cbind(AFR=AF_ref_afr[,1]-refdat[,1],EUR=AF_ref_eur[,1]-refdat[,2]))

DF5 <- reshape(Diff, times = c("AFR","EUR"), timevar = "Anc_Grp",
               varying = list(names(Diff)), direction = "long") # convert wide to long. 
DF5$size <- rep("N=25", 2000)

##########################################################################

N = 50 # number of people (--> allele number = 2 * N)
rep = 1000 # number of replicates 
AF_ref_afr <- t(sapply(refdat[,1], function(x){rbinom(rep, 2 * N, x)/ (2 * N)})) 
AF_ref_eur <- t(sapply(refdat[,2], function(x){rbinom(rep, 2 * N, x)/ (2 * N)}))

##  Creating Storage for Simulation Outputs ## 
Res_50a <- matrix(data = NA, ncol = 5, nrow = rep)
DAT_arr_AFR <- array(data=NA,dim = c(samp,3,rep))

for (i in 1:rep){
  
  DAT_arr_AFR[ , ,i] = cbind(AF_ref_afr[,i],AF_ref_eur[,i],refdat[,4])

  Res_50a[i,] = ancestr(as.data.frame(DAT_arr_AFR[,1:2,i]), as.data.frame(DAT_arr_AFR[,3,i])) # AFR

}

Res50_afr <- as.data.frame(cbind(Res_50a[,1:2]))
colnames(Res50_afr) <- c("Anc_Prop_AFR","Anc_Prop_EUR")

LSEs_50a_1000 <- as.data.frame(cbind(LS = Res_50a[,3], N = rep("N=50", 1000))) # extract least square values for each replicate

df50 <- reshape(Res50_afr, times = c("AFR","EUR"), timevar = "Anc_Grp",
              varying = list(names(Res50_afr)), direction = "long") # convert wide to long. 
df50$size <- rep("N=50",2000)


#######################################################################
## Data Frame for AF* - AFobs (difference between simulated
##          and observed AF for each ancestry group)
#######################################################################

Diff <- as.data.frame(cbind(AFR=AF_ref_afr[,1]-refdat[,1],EUR=AF_ref_eur[,1]-refdat[,2]))

DF6 <- reshape(Diff, times = c("AFR","EUR"), timevar = "Anc_Grp",
               varying = list(names(Diff)), direction = "long") # convert wide to long. 
DF6$size <- rep("N=50", 2000)


##########################################################################

N = 100 # number of people 
rep = 1000 # number of replicates 
AF_ref_afr <- t(sapply(refdat[,1], function(x){rbinom(rep, 2 * N, x)/ (2 * N)})) # refdat[,1] is p in binomial dist
AF_ref_eur <- t(sapply(refdat[,2], function(x){rbinom(rep, 2 * N, x)/ (2 * N)}))

##  Creating Storage for Simulation Outputs ## 
Res_100a <- matrix(data = NA, ncol = 5, nrow = rep)
DAT_arr_AFR <- array(data=NA,dim = c(samp,3,rep))

for (i in 1:rep){
  
  DAT_arr_AFR[ , ,i] = cbind(AF_ref_afr[,i],AF_ref_eur[,i],refdat[,4])

  Res_100a[i,] = ancestr(as.data.frame(DAT_arr_AFR[,1:2,i]), as.data.frame(DAT_arr_AFR[,3,i])) # AFR

}

Res100_afr <- as.data.frame(cbind(Res_100a[,1:2]))
colnames(Res100_afr) <- c("Anc_Prop_AFR","Anc_Prop_EUR")

LSEs_100a_1000 <- as.data.frame(cbind(LS = Res_100a[,3], N = rep("N=100", 1000))) # extract least square values for each replicate

df100 <- reshape(Res100_afr, times = c("AFR","EUR"), timevar = "Anc_Grp",
                varying = list(names(Res100_afr)), direction = "long") # convert wide to long. 
df100$size <- rep("N=100", 2000)


#######################################################################
## Data Frame for AF* - AFobs (difference between simulated
##          and observed AF for each ancestry group)
#######################################################################


Diff <- as.data.frame(cbind(AFR=AF_ref_afr[,1]-refdat[,1],EUR=AF_ref_eur[,1]-refdat[,2]))

DF7 <- reshape(Diff, times = c("AFR","EUR"), timevar = "Anc_Grp",
               varying = list(names(Diff)), direction = "long") # convert wide to long. 
DF7$size <- rep("N=100", 2000)


###################################################################

N = 500 # number of people
rep = 1000 # number of replicates 
AF_ref_afr <- t(sapply(refdat[,1], function(x){rbinom(rep, 2 * N, x)/ (2 * N)})) # refdat[,1] is p 
AF_ref_eur <- t(sapply(refdat[,2], function(x){rbinom(rep, 2 * N, x)/ (2 * N)}))

##  Creating Storage for Simulation Outputs ## 
Res_500a <- matrix(data = NA, ncol = 5, nrow = rep)
DAT_arr_AFR <- array(data=NA,dim = c(samp,3,rep))

for (i in 1:rep){
  
  DAT_arr_AFR[ , ,i] = cbind(AF_ref_afr[,i],AF_ref_eur[,i],refdat[,4])

  Res_500a[i,] = ancestr(as.data.frame(DAT_arr_AFR[,1:2,i]), as.data.frame(DAT_arr_AFR[,3,i])) # AFR

}

Res500_afr <- as.data.frame(cbind(Res_500a[,1:2]))
colnames(Res500_afr) <- c("Anc_Prop_AFR","Anc_Prop_EUR")

LSEs_500a_1000 <- as.data.frame(cbind(LS = Res_500a[,3], N = rep("N=500", 1000))) # extract least square values for each replicate

df500 <- reshape(Res500_afr, times = c("AFR","EUR"), timevar = "Anc_Grp",
                 varying = list(names(Res500_afr)), direction = "long") # convert wide to long. 
df500$size <- rep("N=500",2000)


#######################################################################
## Data Frame for AF* - AFobs (difference between simulated
##          and observed AF for each ancestry group)
#######################################################################

Diff <- as.data.frame(cbind(AFR=AF_ref_afr[,1]-refdat[,1],EUR=AF_ref_eur[,1]-refdat[,2]))

DF8 <- reshape(Diff, times = c("AFR","EUR"), timevar = "Anc_Grp",
               varying = list(names(Diff)), direction = "long") # convert wide to long. 
DF8$size <- rep("N=500", 2000)



####################################################################

N = 1000 # number of people
rep = 1000 # number of replicates 
AF_ref_afr <- t(sapply(refdat[,1], function(x){rbinom(rep, 2 * N, x)/ (2 * N)})) # refdat[,1] is p 
AF_ref_eur <- t(sapply(refdat[,2], function(x){rbinom(rep, 2 * N, x)/ (2 * N)}))

##  Creating Storage for Simulation Outputs ## 
Res_1000a <- matrix(data = NA, ncol = 5, nrow = rep)
DAT_arr_AFR <- array(data=NA,dim = c(samp,3,rep))

for (i in 1:rep){
  
  DAT_arr_AFR[ , ,i] = cbind(AF_ref_afr[,i],AF_ref_eur[,i],refdat[,4])

  Res_1000a[i,] = ancestr(as.data.frame(DAT_arr_AFR[,1:2,i]), as.data.frame(DAT_arr_AFR[,3,i])) # AFR

}

Res1000_afr <- as.data.frame(cbind(Res_1000a[,1:2]))
colnames(Res1000_afr) <-  c("Anc_Prop_AFR","Anc_Prop_EUR")

LSEs_1000a_1000 <- as.data.frame(cbind(LS = Res_1000a[,3], N = rep("N=1000", 1000))) # extract least square values for each replicate

df1000 <- reshape(Res1000_afr, times = c("AFR","EUR"), timevar = "Anc_Grp",
                 varying = list(names(Res1000_afr)), direction = "long") # convert wide to long. 
df1000$size <- rep("N=1000",2000)


#######################################################################
## Data Frame for AF* - AFobs (difference between simulated
##          and observed AF for each ancestry group)
#######################################################################

Diff <- as.data.frame(cbind(AFR=AF_ref_afr[,1]-refdat[,1],EUR=AF_ref_eur[,1]-refdat[,2]))

DF9 <- reshape(Diff, times = c("AFR","EUR"), timevar = "Anc_Grp",
               varying = list(names(Diff)), direction = "long") # convert wide to long. 
DF9$size <- rep("N=1000", 2000)


###################################################
## Combine all data frames (outputs) for each N 
## Both ancestry proportions and difference in AFs
###################################################

df.combined <- rbind(df5,df10, df15, df20, df25, df50, df100, df500, df1000) # Combining all ancestry proportions estimates for each N

DF.combine <- rbind(DF1, DF2, DF3, DF4, DF5, DF6, DF7, DF8, DF9) # Combining difference in AFs for each N

LS.combine <- rbind(LSEs_5a_1000, LSEs_10a_1000, LSEs_15a_1000, LSEs_20a_1000, LSEs_25a_1000, LSEs_50a_1000,
                    LSEs_100a_1000, LSEs_500a_1000, LSEs_1000a_1000) # combining all least square values for each N

## Create plot for ancestry proportion estimates for each N

pdf('Rplot_SNP1000_fac.pdf')
ggplot(df.combined, aes(x=size, y=Anc_Prop_AFR, color=Anc_Grp)) + #x = N, y = ancestry prop, fill = ancestry groups
  geom_boxplot() + labs(x = "Allele Number (2 * N)", y = "Ancestry Proportion")+
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90, vjust = 0.5))+
  ggtitle("Variation of Summix Ancestry Proportion Estimate (SNP = 1,000) \n African/African American GnomAD")+
  scale_x_discrete(limits=c("N=5","N=10","N=15","N=20","N=25","N=50","N=100","N=500","N=1000"))+
  facet_wrap(~Anc_Grp, scales = "free")
dev.off()


## Create plot for difference in simulated and observed AFs for each N

pdf('Rplot_SNP1000_diff.pdf')
ggplot(DF.combine, aes(x=size, y=AFR, color=Anc_Grp)) + #x = N, y = ancestry prop, fill = ancestry groups
  geom_boxplot() + labs(x = "Allele Number (2 * N)", y = "Difference in Allele Frequencies")+
  theme(plot.title = element_text(hjust = 0.5),  axis.text.x = element_text(angle = 90, vjust = 0.5))+
  ggtitle("Distribution of Difference in Allele Frequencies (SNP = 1,000) \n African/African American GnomAD")+
  scale_x_discrete(limits=c("N=5","N=10","N=15","N=20","N=25","N=50","N=100","N=500","N=1000"))+
  facet_wrap(~Anc_Grp, scales = "free")
dev.off()

## Create plot for least square values for each N

pdf('Rplot_SNP1000_LS.pdf')
ggplot(LS.combine, aes(x=N, y=as.numeric(LS))) + #x = N, y = ancestry prop, fill = ancestry groups
  geom_boxplot() + labs(x = "Allele Number (2 * N)", y = "Least Square Values")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Distribution of Least Square Values (SNP = 1,000) \n African/African American GnomAD")+
  scale_x_discrete(limits=c("N=5","N=10","N=15","N=20","N=25","N=50","N=100","N=500","N=1000"))
dev.off()


## Replace 'N=5' with 5 (and so on) so that results can be grouped in order of magnitude of N 

df.combined <- df.combined %>% 
  mutate(N = ifelse(size=="N=5",5,
                    ifelse(size == "N=10", 10,
                           ifelse(size == "N=15", 15,
                                  ifelse(size == "N=20", 20,
                                         ifelse(size == "N=25", 25,
                                                ifelse(size == "N=50", 50,
                                                       ifelse(size == "N=100", 100,
                                                              ifelse(size == "N=500", 500, 1000)))))))))

## Generate summary statistics on ancestry proportion estimates for each N

dfA <- df.combined %>%
  group_by(Anc_Grp,N) %>%
  summarize(mean = mean(Anc_Prop_AFR),
            std = sd(Anc_Prop_AFR),
            minimum =min(Anc_Prop_AFR), 
            Q1=quantile(Anc_Prop_AFR, .25), 
            median = median(Anc_Prop_AFR),
            Q3=quantile(Anc_Prop_AFR, .75),
            maximum = max(Anc_Prop_AFR), .groups = 'keep')

write.table(dfA, file = "SNPs_1000_Anc.txt", sep = ' ') # export as a text file


## Replace 'N=5' with 5 (and so on) so that results can be grouped in order of magnitude of N 

DF.combine <- DF.combine %>% 
  mutate(N = ifelse(size=="N=5",5,
                    ifelse(size == "N=10", 10,
                           ifelse(size == "N=15", 15,
                                  ifelse(size == "N=20", 20,
                                         ifelse(size == "N=25", 25,
                                                ifelse(size == "N=50", 50,
                                                       ifelse(size == "N=100", 100,
                                                              ifelse(size == "N=500", 500, 1000)))))))))

## Generate summary statistics on AF* - AFobs for each N

dfA1 <- DF.combine %>%
  group_by(Anc_Grp,N) %>%
  summarize(mean = mean(AFR),
            std = sd(AFR),
            minimum =min(AFR), 
            Q1=quantile(AFR, .25), 
            median = median(AFR),
            Q3=quantile(AFR, .75),
            maximum = max(AFR), .groups = 'keep')

write.table(dfA1, file = "SNPs_1000_Dif.txt", sep = ' ') # export as a text file

## Create plot for variability in simulated and observed AFs for each N and N^2

pdf('Rplot_SNP1000_Vardiff_N.pdf')
ggplot(dfA1, aes(x=as.numeric(N), y=as.numeric(std^2), color=Anc_Grp)) + 
  geom_point() + labs(x = "Number of People (N)", y = "Variance")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Variability of Difference in Allele Frequencies (SNP = 1,000) \n African/African American GnomAD")+
  facet_wrap(~Anc_Grp, scales = "free")
dev.off()

pdf('Rplot_SNP1000_Vardiff_Nsq.pdf')
ggplot(dfA1, aes(x=as.numeric(N^2), y=as.numeric(std^2), color=Anc_Grp)) + 
  geom_point() + labs(x = TeX("$N^2$"), y = "Variance")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Variability of Difference in Allele Frequencies (SNP = 1,000) \n African/African American GnomAD")+
  facet_wrap(~Anc_Grp, scales = "free")
dev.off()

pdf('Rplot_SNP1000_Vardiff_Ninv.pdf')
ggplot(dfA1, aes(x=as.numeric(1/N), y=as.numeric(std^2), color=Anc_Grp)) + 
  geom_point() + labs(x = TeX("$\\frac{1}{N}$"), y = "Variance")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Variability of Difference in Allele Frequencies (SNP = 1,000) \n African/African American GnomAD")+
  facet_wrap(~Anc_Grp, scales = "free")
dev.off()

pdf('Rplot_SNP1000_Vardiff_Nsqinv.pdf')
ggplot(dfA1, aes(x=as.numeric(1/(N^2)), y=as.numeric(std^2), color=Anc_Grp)) + 
  geom_point() + labs(x = TeX("$\\frac{1}{N^2}$"), y = "Variance")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Variability of Difference in Allele Frequencies (SNP = 1,000) \n African/African American GnomAD")+
  facet_wrap(~Anc_Grp, scales = "free")
dev.off()

######################
### SNP = 10,000 ####
#####################


## SIMULATION ##
# Set Seed
set.seed(122356)
samp1 = 10000

# Select 10K SNPS from reference data 
refdat1 = data_red %>% 
  sample_n(samp1) %>% 
  select(ref_AF_afr_1000G, ref_AF_eur_1000G, gnomad_AF_nfe, gnomad_AF_afr) %>% 
  rename(ref_AFR = ref_AF_afr_1000G, ref_EUR = ref_AF_eur_1000G, obs_EUR = gnomad_AF_nfe, obs_AFR = gnomad_AF_afr)


#######################################################################
##### For N = 5, 10, 15, 20, 25, 10, 50, 100, 500, 1000
#######################################################################


# Generate allele counts for population using binom for N = 5 
N1 = 5 # number of people (--> allele number = 2 * N)
rep = 1000 # number of replicates 
AF_ref_afr1 <- t(sapply(refdat1[,1], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})) # AF_ref_afr = AF_ref_afr*, refdat[,1] is allele freq observed
AF_ref_eur1 <- t(sapply(refdat1[,2], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)}))

##  Creating Storage for Simulation Outputs ## 
Res_5a1 <- matrix(data = NA, ncol = 5, nrow = rep) #ncol is the same as the output from ancestr prop function
DAT_arr_AFR1 <- array(data=NA,dim = c(samp1,3,rep))

for (i in 1:rep){
  
  DAT_arr_AFR1[ , ,i] = cbind(AF_ref_afr1[,i],AF_ref_eur1[,i],refdat1[,4])

  Res_5a1[i,] = ancestr(as.data.frame(DAT_arr_AFR1[,1:2,i]), as.data.frame(DAT_arr_AFR1[,3,i])) # AFR
  
}

Res5_afr1 <- as.data.frame(cbind(Res_5a1[,1:2]))
colnames(Res5_afr1) <- c("Anc_Prop_AFR","Anc_Prop_EUR") ### 

LSEs_5a1_10000 <- as.data.frame(cbind(LS = Res_5a1[,3], N = rep("N=5", 10000))) # extract least square values for each replicate

df5_1 <- reshape(Res5_afr1, times = c("AFR","EUR"), timevar = "Anc_Grp",
                 varying = list(names(Res5_afr1)), direction = "long") # convert wide to long. 
df5_1$size <- rep("N=5", 2000)


#######################################################################
## Data Frame for AF* - AFobs (difference between simulated
##          and observed AF for each ancestry group)
#######################################################################



Diff <- as.data.frame(cbind(AFR=AF_ref_afr1[,1]-refdat1[,1],EUR=AF_ref_eur1[,1]-refdat1[,2]))

DF10 <- reshape(Diff, times = c("AFR","EUR"), timevar = "Anc_Grp",
               varying = list(names(Diff)), direction = "long") # convert wide to long. 
DF10$size <- rep("N=5", 20000)

##########################################################################


# Generate allele counts for population using binom for N = 10 
N1 = 10 # number of people (--> allele number = 2 * N)
rep = 1000 # number of replicates 
AF_ref_afr1 <- t(sapply(refdat1[,1], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})) # AF_ref_afr = AF_ref_afr*, refdat[,1] is allele freq observed
AF_ref_eur1 <- t(sapply(refdat1[,2], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)}))

##  Creating Storage for Simulation Outputs ## 
Res_10a1 <- matrix(data = NA, ncol = 5, nrow = rep) #ncol is the same as the output from ancestr prop function
DAT_arr_AFR1 <- array(data=NA,dim = c(samp1,3,rep))

for (i in 1:rep){
  
  DAT_arr_AFR1[ , ,i] = cbind(AF_ref_afr1[,i],AF_ref_eur1[,i],refdat1[,4])

  Res_10a1[i,] = ancestr(as.data.frame(DAT_arr_AFR1[,1:2,i]), as.data.frame(DAT_arr_AFR1[,3,i])) # AFR
  
}

Res10_afr1 <- as.data.frame(cbind(Res_10a1[,1:2]))
colnames(Res10_afr1) <- c("Anc_Prop_AFR","Anc_Prop_EUR") ### 

LSEs_10a1_10000 <- as.data.frame(cbind(LS = Res_10a1[,3], N = rep("N=10", 10000))) # extract least square values for each replicate

df10_1 <- reshape(Res10_afr1, times = c("AFR","EUR"), timevar = "Anc_Grp",
                  varying = list(names(Res10_afr1)), direction = "long") # convert wide to long. 
df10_1$size <- rep("N=10", 2000)


#######################################################################
## Data Frame for AF* - AFobs (difference between simulated
##          and observed AF for each ancestry group)
#######################################################################



Diff <- as.data.frame(cbind(AFR=AF_ref_afr1[,1]-refdat1[,1],EUR=AF_ref_eur1[,1]-refdat1[,2]))

DF11 <- reshape(Diff, times = c("AFR","EUR"), timevar = "Anc_Grp",
               varying = list(names(Diff)), direction = "long") # convert wide to long. 
DF11$size <- rep("N=10", 20000)



##########################################################################


# Generate allele counts for population using binom for N = 15 
N1 = 15 # number of people (--> allele number = 2 * N)
rep = 1000 # number of replicates 
AF_ref_afr1 <- t(sapply(refdat1[,1], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})) # AF_ref_afr = AF_ref_afr*, refdat[,1] is allele freq observed
AF_ref_eur1 <- t(sapply(refdat1[,2], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)}))

##  Creating Storage for Simulation Outputs ## 
Res_15a1 <- matrix(data = NA, ncol = 5, nrow = rep) #ncol is the same as the output from ancestr prop function
DAT_arr_AFR1 <- array(data=NA,dim = c(samp1,3,rep))

for (i in 1:rep){
  
  DAT_arr_AFR1[ , ,i] = cbind(AF_ref_afr1[,i],AF_ref_eur1[,i],refdat1[,4])

  Res_15a1[i,] = ancestr(as.data.frame(DAT_arr_AFR1[,1:2,i]), as.data.frame(DAT_arr_AFR1[,3,i])) # AFR

}

Res15_afr1 <- as.data.frame(cbind(Res_15a1[,1:2]))
colnames(Res15_afr1) <- c("Anc_Prop_AFR","Anc_Prop_EUR") ### 

LSEs_15a1_10000 <- as.data.frame(cbind(LS = Res_15a1[,3], N = rep("N=15", 10000))) # extract least square values for each replicate

df15_1 <- reshape(Res15_afr1, times = c("AFR","EUR"), timevar = "Anc_Grp",
                  varying = list(names(Res15_afr1)), direction = "long") # convert wide to long. 
df15_1$size <- rep("N=15", 2000)


#######################################################################
## Data Frame for AF* - AFobs (difference between simulated
##          and observed AF for each ancestry group)
#######################################################################


Diff <- as.data.frame(cbind(AFR=AF_ref_afr1[,1]-refdat1[,1],EUR=AF_ref_eur1[,1]-refdat1[,2]))

DF12 <- reshape(Diff, times = c("AFR","EUR"), timevar = "Anc_Grp",
               varying = list(names(Diff)), direction = "long") # convert wide to long. 
DF12$size <- rep("N=15", 20000)



##########################################################################


# Generate allele counts for population using binom for N = 20 
N1 = 20 # number of people (--> allele number = 2 * N)
rep = 1000 # number of replicates 
AF_ref_afr1 <- t(sapply(refdat1[,1], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})) # AF_ref_afr = AF_ref_afr*, refdat[,1] is allele freq observed
AF_ref_eur1 <- t(sapply(refdat1[,2], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)}))

##  Creating Storage for Simulation Outputs ## 
Res_20a1 <- matrix(data = NA, ncol = 5, nrow = rep) #ncol is the same as the output from ancestr prop function
Res_20e1 <- matrix(data = NA, ncol = 5, nrow = rep)
DAT_arr_AFR1 <- array(data=NA,dim = c(samp1,3,rep))

for (i in 1:rep){
  
  DAT_arr_AFR1[ , ,i] = cbind(AF_ref_afr1[,i],AF_ref_eur1[,i],refdat1[,4])

  Res_20a1[i,] = ancestr(as.data.frame(DAT_arr_AFR1[,1:2,i]), as.data.frame(DAT_arr_AFR1[,3,i])) # AFR

}

Res20_afr1 <- as.data.frame(cbind(Res_20a1[,1:2]))
colnames(Res20_afr1) <- c("Anc_Prop_AFR","Anc_Prop_EUR") ### 

LSEs_20a1_10000 <- as.data.frame(cbind(LS = Res_20a1[,3], N = rep("N=20", 10000))) # extract least square values for each replicate

df20_1 <- reshape(Res20_afr1, times = c("AFR","EUR"), timevar = "Anc_Grp",
                  varying = list(names(Res20_afr1)), direction = "long") # convert wide to long. 
df20_1$size <- rep("N=20", 2000)


#######################################################################
## Data Frame for AF* - AFobs (difference between simulated
##          and observed AF for each ancestry group)
#######################################################################


Diff <- as.data.frame(cbind(AFR=AF_ref_afr1[,1]-refdat1[,1],EUR=AF_ref_eur1[,1]-refdat1[,2]))

DF13 <- reshape(Diff, times = c("AFR","EUR"), timevar = "Anc_Grp",
               varying = list(names(Diff)), direction = "long") # convert wide to long. 
DF13$size <- rep("N=20", 20000)


##########################################################################


# Generate allele counts for population using binom for N = 25 
N1 = 25 # number of people (--> allele number = 2 * N)
rep = 1000 # number of replicates 
AF_ref_afr1 <- t(sapply(refdat1[,1], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})) # AF_ref_afr = AF_ref_afr*, refdat[,1] is allele freq observed
AF_ref_eur1 <- t(sapply(refdat1[,2], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)}))

##  Creating Storage for Simulation Outputs ## 
Res_25a1 <- matrix(data = NA, ncol = 5, nrow = rep) #ncol is the same as the output from ancestr prop function
DAT_arr_AFR1 <- array(data=NA,dim = c(samp1,3,rep))

for (i in 1:rep){
  
  DAT_arr_AFR1[ , ,i] = cbind(AF_ref_afr1[,i],AF_ref_eur1[,i],refdat1[,4])

  Res_25a1[i,] = ancestr(as.data.frame(DAT_arr_AFR1[,1:2,i]), as.data.frame(DAT_arr_AFR1[,3,i])) # AFR

}

Res25_afr1 <- as.data.frame(cbind(Res_25a1[,1:2]))
colnames(Res25_afr1) <- c("Anc_Prop_AFR","Anc_Prop_EUR") ### 

LSEs_25a1_10000 <- as.data.frame(cbind(LS = Res_25a1[,3], N = rep("N=25", 10000))) # extract least square values for each replicate

df25_1 <- reshape(Res25_afr1, times = c("AFR","EUR"), timevar = "Anc_Grp",
                  varying = list(names(Res25_afr1)), direction = "long") # convert wide to long. 
df25_1$size <- rep("N=25", 2000)


#######################################################################
## Data Frame for AF* - AFobs (difference between simulated
##          and observed AF for each ancestry group)
#######################################################################


Diff <- as.data.frame(cbind(AFR=AF_ref_afr1[,1]-refdat1[,1],EUR=AF_ref_eur1[,1]-refdat1[,2]))

DF14 <- reshape(Diff, times = c("AFR","EUR"), timevar = "Anc_Grp",
                varying = list(names(Diff)), direction = "long") # convert wide to long. 
DF14$size <- rep("N=25", 20000)



##########################################################################


# Generate allele counts for population using binom for N = 50
N1 = 50 # number of people (--> allele number = 2 * N)
rep = 1000 # number of replicates 
AF_ref_afr1 <- t(sapply(refdat1[,1], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})) # AF_ref_afr = AF_ref_afr*, refdat[,1] is allele freq observed 
AF_ref_eur1 <- t(sapply(refdat1[,2], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)}))

##  Creating Storage for Simulation Outputs ## 
Res_50a1 <- matrix(data = NA, ncol = 5, nrow = rep)
DAT_arr_AFR1 <- array(data=NA,dim = c(samp1,3,rep))

for (i in 1:rep){
  
  DAT_arr_AFR1[ , ,i] = cbind(AF_ref_afr1[,i],AF_ref_eur1[,i],refdat1[,4])

  Res_50a1[i,] = ancestr(as.data.frame(DAT_arr_AFR1[,1:2,i]), as.data.frame(DAT_arr_AFR1[,3,i])) # AFR

}

Res50_afr1 <- as.data.frame(cbind(Res_50a1[,1:2]))
colnames(Res50_afr1) <- c("Anc_Prop_AFR","Anc_Prop_EUR")

LSEs_50a1_10000 <- as.data.frame(cbind(LS = Res_50a1[,3], N = rep("N=50", 10000))) # extract least square values for each replicate

df50_1 <- reshape(Res50_afr1, times = c("AFR","EUR"), timevar = "Anc_Grp",
                varying = list(names(Res50_afr1)), direction = "long") # convert wide to long. 
df50_1$size <- rep("N=50",2000)


#######################################################################
## Data Frame for AF* - AFobs (difference between simulated
##          and observed AF for each ancestry group)
#######################################################################


Diff <- as.data.frame(cbind(AFR=AF_ref_afr1[,1]-refdat1[,1],EUR=AF_ref_eur1[,1]-refdat1[,2]))

DF15 <- reshape(Diff, times = c("AFR","EUR"), timevar = "Anc_Grp",
               varying = list(names(Diff)), direction = "long") # convert wide to long. 
DF15$size <- rep("N=50", 20000)


##########################################################################

# Generate allele counts for population using binom for N = 100

N1 = 100 # number of people 
rep = 1000 # number of replicates 
AF_ref_afr1 <- t(sapply(refdat1[,1], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})) # refdat[,1] is p 
AF_ref_eur1 <- t(sapply(refdat1[,2], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)}))

##  Creating Storage for Simulation Outputs ## 
Res_100a1 <- matrix(data = NA, ncol = 5, nrow = rep)
DAT_arr_AFR1 <- array(data=NA,dim = c(samp1,3,rep))

for (i in 1:rep){
  
  DAT_arr_AFR1[ , ,i] = cbind(AF_ref_afr1[,i],AF_ref_eur1[,i],refdat1[,4])

  Res_100a1[i,] = ancestr(as.data.frame(DAT_arr_AFR1[,1:2,i]), as.data.frame(DAT_arr_AFR1[,3,i])) # AFR
  
}

Res100_afr1 <- as.data.frame(cbind(Res_100a1[,1:2]))
colnames(Res100_afr1) <- c("Anc_Prop_AFR","Anc_Prop_EUR")

LSEs_100a1_10000 <- as.data.frame(cbind(LS = Res_100a1[,3], N = rep("N=100", 10000))) # extract least square values for each replicate

df100_1 <- reshape(Res100_afr1, times = c("AFR","EUR"), timevar = "Anc_Grp",
                 varying = list(names(Res100_afr1)), direction = "long") # convert wide to long. 
df100_1$size <- rep("N=100", 2000)


#######################################################################
## Data Frame for AF* - AFobs (difference between simulated
##          and observed AF for each ancestry group)
#######################################################################


Diff <- as.data.frame(cbind(AFR=AF_ref_afr1[,1]-refdat1[,1],EUR=AF_ref_eur1[,1]-refdat1[,2]))

DF16 <- reshape(Diff, times = c("AFR","EUR"), timevar = "Anc_Grp",
               varying = list(names(Diff)), direction = "long") # convert wide to long. 
DF16$size <- rep("N=100", 20000)


###################################################################

N1 = 500 # number of people
rep = 1000 # number of replicates 
AF_ref_afr1 <- t(sapply(refdat1[,1], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})) # refdat[,1] is p 
AF_ref_eur1 <- t(sapply(refdat1[,2], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)}))

##  Creating Storage for Simulation Outputs ## 
Res_500a1 <- matrix(data = NA, ncol = 5, nrow = rep)
DAT_arr_AFR1 <- array(data=NA,dim = c(samp1,3,rep))

for (i in 1:rep){
  
  DAT_arr_AFR1[ , ,i] = cbind(AF_ref_afr1[,i],AF_ref_eur1[,i],refdat1[,4])
  
  Res_500a1[i,] = ancestr(as.data.frame(DAT_arr_AFR1[,1:2,i]), as.data.frame(DAT_arr_AFR1[,3,i])) # AFR
  
}

Res500_afr1 <- as.data.frame(cbind(Res_500a1[,1:2]))
colnames(Res500_afr1) <- c("Anc_Prop_AFR","Anc_Prop_EUR")

LSEs_500a1_10000 <- as.data.frame(cbind(LS = Res_500a1[,3], N = rep("N=500", 10000))) # extract least square values for each replicate

df500_1 <- reshape(Res500_afr1, times = c("AFR","EUR"), timevar = "Anc_Grp",
                 varying = list(names(Res500_afr1)), direction = "long") # convert wide to long. 
df500_1$size <- rep("N=500",2000)



#######################################################################
## Data Frame for AF* - AFobs (difference between simulated
##          and observed AF for each ancestry group)
#######################################################################


Diff <- as.data.frame(cbind(AFR=AF_ref_afr1[,1]-refdat1[,1],EUR=AF_ref_eur1[,1]-refdat1[,2]))

DF17 <- reshape(Diff, times = c("AFR","EUR"), timevar = "Anc_Grp",
               varying = list(names(Diff)), direction = "long") # convert wide to long. 
DF17$size <- rep("N=500", 20000)


####################################################################

N1 = 1000 # number of people
rep = 1000 # number of replicates 
AF_ref_afr1 <- t(sapply(refdat1[,1], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)})) # refdat[,1] is p 
AF_ref_eur1 <- t(sapply(refdat1[,2], function(x){rbinom(rep, 2 * N1, x)/ (2 * N1)}))

##  Creating Storage for Simulation Outputs ## 
Res_1000a1 <- matrix(data = NA, ncol = 5, nrow = rep)
DAT_arr_AFR1 <- array(data=NA,dim = c(samp1,3,rep))


for (i in 1:rep){
  
  DAT_arr_AFR1[ , ,i] = cbind(AF_ref_afr1[,i],AF_ref_eur1[,i],refdat1[,4])
  
  Res_1000a1[i,] = ancestr(as.data.frame(DAT_arr_AFR1[,1:2,i]), as.data.frame(DAT_arr_AFR1[,3,i])) # AFR
  
}

Res1000_afr1 <- as.data.frame(cbind(Res_1000a1[,1:2]))
colnames(Res1000_afr1) <-  c("Anc_Prop_AFR","Anc_Prop_EUR")

LSEs_1000a1_10000 <- as.data.frame(cbind(LS = Res_1000a1[,3], N = rep("N=1000", 10000))) # extract least square values for each replicate

df1000_1 <- reshape(Res1000_afr1, times = c("AFR","EUR"), timevar = "Anc_Grp",
                  varying = list(names(Res1000_afr1)), direction = "long") # convert wide to long. 
df1000_1$size <- rep("N=1000",2000)



#######################################################################
## Data Frame for AF* - AFobs (difference between simulated
##          and observed AF for each ancestry group)
#######################################################################


Diff <- as.data.frame(cbind(AFR=AF_ref_afr1[,1]-refdat1[,1],EUR=AF_ref_eur1[,1]-refdat1[,2]))

DF18 <- reshape(Diff, times = c("AFR","EUR"), timevar = "Anc_Grp",
               varying = list(names(Diff)), direction = "long") # convert wide to long. 
DF18$size <- rep("N=1000", 20000)



##################################################
####Combined data
##################################################

df.combined1 <- rbind(df5_1,df10_1, df15_1, df20_1, df25_1, df50_1, df100_1, df500_1, df1000_1)

DF.combine1 <- rbind(DF10, DF11, DF12, DF13, DF14, DF15, DF16, DF17, DF18)

LS.combine1 <- rbind(LSEs_5a1_10000, LSEs_10a1_10000, LSEs_15a1_10000, LSEs_20a1_10000, LSEs_25a1_10000, LSEs_50a1_10000,
                     LSEs_100a1_10000, LSEs_500a1_10000, LSEs_1000a1_10000)
  
pdf('Rplot_SNP10000_fac.pdf')
ggplot(df.combined1, aes(x=size, y=Anc_Prop_AFR, color=Anc_Grp)) + #x = N, y = ancestry prop, fill = ancestry groups
  geom_boxplot() + labs(x = "Allele Number (2 * N)", y = "Ancestry Proportion")+
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90, vjust = 0.5))+
  ggtitle("Variation of Summix Ancestry Proportion Estimate (SNP = 10,000) \n African/African American GnomAD")+
  scale_x_discrete(limits=c("N=5","N=10","N=15","N=20","N=25","N=50","N=100","N=500","N=1000"))+
  facet_wrap(~Anc_Grp, scales = "free")
dev.off()


pdf('Rplot_SNP10000_diff.pdf')
ggplot(DF.combine1, aes(x=size, y=AFR, color=Anc_Grp)) + #x = N, y = ancestry prop, fill = ancestry groups
  geom_boxplot() + labs(x = "Allele Number (2 * N)", y = "Difference in Allele Frequencies")+
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90, vjust = 0.5))+
  ggtitle("Distribution of Difference in Allele Frequencies (SNP = 10,000) \n African/African American GnomAD")+
  scale_x_discrete(limits=c("N=5","N=10","N=15","N=20","N=25","N=50","N=100","N=500","N=1000"))+
  facet_wrap(~Anc_Grp, scales = "free")
dev.off()


## Create plot for least square values for each N

pdf('Rplot_SNP10000_LS.pdf')
ggplot(LS.combine1, aes(x=N, y=as.numeric(LS))) + #x = N, y = ancestry prop, fill = ancestry groups
  geom_boxplot() + labs(x = "Allele Number (2 * N)", y = "Least Square Values")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Distribution of Least Square Values (SNP = 10,000) \n African/African American GnomAD")+
  scale_x_discrete(limits=c("N=5","N=10","N=15","N=20","N=25","N=50","N=100","N=500","N=1000"))
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
  summarize(mean = mean(Anc_Prop_AFR),
            std = sd(Anc_Prop_AFR),
            minimum =min(Anc_Prop_AFR), 
            Q1=quantile(Anc_Prop_AFR, .25), 
            median = median(Anc_Prop_AFR),
            Q3=quantile(Anc_Prop_AFR, .75),
            maximum = max(Anc_Prop_AFR), .groups = 'keep')

write.table(dfB, file = "SNPs_10000_Anc.txt", sep = ' ')


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
  summarize(mean = mean(AFR),
            std = sd(AFR),
            minimum =min(AFR), 
            Q1=quantile(AFR, .25), 
            median = median(AFR),
            Q3=quantile(AFR, .75),
            maximum = max(AFR), .groups = 'keep')

write.table(dfB1, file = "SNPs_10000_Dif.txt", sep = ' ')


## Create plot for variability in simulated and observed AFs for each N and N^2

pdf('Rplot_SNP10000_Vardiff_N.pdf')
ggplot(dfB1, aes(x=as.numeric(N), y=as.numeric(std^2), color=Anc_Grp)) + #x = N, y = ancestry prop, fill = ancestry groups
  geom_point() + labs(x = "Number of People (N)", y = "Variance")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Variability of Difference in Allele Frequencies (SNP = 10,000) \n African/African American GnomAD")+
  facet_wrap(~Anc_Grp, scales = "free")
dev.off()


pdf('Rplot_SNP10000_Vardiff_Nsq.pdf')
ggplot(dfB1, aes(x=as.numeric(N^2), y=as.numeric(std^2), color=Anc_Grp)) +
  geom_point() + labs(x = TeX("$N^2$"), y = "Variance")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Variability of Difference in Allele Frequencies (SNP = 10,000) \n African/African American GnomAD")+
  facet_wrap(~Anc_Grp, scales = "free")
dev.off()

pdf('Rplot_SNP10000_Vardiff_Ninv.pdf')
ggplot(dfB1, aes(x=as.numeric(1/N), y=as.numeric(std^2), color=Anc_Grp)) + 
  geom_point() + labs(x = TeX("$\\frac{1}{N}$"), y = "Variance")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Variability of Difference in Allele Frequencies (SNP = 10,000) \n African/African American GnomAD")+
  facet_wrap(~Anc_Grp, scales = "free")
dev.off()

pdf('Rplot_SNP10000_Vardiff_Nsqinv.pdf')
ggplot(dfB1, aes(x=as.numeric(1/(N^2)), y=as.numeric(std^2), color=Anc_Grp)) + 
  geom_point() + labs(x = TeX("$\\frac{1}{N^2}$"), y = "Variance")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Variability of Difference in Allele Frequencies (SNP = 10,000) \n African/African American GnomAD")+
  facet_wrap(~Anc_Grp, scales = "free")
dev.off()

######################
### SNP = 100,000 ####
#####################


## SIMULATION ##
# Set Seed
set.seed(2656256)
samp2 = 100000

# Select 100K SNPS from reference data 
refdat2 = data_red %>% 
  sample_n(samp2) %>% 
  select(ref_AF_afr_1000G, ref_AF_eur_1000G, gnomad_AF_nfe, gnomad_AF_afr) %>% 
  rename(ref_AFR = ref_AF_afr_1000G, ref_EUR = ref_AF_eur_1000G, obs_EUR = gnomad_AF_nfe, obs_AFR = gnomad_AF_afr)



#######################################################################
##### For N = 5, 10, 15, 20, 25, 10, 50, 100, 500, 1000
#######################################################################

# Generate allele counts for population using binom for N = 5
N2 = 5 # number of people (--> allele number = 2 * N)
rep = 1000 # number of replicates 
AF_ref_afr2 <- t(sapply(refdat2[,1], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})) # AF_ref_afr = AF_ref_afr*, refdat[,1] is allele freq observed
AF_ref_eur2 <- t(sapply(refdat2[,2], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)}))

##  Creating Storage for Simulation Outputs ## 
Res_5a2 <- matrix(data = NA, ncol = 5, nrow = rep) #ncol is the same as the output from ancestr prop function
DAT_arr_AFR2 <- array(data=NA,dim = c(samp2,3,rep))

for (i in 1:rep){
  
  DAT_arr_AFR2[ , ,i] = cbind(AF_ref_afr2[,i],AF_ref_eur2[,i],refdat2[,4])
  
  Res_5a2[i,] = ancestr(as.data.frame(DAT_arr_AFR2[,1:2,i]), as.data.frame(DAT_arr_AFR2[,3,i])) # AFR
  
}

Res5_afr2 <- as.data.frame(cbind(Res_5a2[,1:2]))
colnames(Res5_afr2) <- c("Anc_Prop_AFR","Anc_Prop_EUR") ### 

LSEs_5a2_100000 <- as.data.frame(cbind(LS = Res_5a2[,3], N = rep("N=5", 100000))) # extract least square values for each replicate

df5_2 <- reshape(Res5_afr2, times = c("AFR","EUR"), timevar = "Anc_Grp",
                 varying = list(names(Res5_afr2)), direction = "long") # convert wide to long. 
df5_2$size <- rep("N=5", 2000)


#######################################################################
## Data Frame for AF* - AFobs (difference between simulated
##          and observed AF for each ancestry group)
#######################################################################


Diff <- as.data.frame(cbind(AFR=AF_ref_afr2[,1]-refdat2[,1],EUR=AF_ref_eur2[,1]-refdat2[,2]))

DF19 <- reshape(Diff, times = c("AFR","EUR"), timevar = "Anc_Grp",
                varying = list(names(Diff)), direction = "long") # convert wide to long. 
DF19$size <- rep("N=5", 200000)



###########################################################################################

# Generate allele counts for population using binom for N = 10 
N2 = 10 # number of people (--> allele number = 2 * N)
rep = 1000 # number of replicates 
AF_ref_afr2 <- t(sapply(refdat2[,1], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})) # AF_ref_afr = AF_ref_afr*, refdat[,1] is allele freq observed
AF_ref_eur2 <- t(sapply(refdat2[,2], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)}))

##  Creating Storage for Simulation Outputs ## 
Res_10a2 <- matrix(data = NA, ncol = 5, nrow = rep) #ncol is the same as the output from ancestr prop function
DAT_arr_AFR2 <- array(data=NA,dim = c(samp2,3,rep))

for (i in 1:rep){
  
  DAT_arr_AFR2[ , ,i] = cbind(AF_ref_afr2[,i],AF_ref_eur2[,i],refdat2[,4])
  
  Res_10a2[i,] = ancestr(as.data.frame(DAT_arr_AFR2[,1:2,i]), as.data.frame(DAT_arr_AFR2[,3,i])) # AFR
  
}

Res10_afr2 <- as.data.frame(cbind(Res_10a2[,1:2]))
colnames(Res10_afr2) <- c("Anc_Prop_AFR","Anc_Prop_EUR") ### 

LSEs_10a2_100000 <- as.data.frame(cbind(LS = Res_10a2[,3], N = rep("N=10", 100000))) # extract least square values for each replicate

df10_2 <- reshape(Res10_afr2, times = c("AFR","EUR"), timevar = "Anc_Grp",
                  varying = list(names(Res10_afr2)), direction = "long") # convert wide to long. 
df10_2$size <- rep("N=10", 2000)


#######################################################################
## Data Frame for AF* - AFobs (difference between simulated
##          and observed AF for each ancestry group)
#######################################################################


Diff <- as.data.frame(cbind(AFR=AF_ref_afr2[,1]-refdat2[,1],EUR=AF_ref_eur2[,1]-refdat2[,2]))

DF20 <- reshape(Diff, times = c("AFR","EUR"), timevar = "Anc_Grp",
                varying = list(names(Diff)), direction = "long") # convert wide to long. 
DF20$size <- rep("N=10", 200000)

###########################################################################################

# Generate allele counts for population using binom for N = 15 
N2 = 15 # number of people (--> allele number = 2 * N)
rep = 1000 # number of replicates 
AF_ref_afr2 <- t(sapply(refdat2[,1], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})) # AF_ref_afr = AF_ref_afr*, refdat[,1] is allele freq observed
AF_ref_eur2 <- t(sapply(refdat2[,2], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)}))

##  Creating Storage for Simulation Outputs ## 
Res_15a2 <- matrix(data = NA, ncol = 5, nrow = rep) #ncol is the same as the output from ancestr prop function
DAT_arr_AFR2 <- array(data=NA,dim = c(samp2,3,rep))

for (i in 1:rep){
  
  DAT_arr_AFR2[ , ,i] = cbind(AF_ref_afr2[,i],AF_ref_eur2[,i],refdat2[,4])
  
  Res_15a2[i,] = ancestr(as.data.frame(DAT_arr_AFR2[,1:2,i]), as.data.frame(DAT_arr_AFR2[,3,i])) # AFR
  
}

Res15_afr2 <- as.data.frame(cbind(Res_15a2[,1:2]))
colnames(Res15_afr2) <- c("Anc_Prop_AFR","Anc_Prop_EUR") ### 

LSEs_15a2_100000 <- as.data.frame(cbind(LS = Res_15a2[,3], N = rep("N=15", 100000))) # extract least square values for each replicate

df15_2 <- reshape(Res15_afr2, times = c("AFR","EUR"), timevar = "Anc_Grp",
                  varying = list(names(Res15_afr2)), direction = "long") # convert wide to long. 
df15_2$size <- rep("N=15", 2000)


#######################################################################
## Data Frame for AF* - AFobs (difference between simulated
##          and observed AF for each ancestry group)
#######################################################################


Diff <- as.data.frame(cbind(AFR=AF_ref_afr2[,1]-refdat2[,1],EUR=AF_ref_eur2[,1]-refdat2[,2]))

DF21 <- reshape(Diff, times = c("AFR","EUR"), timevar = "Anc_Grp",
                varying = list(names(Diff)), direction = "long") # convert wide to long. 
DF21$size <- rep("N=15", 200000)



###########################################################################################

# Generate allele counts for population using binom for N = 20 
N2 = 20 # number of people (--> allele number = 2 * N)
rep = 1000 # number of replicates 
AF_ref_afr2 <- t(sapply(refdat2[,1], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})) # AF_ref_afr = AF_ref_afr*, refdat[,1] is allele freq observed
AF_ref_eur2 <- t(sapply(refdat2[,2], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)}))

##  Creating Storage for Simulation Outputs ## 
Res_20a2 <- matrix(data = NA, ncol = 5, nrow = rep) #ncol is the same as the output from ancestr prop function
DAT_arr_AFR2 <- array(data=NA,dim = c(samp2,3,rep))

for (i in 1:rep){
  
  DAT_arr_AFR2[ , ,i] = cbind(AF_ref_afr2[,i],AF_ref_eur2[,i],refdat2[,4])
  
  Res_20a2[i,] = ancestr(as.data.frame(DAT_arr_AFR2[,1:2,i]), as.data.frame(DAT_arr_AFR2[,3,i])) # AFR
  
}

Res20_afr2 <- as.data.frame(cbind(Res_20a2[,1:2]))
colnames(Res20_afr2) <- c("Anc_Prop_AFR","Anc_Prop_EUR") ### 

LSEs_20a2_100000 <- as.data.frame(cbind(LS = Res_20a2[,3], N = rep("N=20", 100000))) # extract least square values for each replicate

df20_2 <- reshape(Res20_afr2, times = c("AFR","EUR"), timevar = "Anc_Grp",
                  varying = list(names(Res20_afr2)), direction = "long") # convert wide to long. 
df20_2$size <- rep("N=20", 2000)


#######################################################################
## Data Frame for AF* - AFobs (difference between simulated
##          and observed AF for each ancestry group)
#######################################################################


Diff <- as.data.frame(cbind(AFR=AF_ref_afr2[,1]-refdat2[,1],EUR=AF_ref_eur2[,1]-refdat2[,2]))

DF22 <- reshape(Diff, times = c("AFR","EUR"), timevar = "Anc_Grp",
                varying = list(names(Diff)), direction = "long") # convert wide to long. 
DF22$size <- rep("N=20", 200000)


###########################################################################################

# Generate allele counts for population using binom for N = 25
N2 = 25 # number of people (--> allele number = 2 * N)
rep = 1000 # number of replicates 
AF_ref_afr2 <- t(sapply(refdat2[,1], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})) # AF_ref_afr = AF_ref_afr*, refdat[,1] is allele freq observed
AF_ref_eur2 <- t(sapply(refdat2[,2], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)}))

##  Creating Storage for Simulation Outputs ## 
Res_25a2 <- matrix(data = NA, ncol = 5, nrow = rep) #ncol is the same as the output from ancestr prop function
DAT_arr_AFR2 <- array(data=NA,dim = c(samp2,3,rep))

for (i in 1:rep){
  
  DAT_arr_AFR2[ , ,i] = cbind(AF_ref_afr2[,i],AF_ref_eur2[,i],refdat2[,4])
  
  Res_25a2[i,] = ancestr(as.data.frame(DAT_arr_AFR2[,1:2,i]), as.data.frame(DAT_arr_AFR2[,3,i])) # AFR

  
}

Res25_afr2 <- as.data.frame(cbind(Res_25a2[,1:2]))
colnames(Res25_afr2) <- c("Anc_Prop_AFR","Anc_Prop_EUR") ### 

LSEs_25a2_100000 <- as.data.frame(cbind(LS = Res_25a2[,3], N = rep("N=25", 100000))) # extract least square values for each replicate

df25_2 <- reshape(Res25_afr2, times = c("AFR","EUR"), timevar = "Anc_Grp",
                  varying = list(names(Res25_afr2)), direction = "long") # convert wide to long. 
df25_2$size <- rep("N=25", 2000)


#######################################################################
## Data Frame for AF* - AFobs (difference between simulated
##          and observed AF for each ancestry group)
#######################################################################


Diff <- as.data.frame(cbind(AFR=AF_ref_afr2[,1]-refdat2[,1],EUR=AF_ref_eur2[,1]-refdat2[,2]))

DF23 <- reshape(Diff, times = c("AFR","EUR"), timevar = "Anc_Grp",
                varying = list(names(Diff)), direction = "long") # convert wide to long. 
DF23$size <- rep("N=25", 200000)




###########################################################################################


# Generate allele counts for population using binom for N = 50
N2 = 50 # number of people (--> allele number = 2 * N)
rep = 1000 # number of replicates 
AF_ref_afr2 <- t(sapply(refdat2[,1], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})) # AF_ref_afr = AF_ref_afr*, refdat[,1] is allele freq observed 
AF_ref_eur2 <- t(sapply(refdat2[,2], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)}))

##  Creating Storage for Simulation Outputs ## 
Res_50a2 <- matrix(data = NA, ncol = 5, nrow = rep)
DAT_arr_AFR2 <- array(data=NA,dim = c(samp2,3,rep))

for (i in 1:rep){
  
  DAT_arr_AFR2[ , ,i] = cbind(AF_ref_afr2[,i],AF_ref_eur2[,i],refdat2[,4])
  
  Res_50a2[i,] = ancestr(as.data.frame(DAT_arr_AFR2[,1:2,i]), as.data.frame(DAT_arr_AFR2[,3,i])) # AFR
  
}

Res50_afr2 <- as.data.frame(cbind(Res_50a2[,1:2]))
colnames(Res50_afr2) <- c("Anc_Prop_AFR","Anc_Prop_EUR")

LSEs_50a2_100000 <- as.data.frame(cbind(LS = Res_50a2[,3], N = rep("N=50", 100000))) # extract least square values for each replicate

df50_2 <- reshape(Res50_afr2, times = c("AFR","EUR"), timevar = "Anc_Grp",
                  varying = list(names(Res50_afr2)), direction = "long") # convert wide to long. 
df50_2$size <- rep("N=50",2000)



#######################################################################
## Data Frame for AF* - AFobs (difference between simulated
##          and observed AF for each ancestry group)
#######################################################################


Diff <- as.data.frame(cbind(AFR=AF_ref_afr2[,1]-refdat2[,1],EUR=AF_ref_eur2[,1]-refdat2[,2]))

DF24 <- reshape(Diff, times = c("AFR","EUR"), timevar = "Anc_Grp",
                varying = list(names(Diff)), direction = "long") # convert wide to long. 
DF24$size <- rep("N=50", 200000)


##########################################################################

# Generate allele counts for population using binom for N = 100

N2 = 100 # number of people 
rep = 1000 # number of replicates 
AF_ref_afr2 <- t(sapply(refdat2[,1], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})) # refdat[,1] is p 
AF_ref_eur2 <- t(sapply(refdat2[,2], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)}))

##  Creating Storage for Simulation Outputs ## 
Res_100a2 <- matrix(data = NA, ncol = 5, nrow = rep)
DAT_arr_AFR2 <- array(data=NA,dim = c(samp2,3,rep))

for (i in 1:rep){
  
  DAT_arr_AFR2[ , ,i] = cbind(AF_ref_afr2[,i],AF_ref_eur2[,i],refdat2[,4])
  
  Res_100a2[i,] = ancestr(as.data.frame(DAT_arr_AFR2[,1:2,i]), as.data.frame(DAT_arr_AFR2[,3,i])) # AFR

  
}

Res100_afr2 <- as.data.frame(cbind(Res_100a2[,1:2]))
colnames(Res100_afr2) <- c("Anc_Prop_AFR","Anc_Prop_EUR")

LSEs_100a2_100000 <- as.data.frame(cbind(LS = Res_100a2[,3], N = rep("N=100", 100000))) # extract least square values for each replicate

df100_2 <- reshape(Res100_afr2, times = c("AFR","EUR"), timevar = "Anc_Grp",
                   varying = list(names(Res100_afr2)), direction = "long") # convert wide to long. 
df100_2$size <- rep("N=100", 2000)



#######################################################################
## Data Frame for AF* - AFobs (difference between simulated
##          and observed AF for each ancestry group)
#######################################################################


Diff <- as.data.frame(cbind(AFR=AF_ref_afr2[,1]-refdat2[,1],EUR=AF_ref_eur2[,1]-refdat2[,2]))

DF25 <- reshape(Diff, times = c("AFR","EUR"), timevar = "Anc_Grp",
                varying = list(names(Diff)), direction = "long") # convert wide to long. 
DF25$size <- rep("N=100", 200000)



###################################################################

N2 = 500 # number of people
rep = 1000 # number of replicates 
AF_ref_afr2 <- t(sapply(refdat2[,1], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})) # refdat[,1] is p 
AF_ref_eur2 <- t(sapply(refdat2[,2], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)}))

##  Creating Storage for Simulation Outputs ## 
Res_500a2 <- matrix(data = NA, ncol = 5, nrow = rep)
DAT_arr_AFR2 <- array(data=NA,dim = c(samp2,3,rep))

for (i in 1:rep){
  
  DAT_arr_AFR2[ , ,i] = cbind(AF_ref_afr2[,i],AF_ref_eur2[,i],refdat2[,4])
  
  Res_500a2[i,] = ancestr(as.data.frame(DAT_arr_AFR2[,1:2,i]), as.data.frame(DAT_arr_AFR2[,3,i])) # AFR
  
}

Res500_afr2 <- as.data.frame(cbind(Res_500a2[,1:2]))
colnames(Res500_afr2) <- c("Anc_Prop_AFR","Anc_Prop_EUR")

LSEs_500a2_100000 <- as.data.frame(cbind(LS = Res_500a2[,3], N = rep("N=500", 100000))) # extract least square values for each replicate

df500_2 <- reshape(Res500_afr2, times = c("AFR","EUR"), timevar = "Anc_Grp",
                   varying = list(names(Res500_afr2)), direction = "long") # convert wide to long. 
df500_2$size <- rep("N=500",2000)



#######################################################################
## Data Frame for AF* - AFobs (difference between simulated
##          and observed AF for each ancestry group)
#######################################################################


Diff <- as.data.frame(cbind(AFR=AF_ref_afr2[,1]-refdat2[,1],EUR=AF_ref_eur2[,1]-refdat2[,2]))

DF26 <- reshape(Diff, times = c("AFR","EUR"), timevar = "Anc_Grp",
                varying = list(names(Diff)), direction = "long") # convert wide to long. 
DF26$size <- rep("N=500", 200000)


####################################################################

N2 = 1000 # number of people
rep = 1000 # number of replicates 
AF_ref_afr2 <- t(sapply(refdat2[,1], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)})) # refdat[,1] is p 
AF_ref_eur2 <- t(sapply(refdat2[,2], function(x){rbinom(rep, 2 * N2, x)/ (2 * N2)}))

##  Creating Storage for Simulation Outputs ## 
Res_1000a2 <- matrix(data = NA, ncol = 5, nrow = rep)
DAT_arr_AFR2 <- array(data=NA,dim = c(samp2,3,rep))

for (i in 1:rep){
  
  DAT_arr_AFR2[ , ,i] = cbind(AF_ref_afr2[,i],AF_ref_eur2[,i],refdat2[,4])
  
  Res_1000a2[i,] = ancestr(as.data.frame(DAT_arr_AFR2[,1:2,i]), as.data.frame(DAT_arr_AFR2[,3,i])) # AFR
  
}

Res1000_afr2 <- as.data.frame(cbind(Res_1000a2[,1:2]))
colnames(Res1000_afr2) <-  c("Anc_Prop_AFR","Anc_Prop_EUR")

LSEs_1000a2_100000 <- as.data.frame(cbind(LS = Res_1000a2[,3], N = rep("N=1000", 100000))) # extract least square values for each replicate

df1000_2 <- reshape(Res1000_afr2, times = c("AFR","EUR"), timevar = "Anc_Grp",
                    varying = list(names(Res1000_afr2)), direction = "long") # convert wide to long. 
df1000_2$size <- rep("N=1000",2000)




#######################################################################
## Data Frame for AF* - AFobs (difference between simulated
##          and observed AF for each ancestry group)
#######################################################################


Diff <- as.data.frame(cbind(AFR=AF_ref_afr2[,1]-refdat2[,1],EUR=AF_ref_eur2[,1]-refdat2[,2]))

DF27 <- reshape(Diff, times = c("AFR","EUR"), timevar = "Anc_Grp",
                varying = list(names(Diff)), direction = "long") # convert wide to long. 
DF27$size <- rep("N=1000", 200000)




##################################################
####Combined data
##################################################

df.combined2 <- rbind(df5_2,df10_2, df15_2, df20_2, df25_2, df50_2, df100_2, df500_2, df1000_2)

DF.combine2 <- rbind(DF19, DF20, DF21, DF22, DF23, DF24, DF25, DF26, DF27)

LS.combine2 <- rbind(LSEs_5a2_100000, LSEs_10a2_100000, LSEs_15a2_100000, LSEs_20a2_100000, LSEs_25a2_100000,LSEs_50a2_100000,
                     LSEs_100a2_100000, LSEs_500a2_100000, LSEs_1000a2_100000)

pdf('Rplot_SNP100000_fac.pdf')
ggplot(df.combined2, aes(x=size, y=Anc_Prop_AFR, color=Anc_Grp)) + #x = N, y = ancestry prop, fill = ancestry groups
   geom_boxplot() + labs(x = "Allele Number (2 * N)", y = "Ancestry Proportion")+
   theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90, vjust = 0.5))+
   ggtitle("Variation of Summix Ancestry Proportion Estimate (SNP = 100,000) \n African/African American GnomAD")+
   scale_x_discrete(limits=c("N=5","N=10","N=15","N=20","N=25","N=50","N=100","N=500","N=1000"))+
   facet_wrap(~Anc_Grp, scales = "free")
dev.off()


pdf('Rplot_SNP100000_diff.pdf')
ggplot(DF.combine2, aes(x=size, y=AFR, color=Anc_Grp)) + #x = N, y = ancestry prop, fill = ancestry groups
  geom_boxplot() + labs(x = "Allele Number (2 * N)", y = "Difference in Allele Frequencies")+
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90, vjust = 0.5))+
  ggtitle("Distribution of Difference in Allele Frequencies (SNP = 100,000) \n African/African American GnomAD")+
  scale_x_discrete(limits=c("N=5","N=10","N=15","N=20","N=25","N=50","N=100","N=500","N=1000"))+
  facet_wrap(~Anc_Grp, scales = "free")
dev.off()

## Create plot for least square values for each N

pdf('Rplot_SNP100000_LS.pdf')
ggplot(LS.combine2, aes(x=N, y=as.numeric(LS))) + #x = N, y = ancestry prop, fill = ancestry groups
  geom_boxplot() + labs(x = "Allele Number (2 * N)", y = "Least Square Values")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Distribution of Least Square Values (SNP = 100,000) \n African/African American GnomAD")+
  scale_x_discrete(limits=c("N=5","N=10","N=15","N=20","N=25","N=50","N=100","N=500","N=1000"))
dev.off()


df.combined2 <- df.combined2 %>% 
  mutate(N = ifelse(size=="N=5",5,
                    ifelse(size == "N=10", 10,
                           ifelse(size == "N=15", 15,
                                  ifelse(size == "N=20", 20,
                                         ifelse(size == "N=25", 25,
                                                ifelse(size == "N=50", 50,
                                                       ifelse(size == "N=100", 100,
                                                              ifelse(size == "N=500", 500, 1000)))))))))

dfC <- df.combined2 %>%
  group_by(Anc_Grp,N) %>%
  summarize(mean = mean(Anc_Prop_AFR),
            std = sd(Anc_Prop_AFR),
            minimum =min(Anc_Prop_AFR), 
            Q1=quantile(Anc_Prop_AFR, .25), 
            median = median(Anc_Prop_AFR),
            Q3=quantile(Anc_Prop_AFR, .75),
            maximum = max(Anc_Prop_AFR), .groups = 'keep')

write.table(dfC, file = "SNPs_100000_Anc.txt", sep = ' ')


DF.combine2 <- DF.combine2 %>% 
  mutate(N = ifelse(size=="N=5",5,
                    ifelse(size == "N=10", 10,
                           ifelse(size == "N=15", 15,
                                  ifelse(size == "N=20", 20,
                                         ifelse(size == "N=25", 25,
                                                ifelse(size == "N=50", 50,
                                                       ifelse(size == "N=100", 100,
                                                              ifelse(size == "N=500", 500, 1000)))))))))

dfC1 <- DF.combine2 %>%
  group_by(Anc_Grp,N) %>%
  summarize(mean = mean(AFR),
            std = sd(AFR),
            minimum =min(AFR), 
            Q1=quantile(AFR, .25), 
            median = median(AFR),
            Q3=quantile(AFR, .75),
            maximum = max(AFR), .groups = 'keep')

write.table(dfC1, file = "SNPs_100000_Diff.txt", sep = ' ')

## Create plot for variability in simulated and observed AFs for each N and N^2

pdf('Rplot_SNP100000_Vardiff_N.pdf')
ggplot(dfC1, aes(x=as.numeric(N), y=as.numeric(std^2), color=Anc_Grp)) + 
  geom_point() + labs(x = "Number of People (N)", y = "Variance")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Variability of Difference in Allele Frequencies (SNP = 100,000) \n African/African American GnomAD")+
  facet_wrap(~Anc_Grp, scales = "free")
dev.off()

pdf('Rplot_SNP100000_Vardiff_Nsq.pdf')
ggplot(dfC1, aes(x=as.numeric(N^2), y=as.numeric(std^2), color=Anc_Grp)) +
  geom_point() + labs(x = TeX("$N^2$"), y = "Variance")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Variability of Difference in Allele Frequencies (SNP = 100,000) \n African/African American GnomAD")+
  facet_wrap(~Anc_Grp, scales = "free")
dev.off()

pdf('Rplot_SNP100000_Vardiff_Ninv.pdf')
ggplot(dfC1, aes(x=as.numeric(1/N), y=as.numeric(std^2), color=Anc_Grp)) + 
  geom_point() + labs(x = TeX("$\\frac{1}{N}$"), y = "Variance")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Variability of Difference in Allele Frequencies (SNP = 100,000) \n African/African American GnomAD")+
  facet_wrap(~Anc_Grp, scales = "free")
dev.off()

pdf('Rplot_SNP100000_Vardiff_Nsqinv.pdf')
ggplot(dfC1, aes(x=as.numeric(1/(N^2)), y=as.numeric(std^2), color=Anc_Grp)) + 
  geom_point() + labs(x = TeX("$\\frac{1}{N^2}$"), y = "Variance")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Variability of Difference in Allele Frequencies (SNP = 100,000) \n African/African American GnomAD")+
  facet_wrap(~Anc_Grp, scales = "free")
dev.off()
