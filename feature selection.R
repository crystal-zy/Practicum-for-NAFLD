library(tidyverse)
library(MatchIt)
library(twang)
library(glmnet)
library(ggplot2)
library(dplyr)
library(lsr)
library(survey)
library(rbounds)
library(naniar)
library(tableone)
library(optmatch)
library(mice)
library(plotmo)
library(corrplot)
library(ggplot2)
library(forestplot)
library(ggcorrplot)
library(MatchThem)
library(gtsummary)
library(rstatix)
library(future.apply)
library(parallel)
setwd("/Users/Crystal/Desktop/Practicum/Data")
NAFLD <-read_csv("data_NAFLD_PS.csv",show_col_types = FALSE)
NAFLD<-NAFLD %>%mutate(Sex=as.factor(ifelse(Sex=="Male", 0, 1)))
NAFLDtotal<-NAFLD[,c(1:27,41)]
NAFLD<-NAFLD[,c(3:4,7,8,10:12,14,27,41)]
variables<-c("Sex","LDL","Triglycerides","Glycated_haemoglobin__HbA1c_",
             "HDL_cholesterol","Waist_circumference","Systolic_Blood_Pressure","BMIdich")
omitNAFLDtotal<-na.omit(NAFLDtotal)
omitNAFLD<- na.omit(NAFLD)

## Create a matrix of features for each model

perflistcreator = function(model_list = fitlist){
  # Extract Features selection for each model
  ranklist = lapply(1:length(model_list), function(x){
    # Extract coefficients of all the features including intercept term
    df = coef(model_list[[x]])
    # Create a data frame
    selfeat = data.frame(variable = rownames(df)[-1], importance=df[-1,1], stringsAsFactors = F)
    names(selfeat) = c("variable", "importance")
    selfeat$importance[is.na(selfeat$importance) | is.nan(selfeat$importance)] = 0
    selfeat$importance[selfeat$importance != 0] = 1
    selfeat
  })
  # Combine the list of selected features into a single dataframe
  res = ranklist
  res = Reduce(function(...) merge(..., by="variable", all=TRUE), res)
  library(naturalsort)
  var_organise=function(inputlist, symbol=":"){
    org_var=lapply(inputlist, function(x) {
      a=strsplit(x, symbol);
      ifelse(length(a[[1]])>1, paste0(naturalsort::naturalsort(a[[1]]),collapse = symbol), a[[1]])})
    unlist(org_var)
    return(unlist(org_var))
  } 
  tres = t(res[,-1])
  #   colnames(tres)=var_organise(inputlist = res[,1], symbol = ":")
  #   ## OUTCOME: c("a.a_b", "a.a_c", "a.a_b")
  #   tres = data.frame(tres)
  #   # remove the ".1" so the column become duplicated
  #   names(tres)<-gsub(".1", "", names(tres), fixed = TRUE)
  #   # get the dataframe of not duplicated
  #   nondups <- tres[,!duplicated(colnames(tres))]
  #   # get the dataframe of duplicated
  #   dups <- tres[,duplicated(colnames(tres))]
  #   name <- intersect(colnames(nondups), colnames(dups))
  #   for (i in 1:length(name)){
  #     nondups[,name[i]]<-ifelse(nondups[,name[i]] == 1 | dups[,name[i]] == 1, 1, 0)
  #   }
  #   res = nondups
  #   return(res)
  
  ## FIX
  tres = data.frame(tres)              
  names(tres) = res[,1]
  tres$boot = 1:nrow(tres)
  df1 = reshape2::melt(data = tres, id.vars = "boot")
  df1 = df1[complete.cases(df1),]
  df1$variable = as.character(df1$variable)
  # str(df1)
  df1$variable = var_organise(inputlist = df1$variable, symbol = ":")
  res = reshape2::dcast(df1,boot~variable, value.var = "value")
  res$boot = NULL
  return(res)             
}

# Perform Boot straps with LASSO contain interaction for NAFLD
# Perform boot straps and feature sampling for LASSO with interaction
boot_rep = 500 # Number of Boot straps
numpred = 8 # Total number of predictors or independent variables used in the model.
# PARALLEL RUN
plan(multisession(workers = 8))
## Prepare the Models
res = lapply(1:boot_rep, function(seed) {
  # Generate Bootstrap rows
  set.seed(seed)
  bootrows = sample(1:nrow(omitNAFLD), nrow(omitNAFLD), replace = T)
  splitdf = omitNAFLD[bootrows,]
  # Generate Feature sample
  featname = names(omitNAFLD)[!names(omitNAFLD) %in% c("NAFLD_status","T2D")]# Get the names of all input features. "y" is used because it is the outcome variable. Hence it should not be selected.
  set.seed(seed)
  featnum = sample(2:numpred, 1, replace = F)
  set.seed(seed)
  selfeat = sample(featname, featnum, replace = F)
  
  # Bootstrapped data
  traindf = splitdf[, c(selfeat, "NAFLD_status")]
  # Run LASSO with no covariate control: Continuous Outcome
  ## Convert dataframe into matrix
  X_cont = model.matrix(NAFLD_status~.*.,traindf[,!names(traindf) %in% "T2D"])[,-1]
  Y_cont = traindf$NAFLD_status
  
  # Find best lambda
  cvfit_cont = glmnet::cv.glmnet(x= X_cont, y = Y_cont, nfolds = 5,family="binomial", alpha = 1, standardize = F)
  lambda = cvfit_cont$lambda.1se
  
  ## Run Lasso
  fit_cont = glmnet::glmnet(x= X_cont, y = Y_cont,nfolds = 5,family="binomial", alpha = 1, standardize = F, lambda = lambda)
  # ## View selected features and coefficients
  # coef_cont = coefextract(model = fit_cont)
  # p <- ggplot(coef_cont[-1,], aes(x = variable, y = coefficient))+
  #   geom_bar(stat = "identity")+
  #   scale_y_continuous(limits = c(-0.5,1))+
  #   theme_bw()
  # p
  # print(p)
})

freqdf = perflistcreator(model_list = res)
nboots<-colSums(freqdf, na.rm = T)
totalboots<-sapply(freqdf, function(x) length(x[!is.na(x)]))
f1<- data.frame(t(rbind(nboots,totalboots)))
## Calculate frequency
freq = colSums(freqdf, na.rm = T)/ sapply(freqdf, function(x) length(x[!is.na(x)]))
freq =ifelse(is.na(freq),0,freq)
freqdf = as.data.frame(t(freqdf))
freqdf = freqdf %>% mutate(freq=c(t(freq)),feature=rownames(freqdf))
print(t(rbind(freqdf$feature,freqdf$freq)))
print(freqdf$feature[freqdf$freq>0.7])



boot_rep = 500 # Number of Boot straps
numpred = 8 # Total number of predictors or independent variables used in the model.
## Prepare the Models
# PARALLEL RUN
plan(multisession(workers = 8))

res2 = lapply(1:boot_rep, function(seed) {
  # Generate Bootstrap rows
  set.seed(seed)
  bootrows = sample(1:nrow(omitNAFLD), nrow(omitNAFLD), replace = T)
  splitdf = omitNAFLD[bootrows,]
  # Generate Feature sample
  featname = names(omitNAFLD)[!names(omitNAFLD) %in% c("NAFLD_status","T2D")]# Get the names of all input features. "y" is used because it is the outcome variable. Hence it should not be selected.
  set.seed(seed)
  featnum = sample(2:numpred, 1, replace = F)
  set.seed(seed)
  selfeat = sample(featname, featnum, replace = F)
  
  # Bootstrapped data
  traindf = splitdf[, c(selfeat, "T2D")]
  # Run LASSO with no covariate control: Continuous Outcome
  ## Convert dataframe into matrix
  X_cont = model.matrix(T2D~.*.,traindf[,!names(traindf) %in% "NAFLD_status"])[,-1]
  Y_cont = traindf$T2D
  
  # Find best lambda
  cvfit_cont = glmnet::cv.glmnet(x= X_cont, y = Y_cont, nfolds = 5,family="binomial", alpha = 1, standardize = F)
  lambda = cvfit_cont$lambda.1se
  
  ## Run Lasso
  fit_cont = glmnet::glmnet(x= X_cont, y = Y_cont,nfolds = 5,family="binomial", alpha = 1, standardize = F, lambda = lambda)
  # ## View selected features and coefficients
  # coef_cont = coefextract(model = fit_cont)
  # p <- ggplot(coef_cont[-1,], aes(x = variable, y = coefficient))+
  #   geom_bar(stat = "identity")+
  #   scale_y_continuous(limits = c(-0.5,1))+
  #   theme_bw()
  # p
  # print(p)
})


## Create a matrix of features for each model

freqdf2 = perflistcreator(model_list = res2)
nboots2<-colSums(freqdf2, na.rm = T)
totalboots2<-sapply(freqdf2, function(x) length(x[!is.na(x)]))
f2<- data.frame(t(rbind(nboots2,totalboots2)))
## Calculate frequency
freq2= colSums(freqdf2, na.rm = T)/ sapply(freqdf2, function(x) length(x[!is.na(x)]))
freq2 =ifelse(is.na(freq2),0,freq2)
freqdf2 = as.data.frame(t(freqdf2))
freqdf2 = freqdf2 %>% mutate(freq=c(t(freq2)),feature=rownames(freqdf2))
print(t(rbind(freqdf2$feature,freqdf2$freq)))
print(freqdf2$feature[freqdf2$freq>0.7])