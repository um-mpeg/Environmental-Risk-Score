###Sample code of construction of the Environmental Risk Score using Adaptive Elastic-NET, Xin Wang, University of Michigan###

###metal denotes all the (log-10 transformed) metals in the training set##  
###Y denotes the continuous dependent variable,that is log-10 waist circumference (WC) in this study, in the training set.#  
###covs denotes the covariates we need to force in the selected model, in the training set.# 

#n denotes number of participants in the training set
n<- nrow(metal)

#Standardize all (log10-transformed) metal concentrations
metal<-scale(metal, center = TRUE, scale = TRUE)

#Create square terms of heavy metals 
square_metal<-data.frame(matrix(0,n,ncol(metal)))

for(i in 1:ncol(metal)) {
  square_metal[, i] <- metal[, i]^2
  colnames(square_metal)[i] <- paste(colnames(metal)[i], "^2", sep = "")
}


#Create pair-wise interaction 
pairs_metal <- data.frame(matrix(0, n, choose(ncol(metal),2)))
k = 1
for(i in 1:(ncol(metal)-1)) {
  for(j in (i+1):ncol(metal)) {
    pairs_metal[, k] <- metal[, i]*metal[, j]
    colnames(pairs_metal)[k] <- paste(colnames(metal)[i], ".", colnames(metal)[j], sep = "")
    k = k + 1
  }
}

#Combine covariates and metal variables
X<-cbind(covs,metal,square_metal,pairs_metal)



############################################################## 
######################elastic net (ENET)###################### 
##############################################################

library(gcdnet)

set.seed(0)
foldid_mat <- cbind(sample(n), c(rep(1:5, each = n / 5))) #for cross-validation
weight_enet<-c(rep(0,ncol(covs)),rep(1,ncol(X)-ncol(covs))) #zero is assigned to weights of covariates, to force them in the model selection

#lambda2 is chosen based on the minimum cv error
lambda2_enet <- seq(0.0001, 0.01, by=0.0001) 
min_cv_enet_lambda2=numeric(100)
for(i in 1:length(lambda2_enet)) {
  cv_enet_lambda2 <- cv.gcdnet(X, Y, lambda2 = lambda2_enet[i], foldid = foldid_mat[, 2], method = "ls",pf=weight_enet) 
  min_cv_enet_lambda2[i] <- min(cv_enet_lambda2$cvm)
  print(i)
}
lambda2_enet_opt<-lambda2_enet[which.min(min_cv_enet_lambda2)]

##lambda 1
cv_enet_lambda1 <- cv.gcdnet(X, Y, foldid = foldid_mat[, 2], lambda2 = lambda2_enet_opt,method = "ls",pf=weight_enet)  
lambda1_enet_opt <- cv_enet_lambda1$lambda.min

#Regression
fit_enet <- gcdnet(X, Y, lambda=lambda1_enet_opt, lambda2 = lambda2_enet_opt, method="ls",pf=weight_enet)
beta_enet <- coef(fit_enet)
beta_enet



############################################################## 
################adaptive elastic net (AENET)################## 
############################################################## 
v_star <- log(sum(beta_enet!=0)) / log(n)
gamma <- ceiling(2*v_star/(1 - v_star)) + 1
X_sd <- apply(X, 2, sd)
weights_ad <- (abs(beta_enet[-1]*X_sd) + 1/n)^(-gamma)
X_non0 <- X[, which(beta_enet[-1]!=0)]
weights_aenet <- weights_ad[which(beta_enet[-1] != 0)] 
weights_aenet[1:ncol(covs)]=0 #force in confounders


#lambda2
lambda2_list_aenet <- seq(0.001, 0.025, by = 0.001)
min_cv_aenet <- numeric(length(lambda2_list_aenet))
for(i in 1:length(lambda2_list_aenet)) {
  cv_aenet <- cv.gcdnet(X_non0, Y, foldid = foldid_mat[, 2],
                        lambda2 = lambda2_list_aenet[i], pf = weights_aenet, method="ls") 
  min_cv_aenet[i] <- min(cv_aenet$cvm)
  print(i)
}
lambda2_aenet_opt<-lambda2_list_aenet[which.min(min_cv_aenet)]


#lambda1
cv_aenet <- cv.gcdnet(X_non0, Y, foldid = foldid_mat[, 2],  
                      lambda2 = lambda2_aenet_opt,method = "ls", pf = weights_aenet) 
lambda1_opt_aenet <- cv_aenet$lambda.min


#Regression
fit_aenet <- gcdnet(X_non0, Y, lambda = lambda1_opt_aenet, lambda2 = lambda2_aenet_opt, 
                    method="ls", pf = weights_aenet)
beta_aenet <- coef(fit_aenet)
aenet_tab <- matrix(0, sum(beta_aenet!=0), 1)
rownames(aenet_tab) <- rownames(beta_aenet)[which(beta_aenet!=0)]
aenet_tab[, 1] <- beta_aenet[which(beta_aenet!=0), ] #beta
aenet_tab_metals<-subset(aenet_tab,!row.names(aenet_tab)%in%colnames(training)[3:22])
aenet_tab_metals



############################################################## 
##############Construct ERS in training set################### 
############################################################## 

# Create dataset include standardized log-10 transformed metal variables
metals<-cbind(metal,square_metal,pairs_metal)

#ERS in the training set
ers_train<-as.matrix(metals[,rownames(aenet_tab_metals)[2:ncol(aenet_tab_metals)]])%*%as.vector(aenet_tab_metals[2:ncol(aenet_tab_metals),1])

## ERS in the testing set can be calculated in a similar manner
