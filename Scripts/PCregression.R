library(pls)
library(readr)

CPMdf<-read_csv("cpmTMM.csv")
fitness <- read_csv("datasheets/fitness.csv")

df<-as.data.frame(CPMdf)
df$Seednum<-fitness$Seednum
df<-df[df$Seednum,50,]
SDdf<-apply(df[,2:2754], 2, FUN=scale)
rownames(SDdf)<-df$ID


SDdf<-data.frame(test)
fit<-df$Seednum/mean(df$Seednum)

#make this example reproducible
set.seed(5)

#fit PCR model
model <- prcomp(SDdf)
PC55<-pca$x[,1:55]
PC10<-pca$x[,1:10] 
#Scree plot drops off at 10
#First 55 have over 0.5% of variation

FitnessPCdf<-cbind(PC10, fit)
pcrmodel<-lm(fit~., data =as.data.frame(FitnessPCdf))

# coefficients for PCA linear regression

beta0 <- pcrmodel$coefficients[1]

betas <- pcrmodel$coefficients[2:11]

beta7 <- pcrmodel$coefficients[8]

beta0 #intercept

(Intercept) 

    905 

# transformation

alphas <- pca$rotation[,1:10] %*% betas

alpha7<-pca$rotation[,7] * beta7

alpha7df<-data.frame(alpha = alpha7, gene = names(alpha7))


t(alphas) # scaled data coefficients

#To obtain the coefficients of the unscaled data:

alphas_unscaled <- alphas / sapply(df[,2:2754],sd)

beta0_unscaled <- beta0 - sum(alphas * sapply(df[,2:2754],mean)/sapply(df[,2:2754],sd))

# coefficients for unscaled data

t(alphas_unscaled)

beta0_unscaled

#We can obtain the estimate of our unscaled data

est <- as.matrix(df[,2:2754]) %*% alphas_unscaled + beta0_unscaled

est

SSE <- sum((est - fit)^2)
SST <- sum((fit- mean(fit))^2)

Rsquared <- 1 - SSE/SST
#0.110133