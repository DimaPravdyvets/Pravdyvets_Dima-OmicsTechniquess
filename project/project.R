pr_mat <- GSE21667_series_matrix
#########################################################

boxplot(pr_mat)
#########################################################

means <- colMeans(pr_mat)
#########################################################

par(mfrow=c(1,4))
names(pr_mat)
apply(pr_mat,2,hist)
#########################################################

library(ggpubr)
par(mfrow=c(1,4))
apply(pr_mat,2,ggqqplot)
apply(pr_mat,2,ggdensity)
########################################################
t.test(pr_mat)
#########################################################

apply(pr_mat,2,shapiro.test)
#########################################################


######################################################### ######################################################### 

morvar <- GSE21667_series_matrix
attach(morvar)
names(morvar)

first<-glm(ID_REF~GSM534572,family=binomial(link="logit"))
first<-glm(ID_REF~GSM534572,family=poisson(link="logit"))
first<-glm(ID_REF~GSM534572,family=quasipoisson(link="logit"))
summary(first)

second<-glm(ID_REF~GSM536096,family=binomial(link="logit"))
second<-glm(ID_REF~GSM536096,family=poisson(link="logit"))
second<-glm(ID_REF~GSM536096,family=quasipoisson(link="logit"))
summary(second)

third<-glm(ID_REF~GSM536097,family=binomial(link="logit"))
third<-glm(ID_REF~GSM536097,family=poisson(link="logit"))
third<-glm(ID_REF~GSM536097,family=quasipoisson(link="logit"))
summary(third)

fourth<-glm(ID_REF~GSM536098,family=binomial(link="logit"))
fourth<-glm(ID_REF~GSM536098,family=poisson(link="logit"))
fourth<-glm(ID_REF~GSM536098,family=quasipoisson(link="logit"))
summary(fourth)

BigCheck<-glm(ID_REF~GSM534572+GSM536096+GSM536097+GSM536098,family=poisson(link="logit"))
summary(BigCheck)
#########################################################

#########################################################

#########################################################