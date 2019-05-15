matt <- matrix


means <- colMeans(matt)

par(mfrow=c(1,4))
apply(matt,2,hist)

attach(matt)
matt

hist(matt$GSM534572)

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