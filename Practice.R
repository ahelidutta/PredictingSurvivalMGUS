#Libraries
library(dplyr)
library("survival")
library("survminer")

#Data
train <- survival::mgus2
covariates <- c("age", "sex", "hgb", "creat", "mspike", "ptime", "pstat")

#Data censoring
train$upperBound = ifelse(train$death == 1, train$futime, Inf)
train

#univariate function for multiple covariates
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(upperBound, death)~', x)))
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = a)})
signRes <- lapply(univ_models,
                  
                  
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta <-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- t(as.data.frame(signRes, check.names = FALSE))
as.data.frame(res)

#Multivariable
significance.cox <- coxph(Surv(upperBound, death) ~ age + sex + hgb + creat + mspike + ptime + pstat, data =  a)
summary(significance.cox)

#Graphing for hgb
meanAge <- mean(train$age)
maxHgb <- max(train$hgb, na.rm=TRUE)
minHgb <- min(train$hgb, na.rm=TRUE)
significance.cox <- coxph(Surv(upperBound, death) ~ meanAge + sex + hgb, data =  a)
ggsurvplot(survfit(significance.cox), data=significance.cox, color = "#2E9FDF",
           ggtheme = theme_minimal())
minHgb
maxHgb


