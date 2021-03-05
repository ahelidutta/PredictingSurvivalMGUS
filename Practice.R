#Libraries
library(dplyr)
library("survival")
library("survminer")

#Data
train <- survival::mgus2
train$age
covariates <- c("age", "sex", "hgb", "creat", "mspike", "ptime", "pstat")

#Data censoring

#univariate function for multiple covariates
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(futime, death)~', x)))
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = a)})
signRes <- lapply(univ_models,
                  
                  
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta <-signif(x$coef[1], digits=2);#coefficient beta
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
significance.cox <- coxph(Surv(futime, death) ~ age + sex + hgb + creat + mspike + ptime + pstat, data =  a)
summary(significance.cox)


#bin hemoglobin
train$hgb_bin= cut(train$hgb, breaks=3, labels=FALSE)
train
significance.cox <- coxph(Surv(futime, death) ~ age + sex + hgb_bin + creat + mspike + ptime + pstat, data =  train)


hgb_df <- with(train,
               data.frame(hgb_bin = c(1,2,3), 
                          age = rep(mean(age, na.rm = TRUE), 3),
                          sex = c("M", "M","M"),
                          creat = rep(mean(creat, na.rm = TRUE), 3),
                          mspike = rep(mean(mspike, na.rm = TRUE), 3),
                          ptime = rep(mean(ptime, na.rm = TRUE), 3),
                          pstat = c(1, 1,1)
               )
)

hgb_df_F <- with(train,
               data.frame(hgb_bin = c(1,2,3), 
                          age = rep(mean(age, na.rm = TRUE), 3),
                          sex = c("F", "F","F"),
                          creat = rep(mean(creat, na.rm = TRUE), 3),
                          mspike = rep(mean(mspike, na.rm = TRUE), 3),
                          ptime = rep(mean(ptime, na.rm = TRUE), 3),
                          pstat = c(1, 1,1)
               )
)
hgb_df

fit <- survfit(significance.cox, newdata = hgb_df)
ggsurvplot(fit, conf.int = TRUE, data = hgb_df, legend.labs=c("Low Hemoglobin","Middle Hemoglobin", "High Hemoglobin"),
           ggtheme = theme_minimal())+ ggtitle("Men")
ggsurvplot(fit, conf.int = TRUE, data = hgb_df_F, legend.labs=c("Low Hemoglobin","Middle Hemoglobin", "High Hemoglobin"),
           ggtheme = theme_minimal())+ ggtitle("Women")

#bin age
train$age_bin= cut(train$age, breaks=3, labels=FALSE)
train 
significance.cox <- coxph(Surv(futime, death) ~ age_bin + sex + hgb + creat + mspike + ptime + pstat, data =  train)
age_df <- with(train,
                 data.frame(age_bin = c(1,2,3), 
                            hgb = rep(mean(hgb, na.rm = TRUE), 3),
                            sex = c("F", "F","F"),
                            creat = rep(mean(creat, na.rm = TRUE), 3),
                            mspike = rep(mean(mspike, na.rm = TRUE), 3),
                            ptime = rep(mean(ptime, na.rm = TRUE), 3),
                            pstat = c(1, 1,1)
                 )
)

age_df_M <- with(train,
               data.frame(age_bin = c(1,2,3), 
                          hgb = rep(mean(hgb, na.rm = TRUE), 3),
                          sex = c("M", "M","M"),
                          creat = rep(mean(creat, na.rm = TRUE), 3),
                          mspike = rep(mean(mspike, na.rm = TRUE), 3),
                          ptime = rep(mean(ptime, na.rm = TRUE), 3),
                          pstat = c(1, 1,1)
               )
)
age_df_M
fit <- survfit(significance.cox, newdata = age_df)
fit_M <- survfit(significance.cox, newdata = age_df_M)

ggsurvplot(fit, conf.int = TRUE, data = age_df, legend.labs=c("Low Age","Middle Age", "High Age"),
           ggtheme = theme_minimal())+ ggtitle("Female")

ggsurvplot(fit_M, conf.int = TRUE, data = age_df_M, legend.labs=c("Low Age","Middle Age", "High Age"),
           ggtheme = theme_minimal())+ ggtitle("Male")

#bin mspike
train$mspike_bin= cut(train$mspike, breaks=3, labels=FALSE)
train 
significance.cox <- coxph(Surv(futime, death) ~ age + sex + hgb + creat + mspike_bin + ptime + pstat, data =  train)
mspike_df <- with(train,
               data.frame(age = rep(mean(age, na.rm = TRUE), 3), 
                          hgb = rep(mean(hgb, na.rm = TRUE), 3),
                          sex = c("F", "F","F"),
                          creat = rep(mean(creat, na.rm = TRUE), 3),
                          mspike_bin = c(1,2,3),
                          ptime = rep(mean(ptime, na.rm = TRUE), 3),
                          pstat = c(1, 1,1)
               )
)
fit <- survfit(significance.cox, newdata = mspike_df)
ggsurvplot(fit, conf.int = TRUE, data = mspike_df, legend.labs=c("Low mspike","Middle mspike", "High mspike"),
           ggtheme = theme_minimal())+ ggtitle("Female")

#bin Creat
train$creat_bin= cut(train$creat, breaks=3, labels=FALSE)
train 
significance.cox <- coxph(Surv(futime, death) ~ age + sex + hgb + creat_bin + mspike + ptime + pstat, data =  train)
creat_df <- with(train,
                  data.frame(age = rep(mean(age, na.rm = TRUE), 3), 
                             hgb = rep(mean(hgb, na.rm = TRUE), 3),
                             sex = c("M", "M","M"),
                             creat_bin = c(1,2,3),
                             mspike = rep(mean(mspike, na.rm = TRUE), 3),
                             ptime = rep(mean(ptime, na.rm = TRUE), 3),
                             pstat = c(1, 1,1)
                  )
)
fit <- survfit(significance.cox, newdata = creat_df)
ggsurvplot(fit, conf.int = TRUE, data = creat_df, legend.labs=c("Low creat","Middle creat", "High creat"),
           ggtheme = theme_minimal())+ ggtitle("Male")


require("survival")
fit <- survfit(Surv(futime, death) ~ age + sex + hgb + creat + mspike + ptime + pstat, data = train)
ggsurv <- ggsurvplot(fit, data= train, fun = "event", conf.int = TRUE,
                     ggtheme = theme_bw())
ggsurv$plot +theme_bw() + 
  theme (legend.position = "right")+
  facet_grid(age ~ hgb)

a$pstat
mgus
