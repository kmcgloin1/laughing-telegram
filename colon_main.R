library(survival)
library(survminer)
library(ggpubr)

# import the data
colonMod <- read.csv("colonModified.csv",header=TRUE)

# make sure the data looks correct
head(colonMod)

# get summary statistics of the data
summary(colonMod)

# ----- Prelimiary Analysis -----

# this code will look at one-variable Cox PH models, and report
# the resulting beta coefficients, test statistics, and p-values
testCovariates <- c("rx","sex","age","obstruct","perfor","adhere","extent",
                    "surg","node4")
univ_formulas <- sapply(testCovariates,
                        function(x) as.formula(paste('Surv(time, status)~', x)))
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = colonMod)})
# Extract data
univ_results <- lapply(univ_models,
                       function(x){
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=3)
                         wald.test<-signif(x$wald["test"], digits=3)
                         beta<-signif(x$coef[1], digits=3);#coeficient beta
                         res<-c(beta, wald.test, p.value)
                         names(res)<-c("beta", "wald.test",
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- t(as.data.frame(univ_results))
format(as.data.frame(res),scientific=FALSE)

# From our results of the above function, we will consider the predictors rx, 
# obstruct, adhere, extent, surg, and node4 since they all had p-values below 
# 0.05. Note that these were the same variables I had gotten before using this
# above function. 

# ----- Numerical Variables: Extent and Age -----

# ----- Martingale Residual Plots -----

ggcoxfunctional(Surv(time, status) ~ extent, data = colonMod)
# Looks to be linear

ggcoxfunctional(Surv(time, status) ~ age, data = colonMod)
# Not linear, could use a transfomration function

fit <- coxph(Surv(time,status)~age,data=colonMod)
summary(fit) # age here is not statistically significant

fit2 <- coxph(Surv(time,status)~I(age>50),data=colonMod)
summary(fit2) # age here is not statistically significant

fit3 <- coxph(Surv(time,status)~age*I(age>50),data=colonMod) 
summary(fit3) # age here is statistically significant

fit4 <- coxph(Surv(time,status)~I(age>45),data=colonMod)
summary(fit4) # age here is statistically significant

fit5 <- coxph(Surv(time,status)~age*I(age>45),data=colonMod)
summary(fit5) # age here is not statistically significant

# ----- MODEL FITTING -----

fit1 <- coxph(Surv(time, status) ~ rx + age*I(age>50) + adhere + node4 + extent 
              + obstruct + surg, data=colonMod)
summary(fit1)

fit2 <- coxph(Surv(time, status) ~ rx + age*I(age>45) + adhere + node4 + extent 
              + obstruct + surg, data=colonMod)
summary(fit2)

fit3 <- coxph(Surv(time, status) ~ rx + adhere + node4 + extent + obstruct + surg,
              data=colonMod)
summary(fit3)

AIC(fit1, fit2, fit3)

# ----- Testing Cox PH Assumption of our Models -----

test.fit1 <- cox.zph(fit1)
test.fit1
test.fit2 <- cox.zph(fit2)
test.fit2
test.fit3 <- cox.zph(fit3)
test.fit3
# We reject our global model for all three fits, which is a problem. Thus, we 
# need to investigate node4 and obstruct.

# perform Global Schoenfeld Test for node4 and obstruct to determine if 
# proportional hazards assumption is violated
ggcoxzph(test.fit1,var="node4")
ggcoxzph(test.fit1,var="obstruct")

ggcoxzph(test.fit2,var="node4")
ggcoxzph(test.fit2,var="obstruct")

ggcoxzph(test.fit3,var="node4")
ggcoxzph(test.fit3,var="obstruct")

# new models after removing node4 and obstruct since violated proportional 
# hazards assumption
fit1 <- coxph(Surv(time, status) ~ rx + age*I(age>50) + adhere + extent + surg, data=colonMod)
summary(fit1)

fit2 <- coxph(Surv(time, status) ~ rx + I(age>45) + adhere + extent + surg, data=colonMod)
summary(fit2)

fit3 <- coxph(Surv(time, status) ~ rx + adhere + extent + surg, data=colonMod)
summary(fit3)

AIC(fit1, fit2, fit3) # fit1 is only barely better now, but in reality they are
# all essentially the same in performance

test.fit1 <- cox.zph(fit1)
test.fit1
test.fit2 <- cox.zph(fit2)
test.fit2
test.fit3 <- cox.zph(fit3)
test.fit3
# All three models satisify the proportional hazards assumption. 
# Removing variables does not help our above models. This was checked.


# ----- Model Fitting with Interactions -----

fit5 <- coxph(Surv(time, status) ~ rx + I(age>45) + sex + surg + extent + 
                adhere +  rx:extent + rx:sex + rx:I(age>45) + surg:extent + 
                extent:adhere, data=colonMod)
summary(fit5)
AIC(fit5) # 8033.14
test.fit5 <- cox.zph(fit5)
test.fit5
# node4 fails proportional hazards assumption

fit6 <- coxph(Surv(time, status) ~ rx + I(age>45) + sex + surg + extent + 
                adhere +  rx:extent + rx:sex + rx:I(age>45) + surg:extent,
                data=colonMod)
summary(fit6)
AIC(fit6) # 8031.247
test.fit6 <- cox.zph(fit6)
test.fit6 # satisfies proportional hazards assumption

fit7 <- coxph(Surv(time, status) ~ rx + I(age>45) + sex + surg + extent + 
                adhere +  rx:extent + rx:sex + rx:I(age>45), data=colonMod)
summary(fit7)
AIC(fit7) # 8029.443
test.fit7 <- cox.zph(fit7)
test.fit7 # satisfies proportional hazards assumption

fit8 <- coxph(Surv(time, status) ~ rx + I(age>45) + sex + surg + extent + 
                adhere + rx:sex + rx:I(age>45), data=colonMod)
summary(fit8)
AIC(fit8) # 8026.284
test.fit8 <- cox.zph(fit8)
test.fit8 # satisfies proportional hazards assumption
ggcoxzph(test.fit8) # view Schoenfeld Indvidual Plots
# fit8 appears to be our best model


# ----- MY REALLY GOOD MODEL V.> OTHER STUDENT'S MODELS -----

# this was my model; it has an AIC of 7925.409 (the second lowest of all these 
# models) and it satisfies the Cox PH assumption
fit9 <- coxph(Surv(time, status) ~ rx:extent + rx:node4 + rx:sex + rx:I(age>45)
              + surg:extent + extent:adhere, data=colonMod)
AIC(fit9) # 7925.409
test.fit9 <- cox.zph(fit9)
test.fit9 # satisfies proportional hazards assumption


fit10 <- coxph(Surv(time, status) ~ rx + obstruct:adhere + extent:surg +
                 node4:rx, data=colonMod)
AIC(fit10) # 7959.903
test.fit10 <- cox.zph(fit10)
test.fit10 # satisfies proportional hazards assumption


fit11 <- coxph(Surv(time, status) ~ extent + rx + surg:node4 + 
                 adhere:obstruct + rx:node4 + sex:rx + 
                 I(age < 45):extent, data=colonMod)
AIC(fit11) # 7970.401
test.fit11 <- cox.zph(fit11)
test.fit11 # satisfies proportional hazards assumption

fit12 <- coxph(Surv(time, status) ~ rx:sex + age*I(age>49) + age*perfor + 
                 adhere + extent + surg:sex + rx:node4, data=colonMod)
AIC(fit12) # 7921.122
test.fit12 <- cox.zph(fit12)
test.fit12 # satisfies proportional hazards assumption

fit13 <- coxph(Surv(time, status) ~ rx + age*I(age>50) + extent + surg + 
                 obstruct:node4, data=colonMod)
AIC(fit13) # 8034.863
test.fit13 <- cox.zph(fit13)
test.fit13 # satisfies proportional hazards assumption


# models that do not appear to satisfy the Cox PH assumption unless you are 
# working with a signficance level of 0.10, but since this is medical data 
# you should really use a smaller significance level like 0.05
fit14 <- coxph(Surv(time, status) ~ I(age>47) + rx + adhere + extent + surg + 
                 node4, data=colonMod)
AIC(fit14) # 7939.851
test.fit14 <- cox.zph(fit14)
test.fit14 # does not satisfy proportional hazards assumption

fit15 <- coxph(Surv(time, status) ~ rx + adhere + extent + surg + node4, 
               data=colonMod)
AIC(fit15) # 7938.981
test.fit15 <- cox.zph(fit15)
test.fit15 # does not satisfy proportional hazards assumption
