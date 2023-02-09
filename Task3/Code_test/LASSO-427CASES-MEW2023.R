rm(list = ls())
install.packages('glmnet')
install.packages('readxl')
# LASSO and a ten-fold cross validation

library(Matrix)
library(glmnet)
library(survival)
library(readxl)
# LASSO and a ten-fold cross validation
data <- read.csv("C:/Users/zuoya/Desktop/fs_pre_425-fuscc-os.csv")


str(data)
# Change the age to a numeric
fix(data) # Select the numeric after clicking on age
colnames(data)
x <- as.matrix(data[,2:444])
x <- scale(x, center = TRUE, scale = TRUE)
y <- as.matrix(Surv(data$time,data$status))
lasso <- glmnet(x,y,alpha=1,family="cox")
print(lasso)
plot(lasso,xvar="lambda",label = "TRUE")
set.seed(2)
cv.lasso <- cv.glmnet(x,y,family="cox",alpha=1,nfolds=10)
plot(cv.lasso)
coef <- coef(cv.lasso,s="lambda.min")
lasso <- cv.lasso$lambda.min
lasso  #0.04810934
log(lasso)  # -3.034279

index <- which(coef != 0)
actCoef <- coef[index]
lassoGene=row.names(coef)[index]
geneCoef=cbind(Gene=lassoGene,Coef=actCoef)
write.csv(geneCoef,file="geneCoef.csv",quote=F,row.names=T)
#最终获得了12个m6A-LPS

riskScore=predict(cv.lasso, newx = x, s = "lambda.min",type="response")
median(riskScore) # 0.8921932
outCol=c("time","status",lassoGene)
risk=as.vector(ifelse(riskScore>median(riskScore),"high","low"))
View(risk)
table(risk) #高风险208，低风险209
outTab=cbind(data[,outCol],riskScore=as.vector(riskScore),risk)

write.csv(cbind(id=rownames(outTab),outTab),
          file="lassoRisk.csv",
          quote=F,
          row.names=T)


# Univariate Cox regression
library(survival)
library(plyr)
y<- Surv(time = data$time,event = data$status==1)
Uni_cox_model<-
  function(x){
    FML <- as.formula(paste0 ("y~",x))
    cox<- coxph(FML,data=data)
    cox1<-summary(cox)
    HR <- round(cox1$coefficients[,2],3)
    PValue <- round(cox1$coefficients[,5],3)
    CI5 <-round(cox1$conf.int[,3],3)
    CI95 <-round(cox1$conf.int[,4],3)
    Uni_cox_model<- data.frame(
      names <-rownames(cox1$conf.int),
      'HR' = HR,
      'CI5' = CI5,
      'CI95' = CI95,
      'P' = PValue)
    return(Uni_cox_model)
  }  













# variable.names<- colnames(data)[c(4,6,11,21,27,28,34,37,38,44,50,56,57,60,71)];variable.names
variable.names<- colnames(data)[c(445,446,447,448,449,450,451)];variable.names
Uni_cox <- lapply(variable.names, Uni_cox_model)
Uni_cox<- ldply(Uni_cox,data.frame)
Uni_cox$HR.CI95<-paste0(Uni_cox$HR," (",Uni_cox$CI5,'-',Uni_cox$CI95,")");Uni_cox



# Multivariable  Cox regression
library(survival)
mul_cox<-coxph(Surv(time,status==1)~
                 # age+RA+Immunosuppressive.agents+DLCO+RVD+RAA+PASP+LVEF+CRP+BNP+AST+GGT+ALB+Honeycombing,
                 Gender+Grades+TNM+LD+SUVmax,
               data=data);summary(mul_cox)
cox<-summary(mul_cox) 
cox$coefficients    
cox$conf.int  
mul_HR<- round(cox$coefficients[,2],3) 
mul_PValue<- round(cox$coefficients[,5],3) 
mul_CI1<-round(cox$conf.int[,3],3)
mul_CI2<-round(cox$conf.int[,4],3)
mul_CI95<-paste(mul_CI1,'-',mul_CI2)
mul_cox1 <- data.frame("HR" =mul_HR,
                       "CI95" =mul_CI95,
                       "P"=mul_PValue);mul_cox1

library(ggplot2)
library(ggpubr)
library(survminer)
library(survival)

ggsurvplot(mul_cox, data = data,
           # surv.median.line = "hv", # 添加中位数生存时间线
           # font.title=20,
           # break.time.by =1,
           # title='Survival analysis',
           # size=1,
           # legend.labs=c('High','Low'),
           # legend.title="Risk status",
           # xlab="Survival time(years)",
           # pval = T,
           # conf.int = F, # 设置添加置信区间
           # pval.size=8, # 指定p值文本大小的数字，默认为 5。
           # pval.coord=c(0.1,0.2),
           # #palette ="aaas", # 设置颜色画板
           # palette ="hue", # 设置颜色画板
           # ggtheme = theme_classic()+
           #   theme(plot.title = element_text(size=18.8,hjust = 0.5),
           #         legend.title = element_text(size = 18.8),
           #         legend.text = element_text(size = 18.8))+
           #   theme(
           #     axis.title.x = element_text(size = 18.8),
           #     axis.title.y = element_text(size = 18.8),
           #     axis.text.x = element_text(size = 18.8),
           #     axis.text.y = element_text(size = 18.8)),#调节坐标轴刻度大小 
           # risk.table = T,risk.table.fontsize=7,
           # tables.height = 0.3,
           # risk.table.pos='out',
           # risk.table.title='',
           # tables.theme =  theme_cleantable()+theme(plot.title = element_text(hjust = 0.5))
           surv.median.line = "hv", # 添加中位数生存时间线
           font.title=18,
           break.time.by =1,
           title='Survival analysis',
           #ylim=c(0.6,1),#设置Y轴显示范围
           xlim=c(0,8),
           legend.labs=c("High risk", "Low risk"),
           legend.title="Risk status",
           xlab="Survival time(years)",
           ncensor.plot = F,#绘制在时间t的审查主题的数量
           pval = T,
           pval.method = T, #设置添加P值计算方法
           pval.coord=c(1,0.1),
           conf.int = TRUE, # 设置添加置信区间
           pval.size=4, # 指定p值文本大小的数字，默认为 5。
           pval.method.size=4, # 指定检验方法 log.rank 文本的大小
           pval.method.coord=c(0.1,0.1), # 指定检验方法 log.rank 文本的坐标
           risk.table = TRUE, # 设置添加风险因子表
           tables.height = 0.3, # 设置风险表的高度
           tables.theme =  theme_survminer()+theme(plot.title = element_text(hjust = 0.5))+
             theme(panel.grid=element_blank()),#theme_cleantable(), # 设置风险表的主题
           palette ="hue", # 设置颜色画板
           #"hue","grey","npg","aaas","lancet","jco", "ucscgb","uchicago","simpsons"
           ggtheme = theme_bw()+theme(plot.title = element_text(size=22,hjust = 0.5))+
             theme(panel.grid=element_blank()),#使标题居中
           fontsize=4,#风险表的字体大小
           risk.table.fontsize=4,
           font.x=14,#指定x轴字体大小
           font.y=14#指定y轴字体大小
)


# Nomogram
data  <- read.csv("C:/Users/zuoya/Desktop/fs_pre_425-fuscc-os-riskscore.csv")

data$time<-data$time/365
colnames(data)
suppressMessages(library(rms))
dd <- datadist(data)
options(datadist="dd")
cox_nomo1 <-  cph(Surv(time,status)~Gender+Grades+TNM+LD+SUVmax+RadScore,data=data,x=T,y=T,surv=T)
surv <- Survival(cox_nomo1)
surv1 <- function(x)surv(3,lp=1-x)#"lp=1-x)
surv2 <- function(x)surv(5,lp=1-x)
nomo_2a<-nomogram(cox_nomo1, fun=list(surv1,surv2), lp=F,funlabel =c("3-year Death", "5-year Death"),fun.at =c(0.05, seq(0.1,0.9, by=0.1), 0.95))
plot(nomo_2a, col.grid=c("pink","cyan"),xfrac = 0.3,cex.var = 1,cex.axis = 1,lmgp = 0.3)


# Harrell¡¯s C index
data$p_lp <- predict(cox_nomo1, data, type="lp")
c_harrell_nomogram <- (cph(Surv(time,status)~p_lp, data=data,x=TRUE,y=TRUE)$stats["Dxy"]+1)/2
c_harrell_nomogram

cox_GAP <-  cph(Surv(time,status)~GAP,data=data,x=T,y=T,surv=T)
data$p_lp <- predict(cox_GAP, data, type="lp")
c_harrell_GAP <- (cph(Surv(time,status)~p_lp, data=data,x=TRUE,y=TRUE)$stats["Dxy"]+1)/2
c_harrell_GAP


# internal validation
library(MASS)
N_bootstrap <- 1000      
c_harrell_resample <- 0      
c_harrell_original <- 0    
for (i in 1:N_bootstrap){
  data.train <- data[sample(1:nrow(data), replace=TRUE),]
  Outcome <- "Surv(time,status)"
  CandidateVariables <- c("Gender","rades","TNM","LD","SUVmax","RadScore")
  Formula <- formula(paste(paste(Outcome,"~", collapse=" "), 
                           paste(CandidateVariables, collapse=" + ")))
  
  model.full <- coxph(Formula, data=data.train,x=TRUE)
  
  model.train <- stepAIC(model.full, direction="both")
  
  
  c_harrell_resample[i] <- (cph(Surv(time,status)~p_lp, data=data.train,x=TRUE,y=TRUE)$stats["Dxy"]+1)/2
  data.test <- data
  data.test$p_lp <- predict(model.train, data.test, type="lp")
  c_harrell_original[i] <- (cph(Surv(time,status)~p_lp, data=data.test,x=TRUE,y=TRUE)$stats["Dxy"]+1)/2
}

c_harrell_optimism <- mean(c_harrell_resample - c_harrell_original)
c_harrell_optimism
c_harrell_nomogram - c_harrell_optimism



# Calibration plots
nomogram_3y <- cph(Surv(time,status)~Gender+rades+TNM+LD+SUVmax+RadScore,x=T, y=T, surv=T, time.inc = 3, data=data)
cal3_nomogram <- calibrate(nomogram_3y, cmethod="KM", method="boot", u=3, m= 168, B= 1000,conf.int="TRUE")
plot(cal3_nomogram,xlab="Predicted 3 Years Survival",ylab="Fraction Surviving 3 years")

nomogram_5y <- cph(Surv(time,status)~Gender+rades+TNM+LD+SUVmax+RadScore,x=T, y=T, surv=T, time.inc = 5, data=data)
cal5_nomogram <- calibrate(nomogram_5y, cmethod="KM", method="boot", u=5, m= 168, B= 1000,conf.int="TRUE")
plot(cal5_nomogram,xlab="Predicted 5 Years Survival",ylab="Fraction Surviving 5 years")

# GAP_3y <- cph(Surv(time,status)~GAP,x=T, y=T, surv=T, time.inc = 3, data=data)
# cal3_GAP <- calibrate(GAP_3y, cmethod="KM", method="boot", u=3, m= 168, B= 1000,conf.int="TRUE")
# plot(cal3_GAP,xlab="Predicted 3 Years Survival",ylab="Fraction Surviving 3 years")
# 
# GAP_5y <- cph(Surv(time,status)~GAP,x=T, y=T, surv=T, time.inc = 5, data=data)
# cal5_GAP <- calibrate(GAP_5y, cmethod="KM", method="boot", u=5, m= 168, B= 1000,conf.int="TRUE")
# plot(cal5_GAP,xlab="Predicted 5 Years Survival",ylab="Fraction Surviving 5 years")
abline(0,1, lwd = 2, lty = 3, col = c("#666666"))

legend("bottomright", #图例的位置
       legend = c("1-year","3-year"), #图例文字
       col =c("#BC3C29FF","#0072B5FF"), #图例线的颜色，与文字对应
       lwd = 2,#图例中线的粗细
       bty = "n",seg.len=2,cex=1)

# net reclassification improvement (NRI), integrated discrimination improvement (IDI) M1 is IDI, M2 is NRI
library(survival)
install.packages(survC1)
library(survC1)
library(survIDINRI)
# 3year
model.1 <- coxph(Surv(time,status)~Gender+rades+TNM+LD+SUVmax+RadScore,data=data,x=TRUE)
m11 <- predict(model.1, data=data, type="lp")
model.2 <- coxph(Surv(time,status)~GAP,data=data,x=TRUE)
m22 <- predict(model.2, data=data, type="lp")
IDI<-IDI.INF(data[,c("time","status")],m22, m11, 36, npert = 300, npert.rand = NULL, seed1 = NULL, alpha = 0.05)
IDI.INF.OUT(IDI)
# 5year
IDI<-IDI.INF(data[,c("time","status")],m22, m11, 60, npert = 300, npert.rand = NULL, seed1 = NULL, alpha = 0.05)
IDI.INF.OUT(IDI)


# likelihood-ratio test
TN.GAP <- cph(Surv(time,status)~ GAP, 
              data=data, na.action=na.omit )
TNC.nomogram <- cph(Surv(time,status)~ Gender+rades+TNM+LD+SUVmax+RadScore, 
                    data=data, na.action=na.omit )
TNC.combine <- cph(Surv(time,status)~ GAP+Gender+rades+TNM+LD+SUVmax+RadScore, 
                   data=data, na.action=na.omit )
TN.GAP
TNC.nomogram
TNC.combine
TNC1 <- lrtest(TN.GAP, TNC.combine)
TNC1
TNC2 <- lrtest(TNC.nomogram, TNC.combine)
TNC2


# decision curve analysis (DCA)
library(rms)
library(ggDCA)
setwd("C:/Users/zuoya/Documents/fuscc_os_nomogram/lasso_nomogram/lasso-cox-survival/lasso-cox-nomogram/CTD-ILD")
source("stdca.R")
# 3 year
coxmod <- cph(Surv(time,status)~Gender+rades+TNM+LD+SUVmax+RadScore,x=T, y=T, surv=T, time.inc = 3,data=data)
data$our_model <- c(1 - (summary(survfit(coxmod,newdata=data), times=3)$surv))
coxmod1 <- cph(Surv(time,status)~GAP,x=T, y=T, surv=T, time.inc = 3,data=data)
data$GAP <- c(1 - (summary(survfit(coxmod1,newdata=data), times=3)$surv))
stdca(data=data, outcome="status", ttoutcome="time", timepoint=3, predictors=c("our_model","GAP"), xstop=0.5, smooth=TRUE)
# 5 year
coxmod <- cph(Surv(time,status)~Gender+rades+TNM+LD+SUVmax+RadScore,x=T, y=T, surv=T, time.inc = 5,data=data)
data$our_model <- c(1 - (summary(survfit(coxmod,newdata=data), times=5)$surv))
coxmod1 <- cph(Surv(time,status)~GAP,x=T, y=T, surv=T, time.inc = 5,data=data)
data$GAP <- c(1 - (summary(survfit(coxmod1,newdata=data), times=5)$surv))
stdca(data=data, outcome="status", ttoutcome="time", timepoint=5, predictors=c("our_model","GAP"), xstop=0.5, smooth=TRUE)


# All-cause mortality among 675 CTD-ILD patients
library(Matrix)
library(glmnet)
library(survival)
library(ggplot2)
library(ggpubr)
library(survminer)
A <- read_excel("C:/Users/zuoya/Documents/fuscc_os_nomogram/lasso_nomogram/lasso-cox-survival/lasso-cox-nomogram/CTD-ILD/data/B.xlsx")
fix(A)
fit <- survfit(Surv(time,status)~1,data=A)
res.sum <- surv_summary((fit))
res.sum
ggsurvplot(fit,
           pval=T,conf.int = T,
           risk.table=T,
           risik.table.col="strata",
           linetype="strata",
           surv.median.line="hv",
           ggtheme=theme_bw(),
           palette=c("#E7B800","#2E9FDF"),
           fun = "event")

