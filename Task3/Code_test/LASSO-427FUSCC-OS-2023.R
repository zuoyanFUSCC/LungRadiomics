rm(list = ls())

data<-read.csv("C:/Users/zuoya/Documents/fuscc_os_nomogram/lasso_nomogram/lassoRisk.csv",header = T)
colnames(data)
#归一化
data[,c(9:30)]<-scale(data[,c(9:30)])

###
data_CT<-data[,c(9:30)]

###
data_surv<-subset(data,select = c(`status`,`time`))

colnames(data_surv)<-c('status','time')

# data_surv$OS<-ifelse(data_surv$OS=='Dead',1,0)

###

data_clinic<-data[,c(2,3,4,5,6)]
colnames(data_clinic)<-c("Gender","Grades","TNM","LD","SUVmax")

###
data_se<-cbind(data_clinic,data_surv,data_CT)
str(data_se)

###
#去掉NA
data_surv.raw<-na.omit(data_se)
str(data_surv.raw)

data_surv.raw$time<-data_surv.raw$time/365
colnames(data_surv.raw)


PCT_fea<-colnames(data_surv.raw)[8:29]
PCT_fea
formula_multicox<-as.formula(paste("Surv(time,status)~",paste(PCT_fea,collapse = "+")))
library(ggplot2)
library(ggpubr)
library(survminer)
library(survival)
multicox<-coxph(formula=formula_multicox,data=data_surv.raw)
ggforest(
  multicox,  #coxph得到的Cox回归结果
  # data = data_surv.raw,  #数据集
  data = data_surv.raw,  #数据集
  main = 'Hazard ratio of multicox PET/CT radiomics features',  #标题
  cpositions = c(0, 0.3, 0.4),  #前三列距离
  fontsize = 0.9, #字体大小
  refLabel = 1, #相对变量的数值标签，也可改为1
  noDigits =2 #保留HR值以及95%CI的小数位数
)


data_surv.raw$PETCTRad_score<-predict(multicox,type = 'risk',newdata =data_surv.raw)
colnames(data_surv.raw)

data_surv.raw$PETCTRad_score.status<-ifelse(data_surv.raw$PETCTRad_score<median(data_surv.raw$PETCTRad_score),'low','high')
dim(data_surv.raw) #417  31

a <- median(data_surv.raw$PETCTRad_score)*365
a  #354.8552          
fit.PCT <- survfit(Surv(time,status) ~PETCTRad_score.status,
                   data = data_surv.raw)


ggsurvplot(fit.PCT, data = data_surv.raw,
           surv.median.line = "hv", # 添加中位数生存时间线
           font.title=18,
           break.time.by =1,
           title='Survival analysis',
           #ylim=c(0.6,1),#设置Y轴显示范围
           xlim=c(0,8),
           legend.labs=c("High risk", "Low risk"),
           # legend.labs=417,
           legend.title="Risk status",
           xlab="Follow up time(years)",
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

data.compare<-data_surv.raw
write.csv(data.compare,"./fuscc425-os_clinic.csv")
colnames(data.compare)
# data.compare<-subset(data.compare,Smoking.status!=''&TNM!='')
# data.compare<-subset(data.compare,Grades!=''&TNM!='')
str(data.compare)
colnames(data.compare)

# for(i in names(data.compare)[c(1,3,4,5,6,34)]){data.compare[,i]<-as.factor(data.compare[,i])}

# colnames(data.compare)
library(rms)
library(Hmisc)
library(lattice)
library(Formula)
bb<-datadist(data.compare)
options(datadist='bb')



model1 <- coxph(Surv(time,status)~Grades,
                data=data.compare)
model2 <- coxph(Surv(time,status)~Gender,
                data=data.compare)
model3 <- coxph(Surv(time,status)~TNM,
                data=data.compare)
model4 <- coxph(Surv(time,status)~PETCTRad_score,
                data=data.compare)
model5 <- coxph(Surv(time,status)~LD,
                data=data.compare)
model6<-coxph(Surv(time,status)~SUVmax,
              data=data.compare)

library(ggDCA)

dca2<- dca(model1,model2,model3,model4,model5,model6,times=c(1,3,5))#c(1,5,8)
AUDC(dca2)

# name1<-c('Grades','age','gender')
name1<-c('Grades','Gender','TNM','PETCTRad_score','LD','SUVmax')
library(ggsci)
ggplot(dca2,
       linetype =1,
       lwd = 1)+
  theme_bw()+
  theme(legend.position='top')+
  scale_colour_d3(
    label = c(name1,"ALL","None")
  )+
  labs(title = "DCA curve of PETCTRad_score and clinical factors")+
  theme(plot.title = element_text(hjust = 0.5))


f<-psm(formula = Surv(time,status) ~Gender+Grades+TNM+LD+SUVmax+PETCTRad_score,
       data=data.compare,
       x=T,y=T,
       #surv = T
)  #,time.inc =2920
# f<-psm(formula = Surv(OS.time,OS) ~TNM+Grades+Age+Gender+PCTscore.status+Smoking.status, +PCTscore.status+Grades
#        data=data.compare,
#        x=T,y=T,
#        #surv = T
# )  #,time.inc =2920
summary(f)
c_index1<-  rcorrcens(Surv(time,status) ~ predict(f), 
                      data = data.compare)[1]
c_index1
#####
####普通画图函数########
surv<- Survival(f)
surv1<-function(x) surv(1,x)
surv2<-function(x) surv(3,x)
surv3<-function(x) surv(5,x)
data.compare$Gender <- factor(data.compare$Gender,levels = c(0,1),
                              labels = c('Female','Male'))
nom<- nomogram(f, fun = list(surv1,surv2,surv3), #感兴趣的预测时间
               lp=T, #是否显示回归系数轴
               maxscale = 100,#按百分制
               funlabel=c("1-year survival", "3-year survival", "5-year survival")#,
               #fun.at=c(0.99,0.95, 0.85,0.7, 0.5, 0.1)#设置预测复发率的范围
)
plot(nom,
     xfrac=0.5,#图形与变量占比
     cex.var=1.1,#变量字体加粗
     cex.axis=0.9,#数轴 字体的大小
     lmgp=0.2,#文字与刻度的距离
     # lplabel='Linear Predictorlp',#线性预测轴名字
     points.label='Risk points',#变量分数名字
     total.points.label='Total Points',#总分名字
     ia.space=.3,
     cap.labels=FALSE,
     col.grid = gray(c(0.8, 0.95)),
     label.every = 1,
     force.label=F,#强制标记的每个刻度线都绘制标签
)


library(regplot)

nom1<-regplot(f,clickable = F,
              plots = c('density','no plot'),
              title = 'Nomogram',
              # observation=someone1,
              showP = T,
              subticks = T,
              dencol='#A6CEE3',
              boxcol='grey60',
              obscol='red',
              spkcol='red',
              leftlabel=T,
              points = T,
              odds = F,
              rank = 'range',
              failtime = c(1,3,5),
              prafail=T,
              droplines=F,
)



#+Grades

f3<-psm(formula= Surv(time,status) ~ Gender+Grades+TNM+LD+SUVmax+PETCTRad_score,
        data=data.compare, x=T, y=T,#surv = T, 
        time.inc = 1) 

cal3<-calibrate(f3, cmethod="KM", method="boot",u=1,m=35,B=1000) 


f5<-psm(formula= Surv(time,status) ~Gender+Grades+TNM+LD+SUVmax+PETCTRad_score,
        data=data.compare, x=T, y=T,#surv = T, 
        time.inc = 3) 

cal5<-calibrate(f5, cmethod="KM", method="boot",u=3,m=35,B=1000) 

f8<-psm(formula= Surv(time,status) ~Gender+Grades+TNM+LD+SUVmax+PETCTRad_score,
        data=data.compare, x=T, y=T,#surv = T, 
        time.inc = 5) 

cal8<-calibrate(f8, cmethod="KM", method="boot",u=5,m=35,B=1000) 
##
par(mar=c(8,5,3,2),cex = 1.0)
plot(cal3,lwd = 2,lty = 1,errbar.col = c("#BC3C29FF"),
     bty = "o", #只画左边和下边框
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col = c("#BC3C29FF"),
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal3[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#BC3C29FF"), pch = 16)
mtext("")

plot(cal5,lwd = 2,lty = 1,errbar.col = c("#0072B5FF"),
     xlim = c(0,1),ylim= c(0,1),col = c("#0072B5FF"),add = T)
lines(cal5[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#0072B5FF"), pch = 16)
mtext("")
plot(cal8,lwd = 2,lty = 1,errbar.col = c("#E18727FF"),
     bty = "l", #只画左边和下边框
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col = c("#E18727FF"),
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6,add = T)
lines(cal8[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#E18727FF"), pch = 16)

abline(0,1, lwd = 2, lty = 3, col = c("#666666"))

legend("bottomright", #图例的位置
       legend = c("1-year","3-year","5-year"), #图例文字
       col =c("#BC3C29FF","#0072B5FF","#E18727FF"), #图例线的颜色，与文字对应
       lwd = 2,#图例中线的粗细
       bty = "n",seg.len=2,cex=1)
#####c-index#####
c_index1<-  rcorrcens(Surv(time,status) ~ predict(f), 
                      data = data.compare)[1]
c_index1

