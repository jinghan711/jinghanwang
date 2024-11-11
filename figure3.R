rm(list = ls())
library(ggplot2)
setwd("D:\\002001\\1_users\\0_manuscript\\027_wangjinghan_Metabolome Signature in Gestational Diabetes Mellitus is Associated with Adverse Birth Outcomes")

glucose1 <- subset(metab_final,metab_final$HIP_met == 0)
glucose1 <- rbind(FBG_median = paste0(sprintf("%.2f",median(glucose1$FBG,na.rm = TRUE))," (",sprintf("%.2f",IQR(glucose1$FBG,na.rm = TRUE)),")"),
                  OGTT1_median = paste0(sprintf("%.2f",median(glucose1$OGTT1,na.rm = TRUE))," (",sprintf("%.2f",IQR(glucose1$OGTT1,na.rm = TRUE)),")"),
                  OGTT2_median = paste0(sprintf("%.2f",median(glucose1$OGTT2,na.rm = TRUE))," (",sprintf("%.2f",IQR(glucose1$OGTT2,na.rm = TRUE)),")"))
glucose2 <- subset(metab_final,metab_final$HIP_met == 1)
glucose2 <- rbind(FBG_median = paste0(sprintf("%.2f",median(glucose2$FBG,na.rm = TRUE))," (",sprintf("%.2f",IQR(glucose2$FBG,na.rm = TRUE)),")"),
                  OGTT1_median = paste0(sprintf("%.2f",median(glucose2$OGTT1,na.rm = TRUE))," (",sprintf("%.2f",IQR(glucose2$OGTT1,na.rm = TRUE)),")"),
                  OGTT2_median = paste0(sprintf("%.2f",median(glucose2$OGTT2,na.rm = TRUE))," (",sprintf("%.2f",IQR(glucose2$OGTT2,na.rm = TRUE)),")"))
glucose3 <- subset(metab_final,metab_final$HIP_met == 2)
glucose3 <- rbind(FBG_median = paste0(sprintf("%.2f",median(glucose3$FBG,na.rm = TRUE))," (",sprintf("%.2f",IQR(glucose3$FBG,na.rm = TRUE)),")"),
                  OGTT1_median = paste0(sprintf("%.2f",median(glucose3$OGTT1,na.rm = TRUE))," (",sprintf("%.2f",IQR(glucose3$OGTT1,na.rm = TRUE)),")"),
                  OGTT2_median = paste0(sprintf("%.2f",median(glucose3$OGTT2,na.rm = TRUE))," (",sprintf("%.2f",IQR(glucose3$OGTT2,na.rm = TRUE)),")"))
glucose4 <- subset(metab_final,metab_final$HIP_met == 3)
glucose4 <- rbind(FBG_median = paste0(sprintf("%.2f",median(glucose4$FBG,na.rm = TRUE))," (",sprintf("%.2f",IQR(glucose4$FBG,na.rm = TRUE)),")"),
                  OGTT1_median = paste0(sprintf("%.2f",median(glucose4$OGTT1,na.rm = TRUE))," (",sprintf("%.2f",IQR(glucose4$OGTT1,na.rm = TRUE)),")"),
                  OGTT2_median = paste0(sprintf("%.2f",median(glucose4$OGTT2,na.rm = TRUE))," (",sprintf("%.2f",IQR(glucose4$OGTT2,na.rm = TRUE)),")"))
glucose <- cbind(glucose1,glucose2,glucose3,glucose4)

BG <- metab_final[,c(35:37,134)]
BG$Group <- factor(ifelse(BG$HIP_met == "0","normoglycemic-non-mGDM",ifelse(BG$HIP_met == "1","hyperglycemic-non-mGDM",ifelse(BG$HIP_met == "2","normoglycemic-mGDM","hyperglycemic-mGDM"))),
                   levels = c("normoglycemic-non-mGDM","hyperglycemic-non-mGDM","normoglycemic-mGDM","hyperglycemic-mGDM"))
colnames(BG)[1:3] <- c("FPG\n(N=2050)","OGTT-1hPG\n(N=2050)","OGTT-2hPG\n(N=2050)")
BG <- BG[,-4]
BG <-  melt(BG,id.vars = "Group",variable.names = "variable",value.name = "value")

box_glu <- ggplot(BG,aes(x = variable,y = value,fill = Group))+
  stat_boxplot(size = 0.6,geom = "errorbar",width = 0.4,position = position_dodge(0.6))+
  geom_boxplot(size = 0.6,width = 0.6,outlier.shape = NA,position = position_dodge(0.6))+
  labs(x = "",y = "",title = "A")+
  scale_fill_manual(values = c("#a6cee3","#00468b","#f99a98","#ed0000"))+
  scale_y_continuous(limits = c(2,13))+
  theme(panel.background = element_blank(),axis.line = element_line(),legend.position = "none",
        axis.text.x = element_text(size = 15,color = "black"),axis.text.y = element_text(size = 15,color = "black"),
        panel.grid = element_blank(),plot.title = element_text(size = 25,face = "bold"))

BGnew <- metab_final[,c(35:37,134)]
colnames(BGnew)[1:3] <- c("FPG","OGTT-1hPG","OGTT-2hPG")
BGnew <-  melt(BGnew,id.vars = "HIP_met",variable.names = "variable",value.name = "value")
pValue <- cbind(wilcox.test(BGnew$value[BGnew$HIP_met==0&BGnew$variable=="FPG"],BGnew$value[BGnew$HIP_met==1&BGnew$variable=="FPG"])$p.value, 
                wilcox.test(BGnew$value[BGnew$HIP_met==0&BGnew$variable=="FPG"],BGnew$value[BGnew$HIP_met==2&BGnew$variable=="FPG"])$p.value, 
                wilcox.test(BGnew$value[BGnew$HIP_met==0&BGnew$variable=="FPG"],BGnew$value[BGnew$HIP_met==3&BGnew$variable=="FPG"])$p.value, 
                
                wilcox.test(BGnew$value[BGnew$HIP_met==0&BGnew$variable=="OGTT-1hPG"],BGnew$value[BGnew$HIP_met==1&BGnew$variable=="OGTT-1hPG"])$p.value, 
                wilcox.test(BGnew$value[BGnew$HIP_met==0&BGnew$variable=="OGTT-1hPG"],BGnew$value[BGnew$HIP_met==2&BGnew$variable=="OGTT-1hPG"])$p.value, 
                wilcox.test(BGnew$value[BGnew$HIP_met==0&BGnew$variable=="OGTT-1hPG"],BGnew$value[BGnew$HIP_met==3&BGnew$variable=="OGTT-1hPG"])$p.value, 
                
                wilcox.test(BGnew$value[BGnew$HIP_met==0&BGnew$variable=="OGTT-2hPG"],BGnew$value[BGnew$HIP_met==1&BGnew$variable=="OGTT-2hPG"])$p.value, 
                wilcox.test(BGnew$value[BGnew$HIP_met==0&BGnew$variable=="OGTT-2hPG"],BGnew$value[BGnew$HIP_met==2&BGnew$variable=="OGTT-2hPG"])$p.value, 
                wilcox.test(BGnew$value[BGnew$HIP_met==0&BGnew$variable=="OGTT-2hPG"],BGnew$value[BGnew$HIP_met==3&BGnew$variable=="OGTT-2hPG"])$p.value)
pValue <- rbind(pValue[,1:3],pValue[,4:6],pValue[,7:9])
colnames(pValue) <- c("0vs1","0vs2","0vs3")
rownames(pValue) <- c("FBG","OGTT1","OGTT2")

write(pValue,"Figure 3A.csv")


bl_t2 <- subset(bl_t2,bl_t2$group == "B"& bl_t2$fetus_n == "1") 
table(duplicated(bl_t2$id)) 

bl_t2 <- merge(bl_t2[,c(4,15,32,36,42:46,48)],metab_final,by = c("enrollment_dt","birth_dt_f","lmp","pre_outcome_date")) #1970
table(duplicated(bl_t2$pid)) #552 
bl_t2$exa_date <- ymd(bl_t2$exa_date)
bl_t2 <- bl_t2 %>%
  group_by(pid) %>%
  filter(exa_date == min(exa_date)) %>%
  ungroup() 

bl_t2_result <- bl_t2 %>%
  group_by(HIP_met) %>%
  summarise(TG_median_t2 = paste0(sprintf("%.2f",median(TG,na.rm = TRUE))," (",sprintf("%.2f",IQR(TG,na.rm = TRUE)),")"),
            TC_median_t2 = paste0(sprintf("%.2f",median(total_cholesterol,na.rm = TRUE))," (",sprintf("%.2f",IQR(total_cholesterol,na.rm = TRUE)),")"),
            HDL_median_t2 = paste0(sprintf("%.2f",median(HDL_C,na.rm = TRUE))," (",sprintf("%.2f",IQR(HDL_C,na.rm = TRUE)),")"),
            LDL_median_t2 = paste0(sprintf("%.2f",median(LDL_C,na.rm = TRUE))," (",sprintf("%.2f",IQR(LDL_C,na.rm = TRUE)),")"))

metab_lipid <- merge(metab_final,bl_t2[,c(8,9,6,7,11)],by = "pid",all.x = TRUE)
colnames(metab_lipid)[135:138] <- c("TG","TC","HDL","LDL")

colnames(pressure)
pressure <- pressure[,c("pid","tgweek2","SBP2","DBP2")]
pressure <- subset(pressure,pressure$pid %in% metab_final$pid) #2047
pressure <- subset(pressure,is.na(pressure$tgweek2) == FALSE) #2013

summary(pressure$SBP2) #80-153
summary(pressure$DBP2) #43-100

metab_pressure <- merge(pressure,metab_final,by = "pid",all.y = TRUE) 


boxplot_cont1 <- metab_lipid %>%
  select(135,136,137,138) %>%
  gather(key = "Variable",value = "Value")
boxplot_cont1$Value_scale <- scale(boxplot_cont1$Value,center = TRUE,scale = TRUE)

boxplot_cont2 <- metab_pressure %>%
  select(3,4) %>%
  gather(key = "Variable",value = "Value")
boxplot_cont2$Value_scale <- scale(boxplot_cont2$Value,center = TRUE,scale = TRUE)

boxplot_cont <- rbind(boxplot_cont1,boxplot_cont2)
boxplot_cont <- cbind(boxplot_cont,HIP_met = rep(metab_lipid$HIP_met,times = 6))

pValue <- cbind(wilcox.test(boxplot_cont$Value_scale[boxplot_cont$HIP_met==0&boxplot_cont$Variable=="TG"],boxplot_cont$Value_scale[boxplot_cont$HIP_met==1&boxplot_cont$Variable=="TG"])$p.value, 
                wilcox.test(boxplot_cont$Value_scale[boxplot_cont$HIP_met==0&boxplot_cont$Variable=="TG"],boxplot_cont$Value_scale[boxplot_cont$HIP_met==2&boxplot_cont$Variable=="TG"])$p.value, 
                wilcox.test(boxplot_cont$Value_scale[boxplot_cont$HIP_met==0&boxplot_cont$Variable=="TG"],boxplot_cont$Value_scale[boxplot_cont$HIP_met==3&boxplot_cont$Variable=="TG"])$p.value, 
                
                wilcox.test(boxplot_cont$Value_scale[boxplot_cont$HIP_met==0&boxplot_cont$Variable=="TC"],boxplot_cont$Value_scale[boxplot_cont$HIP_met==1&boxplot_cont$Variable=="TC"])$p.value, 
                wilcox.test(boxplot_cont$Value_scale[boxplot_cont$HIP_met==0&boxplot_cont$Variable=="TC"],boxplot_cont$Value_scale[boxplot_cont$HIP_met==2&boxplot_cont$Variable=="TC"])$p.value,
                wilcox.test(boxplot_cont$Value_scale[boxplot_cont$HIP_met==0&boxplot_cont$Variable=="TC"],boxplot_cont$Value_scale[boxplot_cont$HIP_met==3&boxplot_cont$Variable=="TC"])$p.value, 
                
                wilcox.test(boxplot_cont$Value_scale[boxplot_cont$HIP_met==0&boxplot_cont$Variable=="LDL"],boxplot_cont$Value_scale[boxplot_cont$HIP_met==1&boxplot_cont$Variable=="LDL"])$p.value, 
                wilcox.test(boxplot_cont$Value_scale[boxplot_cont$HIP_met==0&boxplot_cont$Variable=="LDL"],boxplot_cont$Value_scale[boxplot_cont$HIP_met==2&boxplot_cont$Variable=="LDL"])$p.value, 
                wilcox.test(boxplot_cont$Value_scale[boxplot_cont$HIP_met==0&boxplot_cont$Variable=="LDL"],boxplot_cont$Value_scale[boxplot_cont$HIP_met==3&boxplot_cont$Variable=="LDL"])$p.value, 
                
                wilcox.test(boxplot_cont$Value_scale[boxplot_cont$HIP_met==0&boxplot_cont$Variable=="HDL"],boxplot_cont$Value_scale[boxplot_cont$HIP_met==1&boxplot_cont$Variable=="HDL"])$p.value, 
                wilcox.test(boxplot_cont$Value_scale[boxplot_cont$HIP_met==0&boxplot_cont$Variable=="HDL"],boxplot_cont$Value_scale[boxplot_cont$HIP_met==2&boxplot_cont$Variable=="HDL"])$p.value,
                wilcox.test(boxplot_cont$Value_scale[boxplot_cont$HIP_met==0&boxplot_cont$Variable=="HDL"],boxplot_cont$Value_scale[boxplot_cont$HIP_met==3&boxplot_cont$Variable=="HDL"])$p.value,
                
                wilcox.test(boxplot_cont$Value_scale[boxplot_cont$HIP_met==0&boxplot_cont$Variable=="SBP2"],boxplot_cont$Value_scale[boxplot_cont$HIP_met==1&boxplot_cont$Variable=="SBP2"])$p.value, 
                wilcox.test(boxplot_cont$Value_scale[boxplot_cont$HIP_met==0&boxplot_cont$Variable=="SBP2"],boxplot_cont$Value_scale[boxplot_cont$HIP_met==2&boxplot_cont$Variable=="SBP2"])$p.value, 
                wilcox.test(boxplot_cont$Value_scale[boxplot_cont$HIP_met==0&boxplot_cont$Variable=="SBP2"],boxplot_cont$Value_scale[boxplot_cont$HIP_met==3&boxplot_cont$Variable=="SBP2"])$p.value, 
                
                wilcox.test(boxplot_cont$Value_scale[boxplot_cont$HIP_met==0&boxplot_cont$Variable=="DBP2"],boxplot_cont$Value_scale[boxplot_cont$HIP_met==1&boxplot_cont$Variable=="DBP2"])$p.value, 
                wilcox.test(boxplot_cont$Value_scale[boxplot_cont$HIP_met==0&boxplot_cont$Variable=="DBP2"],boxplot_cont$Value_scale[boxplot_cont$HIP_met==2&boxplot_cont$Variable=="DBP2"])$p.value, 
                wilcox.test(boxplot_cont$Value_scale[boxplot_cont$HIP_met==0&boxplot_cont$Variable=="DBP2"],boxplot_cont$Value_scale[boxplot_cont$HIP_met==3&boxplot_cont$Variable=="DBP2"])$p.value)
pValue <- rbind(pValue[,1:3],pValue[,4:6],pValue[,7:9],pValue[,10:12],pValue[,13:15],pValue[,16:18])
colnames(pValue) <- c("0vs1","0vs2","0vs3")
rownames(pValue) <- c("TG","TC","LDL","HDL","SBP","DBP")
write.csv(pValue,"Figure 3B.csv")

boxplot_cont$Group <- factor(ifelse(boxplot_cont$HIP_met == "0","non-hyperglycemic-non-mGDM",ifelse(boxplot_cont$HIP_met == "1","hyperglycemic-non-mGDM",ifelse(boxplot_cont$HIP_met == "2","normoglycemic-mGDM","hyperglycemic-mGDM"))),
                             levels = c("non-hyperglycemic-non-mGDM","hyperglycemic-non-mGDM","normoglycemic-mGDM","hyperglycemic-mGDM"))
boxplot_cont$Variable <- ifelse(boxplot_cont$Variable == "TG","TG\n(N=1416)",
                                ifelse(boxplot_cont$Variable == "TC","TC\n(N=1415)",
                                       ifelse(boxplot_cont$Variable == "LDL","LDL\n(N=697)",
                                              ifelse(boxplot_cont$Variable == "HDL","HDL\n(N=697)",
                                                     ifelse(boxplot_cont$Variable == "SBP2","SBP\n(N=2013)","DBP\n(N=2013)")))))
boxplot_cont$Variable <- factor(boxplot_cont$Variable,
                                levels = c("TG\n(N=1416)","TC\n(N=1415)","LDL\n(N=697)","HDL\n(N=697)","SBP\n(N=2013)","DBP\n(N=2013)"))
boxplot_cont <- boxplot_cont[,-4]
boxplot_cont <- subset(boxplot_cont,is.na(boxplot_cont$Value) == FALSE) #8251
table(boxplot_cont$Variable)
#TG     TC     LDL     HDL     SBP2     DBP2 
#1416   1415   697     697     2013     2013

box_zb <- ggplot(boxplot_cont,aes(x = Variable,y = Value_scale,fill = Group))+
  stat_boxplot(size = 0.6,geom = "errorbar",width = 0.4,position = position_dodge(0.7))+
  geom_boxplot(size = 0.6,width = 0.7,outlier.shape = NA,notch = FALSE,position = position_dodge(0.7))+
  labs(x = "",y = "",title = "B")+
  scale_fill_manual(values = c("#a6cee3","#00468b","#f99a98","#ed0000"))+
  scale_y_continuous(limits = c(-2,3.5))+
  theme(panel.background = element_blank(),axis.line = element_line(),
        axis.text.x = element_text(size = 15,color = "black"),
        axis.text.y = element_text(size = 15,color = "black"),legend.position = "none",
        panel.grid = element_blank(),plot.title = element_text(size = 25,face = "bold"))


baseline <- metab_final[,c("pid","age","bmi_f","parity","hypertension","PTB","LGA_A","SGA_A","any_Structural_Oneyear","nicu","HIP_met")]
baseline$bmi_f <- ifelse(baseline$bmi_f == "24-27.9"|baseline$bmi_f == ">=28","Yes","No")
baseline$bmi_f <- ifelse(baseline$bmi_f == "<35","No","Yes")
baseline$parity <- ifelse(baseline$parity == "multipara","Yes","No")
baseline <- merge(result_t1_new[,c("pid","gain01")],baseline,by = "pid",all.y = TRUE)

baseline <- baseline[,-1]
colnames(baseline)
baseline <- baseline[,c("bmi_f","bmi_f","parity","hypertension","PTB","LGA_A","SGA_A","any_Structural_Oneyear","nicu","HIP_met")]

mi <- NULL
for (i in 1:9)
{
  mi <- rbind(mi,
              rbind((table(baseline[,i],baseline[,10])[2,1])/1122,
                    (table(baseline[,i],baseline[,10])[2,2])/77,
                    (table(baseline[,i],baseline[,10])[2,3])/611,
                    (table(baseline[,i],baseline[,10])[2,4])/240))
}
mi <- data.frame(mi)
colnames(mi)[1] <- "freq"

mi$exposure <- rep(c("normoglycemic-non-mGDM","hyperglycemic-non-mGDM","normoglycemic-mGDM","hyperglycemic-mGDM"),times = 9)
mi$outcome <- rep(c("Age","BMI","Parity","Hypertension","PTB","LGA","SGA","BD","NICU"),each = 4)  

mi$exposure <- factor(mi$exposure,
                      levels = c("normoglycemic-non-mGDM","hyperglycemic-non-mGDM","normoglycemic-mGDM","hyperglycemic-mGDM"))
mi$outcome <- factor(mi$outcome,
                     levels = c("Age","BMI","Parity","Hypertension","PTB","LGA","SGA","BD","NICU"))
mi$bq <- sprintf("%.1f",mi$freq*100)
mi_tz <- subset(mi,mi$outcome == "Age"|mi$outcome == "BMI"|mi$outcome == "Parity"|mi$outcome == "Hypertension")
mi_jj <- subset(mi,mi$outcome == "PTB"|mi$outcome == "LGA"|mi$outcome == "SGA"|mi$outcome == "BD"|mi$outcome == "NICU")

bar_tz <- ggplot(mi_tz,aes(x = outcome,y = freq,fill = exposure))+
  geom_bar(stat = "identity",position = "dodge")+
  labs(x = "",y = "",title = "C")+
  scale_fill_manual(values = c("#a6cee3","#00468b","#f99a98","#ed0000"))+
  scale_y_continuous(limits = c(0,0.4),breaks = seq(0,0.4,by = 0.2))+
  theme(panel.background = element_blank(),axis.line = element_line(),axis.ticks.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey",linetype = "dashed"),legend.position = "none",
        axis.text.x = element_text(size = 15,color = "black"),
        axis.text.y = element_text(size = 15,color = "black"),panel.grid = element_blank(),
        plot.title = element_text(size = 25,face = "bold"))

bar_jj <- ggplot(mi_jj,aes(x = outcome,y = freq,fill = exposure))+
  geom_bar(stat = "identity",position = "dodge")+
  labs(x = "",y = "",title = "D")+
  scale_fill_manual(values = c("#a6cee3","#00468b","#f99a98","#ed0000"))+
  scale_y_continuous(limits = c(0,0.25),breaks = seq(0,0.25,by = 0.1))+
  theme(panel.background = element_blank(),axis.line = element_line(),axis.ticks.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey",linetype = "dashed"),legend.position = "none",
        axis.text.x = element_text(size = 15,color = "black"),
        axis.text.y = element_text(size = 15,color = "black"),panel.grid = element_blank(),
        plot.title = element_text(size = 25,face = "bold"))

figure3_1 <- ggarrange(box_glu,box_zb,ncol = 2,nrow = 1,widths = c(1.2,2))
figure3_2 <- ggarrange(bar_tz,bar_jj,ncol = 2,nrow = 1,widths = c(1,1.2))
figure3 <- ggarrange(figure3_1,figure3_2,ncol = 1,nrow = 2,heights = c(1.2,1))

pdf("Figure 3.pdf",width = 15,height = 10,family = "sans")
print(figure3)
dev.off()
