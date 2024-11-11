rm(list = ls())
library(ggplot2)
library(berryFunctions)
library(forestploter)
setwd("D:\\002001\\1_users\\0_manuscript\\027_wangjinghan_Metabolome Signature in Gestational Diabetes Mellitus is Associated with Adverse Birth Outcomes")

myVars <- c("PTB","LGA","SGA","any_Structural_Oneyear","nicu")
catVars <- c("PTB","LGA","SGA","any_Structural_Oneyear","nicu")

table <- CreateTableOne(vars = myVars,strata = "HIP_met",data = metab_final,factorVars = catVars)  
table <- print(table,quote= FALSE,noSpaces = TRUE,showAllLevels = TRUE)

tmp1 <- NULL
for (i in c(53,64,65,47,56))
{
  fit <- summary(glm(metab_final[,i]~HIP_met,data = metab_final,family = binomial))
  tmp1 <- rbind(tmp1,cbind(cbind(paste0(sprintf("%.2f",exp(coef(fit)[2,1]))," (",sprintf("%.2f",exp(coef(fit)[2,1]-1.96*coef(fit)[2,2])),",",sprintf("%.2f",exp(coef(fit)[2,1]+1.96*coef(fit)[2,2])),")"),sprintf("%.3f",coef(fit)[2,4])),
                           cbind(paste0(sprintf("%.2f",exp(coef(fit)[3,1]))," (",sprintf("%.2f",exp(coef(fit)[3,1]-1.96*coef(fit)[3,2])),",",sprintf("%.2f",exp(coef(fit)[3,1]+1.96*coef(fit)[3,2])),")"),sprintf("%.3f",coef(fit)[3,4])),
                           cbind(paste0(sprintf("%.2f",exp(coef(fit)[4,1]))," (",sprintf("%.2f",exp(coef(fit)[4,1]-1.96*coef(fit)[4,2])),",",sprintf("%.2f",exp(coef(fit)[4,1]+1.96*coef(fit)[4,2])),")"),sprintf("%.3f",coef(fit)[4,4]))))
}

tmp2 <- NULL
for (i in c(53,64,65,47,56))
{
  fit <- summary(glm(metab_final[,i]~HIP_met+age+bmi_f+education+alcohol+smoke+income+residence+parity+hypertension+folate,data = metab_final,family = binomial))
  tmp2 <- rbind(tmp2,cbind(cbind(paste0(sprintf("%.2f",exp(coef(fit)[2,1]))," (",sprintf("%.2f",exp(coef(fit)[2,1]-1.96*coef(fit)[2,2])),",",sprintf("%.2f",exp(coef(fit)[2,1]+1.96*coef(fit)[2,2])),")"),sprintf("%.3f",coef(fit)[2,4])),
                           cbind(paste0(sprintf("%.2f",exp(coef(fit)[3,1]))," (",sprintf("%.2f",exp(coef(fit)[3,1]-1.96*coef(fit)[3,2])),",",sprintf("%.2f",exp(coef(fit)[3,1]+1.96*coef(fit)[3,2])),")"),sprintf("%.3f",coef(fit)[3,4])),
                           cbind(paste0(sprintf("%.2f",exp(coef(fit)[4,1]))," (",sprintf("%.2f",exp(coef(fit)[4,1]-1.96*coef(fit)[4,2])),",",sprintf("%.2f",exp(coef(fit)[4,1]+1.96*coef(fit)[4,2])),")"),sprintf("%.3f",coef(fit)[4,4]))))
}

mm <- data.frame(cbind(tmp1,tmp2))
colnames(mm) <- c("RRcrude1","Pcrude1","RRcrude2","Pcrude2","RRcrude3","Pcrude3","RRadjusted1","Padjusted1","RRadjusted2","Padjusted2","RRadjusted3","Padjusted3")
rownames(mm) <- colnames(metab_final)[c(53,64,65,47,56)]
write.csv(mm,"Figure 4 (RR 95%CI).csv")


ptrend <- metab_final
ptrend$HIP_met <- as.numeric(ptrend$HIP_met)
table(ptrend$HIP_met)

tmp1 <- NULL
for (i in c(53,64,65,47,56))
{
  fit <- summary(glm(ptrend[,i]~HIP_met,data = ptrend,family = binomial))
  tmp1 <- rbind(tmp1,cbind(paste0(sprintf("%.2f",exp(coef(fit)[2,1]))," (",sprintf("%.2f",exp(coef(fit)[2,1]-1.96*coef(fit)[2,2])),",",sprintf("%.2f",exp(coef(fit)[2,1]+1.96*coef(fit)[2,2])),")"),sprintf("%.3f",coef(fit)[2,4])))
}

tmp2 <- NULL
for (i in c(53,64,65,47,56))
{
  fit <- summary(glm(ptrend[,i]~HIP_met+age+bmi_f+education+alcohol+smoke+income+residence+parity+hypertension+folate,data = ptrend,family = binomial))
  tmp2 <- rbind(tmp2,cbind(paste0(sprintf("%.2f",exp(coef(fit)[2,1]))," (",sprintf("%.2f",exp(coef(fit)[2,1]-1.96*coef(fit)[2,2])),",",sprintf("%.2f",exp(coef(fit)[2,1]+1.96*coef(fit)[2,2])),")"),sprintf("%.3f",coef(fit)[2,4])))
}

mm6 <- data.frame(cbind(tmp1,tmp2))
colnames(mm6) <- c("RRcrude","Pcrude","RRadjusted","Padjusted")
rownames(mm6) <- colnames(ptrend)[c(53,64,65,47,56)]
write.csv(mm6,"Figure 4 (Ptrend).csv")

tmp5 <- NULL
for (i in c(53,64,65,47,56))
{
  fit <- summary(glm(metab_final[,i]~HIP_met,data = metab_final,family = binomial))
  tmp5 <- rbind(tmp5,rbind(cbind(paste0(sprintf("%.2f",exp(coef(fit)[2,1]))," (",sprintf("%.2f",exp(coef(fit)[2,1]-1.96*coef(fit)[2,2])),",",sprintf("%.2f",exp(coef(fit)[2,1]+1.96*coef(fit)[2,2])),")"),sprintf("%.3f",coef(fit)[2,4])),
                           cbind(paste0(sprintf("%.2f",exp(coef(fit)[3,1]))," (",sprintf("%.2f",exp(coef(fit)[3,1]-1.96*coef(fit)[3,2])),",",sprintf("%.2f",exp(coef(fit)[3,1]+1.96*coef(fit)[3,2])),")"),sprintf("%.3f",coef(fit)[3,4])),
                           cbind(paste0(sprintf("%.2f",exp(coef(fit)[4,1]))," (",sprintf("%.2f",exp(coef(fit)[4,1]-1.96*coef(fit)[4,2])),",",sprintf("%.2f",exp(coef(fit)[4,1]+1.96*coef(fit)[4,2])),")"),sprintf("%.3f",coef(fit)[4,4]))))
}

tmp6 <- NULL
for (i in c(53,64,65,47,56))
{
  fit <- summary(glm(metab_final[,i]~HIP_met+age+bmi_f+education+alcohol+smoke+income+residence+parity+hypertension+folate,data = metab_final,family = binomial))
  tmp6 <- rbind(tmp6,rbind(cbind(exp(coef(fit)[2,1]),coef(fit)[2,2],exp(coef(fit)[2,1]-1.96*coef(fit)[2,2]),exp(coef(fit)[2,1]+1.96*coef(fit)[2,2]),paste0(sprintf("%.2f",exp(coef(fit)[2,1]))," (",sprintf("%.2f",exp(coef(fit)[2,1]-1.96*coef(fit)[2,2])),",",sprintf("%.2f",exp(coef(fit)[2,1]+1.96*coef(fit)[2,2])),")"),sprintf("%.3f",coef(fit)[2,4])),
                           cbind(exp(coef(fit)[3,1]),coef(fit)[3,2],exp(coef(fit)[3,1]-1.96*coef(fit)[3,2]),exp(coef(fit)[3,1]+1.96*coef(fit)[3,2]),paste0(sprintf("%.2f",exp(coef(fit)[3,1]))," (",sprintf("%.2f",exp(coef(fit)[3,1]-1.96*coef(fit)[3,2])),",",sprintf("%.2f",exp(coef(fit)[3,1]+1.96*coef(fit)[3,2])),")"),sprintf("%.3f",coef(fit)[3,4])),
                           cbind(exp(coef(fit)[4,1]),coef(fit)[4,2],exp(coef(fit)[4,1]-1.96*coef(fit)[4,2]),exp(coef(fit)[4,1]+1.96*coef(fit)[4,2]),paste0(sprintf("%.2f",exp(coef(fit)[4,1]))," (",sprintf("%.2f",exp(coef(fit)[4,1]-1.96*coef(fit)[4,2])),",",sprintf("%.2f",exp(coef(fit)[4,1]+1.96*coef(fit)[4,2])),")"),sprintf("%.3f",coef(fit)[4,4]))))
}

plot <- as.data.frame(cbind(tmp5,tmp6))
colnames(plot) <- c("cRR (95%CI)","cP","RR","se","LCI","UCI","aRR (95%CI)","aP")
plot$` ` <- paste(rep(" ",25),collapse = "") 
plot$group <- rep("")

for (i in 3:6) 
{
  plot[,i] <- as.numeric(plot[,i])
}

write.csv(plot,"2_databases\\plot.csv",row.names = FALSE)
plot <- read.csv("2_databases\\plot_update.csv")
colnames(plot) <- c("","case (%)","cRR (95%CI)","cP","RR","se","LCI","UCI","aRR (95%CI)","aP","P-trend")
plot[,1] <- str_pad(plot[,1],width = 30,side = "right",pad = " ")
plot$`cRR (95%CI)` <- str_pad(plot$`cRR (95%CI)`,width = 20,side = "both",pad = " ")
plot$cP <- str_pad(plot$cP,width = 10,side = "both",pad = " ")
plot$`aRR (95%CI)` <- str_pad(plot$`aRR (95%CI)`,width = 20,side = "both",pad = " ")
plot$aP <- str_pad(plot$aP,width = 10,side = "both",pad = " ")
plot$` ` <- paste(rep(" ",35),collapse = "") 

tm <- forest_theme(
  base_size = 8,
  ci_pch = 19,                                     
  ci_col = c("#6495ED"),                           
  ci_alpha = 1,                                  
  ci_lty = 1,                                      
  ci_lwd = 2,                                    
  ci_Theight = 0.2,
  refline_lwd = 1,                                 
  refline_lty = "dashed",
  refline_col = "grey20",
  footnote_cex = 1,
  footnote_col = "red"
)

g <- forest(plot[,c(1:4,9:10,12,11)],
            est = plot$RR,
            lower = plot$LCI,
            upper = plot$UCI,
            #size = plot$se,
            ci_column = c(7), 
            ref_line = 1,
            xlim = c(0,4),
            theme = tm)

g <- edit_plot(g,row = c(3,8,13,18,23),col = 7,
               which = "ci",gp = gpar(col = "#00468b"))
g <- edit_plot(g,row = c(4,9,14,19,24),col = 7,
               which = "ci",gp = gpar(col = "#f99a98"))
g <- edit_plot(g,row = c(5,10,15,20,25),col = 7,
               which = "ci",gp = gpar(col = "#ed0000"))
g <- edit_plot(g,row = c(1:25),
               which = "background",gp = gpar(fill = "white"))
g <- edit_plot(g,row = c(1,6,11,16,21),
               which = "background",gp = gpar(fill = "black",alpha = 0.1))
g <- edit_plot(g,row = c(1,6,11,16,21),
               gp = gpar(fontface = "bold"))
g <- insert_text(g,text = c("Crude Model","Adjusted Model"),col = c(3,5),
                 part = "header",gp = gpar(fontface = "bold"))

pdf("Figure 4.pdf",width = 16,height = 15,family = "sans")
g
dev.off()

