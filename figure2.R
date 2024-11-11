rm(list = ls())
library(lattice)
library(caret)
library(dendextend)
library(dendsort)
library(gridBase)
library(ComplexHeatmap)
library(circlize)
setwd("D:\\002001\\1_users\\0_manuscript\\027_wangjinghan_Metabolome Signature in Gestational Diabetes Mellitus is Associated with Adverse Birth Outcomes")

grid <- expand.grid(.alpha = seq(0,1,by = 0.1),
                    .lambda = seq(0,0.1,by = 0.001))
metab_final$BL_HIP_new <- as.factor(metab_final$BL_HIP_new)
colnames(metab_final)
data <- metab_final[,c(46,80:129)]
control <- trainControl(method = "cv",number = 10) 

set.seed(123)
enet.train <- train(BL_HIP_new~.,
                    data = data,
                    method = "glmnet",
                    trControl = control,
                    tuneGrid = grid)
enet.train
#The final values used for the model were alpha = 0.8 and lambda = 0.012
enet.train$bestTune

fit <- glmnet(as.matrix(metab_final[,80:129]),metab_final$BL_HIP_new,family = "binomial",alpha = 0.8,lambda = 0.012)
predict <- predict(fit,newx = as.matrix(metab_final[,80:129]),s = 0.012,type = "response")
roc <- roc(metab_final$BL_HIP_new,as.numeric(predict),ci = TRUE)
auc(roc) #0.7709
threshold <- coords(roc,"best")$threshold

predict_new <- ifelse(predict > threshold,1,0) 
table(predict_new,metab_final$BL_HIP_new)
#            BL_HIP_new
#predict     0       1
#      0     1122    77
#      1     611     240

metab_final$HIP_met <- ifelse(predict_new == "0" & metab_final$BL_HIP_new == "0","0",
                              ifelse(predict_new == "0" & metab_final$BL_HIP_new == "1","1",
                                     ifelse(predict_new == "1" & metab_final$BL_HIP_new == "0","2","3")))
metab_final$HIP_met <- as.factor(metab_final$HIP_met)
table(metab_final$HIP_met)
#0     1     2     3 
#1122  77    611   240 

dt_upset <- metab_final[,c("pid","BL_HIP_new","HIP_met")]
dt_upset$GDM <- as.numeric(ifelse(dt_upset$BL_HIP_new == 1,1,0))
dt_upset$mGDM <- as.numeric(ifelse(dt_upset$HIP_met == 2 | dt_upset$HIP_met == 3,1,0))
dt_upset$HIP_met <- ifelse(dt_upset$HIP_met==0,"non-GDM,non-mGDM",
                           ifelse(dt_upset$HIP_met==1,"GDM,non-mGDM",
                                  ifelse(dt_upset$HIP_met==2,"non-GDM,mGDM","GDM,mGDM")))
dt_upset$HIP_met <- factor(dt_upset$HIP_met,levels = c("non-GDM,non-mGDM","GDM,non-mGDM","non-GDM,mGDM","GDM,mGDM"))

p2 <- ComplexUpset::upset(dt_upset,colnames(dt_upset)[4:5],
                          width_ratio = 0.2,
                          height_ratio = 0.3,
                          min_degree = 0, 
                          sort_intersections_by = "degree", 
                          sort_intersections = "ascending", 
                          sort_sets = "ascending",
                          name = "GDM and Predicted GDM from Metabolome",
                          set_sizes = FALSE,
                          base_annotations = list(
                            "Intersection Size" = intersection_size(
                              counts = FALSE, 
                              mapping = aes(fill = HIP_met))+
                              scale_fill_manual(values = c("#9fbbd5","#4e6691","#d69d98","#b8474d"))+
                              #labs(title = "A")+
                              theme(plot.title = element_text(size = 25,face = "bold"))
                          ))
pdf("Figure 2A_1.pdf",width = 5,height = 4,family = "sans")
print(p2)
dev.off()

p3 <- ggplot(dt_upset,mapping = aes(x = GDM,fill = HIP_met))+
  geom_bar(stat = "count",position = "fill")+
  scale_y_continuous(labels = scales::percent_format())+
  scale_fill_manual(values = c("#9fbbd5","#4e6691","#d69d98","#b8474d"))+
  labs(title = "A")+
  theme_minimal()+
  theme(plot.title = element_text(size = 25,face = "bold"))
pdf("Figure 2A_2.pdf",width = 5,height = 4,family = "sans")
print(p3)
dev.off()  


##PCA
df <- metab_final[,c(76,80:129)]
rownames(df) <- df$tube
df <- df[,-1]

pca <- opls(df)
outcomeFc <- as.factor(metab_final[,"HIP_met"])

dev.new()
pca.plot <- plot(pca,typeVc = "x-score",parAscolFcVn = outcomeFc,parEllipsesl = TRUE)
grid.draw(pca.plot)

##PLS-DA
plsda <- opls(df,outcomeFc,orthoI = 0,predI = 2)
sample_score <- plsda@scoreMN %>%
  as.data.frame() %>%
  mutate(outcome = factor(metab_final$HIP_met))
head(sample_score)
sample_score$outcome <- factor(sample_score$outcome,levels = c(0,1,2,3))

plsda@modelDF[1,"R2X"] #0.254
plsda@modelDF[2,"R2X"] #0.084

plsda.plot <- ggplot(sample_score,aes(p1,p2,color = outcome))+
  geom_hline(yintercept = 0,linetype = "solid", size = 0.5)+
  geom_vline(xintercept = 0,linetype = "solid", size = 0.5)+
  labs(x = "P1 (25.4%)",y = "P2 (8.4%)",title = "B")+
  stat_ellipse(aes(fill = outcome),geom = "polygon",alpha = 0.05,level = 0.95,
               linetype = 1,size = 0.6,show.legend = TRUE)+
  geom_point(size = 0.8)+
  geom_point(aes(-10,-8),color = "white")+ 
  stat_ellipse(aes(fill = outcome),geom = "polygon",alpha = 0,level = 0.95,
               linetype = 1,size = 0.6,show.legend = TRUE)+ 
  scale_color_manual(values = c("#a6cee3","#00468b","#f99a98","#ed0000"))+
  scale_fill_manual(values = c("#a6cee3","#00468b","#f99a98","#ed0000"))+
  theme_bw()+
  theme(panel.background = element_blank(),legend.position = "none",
        axis.title.x = element_text(size = 15,color = "black"),
        axis.title.y = element_text(size = 15,color = "black"),
        plot.title = element_text(size = 25,face = "bold"))

pdf("Figure 2C.pdf",width = 9,height = 6,family = "sans")
print(plsda.plot)
dev.off()

plsda.plot2 <- ggplot(sample_score,aes(p1,p2,color = mGDM))+
  geom_hline(yintercept = 0,linetype = "solid", size = 0.5)+
  geom_vline(xintercept = 0,linetype = "solid", size = 0.5)+
  labs(x = "P1 (25.4%)",y = "P2 (8.4%)",title = "B")+
  stat_ellipse(aes(fill = mGDM),geom = "polygon",alpha = 0.05,level = 0.95,
               linetype = 1,size = 0.6,show.legend = TRUE)+ 
  geom_point(size = 0.8)+
  geom_point(aes(-10,-8),color = "white")+ 
  stat_ellipse(aes(colour = mGDM),geom = "polygon",alpha = 0,level = 0.95,
               linetype = 1,size = 0.6,show.legend = TRUE)+
  scale_color_manual(values = c("#278383","#b2446b"))+
  scale_fill_manual(values = c("#278383","#b2446b"))+
  theme_bw()+
  theme(panel.background = element_blank(),legend.position = "none",
        axis.title.x = element_text(size = 15,color = "black"),
        axis.title.y = element_text(size = 15,color = "black"),
        plot.title = element_text(size = 25,face = "bold"))

pdf("Figure 2B.pdf",width = 9,height = 6,family = "sans")
print(plsda.plot2)
dev.off()

##PERMANOVA
tmp <- metab_final
tmp$x <- as.factor(ifelse(tmp$HIP_met == 0|tmp$HIP_met == 1,0,1))
dataumap <- tmp[,80:129]
for(i in 1:50) 
{
  dataumap[,i] <- 2^(dataumap[,i])
}
dataumap <- vegdist(dataumap)

permanova <- adonis2(dataumap~tmp$x,permutations = 999,distance="euclidean")
summary(permanova)
print(as.data.frame(permanova[1,5])) #0.001

tmp <- subset(metab_final,metab_final$HIP_met == 2|metab_final$HIP_met == 3)
tmp$x <- as.factor(ifelse(tmp$HIP_met == 2,0,1))
dataumap <- tmp[,80:129]
for(i in 1:50) 
{
  dataumap[,i] <- 2^(dataumap[,i])
}
dataumap <- vegdist(dataumap)
permanova <- adonis2(dataumap~tmp$x,permutations = 999,distance="euclidean")
summary(permanova)
print(as.data.frame(permanova[1,5])) #0.01

tmp <- subset(metab_final,metab_final$HIP_met == 0|metab_final$HIP_met == 1)
tmp$x <- as.factor(ifelse(tmp$HIP_met == 0,0,1))
dataumap <- tmp[,80:129]
for(i in 1:50) 
{
  dataumap[,i] <- 2^(dataumap[,i])
}
dataumap <- vegdist(dataumap)
permanova <- adonis2(dataumap~tmp$x,permutations = 999,distance="euclidean")
summary(permanova)
print(as.data.frame(permanova[1,5])) #0.191


tmp1 <- subset(metab_final[,c(80:129,134)],metab_final$HIP_met == "0"|metab_final$HIP_met == "2") #1733
tmp2 <- subset(metab_final[,c(80:129,134)],metab_final$HIP_met == "1"|metab_final$HIP_met == "3") #317
tmp3 <- subset(metab_final[,c(80:129,134)],metab_final$HIP_met == "0"|metab_final$HIP_met == "1") #1199
tmp4 <- subset(metab_final[,c(80:129,134)],metab_final$HIP_met == "2"|metab_final$HIP_met == "3") #851

myVars <- colnames(metab_final)[80:129]
table50 <- CreateTableOne(vars = myVars,strata = "HIP_met",data = metab_final)  
table50 <- print(table50,quote= FALSE,noSpaces = TRUE,showAllLevels = TRUE)

mi_nGDM <- NULL
for (i in 1:50)
{
  test <- t.test(tmp1[,i]~tmp1$HIP_met)
  mi_nGDM <- rbind(mi_nGDM,c(colnames(tmp1)[i],test$estimate,test$p.value))
}

mi_nGDM <- as.data.frame(mi_nGDM)
colnames(mi_nGDM)[1] <- "COMP.ID"
colnames(mi_nGDM)[4] <- "P_nGDM"
mi_nGDM$P_nGDM <- as.numeric(mi_nGDM$P_nGDM)
table(mi_nGDM$P_nGDM < 0.05) #49  
mi_nGDM$FDR_nGDM <- p.adjust(mi_nGDM$P_nGDM,method = "fdr")
table(mi_nGDM$FDR_nGDM < 0.05) #49  
mi_nGDM <- merge(mi_nGDM,annotation[,c(1:4,9)],by = "COMP.ID")
mi_nGDM$bz <- ifelse(mi_nGDM$FDR_nGDM < 0.05,"","nosig")

mi_GDM <- NULL
for (i in 1:50)
{
  test <- t.test(tmp2[,i]~tmp2$HIP_met)
  mi_GDM <- rbind(mi_GDM,c(colnames(tmp2)[i],test$estimate,test$p.value))
}

mi_GDM <- as.data.frame(mi_GDM)
colnames(mi_GDM)[1] <- "COMP.ID"
colnames(mi_GDM)[4] <- "P_GDM"
mi_GDM$P_GDM <- as.numeric(mi_GDM$P_GDM)
table(mi_GDM$P_GDM < 0.05) #44   
mi_GDM$FDR_GDM <- p.adjust(mi_GDM$P_GDM,method = "fdr")
table(mi_GDM$FDR_GDM < 0.05) #44   
mi_GDM <- merge(mi_GDM,annotation[,c(1:4,9)],by = "COMP.ID")
mi_GDM$bz <- ifelse(mi_GDM$FDR_GDM < 0.05,"","nosig")

mi_nmGDM <- NULL
for (i in 1:50)
{
  test <- t.test(tmp3[,i]~tmp3$HIP_met)
  mi_nmGDM <- rbind(mi_nmGDM,c(colnames(tmp3)[i],test$estimate,test$p.value))
}

mi_nmGDM <- as.data.frame(mi_nmGDM)
colnames(mi_nmGDM)[1] <- "COMP.ID"
colnames(mi_nmGDM)[4] <- "P_nmGDM"
mi_nmGDM$P_nmGDM <- as.numeric(mi_nmGDM$P_nmGDM)
table(mi_nmGDM$P_nmGDM < 0.05) #5  
mi_nmGDM$FDR_nmGDM <- p.adjust(mi_nmGDM$P_nmGDM,method = "fdr")
table(mi_nmGDM$FDR_nmGDM < 0.05) #0  
mi_nmGDM <- merge(mi_nmGDM,annotation[,c(1:4,9)],by = "COMP.ID")
mi_nmGDM$bz <- ifelse(mi_nmGDM$FDR_nmGDM < 0.05,"","nosig")

mi_mGDM <- NULL
for (i in 1:50)
{
  test <- t.test(tmp4[,i]~tmp4$HIP_met)
  mi_mGDM <- rbind(mi_mGDM,c(colnames(tmp4)[i],test$estimate,test$p.value))
}

mi_mGDM <- as.data.frame(mi_mGDM)
colnames(mi_mGDM)[1] <- "COMP.ID"
colnames(mi_mGDM)[4] <- "P_mGDM"
mi_mGDM$P_mGDM <- as.numeric(mi_mGDM$P_mGDM)
table(mi_mGDM$P_mGDM < 0.05) #15  
mi_mGDM$FDR_mGDM <- p.adjust(mi_mGDM$P_mGDM,method = "fdr")
table(mi_mGDM$FDR_mGDM < 0.05) #1  
mi_mGDM <- merge(mi_mGDM,annotation[,c(1:4,9)],by = "COMP.ID")
mi_mGDM$bz <- ifelse(mi_mGDM$FDR_mGDM < 0.05,"","nosig")

mi_4group <- cbind(substr(rownames(table50),1,6),table50[,2:5])
mi_4group <- data.frame(mi_4group[-1,])
colnames(mi_4group) <- c("COMP.ID","group0","group1","group2","group3")
mi_4group <- merge(mi_4group,mi_nGDM[,c("COMP.ID","P_nGDM","FDR_nGDM")],by = "COMP.ID")
mi_4group <- merge(mi_4group,mi_GDM[,c("COMP.ID","P_GDM","FDR_GDM")],by = "COMP.ID")
mi_4group <- merge(mi_4group,mi_nmGDM[,c("COMP.ID","P_nmGDM","FDR_nmGDM")],by = "COMP.ID")
mi_4group <- merge(mi_4group,mi_mGDM[,c("COMP.ID","P_mGDM","FDR_mGDM")],by = "COMP.ID")
for (i in 6:13)
{
  mi_4group[,i] <- formatC(mi_4group[,i],format = "e",digits = 1)
}
mi_4group <- merge(mi_4group,annotation[,c(1:4,9)],by = "COMP.ID")
mi_4group <- mi_4group[order(match(mi_4group$COMP.ID,sort50$metabolites)),]  

write.csv(mi_4group,"3_tables\\Table S5.csv")

boxplot <- metab_final[,c(80:129,134)]
boxplot$Group <- factor(ifelse(boxplot$HIP_met == "0","normoglycemic-non-mGDM",ifelse(boxplot$HIP_met == "1","hyperglycemic-non-mGDM",ifelse(boxplot$HIP_met == "2","normoglycemic-mGDM","hyperglycemic-mGDM"))),
                        levels = c("normoglycemic-non-mGDM","normoglycemic-mGDM","hyperglycemic-non-mGDM","hyperglycemic-mGDM"))
boxplot <- boxplot[,-51]

standardize <- function(x) {
  return((x - mean(x)) / sd(x))
}
boxplot[,1:50] <- boxplot[,1:50] %>%
  mutate_all(standardize)

for (i in 1:nrow(annotation))
{
  oldname <- as.character(annotation$COMP.ID[i])
  newname <- as.character(annotation$COMPOUND.Name[i])
  
  if(oldname %in% colnames(boxplot)){
    colnames(boxplot)[colnames(boxplot) == oldname] <- newname
  }
}

boxplot <- boxplot[,order(match(colnames(boxplot),unique(sort50$COMPOUND.Name)))]
pheat0 <- subset(boxplot,boxplot$Group == "normoglycemic-non-mGDM")
pheat1 <- subset(boxplot,boxplot$Group == "hyperglycemic-non-mGDM")
pheat2 <- subset(boxplot,boxplot$Group == "normoglycemic-mGDM")
pheat3 <- subset(boxplot,boxplot$Group == "hyperglycemic-mGDM")
pheat0 <- pheat0[,-51]
pheat1 <- pheat1[,-51]
pheat2 <- pheat2[,-51]
pheat3 <- pheat3[,-51]

pheat <- cbind(apply(pheat0,2,mean),apply(pheat1,2,mean),apply(pheat2,2,mean),apply(pheat3,2,mean))
pheat <- cbind(colMeans(pheat0),colMeans(pheat1),colMeans(pheat2),colMeans(pheat3))
pheat <- as.matrix(pheat)
colnames(pheat) <- c("normoglycemic-non-mGDM","hyperglycemic-non-mGDM","normoglycemic-mGDM","hyperglycemic-mGDM")


mycol1 = circlize::colorRamp2(c(-0.6,0,0.6),c("#42b540","white","#ed0000"))
circos.clear()
pdf("Figure 2D.pdf",width = 25,height = 25,family = "sans")
circos.par(gap.after = c(40)) 
circos.heatmap(pheat,col = mycol1,dend.side = "inside",
               rownames.side = "outside", 
               track.height = 0.13, 
               rownames.col = "black",rownames.cex = 1.8,rownames.font = 0.5,
               cluster = FALSE,
               cell.border = NA,
               dend.track.height = 0.2)
draw(packLegend(Legend(title = "Exp",col_fun = mycol1,direction = c("vertical"))),just = "right")
dev.off()

