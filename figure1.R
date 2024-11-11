rm(list = ls())
setwd("D:\\002001\\1_users\\0_manuscript\\027_wangjinghan_Metabolome Signature in Gestational Diabetes Mellitus is Associated with Adverse Birth Outcomes")

mi <- NULL
for (i in 82:786)
{
  fit <- summary(glm(BL_HIP_new~metab_t2_nj[,i]+age+bmi_f,data = metab_t2_nj,family = binomial))
  mi <- rbind(mi,cbind(colnames(metab_t2_nj)[i],coef(fit)[2,1],coef(fit)[2,2],as.numeric(coef(fit)[2,4])))
}

mi <- as.data.frame(mi)
colnames(mi) <- c("metabolites","beta","se","P")
mi$beta <- as.numeric(mi$beta)
mi$P <- as.numeric(mi$P)

mi$qvalue <- p.adjust(mi$P,method = "fdr")
mi <- merge(mi,annotation[,c(1:4,9)],by.x = "metabolites",by.y = "COMP.ID")
table(mi$qvalue < 0.05)
variable_t2 <- subset(mi,mi$qvalue < 0.05) #93 
rownames(variable_t2) <- variable_t2[,1]
tmp <- metab_t2_nj[,which(colnames(metab_t2_nj) %in% variable_t2$metabolites)]
tmp <- as.matrix(tmp)

set.seed(2000) 
cv_model <- cv.glmnet(tmp,metab_t2_nj$BL_HIP_new,family = "binomial",nlambda = 1000,alpha = 0.5,nfolds = 10)
#         Lambda     Index    Measure    SE          Nonzero
# min     0.00601    353      0.6760     0.03655     53
# 1se     0.05151    120      0.7124     0.02885     14
best_lambda <- cv_model$lambda.min
R2 <- cv_model$glmnet.fit$dev.ratio[match(best_lambda,cv_model$lambda)] 

fit <- glmnet(tmp,metab_t2_nj$BL_HIP_new,family = "binomial",alpha = 0.5,lambda = best_lambda)

coefficient <- data.frame(coef(fit,s = fit$lambda.min)@x)
active_index <- which(as.numeric(coef(fit)) != 0)
active_coefficient <- as.numeric(coef(fit))[active_index]
variable_comb <- data.frame(rownames(coef(fit))[active_index])
variable_comb <- cbind(variable_comb,coefficient)
colnames(variable_comb) <- c("Metabolites","coef")
variable_comb <- variable_comb[-1,]

variable_comb <- merge(annotation[,c(1:4,9)],variable_comb,by.y = "Metabolites",by.x = "COMP.ID")
rownames(variable_comb) <- variable_comb[,1]
rm(coefficient,active_index,active_coefficient) 

variable_t2_t3 <- variable_comb %>%
  arrange(desc(abs(variable_comb$coef))) %>%
  group_by(rank = row_number()) 
variable_t2_t3 <- data.frame(variable_t2_t3)
variable_t2_t3$direct <- ifelse(variable_t2_t3$coef > 0,"↑","↓")
variable_t2_t3$rank <- paste0(variable_t2_t3$direct," (",variable_t2_t3$rank,")")
rownames(variable_t2_t3) <- variable_t2_t3$COMP.ID
dt <- variable_t2_t3 %>%
  arrange(SUPER.META.PATHWAY)

metab_nj_final <- metab_t2_nj[,which(colnames(metab_t2_nj) %in% variable_t2_t3$COMP.ID)]
metab_nj_final <- cbind(metab_t2_nj[,1:79],metab_nj_final)

metab_sz_final <- metab_t2_sz[,which(colnames(metab_t2_sz) %in% variable_t2_t3$COMP.ID)]
metab_sz_final <- cbind(metab_t2_sz[,1:79],metab_sz_final)

mi <- NULL
for (i in 80:132)
{
  fit <- summary(glm(BL_HIP_new~metab_sz_final[,i]+age+bmi_f,data = metab_sz_final,family = binomial))
  mi <- rbind(mi,cbind(colnames(metab_sz_final)[i],coef(fit)[2,1],coef(fit)[2,2],as.numeric(coef(fit)[2,4])))
}
mi <- as.data.frame(mi)
colnames(mi) <- c("metabolites","beta","se","P")
mi$beta <- as.numeric(mi$beta)
mi$P <- as.numeric(mi$P)

mi <- merge(mi,variable_t2,by = "metabolites")
mi$bz <- ifelse(mi$beta.x > 0 & mi$beta.y > 0,NA,ifelse(mi$beta.x < 0 & mi$beta.y < 0,NA,"inconsistent"))
table(mi$bz)
exclude_byz <- subset(mi$metabolites,mi$bz == "inconsistent")

nj <- mi[,c(1,5,6)]
colnames(nj) <- c("id","beta","se")
nj$center <- "Nanjing"
sz <- mi[,c(1,2,3)]
colnames(sz) <- c("id","beta","se")
sz$center <- "Suzhou"

meta <- rbind(nj,sz)
meta$se <- as.numeric(meta$se)

pval <- NULL
for (i in 1:53)
{
  tmp <- meta[c(i,i+53),]
  meta_re <- metagen(tmp$beta,tmp$se,sm = "β",data = tmp,byvar = center,print.byvar = FALSE)
  pval <- rbind(pval,unlist(meta_re["pval.Q.b.random"]))
}

pval <- data.frame(pval)
rownames(pval) <- nj$id
pval$metabolites <- rownames(pval)
mii <- merge(mi,pval,by = "metabolites") 
mii <- subset(mii,is.na(mii$bz) == TRUE) #50

mii <- mii %>%
  arrange(desc(abs(mii$beta.y))) %>%
  group_by(rank = row_number()) 
mii <- mii %>%
  arrange(SUPER.META.PATHWAY)

mii$se.x <- as.numeric(mii$se.x)
mii$se.y <- as.numeric(mii$se.y)
mii$`OR(95%CI)_SZ` <- paste0(sprintf("%.2f",exp(mii$beta.x))," (",sprintf("%.2f",exp(mii$beta.x-1.96*mii$se.x)),",",sprintf("%.2f",exp(mii$beta.x+1.96*mii$se.x)),")")
mii$`OR(95%CI)_NJ` <- paste0(sprintf("%.2f",exp(mii$beta.y))," (",sprintf("%.2f",exp(mii$beta.y-1.96*mii$se.y)),",",sprintf("%.2f",exp(mii$beta.y+1.96*mii$se.y)),")")
mii$direct <- ifelse(mii$beta.y > 0,"↑","↓")
mii$rank <- paste0(mii$direct," (",mii$rank,")")

metab_nj_final <- metab_nj_final[,!(colnames(metab_nj_final) %in% exclude_byz)] #1307
metab_sz_final <- metab_sz_final[,!(colnames(metab_sz_final) %in% exclude_byz)] #743

sort_comb <- mii
sort_comb <- sort_comb %>% 
  arrange(SUPER.META.PATHWAY,desc(beta.y)) %>%
  group_by(SUPER.META.PATHWAY)
sort_comb$rank <- c(1:50)

sort_comb <- mii
sort_comb <- sort_comb %>% 
  arrange(desc(beta.y))
sort_comb$rank <- c(1:50)

pdf("Figure 1.pdf",width = 10,height = 7.5,family = "sans")
ggplot(data = sort_comb,mapping = aes(x = reorder(COMPOUND.Name,-rank),y = beta.y,fill = SUPER.META.PATHWAY))+
  geom_bar(stat = "identity",width = 0.75)+
  coord_flip()+
  scale_x_discrete(limits = rev(levels(sort_comb$COMPOUND.Name)))+
  theme_bw()+
  xlab("Compound Name")+
  ylab("Loading")+
  labs(fill = "Metabolite Type")+
  scale_fill_manual(values = c("#00468b","#42b540","#0099b4","#ed0000","#fdaf91","#925e9f"))
dev.off()

