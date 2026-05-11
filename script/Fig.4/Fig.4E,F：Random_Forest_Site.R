#library package
library(randomForest)
library(ggplot2)
library(RColorBrewer)

# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)


#import
metadata <- read.table("../../data/Metadata.txt", sep = "\t",row.names =1 ,header = T, na.strings = "")
metadata$Site <- factor(metadata$Site,levels=c("Village","Factory"))

#Relative abundance
level5 <- as.data.frame(t(read.csv("../../data/level5.csv",row.name=1,check.names = F)))
level5 <- level5[,1:ncol(level5)] / rowSums(level5[,1:ncol(level5)], na=T) 
rowSums(level5)

#除去在 Family 分类水平上分类名为unculture、unknow_Unknown_Family 分类信息的分类群
level5_sub <- level5[, grepl("f__", colnames(level5)) & !grepl("uncultured|Unknown_Family", colnames(level5))]

#Bacterial families with abundances below 0.01% were removed
family_abund_filtered <- level5_sub[,colSums(level5_sub) >= 0.0001]
family_ra <- family_abund_filtered 
rowSums(family_ra)

###Root Site ####
metadata_Root <- metadata[metadata$Area_Microbiome == "Root",]
metadata_Root$Site <- factor(metadata_Root$Site,levels = c("Factory","Village"))

family_Root = family_ra[rownames(family_ra) %in% rownames(metadata_Root),] 


#Random forest regression analysis
set.seed(123)
rf_family_R_F = randomForest(family_Root, metadata_Root$Site, importance=TRUE, proximity=TRUE, ntree = 1000)  #,mtry = 25
rf_family_R_F

rf_family_R_F_importance  <- as.data.frame(rf_family_R_F$importance)
rf_family_R_F_importanceSD  <- as.data.frame(rf_family_R_F$importanceSD)

imp_family_R_Site <- as.data.frame(cbind(rf_family_R_F$importance,rf_family_R_F$importanceSD[,colnames(rf_family_R_F$importanceSD) == "MeanDecreaseAccuracy"]))
names(imp_family_R_Site)[ncol(imp_family_R_Site)] <- "SD"
imp_family_R_Site <- imp_family_R_Site[order(imp_family_R_Site$MeanDecreaseAccuracy, decreasing = T),]

write.csv(imp_family_R_Site,file = "randomForest_importance_family_Root_Site.csv")


##10-fold cross-validation of five replicates 
set.seed(123)
rf_family_R_F_cv <- replicate(5, rfcv(family_Root, metadata_Root$Site, cv.fold = 10, step = 1.1), simplify = FALSE)

rf_family_R_F_cv<- data.frame(sapply(rf_family_R_F_cv, '[[', 'error.cv'))
rf_family_R_F_cv$Feature <- rownames(rf_family_R_F_cv)
rf_family_R_F_cv<- reshape2::melt(rf_family_R_F_cv, id = 'Feature')
rf_family_R_F_cv$Feature <- as.numeric(as.character(rf_family_R_F_cv$Feature))
rf_family_R_F_cv.mean <- aggregate(rf_family_R_F_cv$value, by = list(rf_family_R_F_cv$Feature), FUN = mean)
colnames(rf_family_R_F_cv.mean) <- c("Feature.num","cv.mean.error")
head(rf_family_R_F_cv.mean, 10)

write.table(rf_family_R_F_cv.mean,"rf_family_Root_Site_cv.mean.txt",sep = "\t",row.names = F,quote = FALSE)


###Soil Site ####
metadata_Soil <- metadata[metadata$Area_Microbiome == "Soil",]
metadata_Soil$Site <- factor(metadata_Soil$Site,levels = c("Factory","Village"))

family_Soil = family_ra[rownames(family_ra) %in% rownames(metadata_Soil),] 


#Random forest regression
set.seed(123)
rf_family_S_F = randomForest(family_Soil, metadata_Soil$Site, importance=TRUE, proximity=TRUE, ntree = 1000)  #,mtry = 25
rf_family_S_F

rf_family_S_F_importance  <- as.data.frame(rf_family_S_F$importance)
rf_family_S_F_importanceSD  <- as.data.frame(rf_family_S_F$importanceSD)

imp_family_S_Site <- as.data.frame(cbind(rf_family_S_F$importance,rf_family_S_F$importanceSD[,colnames(rf_family_S_F$importanceSD) == "MeanDecreaseAccuracy"]))
names(imp_family_S_Site)[ncol(imp_family_S_Site)] <- "SD"
imp_family_S_Site <- imp_family_S_Site[order(imp_family_S_Site$MeanDecreaseAccuracy, decreasing = T),]

write.csv(imp_family_S_Site,file = "randomForest_importance_family_Soil_Site.csv")


##cross-validation
set.seed(123) 
rf_family_S_F_cv <- replicate(5, rfcv(family_Soil, metadata_Soil$Site, cv.fold = 10, step = 1.1), simplify = FALSE)

rf_family_S_F_cv<- data.frame(sapply(rf_family_S_F_cv, '[[', 'error.cv'))
rf_family_S_F_cv$Feature <- rownames(rf_family_S_F_cv)
rf_family_S_F_cv<- reshape2::melt(rf_family_S_F_cv, id = 'Feature')
rf_family_S_F_cv$Feature <- as.numeric(as.character(rf_family_S_F_cv$Feature))
rf_family_S_F_cv.mean <- aggregate(rf_family_S_F_cv$value, by = list(rf_family_S_F_cv$Feature), FUN = mean)
colnames(rf_family_S_F_cv.mean) <- c("Feature.num","cv.mean.error")
head(rf_family_S_F_cv.mean, 10)

write.table(rf_family_S_F_cv.mean,"rf_family_Soil_Site_cv.mean.txt",sep = "\t",row.names = F,quote = FALSE)








