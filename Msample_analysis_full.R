
#Macrophage data series - analysis 

#macro_1 agilent data set 
library(limma)
setwd("/Users/new/Downloads/Labs/Helikar_Lab/Macrophage-CNS/Analyzed data/Macro_1/GSE88768_RAW/")
macro_1 <- readTargets("Macro_1")
macro1_x = read.maimages(macro_1, source="agilent",green.only=TRUE)
macro1_y = backgroundCorrect(macro1_x, method="normexp", offset = 16)
macro1_y = normalizeBetweenArrays(macro1_y, method="quantile")
macro1_neg95 <- apply(macro1_y$E[macro1_y$genes$ControlType==-1,],2,function(macro1_x) quantile(macro1_x,probs = 0.95)) # this denotes the control type that you are using and the quantile within which you are designating. This current iteration uses a negative control and looks at the top five percentile The control types that you can have are [-1,1]#identifying which rows are negative controls and establishes quantile
macro1_cutoff <- matrix(1.4*macro1_neg95,nrow(macro1_y),ncol(macro1_y),byrow=TRUE)
macro1_isexpr = rowSums(macro1_y$E>macro1_cutoff)>=3
table(macro1_isexpr)
macro1_y0 <- macro1_y[macro1_y$genes$ControlType==0 & macro1_isexpr,]
write.table(macro1_y0, file="eset1_PMA.xls", quote=F, col.names = NA, sep="\t")

#affy analysis 2 data set
library(affy)
setwd("/Users/new/Downloads/Labs/Helikar_Lab/Macrophage-CNS/Analyzed data/Macro_2/GSE9874_RAW/")
macrophage_2 <-ReadAffy()
eset_2 <- mas5(macrophage_2)
eset2_PMA <- mas5calls(macrophage_2)
write.exprs(eset_2, file="eset_macro_2.txt")
macro2_x <- data.frame(exprs(eset_2), exprs(eset2_PMA), assayDataElement(eset2_PMA, "se.exprs"))
macro2_x <- macro2_x[,sort(names(macro2_x))];
write.table(macro2_x, file="eset2_PMA.xls", quote=F, col.names = NA, sep="\t")

#affy analysis 3 data set 
setwd("/Users/new/Downloads/Labs/Helikar_Lab/Macrophage-CNS/Analyzed data/Macro_3/GSE11886_RAW/")
macrophage_3 <- ReadAffy() 
eset_3 <- mas5(macrophage_3) 
eset3_PMA <- mas5calls(macrophage_3) 
write.exprs(eset_3, file="eset_macro_3_PMA.txt") 
macro3_x <- data.frame(exprs(eset_3), exprs(eset3_PMA), assayDataElement(eset3_PMA, "se.exprs"))
macro3_x <- macro3_x[,sort(names(macro3_x))];
write.table(macro3_x, file="eset3_PMA.xls", quote=F, col.names = NA, sep="\t") 

#affy analysis 4 data set 
setwd("/Users/new/Downloads/Labs/Helikar_Lab/Macrophage-CNS/Analyzed data/Macro_4/GSE13670_RAW/")
macrophage_4 <- ReadAffy() 
eset_4 <- mas5(macrophage_4)
eset4_PMA <- mas5calls(macrophage_4)
write.exprs(eset_4, file="eset_macro_4.txt")
macro4_x <- data.frame(exprs(eset_4), exprs(eset4_PMA), assayDataElement(eset4_PMA, "se.exprs"))
macro4_x <- macro4_x[,sort(names(macro4_x))];
write.table(macro4_x, file="eset4_PMA.xls", quote=F, col.names = NA, sep="\t") 
#reruning up til here 
#macrophage 5 data 
setwd("/Users/new/Downloads/Labs/Helikar_Lab/Macrophage-CNS/Analyzed data/Macro_5/GSE22528_RAW/")
macro_5 <- ReadAffy() 
eset_5 <- mas5(macro_5)
eset5_PMA <- mas5calls(macro_5)
write.exprs(eset_5, file="eset_macro_5.txt")
macro5_x <- data.frame(exprs(eset_5), exprs(eset5_PMA), assayDataElement(eset5_PMA, "se.exprs"))
macro5_x <- macro5_x[,sort(names(macro5_x))];
write.table(macro5_x, file="eset5_PMA.xls", quote=F, col.names = NA, sep="\t") 

#macrophage 6 agilent data
library(limma)
setwd("/Users/new/Downloads/Labs/Helikar_Lab/Macrophage-CNS/Analyzed data/Macro_6/GSE84622_RAW/")
macro_6 <- readTargets("macrophage_6")
macro6_x = read.maimages(macro_6, source="agilent",green.only=TRUE)
macro6_y = backgroundCorrect(macro6_x, method="normexp", offset = 16)
macro6_y = normalizeBetweenArrays(macro6_y, method="quantile")
macro6_neg95 <- apply(macro6_y$E[macro6_y$genes$ControlType==-1,],2,function(macro6_x) quantile(macro6_x,probs = 0.95)) # this denotes the control type that you are using and the quantile within which you are designating. This current iteration uses a negative control and looks at the top five percentile The control types that you can have are [-1,1]#identifying which rows are negative controls and establishes quantile
macro6_cutoff <- matrix(1.4*macro6_neg95,nrow(macro6_y),ncol(macro6_y),byrow=TRUE)
macro6_isexpr = rowSums(macro6_y$E>macro6_cutoff)>=2
table(macro6_isexpr)
macro6_y0 <- macro6_y[macro6_y$genes$ControlType==0 & macro6_isexpr,]
write.table(macro6_y0, file="eset6_PMA.xls", quote=F, col.names = NA, sep="\t") 

#macro7
library(affy)
#set working directory 
setwd("/Users/new/Downloads/Labs/Helikar_Lab/Macrophage-CNS/Analyzed data/Macro_7/GSE20484_RAW/")
macro_7 <- ReadAffy()
eset_7 <- mas5(macro_7)
eset7_PMA <-mas5calls(macro_7)
write.exprs(eset_7, file="eset_macro_7_PMA.txt") 
macro7_x <- data.frame(exprs(eset_7), exprs(eset7_PMA), assayDataElement(eset7_PMA, "se.exprs"))
macro7_x <- macro7_x[,sort(names(macro7_x))];
write.table(macro7_x, file="eset7_PMA.xls", quote=F, col.names = NA, sep="\t") 

#macro8 #didnt do this one
library(affy)
setwd("/Users/new/Downloads/Labs/Helikar_Lab/Macrophage-CNS/Analyzed data/Macro_8/GSE36995_RAW/")
macro_8 <- ReadAffy()
eset_8 <- mas5(macro_8)
eset8_PMA <-mas5calls(macro_8)
write.exprs(eset_8, file="eset_macro_8_PMA.txt") 
macro8_x <- data.frame(exprs(eset_8), exprs(eset8_PMA), assayDataElement(eset8_PMA, "se.exprs"))
macro8_x <- macro8_x[,sort(names(macro8_x))];
write.table(macro8_x, file="eset8_PMA.xls", quote=F, col.names = NA, sep="\t") 


#macro9
setwd("/Users/new/Downloads/Labs/Helikar_Lab/Macrophage-CNS/Analyzed data/Macro_9/GSE2125_RAW/")
macro_9 <- ReadAffy()
eset_9 <- mas5(macro_9)
eset9_PMA <-mas5calls(macro_9)
write.exprs(eset_9, file="eset_macro_9_PMA.txt") 
macro9_x <- data.frame(exprs(eset_9), exprs(eset9_PMA), assayDataElement(eset9_PMA, "se.exprs"))
macro9_x <- macro9_x[,sort(names(macro9_x))];
write.table(macro9_x, file="eset9_PMA.xls", quote=F, col.names = NA, sep="\t") 

#macro10
library(affy)
#set working directory 
setwd("/Users/new/Downloads/Labs/Helikar_Lab/Macrophage-CNS/Analyzed data/Macro_10/GSE3212_RAW/")
macro_10 <- ReadAffy()
eset_10 <- mas5(macro_10)
eset10_PMA <-mas5calls(macro_10)
write.exprs(eset_10, file="eset_macro_10_PMA.txt") 
macro10_x <- data.frame(exprs(eset_10), exprs(eset10_PMA), assayDataElement(eset10_PMA, "se.exprs"))
macro10_x <- macro10_x[,sort(names(macro10_x))];
write.table(macro10_x, file="eset10_PMA.xls", quote=F, col.names = NA, sep="\t") 

#macro11
library(affy)
#set working directory 
setwd("/Users/new/Downloads/Labs/Helikar_Lab/Macrophage-CNS/Analyzed data/Macro_11/GSE5099_RAW/")
macro_11 <- ReadAffy()
eset_11 <- mas5(macro_11)
eset11_PMA <-mas5calls(macro_11)
write.exprs(eset_11, file="eset_macro_11_PMA.txt") 
macro11_x <- data.frame(exprs(eset_11), exprs(eset11_PMA), assayDataElement(eset11_PMA, "se.exprs"))
macro11_x <- macro11_x[,sort(names(macro11_x))];
write.table(macro11_x, file="eset11_PMA.xls", quote=F, col.names = NA, sep="\t") 

#macro12
library(affy)
#set working directory 
setwd("/Users/new/Downloads/Labs/Helikar_Lab/Macrophage-CNS/Analyzed data/Macro_12/GSE8515_RAW/")
macro_12 <- ReadAffy()
eset_12 <- mas5(macro_12)
eset12_PMA <-mas5calls(macro_12)
write.exprs(eset_12, file="eset_macro_12_PMA.txt") 
macro12_x <- data.frame(exprs(eset_12), exprs(eset12_PMA), assayDataElement(eset12_PMA, "se.exprs"))
macro12_x <- macro12_x[,sort(names(macro12_x))];
write.table(macro12_x, file="eset12_PMA.xls", quote=F, col.names = NA, sep="\t") 

