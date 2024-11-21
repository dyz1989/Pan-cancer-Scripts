library(plyr)
library(reshape2)
library(ggpubr)
inputFile="UCHL5geneExp.txt"      
setwd("File path")     

data=read.table(inputFile, header=T, sep="\t", check.names=F)
gene=colnames(data)[2]
colnames(data)[2]="expression"

data$expression[data$expression>45]=45

p=ggboxplot(data, x="CancerType", y="expression", color="Type",
     xlab="",
     ylab=paste0(gene," expression"),
     width=0.6,
     palette = c("blue","red") )
p=p+rotate_x_text(50)
p=p+stat_compare_means(aes(group=Type),
      method="wilcox.test",
      symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
      label = "p.signif")
pdf(file="diff.pdf", width=8, height=5.5)
print(p)
dev.off()


data=data[(data[,"Type"]=="Tumor"),]
med=ddply(data,"CancerType",summarise,med=median(expression))
data$CancerType=factor(data$CancerType, levels=med[order(med[,"med"],decreasing = T), "CancerType"])

p=ggboxplot(data, x="CancerType", y="expression", fill="CancerType",
       xlab="",
       ylab=paste0(gene," expression"),
       width=0.6,
       #palette=rainbow(length(levels(factor(data$CancerType)))),
       legend="")
pdf(file="boxplot.pdf", width=8, height=5.5)    
p+rotate_x_text(50)
dev.off()