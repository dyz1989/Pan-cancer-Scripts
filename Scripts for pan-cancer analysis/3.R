library(ggpubr)          
expFile="UCHL5geneExp.txt"      
cliFile="clinical.txt"      
setwd("File path")     

exp=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
gene=colnames(exp)[1]
exp=exp[(exp[,"Type"]=="Tumor"),]

exp$UCHL5[exp$UCHL5>35]=35

cli=read.table("fustat.txt", header=T, sep="\t", check.names=F, row.names=1)
cliName=colnames(cli)[1]
colnames(cli)=c("Clinical")

sameSample=intersect(row.names(exp), row.names(cli))
exp=exp[sameSample,,drop=F]
cli=cli[sameSample,,drop=F]
data=cbind(exp, cli)

data[,gene]=as.numeric(data[,gene])
p=ggboxplot(data, x="CancerType", y=gene, color="Clinical",
     xlab="",
     ylab=paste0(gene, " expression"),
     legend.title=cliName)
p=p+rotate_x_text(60)
p=p+stat_compare_means(aes(group=Clinical),
      symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
      label = "p.signif")

pdf(file=paste0(cliName,".pdf"), width=7.6, height=5.5)
print(p)
dev.off()




