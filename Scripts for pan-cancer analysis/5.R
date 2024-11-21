inputFile="expTime.txt"        
setwd("File path")      
rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)   
rt$futime=rt$futime/365
gene=colnames(rt)[3]

outTab=data.frame()
for(i in levels(factor(rt[,"CancerType"]))){
	rt1=rt[(rt[,"CancerType"]==i),]
	
	cox=coxph(Surv(futime, fustat) ~ rt1[,gene], data = rt1)
	coxSummary = summary(cox)
    coxP=coxSummary$coefficients[,"Pr(>|z|)"]
	outTab=rbind(outTab,
	             cbind(cancer=i,
	                   HR=coxSummary$conf.int[,"exp(coef)"],
	                   HR.95L=coxSummary$conf.int[,"lower .95"],
	                   HR.95H=coxSummary$conf.int[,"upper .95"],
			           pvalue=coxP) )
	
	group=ifelse(rt1[,gene]>median(rt1[,gene]), "high", "low")
	diff=survdiff(Surv(futime, fustat) ~ group, data=rt1)
	pValue=1-pchisq(diff$chisq,df=1)
	if(pValue<0.05){
		if(pValue<0.001){
			pValue="p<0.001"
		}else{
			pValue=paste0("p=",sprintf("%.03f",pValue))
		}
		
		fit=survfit(Surv(futime, fustat) ~ group, data = rt1)
		surPlot=ggsurvplot(fit, 
				    data=rt1,
				    title=paste0("Cancer: ",i),
				    pval=pValue,
				    pval.size=6,
				    conf.int=F,
				    legend.title=paste0(gene," levels"),
				    legend.labs=c("high","low"),
				    font.legend=12,
				    fontsize=4,
				    xlab="Time(years)",
				    ylab="Overall survival",
				    break.time.by = 2,
				    palette=c("red","blue"),
				    risk.table=TRUE,
				    risk.table.title="",
				    risk.table.height=.25)
		pdf(file=paste0("survival.",i,".pdf"), width=6, height=5, onefile=FALSE)
		print(surPlot)
		dev.off()
	}
}
write.table(outTab, file="cox.result.txt", sep="\t", row.names=F, quote=F)

bioForest=function(coxFile=null, forestFile=null, forestCol=null){
	rt=read.table(coxFile,header=T,sep="\t",row.names=1,check.names=F)
	data=as.matrix(rt)
	HR=data[,1:3]
	hr=sprintf("%.3f",HR[,"HR"])
	hrLow=sprintf("%.3f",HR[,"HR.95L"])
	hrHigh=sprintf("%.3f",HR[,"HR.95H"])
	pVal=data[,"pvalue"]
	pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
	clrs=fpColors(box=forestCol, line="darkblue", summary="royalblue")
	tabletext <- 
		list(c(NA, rownames(HR)),
		    append("pvalue", pVal),
		    append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) )
	pdf(file=forestFile, width=9, height=6, onefile=FALSE)
	forestplot(tabletext, 
	           rbind(rep(NA, 3), HR),
	           col=clrs,
	           graphwidth=unit(50, "mm"),
	           xlog=T,
	           lwd.ci=4,
	           boxsize=0.6,
	           title="Overall survival",
	           xlab="Hazard ratio",
	           txt_gp=fpTxtGp(ticks=gpar(cex=1.1), xlab=gpar(cex = 1.25))
	           )
	
}

bioForest(coxFile="cox.result.txt", forestFile="forest.pdf", forestCol="red")
dev.off()