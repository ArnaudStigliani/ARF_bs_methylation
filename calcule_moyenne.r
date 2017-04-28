rm(list=ls())

means <- NULL
z <- 0
for (elt in list("ARF5_pos.bed","ARF5_neg.bed","ARF2_pos.bed","ARF2_neg.bed"))
{
    z <- z+1
    T1_CHH <- abs(read.table(paste("T1_CHH_met_",elt,sep=""),header=TRUE,sep="\t")[,4])
    T1_CHG <- abs(read.table(paste("T1_CHG_met_",elt,sep=""),header=TRUE,sep="\t")[,4])
    T1_CG <- abs(read.table(paste("T1_CG_met_",elt,sep=""),header=TRUE,sep="\t")[,4])
                                        #
    G2_CHH <- abs(read.table(paste("G2_CHH_met_",elt,sep=""),header=TRUE,sep="\t")[,4])
    G2_CHG <- abs(read.table(paste("G2_CHG_met_",elt,sep=""),header=TRUE,sep="\t")[,4])
    G2_CG <- abs(read.table(paste("G2_CG_met_",elt,sep=""),header=TRUE,sep="\t")[,4])
                                        #
    T3_CHH <- abs(read.table(paste("T3_CHH_met_",elt,sep=""),header=TRUE,sep="\t")[,4])
    T3_CHG <- abs(read.table(paste("T3_CHG_met_",elt,sep=""),header=TRUE,sep="\t")[,4])
    T3_CG <- abs(read.table(paste("T3_CG_met_",elt,sep=""),header=TRUE,sep="\t")[,4])
                                        #
    C4_CHH <- abs(read.table(paste("C4_CHH_met_",elt,sep=""),header=TRUE,sep="\t")[,4])
    C4_CHG <- abs(read.table(paste("C4_CHG_met_",elt,sep=""),header=TRUE,sep="\t")[,4])
    C4_CG <- abs(read.table(paste("C4_CG_met_",elt,sep=""),header=TRUE,sep="\t")[,4])
                                        #
    G5_CHH <- abs(read.table(paste("G5_CHH_met_",elt,sep=""),header=TRUE,sep="\t")[,4])
    G5_CHG <- abs(read.table(paste("G5_CHG_met_",elt,sep=""),header=TRUE,sep="\t")[,4])
    G5_CG <- abs(read.table(paste("G5_CG_met_",elt,sep=""),header=TRUE,sep="\t")[,4])
                                        #
    G6_CHH <- abs(read.table(paste("G6_CHH_met_",elt,sep=""),header=TRUE,sep="\t")[,4])
    G6_CHG <- abs(read.table(paste("G6_CHG_met_",elt,sep=""),header=TRUE,sep="\t")[,4])
    G6_CG <- abs(read.table(paste("G6_CG_met_",elt,sep=""),header=TRUE,sep="\t")[,4])
                                        #
                                        #
    a <- mean(T1_CHH)
    b <- mean(T1_CHG)
    c <- mean(T1_CG)
                                        #
    d <- mean(G2_CHH)
    e <- mean(G2_CHG)
    f <- mean(G2_CG)
                                        #   
    g <- mean(T3_CHH)
    h <- mean(T3_CHG)
    i <- mean(T3_CG)
                                        #
    j <- mean(C4_CHH)
    k <- mean(C4_CHG)
    l <- mean(C4_CG)
                                        #
    m <- mean(G5_CHH)
    n <- mean(G5_CHG)
    o <- mean(G5_CG)
                                        #
    p <- mean(G6_CHH)
    q <- mean(G6_CHG)
    r <- mean(G6_CG)
                                        #
    all_means_ARF<- c(a,d,g,j,m,p,b,e,h,k,n,q,c,f,i,l,o,r)
    all_means_ARF<- as.table(all_means_ARF)
    names(all_means_ARF) <- rep(c("1(T)","2(G)","3(T)","4(C)","5(G)","6(G)"),3) 
    means <- cbind(means,all_means_ARF)
}

#--------------------------ARF2-------------------------
all_means_ARF5 <- NULL
all_means_ARF2 <- NULL
for (z in 1:18)
{
    all_means_ARF5 <- c(all_means_ARF5,means[z,1],means[z,2])
    all_means_ARF2 <- c(all_means_ARF2,means[z,3],means[z,4])
}
names(all_means_ARF5) <- rep(c(1,1,2,2,3,3,4,4,5,5,6,6),3)
names(all_means_ARF2) <- rep(c(1,1,2,2,3,3,4,4,5,5,6,6),3)

png("ARF2_ARF5_methylation_all_positions.png",width=1600,height=1000)
#
par(lwd=2)
par(mfrow=c(1,2))
barplot(all_means_ARF2,density=c(NA,5),space=c(1,0.1,1,0.1),width=0.3,col=c(rep("cornflowerblue",12),rep("coral3",12),rep("cyan4",12)),lend="butt",main="Study of methylation patern on ARF2 binding sites ",xlab="position in TGTCGG",ylim=c(0,0.35),ylab="percentage of methylated cytosine on average",cex.lab=1.5,cex.main=2)
legend('topleft',legend=c("CHH","CHG","CG","negative"),fill=c("cornflowerblue","coral3","cyan4","black"),pch=c(NA,NA,NA),density=c(NA,NA,NA,5),cex=1.5)
#arrows(x0=(0:20),y0=meanDRt+sdDRt,y1=meanDRt-sdDRt,lwd=1, angle=90,length=0.1,code=3)
barplot(all_means_ARF5,density=c(NA,5),space=c(1,0.1,1,0.1),width=0.3,col=c(rep("cornflowerblue",12),rep("coral3",12),rep("cyan4",12)),lend="butt",main="Study of methylation patern on ARF5 binding sites ",xlab="position in TGTCGG",ylim=c(0,0.35),ylab="percentage of methylated cytosine on average",cex.lab=1.5,cex.main=2)
legend('topleft',legend=c("CHH","CHG","CG","negative"),fill=c("cornflowerblue","coral3","cyan4","black"),pch=c(NA,NA,NA),density=c(NA,NA,NA,5),cex=1.5)
#arrows(x0=(0:20),y0=meanDRt+sdDRt,y1=meanDRt-sdDRt,lwd=1, angle=90,length=0.1,code=3)                 #
#
dev.off()





