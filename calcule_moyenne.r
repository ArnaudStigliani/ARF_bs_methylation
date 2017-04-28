rm(list=ls())
means <- NULL
z <- 0
mean_wo_sites <- NULL
eff <- NULL
for (elt in list("ARF5_pos.bed","ARF5_neg.bed","ARF2_pos.bed","ARF2_neg.bed"))
{
    z <- z+1
    stack_C <- NULL
    stack_G <- NULL
    eff_C <- NULL
    eff_G <- NULL
    for (letter in list("T1_C_all_met_","G2_C_all_met_","T3_C_all_met_","C4_C_all_met_","G5_C_all_met_","G6_C_all_met_"))
    {
        C_all <- read.table(paste(letter,elt,sep=""),header=TRUE,sep="\t")
        C_all <- C_all[,4] * C_all[,8]
        stack_C <- c(stack_C , mean(C_all[C_all >= 0]))
        stack_G <- c(stack_G , abs(mean(C_all[C_all <= 0])))
        eff_C <- c(eff_C , sum(C_all >= 0))
        eff_G <- c(eff_G , sum(C_all <= 0))
    }
                                        #
    reg_wo_sites_C_all <- read.table(paste("reg_wo_sites_C_all_met_",elt,sep=""),header=TRUE,sep="\t")
    len_wo_sites <- reg_wo_sites_C_all[,3] - reg_wo_sites_C_all[,2]
    wo_sites <- rep(reg_wo_sites_C_all[,4] , times = len_wo_sites)
    mean_wo_sites <- c(mean_wo_sites,mean(abs(wo_sites)))
    #
    all_means_ARF <- cbind(stack_C,stack_G)
    rownames(all_means_ARF) <- c("1(T)","2(G)","3(T)","4(C)","5(G)","6(G)")
    means <- cbind(means,all_means_ARF)
    all_eff <- cbind(eff_C,eff_G)
    rownames(all_eff) <- c("1(T)","2(G)","3(T)","4(C)","5(G)","6(G)")
    eff <- cbind(eff,all_eff)
}

colnames(means) <- c("ARF5_pos","ARF5_pos","ARF5_neg","ARF5_neg","ARF2_pos","ARF2_pos","ARF2_neg","ARF2_neg")
colnames(eff)   <- c("ARF5_pos","ARF5_pos","ARF5_neg","ARF5_neg","ARF2_pos","ARF2_pos","ARF2_neg","ARF2_neg")

#--------------------------ARF2-------------------------
all_means_ARF5 <- NULL
all_means_ARF2 <- NULL
for (z in 1:6)
{
    all_means_ARF5 <- c(all_means_ARF5,means[z,1],means[z,3],means[z,2],means[z,4])
    all_means_ARF2 <- c(all_means_ARF2,means[z,5],means[z,7],means[z,6],means[z,8])
}
all_means_ARF5 <- c(all_means_ARF5 , mean_wo_sites[1:2])
all_means_ARF2 <- c(all_means_ARF2 , mean_wo_sites[3:4])
    
names(all_means_ARF5) <- c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,6,6,6,6,"no_site","no_site")
names(all_means_ARF2) <- c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,6,6,6,6,"no_site","no_site")

png("ARF2_ARF5_methylation_all_positions.png",width=1600,height=1000)
#

par(lwd=2)
par(mfrow=c(1,2))
barplot(all_means_ARF2,density=c(NA,5),space=c(2,0.1,0.5,0.1),width=0.3,col=c(rep(c("cornflowerblue","cornflowerblue","cyan4","cyan4"),times=6),"coral3","coral3"),lend="butt",main="Study of methylation patern on ARF2 binding sites ",xlab="position in TGTCGG",ylim=c(0,0.15),ylab="percentage of methylated cytosine on average",cex.lab=1.5,cex.main=2)
legend('topright',legend=c("C","G","Regions with bs removed","negative"),fill=c("cornflowerblue","cyan4","coral3","black"),pch=NA,density=c(NA,NA,NA,10),cex=1.5)
#
barplot(all_means_ARF5,density=c(NA,5),space=c(2,0.1,0.5,0.1),width=0.3,col=c(rep(c("cornflowerblue","cornflowerblue","cyan4","cyan4"),times=6),"coral3","coral3"),lend="butt",main="Study of methylation patern on ARF5 binding sites ",xlab="position in TGTCGG",ylim=c(0,0.15),ylab="percentage of methylated cytosine on average",cex.lab=1.5,cex.main=2)
legend('topright',legend=c("C","G","Regions with bs removed","negative"),fill=c("cornflowerblue","cyan4","coral3","black"),pch=NA,density=c(NA,NA,NA,10),cex=1.5)
#

dev.off()





