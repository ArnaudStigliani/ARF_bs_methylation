rm(list=ls())
library(Biostrings)
nRegion <- 600
library(stringr)


#-------------------------------------read matrices ------------------------------------------


pfm_ARF<- read.table("m_ARF2.txt",header=TRUE,sep="\t",skip=1)
pfm_ARF <- round((t(as.matrix(pfm_ARF)))*nRegion)+1 ;pfm_ARF
pfm_ARF <- cbind(rep(151,4),pfm_ARF,rep(151,4))
maxi_ARF <- apply(pfm_ARF,FUN=max, 2)
maxi_ARF <- matrix(nrow=4, rep(maxi_ARF,4),byrow=TRUE)
pwm_ARF <- log(pfm_ARF/maxi_ARF)
pwm_ARF2_rev <- pwm_ARF - minScore(pwm_ARF)/dim(pwm_ARF)[2] 
pwm_ARF2 <-  reverseComplement(pwm_ARF2_rev) ; pwm_ARF2

pfm_ARF<- read.table("m_ARF5.txt",header=TRUE,sep="\t",skip=1)
pfm_ARF <- round((t(as.matrix(pfm_ARF)))*nRegion)+1 ;pfm_ARF
pfm_ARF <- cbind(rep(151,4),rep(151,4),pfm_ARF)
maxi_ARF <- apply(pfm_ARF,FUN=max, 2)
maxi_ARF <- matrix(nrow=4, rep(maxi_ARF,4),byrow=TRUE)
pwm_ARF <- log(pfm_ARF/maxi_ARF)
pwm_ARF5_rev <- pwm_ARF - minScore(pwm_ARF)/dim(pwm_ARF)[2] 
pwm_ARF5 <-  reverseComplement(pwm_ARF5_rev) ; pwm_ARF5

#-------------------------------------Loop ------------------------------------------

i <- 0
for (elt in list("ARF5_pos","ARF5_neg","ARF2_pos","ARF2_neg"))
{
    i <- i+1
    if (i==1 || i==2)
    {
        pwm_ARF <- pwm_ARF5
        pwm_ARF_rev<- pwm_ARF5_rev
        print("ARF5")
    }
    else
    {
        pwm_ARF <- pwm_ARF2
        pwm_ARF_rev<- pwm_ARF2_rev
        print("ARF2")
    }   
#-------------------------------------read fasta files-----------------------------------------

    ARF_pos <- readDNAStringSet(paste(elt,".fas",sep=""))

    width_pos <- width(ARF_pos)
    seq_pos <- as.character(ARF_pos)



#-------------------------------------Compute Scores-----------------------------------------
#
    threshold <-  maxScore(pwm_ARF) - 12
#
    match_ARF <- sapply(FUN=matchPWM,seq_pos,pwm=pwm_ARF,min.score=threshold)
    nb_sites <- sapply(match_ARF,length)
    start_sites <- unlist(sapply(match_ARF,start)) -1
#
    match_ARF_rev <- sapply(FUN=matchPWM,seq_pos,pwm=pwm_ARF_rev,min.score=threshold)
    nb_sites_rev <- sapply(match_ARF_rev,length)
    start_sites_rev <- unlist(sapply(match_ARF_rev,start)) -1
#
#

    seq_names<- names(seq_pos)
#
    seq_names_plus <- rep(seq_names,times=nb_sites)
    seq_names_rev <- rep(seq_names,times=nb_sites_rev)
#
    seq_names_rep <- c(seq_names_plus,seq_names_rev)
#
    list_seq_names <- unlist(str_split(seq_names_rep,"[:-]"))
    tab_names <- matrix(list_seq_names,byrow=TRUE,ncol=3)
    tab_names <- cbind(tab_names,c(rep(1,length(seq_names_plus)),rep(-1,length(seq_names_rev))))
    chr_reg <- tab_names[,1]
    start_all_sites <- c(start_sites, start_sites_rev)

#-------------------------------------T1 -----------------------------------------

    offset <- 3
    pos <- ifelse(tab_names[,4] ==1 , start_all_sites + offset, start_all_sites +11 - offset)
    pos_C <- as.integer(tab_names[,2])+pos
    T1 <- data.frame(chr_reg,pos_C,pos_C+1,tab_names[,4])

#-------------------------------------G2 -----------------------------------------

    offset <- 4
    pos <- ifelse(tab_names[,4] ==1 , start_all_sites + offset, start_all_sites +11 - offset)
    pos_C <- as.integer(tab_names[,2])+pos
    G2 <- data.frame(chr_reg,pos_C,pos_C+1,tab_names[,4])

#-------------------------------------T3 -----------------------------------------

    offset <- 5
    pos <- ifelse(tab_names[,4] ==1 , start_all_sites + offset, start_all_sites +11 - offset)
    pos_C <- as.integer(tab_names[,2])+pos
    T3 <- data.frame(chr_reg,pos_C,pos_C+1,tab_names[,4])

#-------------------------------------C4 -----------------------------------------

    offset <- 6
    pos <- ifelse(tab_names[,4] ==1 , start_all_sites + offset, start_all_sites +11 - offset)
    pos_C <- as.integer(tab_names[,2]) +pos
    C4 <- data.frame(chr_reg,pos_C,pos_C+1,tab_names[,4])

#-------------------------------------G5 -----------------------------------------

    offset <- 7
    pos <- ifelse(tab_names[,4] ==1 , start_all_sites + offset, start_all_sites +11 - offset)
    pos_C <- as.integer(tab_names[,2]) +pos
    G5 <- data.frame(chr_reg,pos_C,pos_C+1,tab_names[,4])

#-------------------------------------G6 -----------------------------------------

    offset <- 8
    pos <- ifelse(tab_names[,4] ==1 , start_all_sites + offset, start_all_sites +11 - offset)
    pos_C <- as.integer(tab_names[,2]) +pos
    G6 <- data.frame(chr_reg,pos_C,pos_C+1,tab_names[,4])

#-------------------------------------sites -----------------------------------------

    offset <- 3
    pos <-  start_all_sites + offset
    pos_begin <- as.integer(tab_names[,2]) + pos
    sites <- data.frame(chr_reg,pos_begin,pos_begin + 6)

    options(scipen=999)
    write.table(T1,paste("T1_",elt,".bed",sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
    write.table(G2,paste("G2_",elt,".bed",sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
    write.table(T3,paste("T3_",elt,".bed",sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
    write.table(C4,paste("C4_",elt,".bed",sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
    write.table(G5,paste("G5_",elt,".bed",sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
    write.table(G6,paste("G6_",elt,".bed",sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
    write.table(sites,paste("sites_",elt,".bed",sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
}
