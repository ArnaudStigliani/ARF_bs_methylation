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
pwm_ARF_rev <- pwm_ARF - minScore(pwm_ARF)/dim(pwm_ARF)[2] 
pwm_ARF <-  reverseComplement(pwm_ARF_rev) ; pwm_ARF

#-------------------------------------read fasta files-----------------------------------------
tair10dir <- "~/Data/tair10/tair10.fas"
tair10 <- readDNAStringSet(tair10dir)
CG <- read.table("mCG.bedGraph",sep="\t",header=FALSE,stringsAsFactors=FALSE)
CHG <- read.table("mCHG.bedGraph",sep="\t",header=FALSE,stringsAsFactors=FALSE)
CHH <- read.table("mCHH.bedGraph",sep="\t",header=FALSE,stringsAsFactors=FALSE)
ARF_pos <- readDNAStringSet('ARF2.fas')

width_pos <- width(ARF_pos)
seq_pos <- as.character(ARF_pos)



#-------------------------------------Compute Scores-----------------------------------------
#
scores_ARF_pos<- mapply(seq_pos,FUN=PWMscoreStartingAt,SIMPLIFY=TRUE,  starting.at=mapply(seq,1,width_pos-dim(pwm_ARF)[2],SIMPLIFY=FALSE),MoreArgs=list(pwm=pwm_ARF)) - maxScore(pwm_ARF)
#
scores_ARF_rev_pos<- mapply(seq_pos,FUN=PWMscoreStartingAt,SIMPLIFY=TRUE,  starting.at=mapply(seq,1,width_pos-dim(pwm_ARF)[2],SIMPLIFY=FALSE),MoreArgs=list(pwm=pwm_ARF_rev)) - maxScore(pwm_ARF)
#
max <- apply(FUN=max,scores_ARF_pos,2)
max_rev<- apply(FUN=max,scores_ARF_rev_pos,2)
#
pos_max <- apply(FUN=which.max,scores_ARF_pos,2)
pos_max_rev<- apply(FUN=which.max,scores_ARF_rev_pos,2)
#
seq_names<- names(seq_pos)
list_seq_names <- unlist(str_split(seq_names,"[:-]"))
tab_names <- matrix(list_seq_names,byrow=TRUE,ncol=3)
chr_reg <- tab_names[,1]

#-------------------------------------G2 -----------------------------------------

offset <- 5
pos <- ifelse(max > max_rev , pos_max + offset, pos_max_rev +11 - offset)
pos_C <- as.integer(tab_names[,2])+pos
G2 <- data.frame(chr_reg,pos_C,pos_C+1)

#-------------------------------------C4 -----------------------------------------

offset <- 7
pos <- ifelse(max > max_rev , pos_max + offset, pos_max_rev +11 - offset)
pos_C <- as.integer(tab_names[,2]) + pos
C4 <- data.frame(chr_reg,pos_C,pos_C+1)

#-------------------------------------G5 -----------------------------------------

offset <- 8
pos <- ifelse(max > max_rev , pos_max + offset, pos_max_rev +11 - offset)
pos_C <- as.integer(tab_names[,2]) + pos
G5 <- data.frame(chr_reg,pos_C,pos_C+1)

#-------------------------------------G6 -----------------------------------------

offset <- 9
pos <- ifelse(max > max_rev , pos_max + offset, pos_max_rev +11 - offset)
pos_C <- as.integer(tab_names[,2]) + pos
G6 <- data.frame(chr_reg,pos_C,pos_C+1)


write.table(G2,"G2_ARF2.bed",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
write.table(C4,"C4_ARF2.bed",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
write.table(G5,"G5_ARF2.bed",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
write.table(G6,"G6_ARF2.bed",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)

