#setwd("C:/Users/whale/Desktop/Project")


## set up
require(rstan, warn.conflicts=FALSE, quietly=TRUE)
require(rethinking, warn.conflicts=FALSE, quietly=TRUE)
require(bayesplot, warn.conflicts=FALSE, quietly=TRUE)
require(data.table, warn.conflicts=FALSE, quietly=TRUE)
require(tidyverse, warn.conflicts=FALSE, quietly=TRUE)

library(KScorrect)
#library(genetics)
#library(ape)

### Read a FASTA file
read.fasta <- function(fname, as.char=FALSE) {
  a = scan(fname,what=character(0),sep="\n",quiet=TRUE,na.strings="")
  wh = as.vector(sapply(a,substr,1,1))==">"
  labs = substr(as.character(a[wh]),2,1000);
  
  lseqs = a[!wh]
  nlines = length(lseqs)%/%length(labs)
  n = length(lseqs)%/%nlines
  seqs = rep("",n);
  names(seqs) <- labs
  for(i in 1:n) {
    ibeg = (i-1)*nlines+1
    iend = i*nlines
    seqs[i] = paste(lseqs[ibeg:iend],collapse="")
  }
  seqlen = as.numeric(sapply(seqs,nchar))
  if(length(seqlen)>1 & var(seqlen)>0) {
    warning("Sequences have differing lengths");
    mx = max(seqlen)
    for(i in 1:n) seqs[i] = paste(seqs[i],paste(rep("-",mx-seqlen[i]),collapse=""),sep="")
  }
  L = as.numeric(nchar(seqs[1]))
  
  SEQ = array("-",dim=c(n,L))
  for(i in 1:n) SEQ[i,] = unlist(strsplit(seqs[i],""))
  rownames(SEQ) <- labs;
  if(as.char==TRUE) {
    return(SEQ);
  } else {
    fSEQ = apply(toupper(SEQ),2,factor,levels=c("A","G","C","T"));
    return(fSEQ);
  }
}

### Write a FASTA file
write.fasta <- function(DNA,filename) {
  ofile <- file(filename,"w");
  for(n in 1:nrow(DNA)) {
    writeLines(paste(">",rownames(DNA)[n],sep=""),ofile);
    writeLines(paste(DNA[n,],collapse=""),ofile);    
  }
  close(ofile);
}

### General
totriplet = function(x) {
  L = floor(length(x)/3)*3
  paste(x[seq(1,L,by=3)],x[seq(2,L,by=3)],x[seq(3,L,by=3)],sep="")
}
geneticCode = list(
  "TTT"="Phe","TTC"="Phe","TTA"="Leu","TTG"="Leu",
  "TCT"="Ser","TCC"="Ser","TCA"="Ser","TCG"="Ser",
  "TAT"="Tyr","TAC"="Tyr","TAA"="STO","TAG"="STO",
  "TGT"="Cys","TGC"="Cys","TGA"="STO","TGG"="Trp",
  "CTT"="Leu","CTC"="Leu","CTA"="Leu","CTG"="Leu",
  "CCT"="Pro","CCC"="Pro","CCA"="Pro","CCG"="Pro",
  "CAT"="His","CAC"="His","CAA"="Gln","CAG"="Gln",
  "CGT"="Arg","CGC"="Arg","CGA"="Arg","CGG"="Arg",
  "ATT"="Ile","ATC"="Ile","ATA"="Ile","ATG"="Met",
  "ACT"="Thr","ACC"="Thr","ACA"="Thr","ACG"="Thr",
  "AAT"="Asn","AAC"="Asn","AAA"="Lys","AAG"="Lys",
  "AGT"="Ser","AGC"="Ser","AGA"="Arg","AGG"="Arg",
  "GTT"="Val","GTC"="Val","GTA"="Val","GTG"="Val",
  "GCT"="Ala","GCC"="Ala","GCA"="Ala","GCG"="Ala",
  "GAT"="Asp","GAC"="Asp","GAA"="Glu","GAG"="Glu",
  "GGT"="Gly","GGC"="Gly","GGA"="Gly","GGG"="Gly")
oneLetterCodes = unlist(list("Ala"="A","Arg"="R","Asn"="N","Asp"="D","Cys"="C","Glu"="E","Gln"="Q","Gly"="G","His"="H","Ile"="I","Leu"="L","Lys"="K","Met"="M","Phe"="F","Pro"="P","Ser"="S","Thr"="T","Trp"="W","Tyr"="Y","Val"="V","STO"="X","---"="-"))
aminoAcids = names(table(unlist(geneticCode)))
oneLetterAminoAcids = names(table(unlist(oneLetterCodes)))
tripletNames = names(geneticCode)

transcribe = function(x) {
  y = t(sapply(1:nrow(x),function(i) totriplet(x[i,])))
  rownames(y) = rownames(x)
  return(y)
}
translate = function(x,oneLetter=FALSE) {
  x = toupper(x)
  tr = t(apply(x,1,function(y)sapply(y,function(i) {aa=geneticCode[[i]];ifelse(is.null(aa),"---",aa)} )))
  if(oneLetter) tr = t(apply(tr,1,function(y) oneLetterCodes[y]))
  rownames(tr) = rownames(x)
  return(tr)
}
view.nucleotide = function(x) {
  image(0:ncol(x),0:nrow(x),t(matrix(as.numeric(factor(x,levels=c("-","A","G","C","T"))),nrow=nrow(x))),col=c("white","red","green","yellow","blue"))
}
view.codon = function(x) {
  image(0:ncol(x),0:nrow(x),t(matrix(as.numeric(factor(x),levels=tripletNames),nrow=nrow(x))),col=rainbow(20))
}
view.protein = function(x,oneLetter=FALSE) {
  levs = aminoAcids
  if(oneLetter) levs = oneLetterAminoAcids
  cols = rainbow(20)
  if(oneLetter) cols = c("white",cols)
  image(0:ncol(x),0:nrow(x),t(matrix(as.numeric(factor(x,levels=levs)),nrow=nrow(x))),col=cols)
}
# Assumes a fasta file representing a single genome, possibly split across contigs
read.fasta.ref = function(ref_file) {
  r = scan(ref_file,what=character(0),sep="\n")
  beg = substr(r,1,1)
  gd = beg!=">"
  rcat = paste(r[gd],collapse="")
  return(toupper(unlist(strsplit(rcat,""))))
}
# Assumes a fasta file representing a single genome, possibly split across contigs
read.fasta.ref.contig = function(ref_file) {
  r = scan(ref_file,what=character(0),sep="\n")
  beg = substr(r,1,1)
  gd = beg!=">"
  contig = rep(cumsum(!gd)[gd],times=nchar(r[gd]))
  return(contig)
}


## data
data <- read.fasta("porB3.carriage.noindels.txt")
seqs <- transcribe(data)
tri <- translate(seqs)

sort(table(seqs))
tripletNames_noSTO <- tripletNames[-c(11, 12, 15)]

## substitution rate matrix 
PDRM <- function(mu, kappa, omega){
  M <- matrix(0, nrow=61, ncol=61)
  
  M[1,2] <- kappa*mu;
  M[1,3] <- omega*mu;
  M[1,4] <- omega*mu;
  M[1,5] <- kappa*omega*mu;
  M[1,9] <- omega*mu;
  M[1,11] <- omega*mu;
  M[1,14] <- kappa*omega*mu;
  M[1,30] <- omega*mu;
  M[1,46] <- omega*mu;
  M[2,3] <- omega*mu;
  M[2,4] <- omega*mu;
  M[2,6] <- kappa*omega*mu;
  M[2,10] <- omega*mu;
  M[2,12] <- omega*mu;
  M[2,15] <- kappa*omega*mu;
  M[2,31] <- omega*mu;
  M[2,47] <- omega*mu;
  M[3,4] <- kappa*mu;
  M[3,7] <- kappa*omega*mu;
  M[3,16] <- kappa*mu;
  M[3,32] <- omega*mu;
  M[3,48] <- omega*mu;
  M[4,8] <- kappa*omega*mu;
  M[4,13] <- omega*mu;
  M[4,17] <- kappa*mu;
  M[4,33] <- omega*mu;
  M[4,49] <- omega*mu;
  M[5,6] <- kappa*mu;
  M[5,7] <- mu;
  M[5,8] <- mu;
  M[5,9] <- omega*mu;
  M[5,11] <- omega*mu;
  M[5,18] <- kappa*omega*mu;
  M[5,34] <- omega*mu;
  M[5,50] <- omega*mu;
  M[6,7] <- mu;
  M[6,8] <- mu;
  M[6,10] <- omega*mu;
  M[6,12] <- omega*mu;
  M[6,19] <- kappa*omega*mu;
  M[6,35] <- omega*mu;
  M[6,51] <- omega*mu;
  M[7,8] <- kappa*mu;
  M[7,20] <- kappa*omega*mu;
  M[7,36] <- omega*mu;
  M[7,52] <- omega*mu;
  M[8,13] <- omega*mu;
  M[8,21] <- kappa*omega*mu;
  M[8,37] <- omega*mu;
  M[8,53] <- omega*mu;
  M[9,10] <- kappa*mu;
  M[9,11] <- kappa*omega*mu;
  M[9,22] <- kappa*omega*mu;
  M[9,38] <- omega*mu;
  M[9,54] <- omega*mu;
  M[10,12] <- kappa*omega*mu;
  M[10,23] <- kappa*omega*mu;
  M[10,39] <- omega*mu;
  M[10,55] <- omega*mu;
  M[11,12] <- kappa*mu;
  M[11,13] <- omega*mu;
  M[11,26] <- kappa*omega*mu;
  M[11,42] <- omega*mu;
  M[11,58] <- omega*mu;
  M[12,13] <- omega*mu;
  M[12,27] <- kappa*omega*mu;
  M[12,43] <- omega*mu;
  M[12,59] <- omega*mu;
  M[13,29] <- kappa*omega*mu;
  M[13,45] <- omega*mu;
  M[13,61] <- omega*mu;
  M[14,15] <- kappa*mu;
  M[14,16] <- mu;
  M[14,17] <- mu;
  M[14,18] <- kappa*omega*mu;
  M[14,22] <- omega*mu;
  M[14,26] <- omega*mu;
  M[14,30] <- omega*mu;
  M[14,46] <- omega*mu;
  M[15,16] <- mu;
  M[15,17] <- mu;
  M[15,19] <- kappa*omega*mu;
  M[15,23] <- omega*mu;
  M[15,27] <- omega*mu;
  M[15,31] <- omega*mu;
  M[15,47] <- omega*mu;
  M[16,17] <- kappa*mu;
  M[16,20] <- kappa*omega*mu;
  M[16,24] <- omega*mu;
  M[16,28] <- omega*mu;
  M[16,32] <- omega*mu;
  M[16,48] <- omega*mu;
  M[17,21] <- kappa*omega*mu;
  M[17,25] <- omega*mu;
  M[17,29] <- omega*mu;
  M[17,33] <- omega*mu;
  M[17,49] <- omega*mu;
  M[18,19] <- kappa*mu;
  M[18,20] <- mu;
  M[18,21] <- mu;
  M[18,22] <- omega*mu;
  M[18,26] <- omega*mu;
  M[18,34] <- omega*mu;
  M[18,50] <- omega*mu;
  M[19,20] <- mu;
  M[19,21] <- mu;
  M[19,23] <- omega*mu;
  M[19,27] <- omega*mu;
  M[19,35] <- omega*mu;
  M[19,51] <- omega*mu;
  M[20,21] <- kappa*mu;
  M[20,24] <- omega*mu;
  M[20,28] <- omega*mu;
  M[20,36] <- omega*mu;
  M[20,52] <- omega*mu;
  M[21,25] <- omega*mu;
  M[21,29] <- omega*mu;
  M[21,37] <- omega*mu;
  M[21,53] <- omega*mu;
  M[22,23] <- kappa*mu;
  M[22,24] <- omega*mu;
  M[22,25] <- omega*mu;
  M[22,26] <- kappa*omega*mu;
  M[22,38] <- omega*mu;
  M[22,54] <- omega*mu;
  M[23,24] <- omega*mu;
  M[23,25] <- omega*mu;
  M[23,27] <- kappa*omega*mu;
  M[23,39] <- omega*mu;
  M[23,55] <- omega*mu;
  M[24,25] <- kappa*mu;
  M[24,28] <- kappa*omega*mu;
  M[24,40] <- omega*mu;
  M[24,56] <- omega*mu;
  M[25,29] <- kappa*omega*mu;
  M[25,41] <- omega*mu;
  M[25,57] <- omega*mu;
  M[26,27] <- kappa*mu;
  M[26,28] <- mu;
  M[26,29] <- mu;
  M[26,42] <- omega*mu;
  M[26,58] <- omega*mu;
  M[27,28] <- mu;
  M[27,29] <- mu;
  M[27,43] <- omega*mu;
  M[27,59] <- omega*mu;
  M[28,29] <- kappa*mu;
  M[28,44] <- mu;
  M[28,60] <- omega*mu;
  M[29,45] <- mu;
  M[29,61] <- omega*mu;
  M[30,31] <- kappa*mu;
  M[30,32] <- mu;
  M[30,33] <- omega*mu;
  M[30,34] <- kappa*omega*mu;
  M[30,38] <- omega*mu;
  M[30,42] <- omega*mu;
  M[30,46] <- kappa*omega*mu;
  M[31,32] <- mu;
  M[31,33] <- omega*mu;
  M[31,35] <- kappa*omega*mu;
  M[31,39] <- omega*mu;
  M[31,43] <- omega*mu;
  M[31,47] <- kappa*omega*mu;
  M[32,33] <- kappa*omega*mu;
  M[32,36] <- kappa*omega*mu;
  M[32,40] <- omega*mu;
  M[32,44] <- omega*mu;
  M[32,48] <- kappa*omega*mu;
  M[33,37] <- kappa*omega*mu;
  M[33,41] <- omega*mu;
  M[33,45] <- omega*mu;
  M[33,49] <- kappa*omega*mu;
  M[34,35] <- kappa*mu;
  M[34,36] <- mu;
  M[34,37] <- mu;
  M[34,38] <- omega*mu;
  M[34,42] <- omega*mu;
  M[34,50] <- kappa*omega*mu;
  M[35,36] <- mu;
  M[35,37] <- mu;
  M[35,39] <- omega*mu;
  M[35,43] <- omega*mu;
  M[35,51] <- kappa*omega*mu;
  M[36,37] <- kappa*mu;
  M[36,40] <- omega*mu;
  M[36,44] <- omega*mu;
  M[36,52] <- kappa*omega*mu;
  M[37,41] <- omega*mu;
  M[37,45] <- omega*mu;
  M[37,53] <- kappa*omega*mu;
  M[38,39] <- kappa*mu;
  M[38,40] <- omega*mu;
  M[38,41] <- omega*mu;
  M[38,42] <- kappa*omega*mu;
  M[38,54] <- kappa*omega*mu;
  M[39,40] <- omega*mu;
  M[39,41] <- omega*mu;
  M[39,43] <- kappa*omega*mu;
  M[39,55] <- kappa*omega*mu;
  M[40,41] <- kappa*mu;
  M[40,44] <- kappa*omega*mu;
  M[40,56] <- kappa*omega*mu;
  M[41,45] <- kappa*omega*mu;
  M[41,57] <- kappa*omega*mu;
  M[42,43] <- kappa*mu;
  M[42,44] <- omega*mu;
  M[42,45] <- omega*mu;
  M[42,58] <- kappa*omega*mu;
  M[43,44] <- omega*mu;
  M[43,45] <- omega*mu;
  M[43,59] <- kappa*omega*mu;
  M[44,45] <- kappa*mu;
  M[44,60] <- kappa*omega*mu;
  M[45,61] <- kappa*omega*mu;
  M[46,47] <- kappa*mu;
  M[46,48] <- mu;
  M[46,49] <- mu;
  M[46,50] <- kappa*omega*mu;
  M[46,54] <- omega*mu;
  M[46,58] <- omega*mu;
  M[47,48] <- mu;
  M[47,49] <- mu;
  M[47,51] <- kappa*omega*mu;
  M[47,55] <- omega*mu;
  M[47,59] <- omega*mu;
  M[48,49] <- kappa*mu;
  M[48,52] <- kappa*omega*mu;
  M[48,56] <- omega*mu;
  M[48,60] <- omega*mu;
  M[49,53] <- kappa*omega*mu;
  M[49,57] <- omega*mu;
  M[49,61] <- omega*mu;
  M[50,51] <- kappa*mu;
  M[50,52] <- mu;
  M[50,53] <- mu;
  M[50,54] <- omega*mu;
  M[50,58] <- omega*mu;
  M[51,52] <- mu;
  M[51,53] <- mu;
  M[51,55] <- omega*mu;
  M[51,59] <- omega*mu;
  M[52,53] <- kappa*mu;
  M[52,56] <- omega*mu;
  M[52,60] <- omega*mu;
  M[53,57] <- omega*mu;
  M[53,61] <- omega*mu;
  M[54,55] <- kappa*mu;
  M[54,56] <- omega*mu;
  M[54,57] <- omega*mu;
  M[54,58] <- kappa*omega*mu;
  M[55,56] <- omega*mu;
  M[55,57] <- omega*mu;
  M[55,59] <- kappa*omega*mu;
  M[56,57] <- kappa*mu;
  M[56,60] <- kappa*omega*mu;
  M[57,61] <- kappa*omega*mu;
  M[58,59] <- kappa*mu;
  M[58,60] <- mu;
  M[58,61] <- mu;
  M[59,60] <- mu;
  M[59,61] <- mu;
  M[60,61] <- kappa*mu;
  
  #Fill in the lower triangle
  for (i in 1:61){
    for (j in (i+1):61){
      M[j,i]=M[i,j];
    }
  }
  
  #Compute the diagonal
  for (i in 1:61){
    rowsum <- rowSums(M[,1:61])
    M[i,i] <- -rowsum[i]
  }
  
}

## constant dN/dS ratio per gene
m1.1 <- map2stan(
  alist(
    # likelihood
    
    
    
    
    pi <- rep(0.0163934426229508, 61),
    # improper priors for theta and kappa: log uniform distribution
    theta <- dlunif(1,0,Inf),
    kappa <- dlunif(1,0,Inf),
    # prior for omega
    for (i in 1:61){
      omega <- dexp(1)
      
    }
    
  ),
  data=seqs, 
  start=list(theta=0.17, kappa=1), iter=1e4, warmup=1e3, chains=2)
