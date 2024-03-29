library(cmdstanr)
set_cmdstan_path(path = "")
setwd("")

#####################################
## Read a FASTA file
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


########################################################
## data
data <- read.fasta("")
# saveRDS(object=data,file="porB_align.rds")

# check invariant sites
l <- ncol(seqs)
invar <- c()
discard <- c()

for (pos in 1:l){
  count <- as.data.frame(table(seqs[,pos]))
  if (length(count$Var1) == 1){
    invar <- append(invar, pos)
  }
  if (any(count$Var1 == 'NANANA')){
    discard <- append(discard, pos)
  }
}

invar_l <- l-length(invar)-length(discard)
l_new <- c(1:l)[-c(invar,discard)]

X <- matrix(NA, nrow=invar_l, ncol=61)
n <- rep(NA, invar_l)
pi <- rep(0.0163934426229508, 61)

i <- 1
for (pos in l_new){
  count <- as.data.frame(table(seqs[,pos]))
  order <- match(tripletNames_noSTO, count$Var1)
  x <- count$Freq[order]
  x[is.na(x)] <- 0
  X[i,] <- x
  n[i] <- sum(x)
  i <- i+1
}

data_var <- list(X=X, n=n, pi=pi, l=invar_l)

file <- file.path(cmdstan_path(), "porB_var.json")
write_stan_json(data_var, file)