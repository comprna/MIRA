#Read deepbind output file and create significant outputs
library("Biostrings")


# infile <- "/projects_rg/babita/TCGA/mutation/mut_pipeline/subs/subs_new/main_tables/splice_sites/deepbind/3SS.deepbind.out" 
# #seqfile <- "/projects_rg/babita/TCGA/mutation/mut_pipeline/subs/subs_new/main_tables/splice_sites/deepbind/3SS.deepbind.seqs"
# #motiffile <- "/projects_rg/babita/TCGA/mutation/mut_pipeline/subs/subs_new/main_tables/splice_sites/deepbind/deepbind.ids.detail"
# #outfile <- "~/Desktop/tmp"

args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]
seqfile <- args[2]
motiffile <- args[3]
outfile <- args[4]

if(length(args)< 4 )stop("incorrect parameters!! usage: script.R <1.deep.output> <2.deep.seq.input> <3.deep.motif.detail> <4.output_file>")

file.create(outfile)
cat("seq_id\tseq\tdeepbind_ids\tRBP_ids\tdeepbind_scores\n",file=outfile)

#read deepbind output
deep.df <- read.delim(infile)


#read deepbind input file
fastaFile <- readDNAStringSet(seqfile)
seq_name = names(fastaFile)
sequence = paste(fastaFile)
seq.df <- data.frame(seq_name, sequence)


#read deepbind annotation file
motif.df <- read.delim(motiffile)


#run for each row of deep.df, significant value (positive?)
cnt <<- 1

apply(deep.df, 1, function(x){
  
  fa.seq.row <- x[2:length(x)]
  #print(row)
  kmer.seq <- x[1]
  fa.seq.df <- as.data.frame(fa.seq.row)

  fa.seq.df$fa.seq.row <- as.numeric(as.character(fa.seq.df$fa.seq.row))
  sig.df <- subset(fa.seq.df,fa.seq.row > 0)
  sig.df$ids <- rownames(sig.df)
 
  
  #sort significant df by value and pick top three.
  sig.df <- sig.df[rev(order(sig.df$fa.seq.row)),]
  sig.df <-  sig.df[!duplicated(sig.df),]
  
  #add annotation ids
  sig.df$rbps.annot <- sapply(sig.df$ids, function(x){
    rbps.sig.db <- subset(motif.df,ID == as.character(x ))
    return(unique(rbps.sig.db$Protein))
  }) 
  

  #remove columns with same RBP annotation
  sig.df = subset(sig.df, !duplicated(rbps.annot))
  
  #----------------------------------------------------#
  #old
  #rbps.annot <- head(sig.df$rbps.annot, 3) 
  
  #pick top three ids
  #top.rbps.deep <- head(sig.df$ids, 3)
  
  #pick top three scores
  #top.scores <- head(sig.df$fa.seq.row,3)

  #top deep-rbps ids annotated
  #sig.rbps.annot <- paste(unlist(rbps.annot), collapse="|")
  #top.rbps.deep <- paste(top.rbps.deep, collapse="|")
  #top.scores <- paste(top.scores, collapse="|")
  #-----------------------------------------------------#
  
  ##changed 26 april
  #pick top 1 ids now
  rbps.annot <- head(sig.df$rbps.annot, 1) 
  top.rbps.deep <- head(sig.df$ids, 1)
  top.scores <- head(sig.df$fa.seq.row,1)
  
  #top deep-rbps ids annotated
  sig.rbps.annot <- as.character(rbps.annot)
  top.rbps.deep <- as.character(top.rbps.deep)
  top.scores <- as.character(top.scores)
  
  #enriched kmer information
  enr.info.df <- seq.df[cnt,]
  seq.info <- as.character(enr.info.df$seq_name)
  
  line <- cbind(seq.info, as.character(kmer.seq), top.rbps.deep, sig.rbps.annot,top.scores )
  
  #print(head(sig.df,3))
  #print(line)
  
  write.table(line,file=outfile,row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t",append=TRUE)

  
  cnt <<- cnt + 1

})

print(paste("file saved: ", outfile))

