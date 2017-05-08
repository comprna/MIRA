#This script reads the 5-columns Kmer-enrichment files (forward and reverse strand) and extract those 
#enriched kmers that are not present on reverse strand (ie to only extract RNA related motifs)



## All input files with *No header* 
args <- commandArgs(trailingOnly = TRUE)
main_enr_file <- args[1] #columns: kmer count mean_control sd_control zscore
rev_enr_file <- args[2]
outfile <- args[3]
label.option <- as.numeric(args[4])
#optional: If user wants to match kmers with known db, provide kmer \t label \t .... files for the labels available for known kmers
#First two columns of the files are necessary where column1 = kmer, column2 = label (Motif-id), rest of the columns are optional (only first two columns will be used)

if (!file.size(main_enr_file) ==0 && !file.size(rev_enr_file) ==0 ) {
  
  #for kmers eriched in main file
  df.main <- read.delim(main_enr_file, header=FALSE)
  df.rev <- read.delim(rev_enr_file, header=FALSE)
  
} else {
  
  stop(paste("Empty file: ", main_enr_file))
}


if(label.option){
  
  motif_labels_file <- args[5]

  #Filter for RBP motifs that are >= 6kmer lenght and only from humans,mouse and rat.
  #cmd: $cat /projects_rg/babita/TCGA/mutation/attract_motifs/attract_rbp_motifs.txt | egrep 'Homo_sapiens|Mus_musculus|Rattus_norvegicus' \
  # | awk '$5 > 5' | cut -f1,2,3,4,7 | awk 'OFS="\t"{print $4,$2,$3,$1,$7}' >scripts/attract_rbp_motif.filtered

  #motif_labels_file <- file.path("/projects_rg/babita/TCGA/mutation/attract_motifs/attract_rbp_motifs.txt")
  labels.df <- read.delim(motif_labels_file,header=FALSE)
  labels.df <- labels.df[,c("V1","V2")]
  labels.df <- labels.df[!duplicated(labels.df), ]
  names(labels.df) <- c("Seq","Motif")
  #convert all Us to Ts
  labels.df$Seq <- gsub("U", "T", labels.df$Seq)
}




kmer_count_cutoff <- 5
  
df.main <- df.main[,c("V1","V2","V5")]
names(df.main) <- c("kmer","count","z.score")
df.main <- subset(df.main, z.score > 1.96)
df.main <- df.main[is.finite(df.main$z.score), ] #to remove NaN and Inf
df.main <- subset(df.main, count >= kmer_count_cutoff) #FOR INTRON IT IS >3
  


  #print(head(df.main))
  #============= for reverse ============================
  df.rev <- df.rev[, c("V1","V5")]
  names(df.rev) <- c("kmer","z.score")
  df.rev <- subset(df.rev, z.score > 1.96 )
  df.rev <- df.rev[is.finite(df.rev$z.score),]

  df.main$Enriched.in.reverse <- sapply(df.main$kmer, function(x){
    d <- subset(df.rev, kmer == as.character(x))
    if(nrow(d)) return("Yes")
    else return("No")
  })
  
  #remove the kmers, also present in reverse enrichment
  df.main <- subset(df.main, Enriched.in.reverse == "No")
  
  
  #print(head(df.main))

#============= label kmers ============================
#(if label file is provided)
if(label.option == 1){ #label kmers
  
  df.main$annot <- sapply(df.main$kmer, function(x){  
    labels.df.main <- subset(labels.df, grepl(x, labels.df$Motif)) 
    
    if(nrow(labels.df.main)){
      rbp.genes <- unique(labels.df.main$Gene.Name)
      rbp.genes <- paste(rbp.genes, collapse=",")
      rbp.genes <- paste(rbp.genes, as.character(x), sep= "_")
    }
    
    else{
      rbp.genes <- as.character(x)
    }
    
    return(rbp.genes)
  })
  
}

  
  #outfile <- file.path(main_dir,paste0(reg, ".kmers.zscore.enriched"))
  write.table(file=outfile, df.main, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

 print(paste("Done: ",outfile))