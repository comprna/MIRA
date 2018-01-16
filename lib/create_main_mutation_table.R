####
#Script to create one file from all mutation files
#and to filter later the main file based on score and pvalues cutoff, and kmer repeats
####
library(stringr)
library(data.table)

fun_create_main_table <- function(in_dir, outfile_main){
  
  setwd(in_dir)
  
  file_list = list.files(path = ".", pattern = "*.mut.out", full.names=TRUE)
  
  #create outfile:
  file.create(outfile_main)
  cat("chr\tkmer_start\tkmer_end\tgene\tkmer\tstrand\tchr_gene\tgene_start\tgene_end\tkmer_mut\tgene_mut\tscore\tpval\n",file=outfile_main)
  
  not_found = 0
  
  lapply(file_list, function(file){
    
    if(file.exists(file) && file.info(file)$size != 0){
      
      df <- read.delim(file, header=FALSE)
      df <- df[complete.cases(df),]
  
        f <- basename(file) 
      
        f <- str_split_fixed(f,"_",4)[,4]
        f <- gsub(".mut.out", "", f)
       
        df$V4 <- f
        #print(head(df, 2))
       
        
      #write to outfile
      write.table(df,file=outfile_main, sep="\t", quote=FALSE, append=TRUE,row.names=FALSE, col.names=FALSE)

    } else {
      
      not_found = not_found + 1
      print(paste("Not found:", file, file.info(file)$size, sep=" | "))
    }
    
  })
 
  print("Done fun: create main table")
}
    



fun_create_cutoff_table <- function(outfile_main, score.cutoff, pval.cutoff,  outfile_cutoff ){
  
  file.create(outfile_cutoff)
    
  ##df.all <- read.delim(outfile_main,header=TRUE)
  ##df.all <- fread(outfile_main, header=FALSE)
  df.all <- fread(outfile_main, header=TRUE)
  df.all <- as.data.frame(df.all)
  
  df.all <- df.all[complete.cases(df.all),]
  #remove non-numeric columns from data
  #header = chr\tkmer_start\tkmer_end\tgene\tkmer\tstrand\tchr_gene\tgene_start\tgene_end\tkmer_mut\tgene_mut\tscore\tpval

  
  #original:
  df.all$padjust <- p.adjust(as.numeric(df.all$pval), method = "BH", n = length(df.all$pval))
  #subset dataframe based on cut-offs
  sub.df.all <- subset(df.all, score > as.numeric(score.cutoff) & padjust < as.numeric(pval.cutoff))
  
  
  
  #following is only for control

  #original:
  ##df.all$padjust <- p.adjust(as.numeric(df.all$V13), method = "BH", n = length(df.all$V13)) 
  ###sub.df.all <- subset(df.all, V12 > as.numeric(score.cutoff) & padjust < as.numeric(pval.cutoff))
  ####sub.df.all <- subset(df.all, V12 > as.numeric(score.cutoff))

  ##print(head(sub.df.all))
  #subset dataframe based on kmers (remove repeat kmers)
  repeats_list <- c("AAAAAAAA", "TTTTTTTT", "GGGGGGGG", "CCCCCCCC", "AAAAAAA", "TTTTTTT", "GGGGGGG", "CCCCCCC")
  pat_repeats <- paste(repeats_list, collapse="|")

  #original:
  sub.df.all <- subset(sub.df.all, !(grepl(pat_repeats, sub.df.all$kmer)))
  ###sub.df.all <- subset(sub.df.all, !(grepl(pat_repeats, sub.df.all$V5)))
  
  write.table(sub.df.all,file=outfile_cutoff, sep="\t", quote=FALSE, append=TRUE,row.names=FALSE, col.names=FALSE)
  
  print(paste("Done fun: create cutoff table: ",outfile_cutoff))
 
}




#mutation_files_dir = "/projects_rg/babita/TCGA/mutation/mut_pipeline/gene_lib_7/outdir"
#main_tables_dir = "/projects_rg/babita/TCGA/mutation/mut_pipeline/main_tables"
#outfile_main = file.path(main_tables_dir,"tcga_all_genes_mutations_7mer.tab")

##get input gene files directory and outfile name
args <- commandArgs(trailingOnly = TRUE)
mutation_files_dir <- args[1] 
outfile_main <- args[2]
score_cutoff <- args[3]
pval_cutoff <- args[4]


if(length(args) < 4){

print("usage: script.R <mut_dir_path> <main_table_name> <score_cutoff> <pval_cutoff>")
stop("input parameters wrong")
}

print(paste(mutation_files_dir, outfile_main, score_cutoff, pval_cutoff, sep=" | "))

#Following function creates one main file for all mutation scores files

fun_create_main_table(mutation_files_dir, outfile_main)

outfile_cutoff = paste0(outfile_main,".filtered")
fun_create_cutoff_table(outfile_main, score_cutoff, pval_cutoff, outfile_cutoff)







