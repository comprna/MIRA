# MIRA
Mutation identification for RNA alterations

MIRA is a pipeline to identify mutational patterns related to RNA processing in human tumors. MIRA performs unbiased search for significantly mutated regions (SMRs) along gene loci through dividing sequence into short kmer windows of length n (default: 7). The method of calculation is based on binomial test and further corrected for local nucleotide biases. Later, the significant kmers thus obtained are overlapped with genic regions and motif enrichment test is performed for siginificant kmers enriched with mutations. These kmers and then labeled using annotation tool like Deepbind and downstream analysis of functional impact is performed.

Command to run MIRA directly:
## A. Overlap mutation (MAF) file with Gene cordinates 

*Command to overlap gene file with MAF file

```
#1. Use fjoin.py to combine both files
python MIRA/fjoin.py -s both -1 $input_mutation_File -2 $input_annotation_File --columns1=1,2,3 --columns2=1,2,3 >"$input_mutfile.tmp"

#2. remove first column returned by fjoin
cat "$input_mutfile.tmp" | awk 'BEGIN{FS=OFS="\t"}{$1="";sub("\t","")}1' >"$input_mutfile.mut.out"
rm $input_mutfile.tmp
````

## B. Run MIRA mutation analysis command:

Prerequisites:
1. Kmer_length: length of the kmer regions to scan. (Integer, usually between 5-8)
2. ref_genome: Reference genome file (eg. hg19.fa)
3. all_gene_overlap_mutfile: MAF or mutation cordinates overlapped with gene cordinates* (see command below to generate such file).
4. gene_lib_dir: An output directory where all files will be created.
5. gene_annotation.bed: Gene annotation file in bed format (used to create file 3 above).
```
python MIRA/run_count_mutation_analysis.py <$kmer_length> <$ref_genome> <$all_gene_overlap_mutfile> <$gene_lib_dir> <gene_annotation.bed>
```

## C. Create clusters of sigificantly mutated regions:

This command will create all files where all mutated regions were found. Run following command to create the clusters of these siginificantly mutated regions (SMRs).

Prerequisite:
1. score_cutoff=6  #Minimum score cutoff to use 
2. pval_cutoff=0.05 #Pvalue cutoff.
3. gene_lib_dir: An output directory where all files were created by script 'run_count_mutation_analysis.py'.
4. out_mutation_file : output file name for all SMRs.

```
Rscript MIRA/create_main_mutation_table.R <$gene_lib_dir> <$out_mutation_file> <$score_cutoff> <$pval_cutoff>
```

## D. Seperate SMRs into region types
Optional Step:
Divide SMRs into region type, for eg. CDS/UTRs/INTRONs etc (Requires annotation files of different region type in bedfile format). It uses bedtools; IntersectBed.


Pipeline overview:

1) run_pipeline_MIRA.sh: 
  - Takes MAF files and gene cordinates, scan mutations for each gene, creating 7-mer windows and reports only significant kmer windows.
  - Overlaps significant kmer windows to genic regions (UTRs/Exons/CDS and Introns) based on cascade method. 
  - Clusters kmers to SMRs (significantly mutated regions)
  - Performs enrichment analysis of the kmers with significant mutations
  

2) run_pipeline_MIRA_enrichment.sh:
  
  - Take the regulated SMRs (\*.clust) files, create control sequences and performs kmer enrichment
  - label kmers with Deepbind, through deepbind analysis.


Results:

Results from MIRA on 505 WGS samples from TCGA and associated RNA-seq samples are available here
http://comprna.upf.edu/Data/MutationsRBPMotifs/INTRON/SNRNP70/
https://www.biorxiv.org/content/early/2017/10/09/200188
