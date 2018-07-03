![Mira1.jpg](https://user-images.githubusercontent.com/23315833/35625867-2fe09b2a-0694-11e8-98b5-a4c6c2caf483.png)


# MIRA
Mutation Identification for RNA Alterations

MIRA is a pipeline to identify mutational patterns related to RNA processing in human tumors. MIRA performs unbiased search for significantly mutated regions (SMRs) along gene loci through dividing sequence into short kmer windows of length n (default: 7). The method of calculation is based on binomial test and further corrected for local nucleotide biases. Later, the significant kmers thus obtained are overlapped with genic regions and motif enrichment test is performed for siginificant kmers enriched with mutations. These kmers and then labeled using annotation tool like Deepbind and downstream analysis of functional impact is performed.


### Please run script 'lib/run_pipeline_MIRA.sh' to run the pipeline.


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

Input parameters:   <length_kmer> <threshold-mut_count> <ref_genome_fasta_file(hg19.fa)> <file *.mut.out> <output_dir> <gene_annotation_file> <bedtools_path> <perl_path>
1. <Kmer_length>: Length of the kmer regions to scan. (Integer, usually between 5-8)
2. <threshold-mut_count>: Minimum mutation count on the kmers (Integer, usually 3 or more)
3. <Ref_genome>: Reference genome file (file path, eg. hg19.fa)
4. <mut_gene_overlap_mutfile>: MAF or mutation cordinates overlapped with gene cordinates* (file path, see command below to generate such file).
5. <output_gene_lib_dir>: An output directory where all files will be created. (file path)
6. <ggene_annotation_file>: Gene annotation file in bed format (file path, used to create file 3 above, provided in ref_files/)
7. <bedtools_path> : Path to bedtools (path, eg. /bedtools/bin) 
8. <perl_path> : Path to Perl dir (path, eg. /perl/bin)

```
python MIRA/run_count_mutation_analysis.py <$length_kmer> <$threshold-mut_count> <$ref_genome_fasta_file> <$mut_gene_overlap_mutfile> <$output_gene_lib_dir> <$gene_annotation_file> <$bedtools_path> <$perl_path>
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

http://comprna.upf.edu/Data/MutationsRBPMotifs/

http://mcr.aacrjournals.org/content/early/2018/04/18/1541-7786.MCR-17-0601
