# MIRA
Mutation identification for RNA alterations

MIRA is a pipeline to identify mutational patterns related to RNA processing in human tumors. MIRA performs unbiased search for significantly mutated regions (SMRs) along gene loci through dividing sequence into short kmer windows of length n (default: 7). The method of calculation is based on binomial test and further corrected for local nucleotide biases. Later, the significant kmers thus obtained are overlapped with genic regions and motif enrichment test is performed for siginificant kmers enriched with mutations. These kmers and then labeled using annotation tool like Deepbind and downstream analysis of functional impact is performed.


Pipeline:

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
