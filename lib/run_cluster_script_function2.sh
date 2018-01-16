#$ -N mut_subs_7mers
#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -q long,bigmem,normal

## This script is to run on cluster.
### Divide the chunk of genes into 100-200 list and run the process of mutation counting on kmers 

#Usage: qsub script/run_mutation_pipeline.sh  #for gencluster

#bed_line is bedfile coordinates for each gene.

#source /usr/local/sge/default/common/settings.sh

export PATH=/soft/devel/python-2.7/bin/python:$PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/soft/lib
export PATH=/.bashrc:$PATH
export PATH=/soft/devel/perl-5.16.3/bin/:$PATH

#export PATH=/soft/devel/python-2.7/lib/python2.7/site-packages/pybedtools-0.6.2-py2.7-linux-x86_64.egg/pybedtools/:$PATH

kmer_len=7 	#change length of kmers
mut_count=3  	#change minimum mutation required
ref_genome="./ref_genome/hg19.fa"

main_dir="./"
out_dir="./test_7mers"
refpath="MIRA/ref_files"  

mutation_bed="$refpath/mutations_tcga_substitution_formatted.tsv"   #add mutation file *.tsv in ref_files (Not provided in repository)
annotation_bed="$refpath/gencode.v19.annotation.gtf.clean.GENE.bed.clustered.new"

all_gene_overlap_mutfile="$mutation_bed.mut.out"   #path with 'fjoined' mutation file obtained from function 1
gene_lib_dir="$out_dir/outdir"


mkdir $gene_lib_dir

script_path="MIRA/lib"


echo "shell start"


LINES_ANNOT=$(wc -l < $annotation_bed)

cnt=0

for ((i=0; i<=LINES_ANNOT; i+=100)); #no. of gene chunks to be cut from main annotation file. Here, run 100 genes at a time.

do

    echo $i
    ii=$((i+1))
    j=$((i+100))

    #subset_gene_file="genes_$j"  into 100 gene chunks (change to whatever number instead of 100)
    
    cat $annotation_bed | sed -n "$ii","$j"p >genes_"$j".bed

    cnt=$((cnt+1))
    #Run for each job in clusters for 100 genes 
           
    job="/soft/devel/python-2.7/bin/python $script_path/run_count_mutation_analysis.py $kmer_len $mut_count $ref_genome $all_gene_overlap_mutfile $gene_lib_dir genes_"$j".bed /soft/bio/sequence/bedtools-2.25.0/bin /soft/devel/perl-5.16.3/bin"

    #source /usr/local/sge/default/common/settings.sh
    #source ~/.bashrc
    export PATH=/soft/devel/python-2.7/lib/python2.7/site-packages/pybedtools-0.6.2-py2.7-linux-x86_64.egg/pybedtools/:$PATH
 
    qsub -N job_"$cnt" -b y -cwd -V -q long $job  #submit jobs to cluster
 

    echo "Done gene_$j.bed"     
    #sleep 1000
        
done


	
echo "shell done "
echo "==================="

