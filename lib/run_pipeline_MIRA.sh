##Define Filepaths
bedtoolspath="/soft/bio/sequence/bedtools-2.25.0/bin"
pythonpath="/soft/devel/python-2.7/bin"
perlpath="/soft/devel/perl-5.16.3/bin"
rscriptpath="/soft/R/R-3.0.0/bin"

scriptpath="/home/babita/Desktop/MIRA/lib"
refpath="/home/babita/Desktop/MIRA/ref_files"
ref_genome="/home/babita/ref_genome/hg19.fa" #human genome fasta file (hg19.fa)

kmer_len=7
main_outdir="/home/babita/Desktop/MIRA/test$kmer_len"
mkdir $main_outdir
mutpath="$main_outdir/main_tables" 
mkdir $mutpath

main2_outdir="/projects_rg/babita/TCGA/mutation/mut_pipeline/subs/subs_new/test_"$kmer_len"mers"

mutation_bed="$refpath/mutations_tcga_substitution_formatted.tsv"
annotation_bed="$refpath/GENE.gencode.v19.clustered.bed"



#Define Functions
######################################################################################################
#1. create a common mutation file by overlapping gene annotation file and refernece mutation file

fun_fjoin_annotate_mut_file () {
input_mutfile=$1
input_annotfile=$2

#1. create *.mut.out file by overlapping mutation .tsv file and genecode annotation.
cd $scriptpath

$pythonpath/python fjoin.py -s both -1 $input_mutfile -2 $input_annotfile --columns1=1,2,3 --columns2=1,2,3 >"$input_mutfile.tmp"

#2. remove first column returned by fjoin
cat "$input_mutfile.tmp" | awk 'BEGIN{FS=OFS="\t"}{$1="";sub("\t","")}1' >"$input_mutfile.mut.out"

rm $input_mutfile.tmp

echo "Done function 1"
}

#----------------------------#


fun_count_mutation (){

#2. Run script: rscripts/run_count_mutation_analysis.py to count mutations on genes (see cluster scripts)

cd $main_outdir
mkdir outdir
out_files_dir="$main_outdir/outdir"

kmer_len=$kmer_len
mut_count=3 	#mutations threshold on each kmer >= 3

#args: python <script.py> <length_kmer> <threshold-mut_count> <ref_genome_fasta_file(hg19.fa)> <file *.mut.out> <output_dir> <gene_annotation_file> <bedtools_path> <perl_path>
$pythonpath/python $script_path/run_count_mutation_analysis.py $kmer_len $mut_count $ref_genome $mutation_bed.mut.out $out_files_dir $annotation_bed $bedtoolspath $perlpath


echo "Output files stored in ./outdir"
}


##create main output file
fun_create_main_table () {

#3. Function to read all mutaion files in outdir (generated through function 'fun_count_mutation') and create one single mutation file.

mutation_file=$1
mut_outdir=$2

#set cutoff value for adjust-Pvalue and score
score_cutoff=6
pval_cutoff=0.05

echo "============"
echo $mut_outdir
echo $mutation_file
#args: rscript <script.R> <mutation_files_output_dir> <output_file_name> <score_cutoff> <pvalue_cutoff>
$rscriptpath/Rscript $scriptpath/create_main_mutation_table.R $mut_outdir "$mutation_file" $score_cutoff $pval_cutoff

echo "Done function 2"
}

#----------------------------#

##create overlap of region files
fun_create_overlaps () {
 mutation_file=$1

 echo $mutation_file
 echo "$refpath/5SS.gencode.v19.clustered.bed"
 #====================
 #Overlap kmer files with 5'SS and 3'SS first
 #====================

 #5'SS
 $bedtoolspath/intersectBed -s -wo -a $mutation_file -b $refpath/5SS.gencode.v19.clustered.bed >$mutation_file.5SS
 $bedtoolspath/intersectBed -s -v -wo -a $mutation_file -b $refpath/5SS.gencode.v19.clustered.bed >$mutation_file.SS.rest

 echo "Done intersectBed $mutation_file.5SS"

 #3'SS
 $bedtoolspath/intersectBed -s -wo -f 0.4 -a $mutation_file.SS.rest -b $refpath/3SS.gencode.v19.clustered.bed >$mutation_file.3SS
 $bedtoolspath/intersectBed -s -v -wo -f 0.4 -a $mutation_file.SS.rest -b $refpath/3SS.gencode.v19.clustered.bed >$mutation_file.all.rest

 echo "Done intersectBed $mutation_file.3SS"

 #====================
 # Label rest of the kmers back to other genic regions.
 #====================

 #fixing the kmers fraction of minimum 90% overlap of file A to B

 $bedtoolspath/intersectBed -s -wo -f 0.9 -a $mutation_file.all.rest -b $refpath/CDS.gencode.v19.clustered.bed >$mutation_file.CDS
 $bedtoolspath/intersectBed -s -v -wo -f 0.9 -a $mutation_file -b $refpath/CDS.gencode.v19.clustered.bed >$mutation_file.notCDS.rest

 echo "Done intersectBed $mutation_file.CDS"

 $bedtoolspath/intersectBed -s -wo -f 0.9 -a $mutation_file.notCDS.rest -b $refpath/3UTR.gencode.v19.clustered.bed >$mutation_file.3UTR
 $bedtoolspath/intersectBed -s -v -wo -f 0.9 -a $mutation_file.notCDS.rest -b $refpath/3UTR.gencode.v19.clustered.bed >$mutation_file.notCDS-3UTR.rest

 echo "Done intersectBed $mutation_file.3UTR"

 $bedtoolspath/intersectBed -s -wo -f 0.9 -a $mutation_file.notCDS-3UTR.rest -b $refpath/5UTR.gencode.v19.clustered.bed >$mutation_file.5UTR
 $bedtoolspath/intersectBed -s -v -wo -f 0.9 -a $mutation_file.notCDS-3UTR.rest -b $refpath/5UTR.gencode.v19.clustered.bed >$mutation_file.notCDS-UTR.rest

 echo "Done intersectBed $mutation_file.5UTR"

 $bedtoolspath/intersectBed -s -wo -f 0.9 -a $mutation_file.notCDS-UTR.rest -b $refpath/EXON.gencode.v19.clustered.bed >$mutation_file.EXON
 $bedtoolspath/intersectBed -s -v -wo -f 0.9 -a $mutation_file.notCDS-UTR.rest -b $refpath/EXON.gencode.v19.clustered.bed >$mutation_file.notCDS-UTR-EXON.rest
 echo "Done intersectBed $mutation_file.EXON"

 $bedtoolspath/intersectBed -s -wo -f 0.9 -a $mutation_file.notCDS-UTR-EXON.rest -b $refpath/INTRON.gencode.v19.clustered.bed >$mutation_file.INTRON
 $bedtoolspath/intersectBed -s -v -wo -f 0.9 -a $mutation_file.notCDS-UTR-EXON.rest -b $refpath/INTRON.gencode.v19.clustered.bed >$mutation_file.notCDS-UTR-EXON-INTRON.rest

 echo "Done intersectBed $mutation_file.INTRON"


 #remove *.rest files 
 rm $mutation_file*.rest
}



fun_create_clusters () {
### create clusters of regions
	mutation_file=$1
        reg=$2
    
	regfile="$mutation_file.$reg"

        $perlpath/perl $scriptpath/cluster_bed_simes.pl $regfile >$regfile.clust.tmp
	$pythonpath/python $scriptpath/extract_seq.py $regfile.clust.tmp $ref_genome $regfile.clust
	rm $regfile.clust.tmp

	echo "Done clustering $regfile.clust"

}



fun_create_controls () {
#create control sets for each cluster of kmers for each regions (UTRs/Exons/Introns/CDSs)
	reg=$1
	clustfile="$mutation_file.$reg.clust"
	
	annotfile="$refpath/$reg.gencode.v19.clustered.bed"
	
	$pythonpath/python $scriptpath/create_control_set.py $clustfile $annotfile $ref_genome 100

	echo "Done control set for $reg. Outfile saved in: $clustfile.control"
	
}



fun_create_enrichment () {

		reg=$1
		mutation_cluster_file="$mutation_file.$reg.clust"
		control_cluster_file="$mutation_file.$reg.clust.control"
		len_kmers=6
		control_outdir="$mutpath/$reg"
		control_set=100
	 	
	echo "$mutation_cluster_file $control_cluster_file $len_kmers $control_outdir $reg $control_set $reg"

	$pythonpath/python $scriptpath/control_enrichment.py $mutation_cluster_file $control_cluster_file \
			$len_kmers $control_outdir $reg $control_set

	echo "Done shell : enrichment $reg"

}


fun_create_reverse_enrichment (){
		#run scripts/reverse_complement_enrichment.py before the following:

		reg=$1
		file_dir=$2
		mutation_file=$3

		cd $file_dir

		#$pythonpath/python $scriptpath/reverse_complement_enrichment.py $file_dir $mutation_file"

		mutation_cluster_file="$file_dir/$mutation_file.$reg.clust.rev"
		control_cluster_file="$file_dir/$mutation_file.$reg.clust.control.rev"
		len_kmers=6
		control_outdir="$file_dir/$reg"
		control_set=100
	 	
	echo "$mutation_cluster_file $control_cluster_file $len_kmers $control_outdir $reg $control_set"

	$pythonpath/python $scriptpath/control_enrichment.py $mutation_cluster_file $control_cluster_file \
			$len_kmers $control_outdir $reg $control_set

	echo "Done shell : enrichment $reg"
}

######################################################################################################

#function 1: Create a reference mutation file by extracting mutations found in geneic regions
#fun_fjoin_annotate_mut_file $mutation_bed $annotation_bed


#function 2
#Note: Before running this funtion it's better to extract the gene ids that overlapped with mutation file in function 1 -
 	#through 'fjoin', convert it to bedfile format and run the following script for only the overlapping genes.
#Note 2: Even better if the overlapping genes are divided into chunks of 100-200 genes and run parallel processes in the cluster.

#fun_count_mutation 


#function 3.
mut_outdir="$main2_outdir/outdir"
mutation_file1="tcga_all_genes_mutations_"$kmer_len"mer.tab" #change filename
#fun_create_main_table "$mutpath/$mutation_file1" $mut_outdir


#function 4.
mutation_file1="$mutpath/tcga_all_genes_mutations_"$kmer_len"mer.tab.filtered"
#fun_create_overlaps "$mutation_file1"



##call functions per regions
regions=("CDS" "3UTR" "5UTR" "EXON" "INTRON" "3SS" "5SS")


for region in "${regions[@]}" ; do

	echo "$region"
	fun_create_clusters "$mutation_file1" $region
	
	#fun_create_controls $region
	#fun_create_enrichment $region

	
	#run the script "reverse_complement_enrichment.py" before running following
	file_dir="$mutpath/rev_enrichment"
	#mutation_file1="tcga_all_genes_mutations_"$kmer_len"mer.tab.filtered"
	
	#mkdir $file_dir
	#fun_create_reverse_enrichment $region $file_dir $mutation_file1
	
done

