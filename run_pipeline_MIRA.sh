##Define Filepaths
bedtoolspath="/soft/bio/sequence/bedtools-2.25.0/bin"
scriptpath="/projects_rg/babita/TCGA/mutation/mut_pipeline/scripts"
pythonpath="/soft/devel/python-2.7/bin"
refpath="/projects_rg/babita/ref_genome/no_pseudogenes"
fastafile="/projects_rg/babita/ref_genome/hg19.fa"


#change paths here
main_dir="/projects_rg/babita/TCGA/mutation/mut_pipeline/subs"
source_files="$main_dir/src_files"
mut_files_dir="$main_dir/outdir"
main_dir="/projects_rg/babita/TCGA/mutation/mut_pipeline/subs/subs_new"
mutpath="$main_dir/main_tables"

mutation_bed="$source_files/mutations_tcga_substitution_formatted.tsv"
annotation_bed="$source_files/gencode.v19.annotation.gtf.clean.GENE.bed.clustered.new"


#Define Functions
######################################################################################################
#1. create All mutation file from annotation and mutation bed file

fun_fjoin_annotate_mut_file () {
input_mutfile=$1
input_annotfile=$2

#1. create *.mut.out file by overlapping mutation .tsv file and genecode annotation.

$pythonpath/python /soft/bio/sequence/fjoin/fjoin.py -s both -1 $input_mutfile -2 $input_annotfile --columns1=1,2,3 --columns2=1,2,3 >"$input_mutfile.tmp"

#2. remove first column returned by fjoin
cat "$input_mutfile.tmp" | awk 'BEGIN{FS=OFS="\t"}{$1="";sub("\t","")}1' >"$input_mutfile.mut.out"

rm $input_mutfile.tmp

echo "Done function 1"
}

#----------------------------#


fun_count_mutation (){

#2. Run script: rscripts/run_count_mutation_analysis.py to count mutations on genes (see cluster scripts)
mkdir outdir/
$pythonpath/python $scriptpath/run_count_mutation_analysis.py 7 $fastafile "$mutation_bed.mut.out" outdir/ $annotation_bed

echo "Output files stored in ./outdir"
}


##create main output file
fun_create_main_table () {

#3. Read all mutaion files in outdir and create one single file
mutation_file=$1

#set cutoff value for adjust-Pvalue and score
score_cutoff=6
pval_cutoff=0.05

/soft/R/R-3.0.0/bin/Rscript $scriptpath/create_main_mutation_table.R $mut_files_dir "$mutation_file" $score_cutoff $pval_cutoff
echo "Done function 2"
}

#----------------------------#

fun_create_overlaps_control () {

cd /projects_rg/babita/TCGA/mutation/mut_pipeline/subs/subs_new/main_tables/controls_1

mutation_file="/projects_rg/babita/TCGA/mutation/mut_pipeline/subs/subs_new/main_tables/controls_1/tcga_subs_control_1_overlap_with_mutated_genes.tab"
#mutation_file="/projects_rg/babita/TCGA/mutation/mut_pipeline/subs/main_tables/tcga_subs_control_1_overlap_with_mutated_genes.tab"
mut_files_dir="./"
score_cutoff=6
pval_cutoff=0.05

#/soft/R/R-3.0.0/bin/Rscript $scriptpath/create_main_mutation_table.R $mut_files_dir "$mutation_file" $score_cutoff $pval_cutoff

mutation_file="$mutation_file.filtered"

#cat $mutation_file >tmp

#$scriptpath/trim_kmer_columns.R tmp $mutation_file 

rm tmp

#start codon
$bedtoolspath/intersectBed -s -wo -a $mutation_file -b $refpath/gencode.v19.annotation.gtf.clean.start_codon.bed.clustered >$mutation_file.start.codon
$bedtoolspath/intersectBed -s -v -wo  -a $mutation_file -b $refpath/gencode.v19.annotation.gtf.clean.start_codon.bed.clustered >$mutation_file.codon1.rest

echo "Done intersectBed $mutation_file.start.codon"

#stop codon
$bedtoolspath/intersectBed -s -wo -a $mutation_file.codon1.rest -b $refpath/gencode.v19.annotation.gtf.clean.stop_codon.bed.clustered >$mutation_file.stop.codon
$bedtoolspath/intersectBed -s -v -wo -a $mutation_file.codon1.rest -b $refpath/gencode.v19.annotation.gtf.clean.stop_codon.bed.clustered >$mutation_file.codon2.rest

echo "Done intersectBed $mutation_file.stop.codon" #none found

#5'SS
$bedtoolspath/intersectBed -s -wo -a $mutation_file.codon2.rest -b $refpath/gencode.v19.annotation.gtf.clean.CDS_5SS.bed.clustered.new >$mutation_file.5SS
$bedtoolspath/intersectBed -s -v -wo -a $mutation_file.codon2.rest -b $refpath/gencode.v19.annotation.gtf.clean.CDS_5SS.bed.clustered.new >$mutation_file.SS.rest

echo "Done intersectBed $mutation_file.5SS"



#3'SS
$bedtoolspath/intersectBed -s -wo -f 0.4 -a $mutation_file.SS.rest -b $refpath/gencode.v19.annotation.gtf.clean.CDS_3SS.bed.clustered.new >$mutation_file.3SS
$bedtoolspath/intersectBed -s -v -wo -f 0.4 -a $mutation_file.SS.rest -b $refpath/gencode.v19.annotation.gtf.clean.CDS_3SS.bed.clustered.new >$mutation_file.all.rest

echo "Done intersectBed $mutation_file.3SS"

rm *.rest
#====================


regions=("5SS" "3SS" "start.codon" "stop.codon")

for reg in "${regions[@]}" ; do

### create clusters of regions

	regfile="$mutation_file.$reg"

        $scriptpath/cluster_bed_simes.pl $regfile >"$regfile".clust.tmp
	$pythonpath/python $scriptpath/extract_seq.py $regfile.clust.tmp $fastafile $regfile.clust
	rm $regfile.clust.tmp

	echo "Done clustering $regfile.clust"

	#cat $regfile.clust | sort -k1,1V -k2,2n -k3,3nr | intersectBed -sorted -wo -a /projects_rg/babita/TCGA/mutation/mut_pipeline/subs/src_files/mutations_tcga_substitution_formatted.tsv.mut.out.sorted -b stdin >"$regfile.clust.ref"

	echo "Done clustering $regfile.clust.ref"

### scan with the mutation file

done
}

#----------------------------#

##create overlap of region files
fun_create_overlaps () {
mutation_file=$1

##change of pipeline : March 13, 2017. 

#Overlap kmer files with 5'SS and 3'SS as well as start and stop codons.
#Remove these kmers and cluster them seperately.

#start codon
$bedtoolspath/intersectBed -s -wo -a $mutation_file -b $refpath/gencode.v19.annotation.gtf.clean.start_codon.bed.clustered >$mutation_file.start.codon
$bedtoolspath/intersectBed -s -v -wo  -a $mutation_file -b $refpath/gencode.v19.annotation.gtf.clean.start_codon.bed.clustered >$mutation_file.codon.rest

echo "Done intersectBed $mutation_file.start.codon"

#stop codon
#$bedtoolspath/intersectBed -s -wo -a $mutation_file.codon.rest -b $refpath/gencode.v19.annotation.gtf.clean.stop_codon.bed.clustered >$mutation_file.stop.codon
#$bedtoolspath/intersectBed -s -v -wo -a $mutation_file.codon.rest -b $refpath/gencode.v19.annotation.gtf.clean.stop_codon.bed.clustered >$mutation_file.codon.rest

echo "Done intersectBed $mutation_file.stop.codon" #none found

#5'SS
$bedtoolspath/intersectBed -s -wo -a $mutation_file.codon.rest -b $refpath/gencode.v19.annotation.gtf.clean.CDS_5SS.bed.clustered >$mutation_file.5SS
$bedtoolspath/intersectBed -s -v -wo -a $mutation_file.codon.rest -b $refpath/gencode.v19.annotation.gtf.clean.CDS_5SS.bed.clustered >$mutation_file.SS.rest

echo "Done intersectBed $mutation_file.5SS"

#3'SS
$bedtoolspath/intersectBed -s -wo -f 0.4 -a $mutation_file.SS.rest -b $refpath/gencode.v19.annotation.gtf.clean.CDS_3SS.bed.clustered >$mutation_file.3SS
$bedtoolspath/intersectBed -s -v -wo -f 0.4 -a $mutation_file.SS.rest -b $refpath/gencode.v19.annotation.gtf.clean.CDS_3SS.bed.clustered >$mutation_file.all.rest

echo "Done intersectBed $mutation_file.3SS"

#====================
# Label rest of the kmers back to genic regions.


#fixing the kmers fraction of minimum 90% overlap of file A to B

$bedtoolspath/intersectBed -s -wo -f 0.9 -a $mutation_file.all.rest -b $refpath/gencode.v19.annotation.gtf.clean.CDS.bed.clustered >$mutation_file.CDS
$bedtoolspath/intersectBed -s -v -wo -f 0.9 -a $mutation_file -b $refpath/gencode.v19.annotation.gtf.clean.CDS.bed.clustered >$mutation_file.notCDS.rest

echo "Done intersectBed $mutation_file.CDS"

$bedtoolspath/intersectBed -s -wo -f 0.9 -a $mutation_file.notCDS.rest -b $refpath/gencode.v19.annotation.gtf.clean.3UTR.bed.clustered >$mutation_file.3UTR
$bedtoolspath/intersectBed -s -v -wo -f 0.9 -a $mutation_file.notCDS.rest -b $refpath/gencode.v19.annotation.gtf.clean.3UTR.bed.clustered >$mutation_file.notCDS-3UTR.rest

echo "Done intersectBed $mutation_file.3UTR"

$bedtoolspath/intersectBed -s -wo -f 0.9 -a $mutation_file.notCDS-3UTR.rest -b $refpath/gencode.v19.annotation.gtf.clean.5UTR.bed.clustered >$mutation_file.5UTR
$bedtoolspath/intersectBed -s -v -wo -f 0.9 -a $mutation_file.notCDS-3UTR.rest -b $refpath/gencode.v19.annotation.gtf.clean.5UTR.bed.clustered >$mutation_file.notCDS-UTR.rest

echo "Done intersectBed $mutation_file.5UTR"

$bedtoolspath/intersectBed -s -wo -f 0.9 -a $mutation_file.notCDS-UTR.rest -b $refpath/gencode.v19.annotation.gtf.clean.EXON.bed.clustered >$mutation_file.EXON
$bedtoolspath/intersectBed -s -v -wo -f 0.9 -a $mutation_file.notCDS-UTR.rest -b $refpath/gencode.v19.annotation.gtf.clean.EXON.bed.clustered >$mutation_file.notCDS-UTR-EXON.rest
echo "Done intersectBed $mutation_file.EXON"

$bedtoolspath/intersectBed -s -wo -f 0.9 -a $mutation_file.notCDS-UTR-EXON.rest -b $refpath/gencode.v19.annotation.gtf.clean.INTRON.bed.clustered >$mutation_file.INTRON
$bedtoolspath/intersectBed -s -v -wo -f 0.9 -a $mutation_file.notCDS-UTR-EXON.rest -b $refpath/gencode.v19.annotation.gtf.clean.INTRON.bed.clustered >$mutation_file.notCDS-UTR-EXON-INTRON.rest

echo "Done intersectBed $mutation_file.INTRON"

}



fun_create_clusters () {
### create clusters of regions
	reg=$1

	regfile="$mutation_file.$reg"

        $scriptpath/cluster_bed_simes.pl $regfile >"$regfile".clust.tmp
	$pythonpath/python $scriptpath/extract_seq.py $regfile.clust.tmp $fastafile $regfile.clust
	rm $regfile.clust.tmp

	echo "Done clustering $regfile.clust"

}



fun_create_controls () {
#create control sets for each cluster of kmers for each regions (UTRs/Exons/Introns/CDSs)
	reg=$1
	clustfile="$mutation_file.$reg.clust"
	
	annotfile="$refpath/gencode.v19.annotation.gtf.clean.$reg.bed.clustered"
	
	$pythonpath/python $scriptpath/create_control_set.py $clustfile $annotfile $fastafile 100

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



mkdir $mutpath

#function 1
fun_fjoin_annotate_mut_file $mutation_bed $annotation_bed



#function 2
#2. better Run script on clusters: #path:  /data/users/babita/2016_11_28_mutation_count/script/run_mutation_pipeline2.sh on cluster
fun_count_mutation 


#function 3.
mutation_file="tcga_all_genes_mutations_7mer.tab"
fun_create_main_table "$mutpath/$mutation_file"



#function 4.
mutation_file="tcga_all_genes_mutations_7mer.tab.filtered"
fun_create_overlaps "$mutpath/$mutation_file"

##call functions per regions
regions=("CDS" "3UTR" "5UTR" "EXON" "INTRON" "3SS" "5SS")


for region in "${regions[@]}" ; do

	echo "$region"
	fun_create_clusters $region


	fun_create_controls $region
	fun_create_enrichment $region

	
	#run the script "reverse_complement_enrichment.py" before running following
	file_dir="$mutpath/rev_enrichment"
	mutation_file="tcga_all_genes_mutations_7mer.tab.filtered"
	
	#mkdir $file_dir
	fun_create_reverse_enrichment $region $file_dir $mutation_file
	
done

