#Pipeline to run splice sites kmer enrichment analysis

#Prerequisite:
#	- *.clust files (SMR files, ie mutated kmers clustered together to create SMRs)
#	- reference file with 5SS/3SS/stop/start cordinates for all genes (to select controls)


#paths fixed
bedtoolspath="/soft/bio/sequence/bedtools-2.25.0/bin"
scriptpath="/projects_rg/babita/TCGA/mutation/mut_pipeline/scripts"
pythonpath="/soft/devel/python-2.7/bin"
refpath="/projects_rg/babita/ref_genome/no_pseudogenes"
fastafile="/projects_rg/babita/ref_genome/hg19.fa"
RscriptPath="/soft/R/R-3.2.3/bin"
motif_filtered_db="/projects_rg/babita/TCGA/mutation/motif_filtered_db"
deepbindpath="~/soft/deepbind"


# Paths to replace
smr_dir="/projects_rg/babita/TCGA/mutation/mut_pipeline/subs/subs_new/main_tables"
ref_dir="/projects_rg/babita/ref_genome/no_pseudogenes/"  
mutation_file="$smr_dir/tcga_all_genes_mutations_7mer.tab.filtered"

#####functions ###

fun_create_controls () {
#create control sets for each cluster of kmers for each regions (UTRs/Exons/Introns/CDSs)
	reg=$1
	
	clustfile="$mutation_file.$reg.clust"
	annotfile="$refpath/gencode.v19.annotation.gtf.clean.$reg.bed.clustered"

	$pythonpath/python $scriptpath/create_control_set.py $clustfile $annotfile $fastafile 100
 	
	#remove coordinates from control file that are already present in regulated set.
	cat "$clustfile.control" >$smr_dir/tmp
	$bedtoolspath/intersectBed -v -a $smr_dir/tmp -b $clustfile >"$clustfile.control" 
}


fun_create_enrichment () {

		reg=$1
		mutation_cluster_file="$mutation_file.$reg.clust"
		control_cluster_file="$mutation_file.$reg.clust.control"
		len_kmers=6
		control_outdir="$smr_dir/$reg"
		control_set=100
	 	min_kmer_count=3 #minimum number of kmer occurence in regulated set
		output_dir=$smr_dir #where enrichment outfile will be stored with .zscore ext.		

	echo "$mutation_cluster_file $control_cluster_file $len_kmers $control_outdir $reg $control_set $reg\
		$min_kmer_count $output_dir"

	$pythonpath/python $scriptpath/control_enrichment.py $mutation_cluster_file $control_cluster_file \
			$len_kmers $control_outdir $reg $control_set $min_kmer_count $output_dir

	echo "Done shell : enrichment $reg"

}


fun_create_reverse_enrichment (){
		#run scripts/reverse_complement_enrichment.py before the following:

		reg=$1
		file_dir="$smr_dir/rev_enrichment"
		input_file="$mutation_file.$reg" #"tcga_all_genes_mutations_7mer.tab.filtered"

		#Following will create reverse complement set for both 'clust' and 'clust.control' files
		$pythonpath/python $scriptpath/reverse_complement_enrichment.py $file_dir $input_file

		#once rev-complement files are created, run enrichment analysis
		input_file=$(basename $input_file)
		mutation_cluster_file="$file_dir/$input_file.clust.rev"
		control_cluster_file="$file_dir/$input_file.clust.control.rev"
		len_kmers=6
		control_outdir="$file_dir/$reg"
		control_set=100
		min_kmer_count=3 	#minimum number of kmer occurence in regulated set
		output_dir=$file_dir 	#where enrichment outfile will be stored with .zscore ext.
	 	
	echo "$mutation_cluster_file $control_cluster_file $len_kmers $control_outdir $reg $control_set \
		$min_kmer_count $output_dir"

	$pythonpath/python $scriptpath/control_enrichment.py $mutation_cluster_file $control_cluster_file \
			$len_kmers $control_outdir $reg $control_set $min_kmer_count $output_dir

	echo "Done shell : enrichment $reg"
}


fun_filter_RNA_motifs () {
	
	reg=$1
	main_enriched_file="$smr_dir/$reg.kmers.zscore"
	rev_enriched_file="$smr_dir/rev_enrichment/$reg.kmers.zscore"
	outdir="$smr_dir/enriched_kmers"
	mkdir $outdir
	outfile="$outdir/$reg.kmers.zscore.enriched"
	label_option=1  				#if kmers annotation is available, else 0 
	labels_file="$motif_filtered_db/attract_rbp_motif.filtered"  #not necessary if label_option = 0

        $RscriptPath/Rscript $scriptpath/filter_enriched_RNA_motifs.R $main_enriched_file $rev_enriched_file \
	$outfile $label_option $labels_file
}


fun_score_with_deepbind() {

	reg=$1
	
	mkdir deepbind
	cd deepbind

	#create fasta sequence files from SMRs as input for deepbind
	mutation_cluster_file="$mutation_file.$reg.clust"
	perl -lane '@l=split; $s=join('_',$l[0],$l[1],$l[2],$l[3],$l[4],$l[5],$l[10]); print ">$s\n$l[4]";' <$mutation_cluster_file >$reg.deepbind.seqs 
	
	deepbind_ids="$motif_filtered_db/deepbind_filtered.ids"

	~/soft/deepbind/deepbind --echo $deepbind_ids < $reg.deepbind.seqs >$reg.deepbind.out
	~/soft/deepbind/deepbind --dump-info $deepbind_ids < $reg.deepbind.seqs >deepbind.ids.detail

	#create script to parse deepbind output
	
}


fun_score_with_deepbind_enriched_kmers() {

	reg=$1
	dir_deepbind="/projects_rg/babita/TCGA/mutation/mut_pipeline/subs/subs_new/main_tables/deepbind"
	
	#for all regions except split-sites	
	dir="/projects_rg/babita/TCGA/mutation/mut_pipeline/subs/subs_new/main_tables/enriched_kmers"
	file_enr="$dir/$reg.kmers.zscore.enriched"

	#for Splice_sites
	#dir="/projects_rg/babita/TCGA/mutation/mut_pipeline/subs/subs_new/main_tables/splice_sites/enriched_kmers"
	#file_enr="$dir/$reg.kmers.zscore.enriched.1"	

	#define files
	deepbind_ids="$motif_filtered_db/deepbind_filtered.ids"
	deep_input="$dir_deepbind/$reg.kmers.deepbind.seqs" 
	deep_output="$dir_deepbind/$reg.kmers.deepbind.out" 	
	deep_detail="$dir_deepbind/$reg.kmers.deepbind.ids.detail"


	#create fasta from kmer file using first three columns
	cat $file_enr | sed '1d' | perl -lane '@l=split; $s=join('_',$l[0],"count",$l[1],"zscore",$l[2]); print ">$s\n$l[0]";' >$deep_input

	
	~/soft/deepbind/deepbind --echo $deepbind_ids < $deep_input  >$deep_output
	~/soft/deepbind/deepbind --dump-info $deepbind_ids < $deep_input >$deep_detail

        #parse deepbind output
	#usage: script.R <1.deep.output> <2.deep.seq.input> <3.deep.motif.detail> <4.output_file>
	
	#outputfile:
	deepout_parsed="$file_enr.deepbind"

	$RscriptPath/Rscript $scriptpath/parse_deepbind.R $deep_output $deep_input $deep_detail $deepout_parsed
	
}



regions=("CDS" "3UTR" "5UTR" "EXON" "INTRON" "3SS" "5SS")


for region in "${regions[@]}" ; do

	echo $region

	#1. Create control sets (starts with SMR files)
	fun_create_controls $region
	fun_create_enrichment $region
	fun_create_reverse_enrichment $region
	fun_filter_RNA_motifs $region
	fun_score_with_deepbind $region

	fun_score_with_deepbind_enriched_kmers $region
	
	
done
