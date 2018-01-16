'''
NOTE: This Script was previously named: "run_pipeline_clusters.py"
Mutation searching Script 
Input: All mutation bedfile
Output: 
'''

## get gene coordinate
## extract sequence
## 


# gene coordinate: BRCA: chr13: 32889611 - 32973805
# - create bed file with gene coordinate
# cmd: nano BRCA.bed

# extract gene sequence
# cmd: bedtools getfasta -s -fi /projects_rg/babita/ref_genome/hg19.fa -bed BRCA.bed -fo BRCA_chr13_32889611_32973805.fa

#extract all the mutations from VCF files falling in above coordinates
# find overlapping regions 

#cmd : cat /projects_rg/TCGA_vcf/BRCA/genome.wustl.edu_BRCA.IlluminaHiSeq_DNASeq_Cont_automated.Level_2.1.3.0/genome.wustl.edu.TCGA-AN-A046.snv.1.vcf | bedtools intersect -u -a stdin -b BRCA.bed >BRCA.snv.vcf
#extracting from all patients:
#cmd2: rm BRCA.snv.vcf ; for file in /projects_rg/TCGA_vcf/BRCA/genome.wustl.edu_BRCA.IlluminaHiSeq_DNASeq_Cont_automated.Level_2.1.3.0/*.snv.*.vcf ; do echo $file ; bedtools intersect -u -a $file -b BRCA.bed >>BRCA.snv.vcf ; done

# Find mutation frequency by creating kmer windows.
#input: gene overlapped VCF files (column 1,2,4,5) + gene coordinate

#read gene coordinate file
#read VCF file

import sys, os
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import random
import pybedtools
from scipy.stats import binom_test
import pandas as pd
import numpy as np
import math

import gene_class
import create_library

def main(argv):

	
	'''
	#read fasta sequence
	genome_path = './ref_genome/hg19.fa'
	bed_file = './testfiles/Refseq_genes_hg19_new.bed'
	tmp_file = './testfiles/tmp'
	kmer_len = 10
	mutation_bed_file = '/projects_rg/TCGA_Fred/tsv/mutations_v2_formatted.tsv'
	output_gene_lib_dir = './testfiles/gene_lib'
	final_outfile = './testfiles/test_final_score.out'
	'''

	##################################################################################
	##################################################################################

	'''
	IMPORTANT: Before running this script create main mutation refrence file, 
		by overlapping reference bedfile (input_gene_file) with MAF file (input_all_mutation_file), as an output file "*.mut.out"
	
	example:
	python create_mutation_table.py <1.mutation_file.bed> <2.Refseq_gene_annotation.bed>
	output: creates file ->> mutation_file.bed.mut.out

	python create_mutation_table.py src_files/SCLC_mutations_sorted.bed src_files/Refseq_genes_hg19_new.bed
	'''
	#Call this script from cluster
	#path:  /data/users/babita/2016_11_15_mutation_count/script/run_mutation_pipeline2.sh
	#usage:
	#python script.py <bedline> <genome.fa> <mutation_overlap> <out_dir>


	kmer_len = sys.argv[1]					#length of kmers
	mut_count = sys.argv[2]					#minimum number of mutation count to scan (usually 3)
	genome_path = sys.argv[3]				#path of genome file hg19.fa
	input_all_mutation_file = sys.argv[4]	#MAF file or mutation file
	output_gene_lib_dir = sys.argv[5]		#output dir to store all temporary output files
	input_gene_file = sys.argv[6]			#input annotation bed file
	#define paths

	path_bedtools = sys.argv[7]				#/soft/bio/sequence/bedtools-2.25.0/bin
	path_perl = sys.argv[8]					#/soft/devel/perl-5.16.3/bin



	check1 = check_files(input_all_mutation_file)
	check2 = check_files(input_gene_file)

	if len(sys.argv) == 9 and (check1 and check2):

		kmer_len = int(kmer_len)
		mut_count = int(mut_count)

		rf = open(input_gene_file, 'r')

		for line in rf:
			#print("{}".format(line))

			'''
			1. extract gene id and coordinates from bed line
			'''
			line=line.rstrip()
			chr_name, chr_start, chr_stop, gname, val, strand = bedline_split(line) #edit to take gene names with full chr ids (chr_234_456_gene)
			##gname = chr_name + "_" + str(chr_start) + "_" + str(chr_stop) + "_" + name
			name = re.split('_', gname)[3] 

			print("Start analysis for : {}".format(line))

			'''
			2.Create gene mutation
			'''
			check = check_files(input_all_mutation_file)

			if check:

				check2 = create_library.create_gene_library(gname, input_all_mutation_file, output_gene_lib_dir, path_perl)

			else:
				print("Mutation overlap file is empty or not found: {}".format(input_all_mutation_file))


			print("Done step : 2, create_gene_library")

			'''
			If mutation file is empty or doesnot exist, skip following calculation
			'''
			if(check2):
				'''
				create mutation bedfile for each gene
				'''
				gname_in = gname + ".mut"
				mutation_bed_file = os.path.join(output_gene_lib_dir, gname_in)

				gname_out = gname + ".mut.out"
				final_outfile = os.path.join(output_gene_lib_dir, gname_out)

				'''
				create temporary bedfile with gene coordinate to overlap with mutation file
				'''

				tmp_kmers_bed_file3, tmp_overlap_bed_file4, kmeroutfile5 = create_random_tmpfiles_names(output_gene_lib_dir)

				#Condition in case if tmp file with same name already exists in the path
				while os.path.isfile(tmp_kmers_bed_file3) or os.path.isfile(tmp_overlap_bed_file4) or os.path.isfile(kmeroutfile5):
					tmp_kmers_bed_file3, tmp_overlap_bed_file4, kmeroutfile5 = create_random_tmpfiles_names(output_gene_lib_dir)

				'''
				check if tmp file path exists, and if does, create new path for tmp files
				'''

				'''
				read overlap mutation file
				'''

				df = pd.read_csv(mutation_bed_file, header=None, sep='\t')
				total_gene_mut = len(df)  #expected value


				check = create_kmers(line, kmer_len, tmp_kmers_bed_file3) #create kmers of window length kmer_len

				print("Done step : 3, create_kmers")

				
				if(check):				
					check = overlap_kmers(mutation_bed_file, tmp_kmers_bed_file3, tmp_overlap_bed_file4, path_bedtools)
					os.remove(tmp_kmers_bed_file3)

				
				print("Done step : 4, overlap_kmers")

				'''
				Extract kmer sequences usin bedtools
				'''
				
				check = extract_kmer_seq(tmp_overlap_bed_file4, genome_path, kmeroutfile5)
				os.remove(tmp_overlap_bed_file4)

				if(check):

					print("Done step : 5, extract_kmer_seq")
					check = statistics(line, kmeroutfile5, kmer_len, mut_count, total_gene_mut, 
					genome_path,  mutation_bed_file, final_outfile)
					
					os.remove(kmeroutfile5)

					if(check):
						print("Done for gene : {} | file created: {}".format(gname, final_outfile))
					else:
						print("No mutation found on gene kmers; condtion >={} mutations not true for gene: {}".format(mut_count,gname))
						os.remove(final_outfile)
						os.remove(mutation_bed_file)
				else:
					print("File : {} was not created. Problem!".format(tmp_kmers_bed_file3))
					os.remove(kmeroutfile5)
				print("Done step : 6, statistics")

			
			else:
				print("No mutations were found on gene: {} | skipping.".format(line))

			print("Done analysis: {}\n".format(gname))

		rf.close()

	else:

		print("Script Usage:\npython run_pipeline_clusters.py <1.kmer_len> <2.minimum mut.count> <3.genome.file> <4.all_mutation_file> <5.out_dir> <6.input_gene.bed>")


def create_kmers(line, N, tmp_bed):
	
	#from start and end position create kmers of length N, as bedfile coordinates 
	chr_name, chr_start, chr_stop, gene, val, strand = bedline_split(line)

	i = int(chr_start)
	N = int(N)

	# write to file
	f = open(tmp_bed, 'w')

	while i <= int(chr_stop):
		end = i + N
		
		f.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(chr_name, i, end, gene, val, strand)) 
		i = i + 1

	f.close()

	file_chk = check_files(tmp_bed)
	return(file_chk) 

	##send kmer file to overlap with variant (mutation file) using bedtools. 


def overlap_kmers(mutation_file, tmp_kmers_file, tmp_overlap_file, path_bedtools):
	'''
	use bedtools to intersect file with variant (mutation) file
	important: Mutation file should be in bed format
	'''

	'''
	kmers = pybedtools.BedTool(tmp_kmers_file)
	mutations = pybedtools.BedTool(mutation_file)

	#find overlap and write to file
	overlap_mutations = kmers.intersect(mutations, wo = True)
	c = overlap_mutations.saveas(tmp_overlap_file)
	'''

	#Directly use os.system command to intersect bedfiles
	tmp_kmers_file_sort = tmp_kmers_file + ".sort"

	os.system("{}/sortBed -i {} >{}".format(
		path_bedtools, tmp_kmers_file, tmp_kmers_file_sort))

	os.system("{}/intersectBed -sorted -wo -a {} -b {} >{}".format(
		path_bedtools, tmp_kmers_file_sort, mutation_file, tmp_overlap_file))

	os.remove(tmp_kmers_file_sort)
	
	file_chk = check_files(tmp_overlap_file)
	return(file_chk) 


def extract_kmer_seq(kmer_bedfile, hg19_fa_file, outfile_kmers_seq):
	'''
	open bed file to write the kmer sequence
	'''	
	with open(outfile_kmers_seq,'w') as of:
		pass

	wf = open(outfile_kmers_seq,'a+')

	# read names and postions from bed file
	records = SeqIO.to_dict(SeqIO.parse(open(hg19_fa_file), 'fasta'))

	with open(kmer_bedfile) as rf:
		for line in rf:
			name,start,stop,gene,val,strand = bedline_split(line)
			long_seq_record = records[name]
			long_seq = long_seq_record.seq
			alphabet = long_seq.alphabet
			short_seq = str(long_seq)[start-1:stop]
			short_seq_record = SeqRecord(Seq(short_seq, alphabet))
			short_seq_record.seq.strip()

			'''
			If the sequence is in reverse strand, rev complement it
			#edit 14-march-2017 (to extract kmers strictly to half open bed format)
			'''
			if strand == "+":
				kmer_seq = short_seq_record.seq
				#kmer created will be one base greater than the kmer size, trim
				kmer_seq = kmer_seq[:1]

			elif strand == "-":
				kmer_seq = short_seq_record.seq.reverse_complement() 
				kmer_seq = kmer_seq[:-1]
			else:
				print("Strand information not found for line:\n {}".format(line))
				quit()

			wf.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(name,start,stop,gene,kmer_seq,strand))

	rf.close()
	wf.close()

	file_chk = check_files(outfile_kmers_seq)
	return(file_chk) 


def statistics(line, tmp_overlap_file, N, mut_count, mut_on_gene, genome_path, mutation_bed_file, outfile):
	'''
	Test1:
	Binomial test for mutations found on kmers
	x = Number of successes ie No. of mutations found on kmer
	n = Number of trials ie. No. of mutations on whole gene
	p = Expected value ie the hypothesized probability of success. (0<=p<=1; defualt 0.5)
		(len of kmer/len of gene)
	''' 

	'''
	calculate x: no. of mutations found in kmer
	step: Open overlap-kmer file
		  for each kmer line
		  	find overlapping regions in overlap file
	'''

	N = int(N) #kmer_len
	name,start,stop,gene,val,strand = bedline_split(line)
	#gene = re.split('_',gene)[3]
	gene = '_'.join(gene.split('_')[3:])  #for geneids separated by '_' eg. Y_RNA

	#rate of A/T/G/C bases on the whole gene

	rateA, rateT, rateG, rateC = gene_class.some_score(name, start, stop, gene, val, strand, genome_path, mutation_bed_file)

	#read overlap_kmer file
	df = pd.read_csv(tmp_overlap_file, header=None, sep='\t')

	#expected value = length of kmer/ length of gene
	expected_2 = N/float(stop - start)

	i = int(start)

	#write in a file
	with open(outfile,'w') as of:
		pass

	while i <= int(stop):   #for each kmer window of len N
		end = i + N

		mut_on_kmer = len(df[(df[1] == i) & (df[2] == end)])  #subset windows with mutations.
		sub_df = df[(df[1] == i) & (df[2] == end)]

		if mut_on_kmer >= mut_count: #(arbitrary cut-off, ideally number of mutations on the given Kmer should be 3 or more & 1 for control)
		
 			'''
 			perform Binomial test
 			b = scipy.stats.binom(x,n,p)
 			'''
			#extract kmer_sequence
			kmer_seq = pd.unique(sub_df[4].values.ravel())
			kmer_seq = str((kmer_seq)[0]).upper()

			'''
			Calculated base count of each kmer-seq
			'''
			kmer_baseA = kmer_seq.count('A')
			kmer_baseT = kmer_seq.count('T')
			kmer_baseG = kmer_seq.count('G')
			kmer_baseC = kmer_seq.count('C')

			'''qstat
			Calculate expected value for each kmer based on rate of bases on the whole gene
			Find expected value of mutation on a gene
			Si = sequqnce
			exp(Si) = Count_As*rateA + Count_Ts*rateT + Count_Gs*rateG + Count_Cs*rateC
			'''

			expected_1 = (kmer_baseA*float(rateA)) +  (kmer_baseT*float(rateT)) +  (kmer_baseG*float(rateG)) + 	(kmer_baseC*float(rateC))

			'''
			Test 1:
			log2(obs mutation on kmers / exp mutation on gene)
			'''

			try:
				score1 = math.log((mut_on_kmer/float(expected_1)), 2) #for log2

			except ZeroDivisionError:
				score1 = 0
				
			'''
			Test 2:
			Normal Binomial
			'''
			pval2 = binomial_test(mut_on_kmer, mut_on_gene, expected_2)

			with open(outfile,'a') as of:
				of.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(name, i, end, gene, kmer_seq,
						strand, name, start, stop, mut_on_kmer, mut_on_gene, score1, pval2))
			of.close()

		i = i + 1


	#of.close()	

	check = check_files(outfile)
	return check


	
def binomial_test(mut_on_kmer, mut_on_gene, e):
	'''
	performs binomial tests for mutations on each kmer, returns pvalue
	'''
	pval = binom_test(mut_on_kmer, mut_on_gene, e)
	return(pval)


def bedline_split(line):
	'''
	Returns first six fields from line in bed format
	'''

	line = line.rstrip()

	name = re.split('[\t,]', line)[0]
	start = re.split('[\t,]', line)[1]
	stop = re.split('[\t,]', line)[2]
	gene = re.split('[\t,]', line)[3]
	val  = re.split('[\t,]', line)[4]
	strand  = re.split('[\t,]', line)[5]

	start = int(start)
	stop = int(stop)

	return(name, start, stop, gene, val, strand)

def check_files(file):
	'''
	checks if file exists and is not empty, returns true
	'''
	if((os.path.isfile(file)) and (os.stat(file).st_size != 0) ):
		return True
	
	return False


def create_random_tmpfiles_names(output_gene_lib_dir):
	rand = random.random()  
	tmp_file1  = os.path.join(output_gene_lib_dir, ("tmp_" + str(rand) + ".kmers.bed3"))		#outfile 3
	tmp_file2 = os.path.join(output_gene_lib_dir, ("tmp_" + str(rand) + ".overlap.bed4"))		#outfile 4
	tmp_file3 = os.path.join(output_gene_lib_dir, ("tmp_" + str(rand) + ".kmerout.bed5")) 		#outfile 5

	return(tmp_file1, tmp_file2, tmp_file3)



if __name__ == "__main__":
	main(sys.argv)
