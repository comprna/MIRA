
'''
Create control sets
	For each kmer, match the length and GC content and create 100 random set files for same
'''
from __future__ import division
import sys, os, re
from Bio import SeqIO
import scipy, scipy.cluster
import numpy as np
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from itertools import islice
import operator
from Bio.SeqUtils import GC
import math
import random
import pandas as pd
import numpy as np


def main(argv):

	'''
	Example input parameters:
	#input_mutated_file = "tcga_all_genes_mutations_7mer_filtered.tab.EXON.clust"
	#input_control_file = "tcga_all_genes_mutations_7mer_filtered.tab.EXON.clust.control"
	#len_kmers = 7
	#control_outdir = "./control"
	#prefix = "EXON"
	#number_control_set = 100
	'''

	input_mutated_file = sys.argv[1]
	input_control_file = sys.argv[2]
	len_kmers = sys.argv[3]
	control_outdir = sys.argv[4]
	prefix = sys.argv[5]
	number_control_set = sys.argv[6]
	min_kmer_occurence = sys.argv[7] #atleast 6 or more kmers should be present in the regulated set
	main_outdir = sys.argv[8]

	chk1 = check_files(input_mutated_file)
	chk2 = check_files(input_control_file)


	if len(sys.argv) == 9 and (chk1 and chk2):
		print("Running enrichment analysis for : {}".format(prefix))
		
		if not os.path.exists(control_outdir):
			os.makedirs(control_outdir)
		
		number_control_set = int(number_control_set)
		len_kmers = int(len_kmers)
		min_kmer_occurence = int(min_kmer_occurence)

	else:
		print("Error in processing, parameters not correct or empty files or files not found!!")
		print("\nUsage: python control_enrichment.py <1.mutation_cluster_file> <2.control_cluster_file>\n<3.kmer_length> <4.control_directory> <5.prefix_identifier> \
				<6.No._randomization> <7.mimimum_kmer_occurence_cut-off> <8.output_directory>\n\n")
		sys.exit()



	##create fasta files
	file_list = [input_mutated_file, input_control_file]
	
	'''
	Function 1: Create fasta files of regulated and control clusters
	'''
	write_fasta(file_list)
	print("1.Done regulated and control : fastafile creation")

	##Read regulated fasta file and create dictonary of kmers
	reg_fafile = input_mutated_file + ".fa"
	cont_fafile = input_control_file + ".fa"
	print(reg_fafile)
	
	'''
	Function 2: Create regulated and control dictonary
	'''	
	reg_dict_fa = get_fasta_seq_dictonary(reg_fafile)
	bg_dict_fa = get_fasta_seq_dictonary(cont_fafile)


	print("2.Done regulated and control : dictonary creation")

	'''
	Function 3: Create cluster fasta files for control and regulated sets
	'''	
	create_control_sets(reg_dict_fa, bg_dict_fa, 
						number_control_set, control_outdir, prefix)

	'''
	Note:
	To create regulated and control files for "INTRONS" only in unix
	cat tcga_all_genes_mutations_7mer_filtered.tab.INTRON.entropy_repeats_filtered.clust | awk '{print ">"$1"_"$2"_"$3"_"$4"_"$6"\n"$5}' >INTRON/INTRON.reg.fa
	for ((i=0;i<100;i++)); do outfile="INTRON.cont.fa.$i"; shuf tcga_all_genes_mutations_7mer_filtered.tab.INTRON.clust.control | head -121913 | awk '{print ">"$1"_"$2"_"$3"_"$4"_"$6"\n"$5}' >INTRON/$outfile; echo "done $outfile" ; done
	'''
	print("3.Done regulated and control : cluster fasta sequence")

	reg_clusterfile =  prefix + ".reg.fa"
	reg_clusterfile = os.path.join(control_outdir, reg_clusterfile)
	
	'''
	Function 4: Create kmer counts file for control and regulated sets
	'''	

	create_kmer_counts(reg_clusterfile, len_kmers)


	for i in range(0, number_control_set):
		cont_clusterfile = prefix + ".cont.fa." + str(i) 
		cont_clusterfile = os.path.join(control_outdir, cont_clusterfile)
		create_kmer_counts(cont_clusterfile, len_kmers)

	
	print("4.Done regulated and control : kmer counts")
	

	'''
	Function 5: Create enrichment by comparing kmer count files:
	perform enrichment by comparing each regulated kmer to 100 control kmers, create zscore 
	by counting frequency of kmers 
	'''
	if not os.path.exists(main_outdir):
		    os.makedirs(main_outdir)

	outfile =  prefix + ".kmers.zscore"
	zscore_outfile = os.path.join(main_outdir, outfile)

	perform_enrichment(control_outdir, prefix, number_control_set, zscore_outfile, min_kmer_occurence)


	print("5.Done regulated and control : Enrichment")
	print("Done for Input parameters:\n{}\t{}\n{}\t{}\t{}".format(input_mutated_file,
			input_control_file, len_kmers, control_outdir, prefix))

	print("\n\nEnrichment file created: {}\n\n".format(outfile))
	#quit()



def write_fasta(file_list):
	'''
	Writes fasta files for cluster sequqnce from bedfile
	'''

	for infile in file_list:

		ifi = open(infile, "r")
		out_fa_file = infile + ".fa"

		with open(out_fa_file, "w") as of:
			of.close()
			pass

		of = open(out_fa_file,'a+')

		for line in ifi:
			f_header, f_seq = fasta_header(line)
			'''
			Create outfile
			'''
			of.write("{}\n{}\n".format(f_header, f_seq))

		print("done {}".format(out_fa_file))
		ifi.close()


def fasta_header(line):
	#sequence is in fifth column
	line = line.rstrip()
	line = line.split('\t')
	faline = ">" + '_'.join([line[0], line[1], line[2], line[3], line[5]])
	faseq = line[4].upper()
	return(faline,faseq)


def create_kmer_counts(fafile, len_kmers):
	'''
	Reads fastafile and creates dictonary of kmers with their count
	'''

	fa_seq = []
	fa_id = []
	matched_motifs_sofar = []
	kmers_count = {}
	kmers_count_outfile = fafile + ".kmer"
	wf = open(kmers_count_outfile, 'w')


	#read fasta file
	for seq_record in SeqIO.parse(fafile, "fasta"):
		#fa_id.append(seq_record.id)
		fa_seq.append(seq_record.seq)
		fa_len = len(fa_seq)

	#zipped_fa = [list(x) for x in zip(fa_id, fa_seq)]	

	#print("\n\nScanning fasta sequences for kmers")
	cnt = 0
	kmer_list = []

	for seq in fa_seq:
		cnt += 1
		#_status_bar(cnt, fa_len)

		kmer_list_per_sequence = []

		kmer_list_per_sequence = (["".join(x) for x in slide_windows(seq, len_kmers)])


		for kmer_found in kmer_list_per_sequence:

			if len(kmer_found) == len_kmers:
				
				if kmers_count.has_key(kmer_found) and kmer_found in kmer_list:
					kmers_count[kmer_found] += 1
							
				else:
					kmers_count[kmer_found] = 1
					kmer_list.append(kmer_found)


	#sort dictonary created above on values (kmer counts)
	sorted_kmers = sorted(kmers_count.items(), key=operator.itemgetter(1))

	for kmer, count in sorted_kmers:
		wf.write("{}\t{}\n".format(kmer, count))
	
	wf.close()
	#print("\nkmers-count file created:\t{}\n\n".format(kmers_count_outfile))


#k-mers sliding window
def slide_windows(seq, leng):
	seq = str(seq)
	seq = seq.rstrip()
	seq = seq.replace("--", "")
	it = iter(seq)

	result = tuple(islice(it, leng))

	if len(result) == leng:
		yield result

	for elem in it:
		result = result[1:] + (elem,)
		yield result



def get_fasta_seq_dictonary(fa_file):
	#returns fasta files dictonary for length and gc content

	dict_fa = {}

	for seq_record in SeqIO.parse(fa_file, "fasta"):
		fa_id = seq_record.id
		faseq = seq_record.seq
		gc_count = GC(faseq)
		seq_len = len(faseq)

		#calculate gc content distribution to nearest 10
		gc_content_decimal_distribution = math.floor(gc_count / 10) * 10 #10-bin window
		#gc_content_decimal_distribution = gc_count/seq_len

		dict_fa[fa_id] = [faseq, seq_len, gc_content_decimal_distribution]
	

	return dict_fa


def create_control_sets(reg_dict_fa, bg_dict_fa, number_control_set, outdir, prefix):

	'''
	foreach key and values in regulated dictonary, 
	matches the control set of same length and gc content distribution
	'''
	
	file_pattern = prefix + ".cont.*"
	remove_files(outdir, file_pattern)

	#create outdir for control files, delete all previous files

	reg_file_name = prefix + ".reg.fa" 
	reg_outfile = os.path.join(outdir, reg_file_name)


	with open(reg_outfile, "w") as rf:
		rf.close()
		pass

	notfound = 0
	notfound_list = []
	found = 0

	for keyr, value in reg_dict_fa.items():

		fa_seq, fa_len, fa_gc = value  #regulated
		fa_gc = math.floor(fa_gc / 10) * 10

		#print("{}: {}\t{}\t{}".format(key, fa_seq, fa_len, fa_gc))


		#subset Keys based on regulated gc content and length (+/- 5)
		match_dict = [ x for x, y in bg_dict_fa.items() if ((fa_len - 5 <= y[1] <= fa_len + 5) and (y[2] == fa_gc))]
		
		#for keys get values
		subset_dict = dict([(i, bg_dict_fa[i]) for i in match_dict if i in bg_dict_fa])

		if len(subset_dict) < 2:
			#print("Not enough control subset found for:")
			#print("{} : {}\t{}\t{}\n".format(key, fa_seq, fa_len, fa_gc))
			notfound += 1
			notfound_list.append(fa_seq)

		else:
			for i in range(0, number_control_set):
				keys = subset_dict.keys() 	#List of keys
				random.shuffle(keys)
				
				for key in keys:
					cont_seq = subset_dict[key][0]

				'''
				Write faseq to control file
				'''
				cont_file_name = prefix + ".cont.fa." + str(i)
				cont_outfile = os.path.join(outdir, cont_file_name)
				append_file(cont_outfile, key, cont_seq)

			
			found += 1
			'''
			Write faseq to regulated file
			'''
			append_file(reg_outfile, keyr, fa_seq)

	if notfound :
		print("Total control-match found clusters: {} for : {}".format(found, prefix))
		print("Total control-match NOT-found clusters: {} for : {}".format(notfound, prefix))
		print(notfound_list)
	




def perform_enrichment(outdir, prefix, cont_set, zscore_outfile, min_kmer_occurence):
	'''
	for each kmer in regulated, count how many times it has appear in the control set
	'''
	cont_kmerlist = []

	#read reg kmer file:
	reg_kmerfile =  prefix + ".reg.fa.kmer"
	reg_kmerfile = os.path.join(outdir, reg_kmerfile)

	rf = open(reg_kmerfile,'r')
	wf = open(zscore_outfile, 'w')
	wf.close()

	for i in range(0, cont_set):
		cont_kmerfile = prefix + ".cont.fa." + str(i) + ".kmer"
		cont_kmerfile = os.path.join(outdir, cont_kmerfile) 
		cont_kmerlist.append(cont_kmerfile)

	#read cont kmer files in dataframes
	cont_df = {cont_kmer: pd.read_csv(cont_kmer, sep="\t", header=None, names=['kmer', 'count']) for cont_kmer in cont_kmerlist} 

	for line in rf:
		line = line.rstrip()
		reg_kmer, reg_count = line.split('\t')
		reg_count = int(reg_count)

		#reg_kmer = 'TTTTTTT'
		#reg_count = 85

		if reg_count >= min_kmer_occurence: #arbitrary cut-off, atleast the regulated kmer counts should be above 5
			count_cont = []

			for cont_kmer, df in cont_df.items():
				#print cont_kmer

				subset_kmer = df[ df['kmer'] == reg_kmer]

				if subset_kmer.empty:
					next

				else:
					count_cont.append(int(subset_kmer['count']))
			
			#if list is not empty calculate zscore
			if count_cont:
				bg_expected_mean = np.mean(count_cont)
				bg_expected_sd = np.std(count_cont)	



				try:
					z_score = float(reg_count - bg_expected_mean) / float(bg_expected_sd)
				
				except ZeroDivisionError:

					#this error is because Standard deviation is 0, so zscore should be 0 as the distribution is not normal
					#https://www.quora.com/While-calculating-a-z-score-what-do-you-do-when-standard-deviation-is-zero
					z_score = 0

				'''
				print(min_kmer_occurence)
				print(reg_count)
				print(bg_expected_mean)
				print(z_score)
				'''

				of = open(zscore_outfile,'a+')
				of.write("{}\t{}\t{}\t{}\t{}\n".format(reg_kmer, reg_count, bg_expected_mean, bg_expected_sd, z_score))
				of.close()

	#print("Enrichment file saved: {}".format(zscore_outfile) )
	rf.close()




def remove_files(dir, pat):
	for f in os.listdir(dir):
		if re.search(pat, f):
			os.remove(os.path.join(dir, f))



def append_file(outfile, header, seq):
	of = open(outfile, "a+")
	of.write(">{}\n{}\n".format(header, seq))
	of.close()


def check_files(file):
	'''
	checks if file exists and is not empty, returns true
	'''
	if((os.path.isfile(file)) and (os.stat(file).st_size != 0) ):
		return True
	
	return False



if __name__ == "__main__":
	main(sys.argv)
