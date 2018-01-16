'''
Script to create control kmer sequences from same gene, same region

#sample input:

mutation file input:
chr7	591532	591545	PRKAR1B	NA	-	chr7	588833	767313	3953	20	6.91135236784	2.17711396042e-35
chr7	597780	597793	PRKAR1B	NA	-	chr7	588833	767313	3953	20	7.86700676711	2.17711396042e-35
chr7	597959	597972	PRKAR1B	NA	-	chr7	588833	767313	3953	15	7.14488926602	4.63033888639e-25

annotation file input:
chr1	70006	70008	OR4F5	.	+
chr1	134901	135802	AL627309.1	.	-
chr1	137621	138532	AL627309.1	.	-

# cutoff adj.pvalue = 
'''
import sys,os
import random
import re

import extract_seq

def main(argv):

	#clust_file = "/projects_rg/babita/TCGA/mutation/mut_pipeline/main_tables/tcga_all_genes_mutations_7mer_filtered.tab.EXON.clust"
	#annot_file = "/projects_rg/babita/ref_genome/gencode.v19.annotation.gtf.clean.EXON.bed"
	#hg19_fa = "/projects_rg/babita/ref_genome/hg19.fa"
	#rand_set = 100 #no. of times control set created for each region coordinate

	clust_file =  sys.argv[1]
	annot_file = sys.argv[2]
	hg19_fa = sys.argv[3]
	rand_set = sys.argv[4]
	
	rand_set  = int(rand_set)

	chk1 = check_files(clust_file)
	chk2 = check_files(annot_file)
	chk3 = check_files(hg19_fa)

	#print("Input parameters:\n{}\n{}\n{}\n{}".format(clust_file,annot_file,hg19_fa,rand_set))

	if len(sys.argv) == 5 and (chk1 and chk2 and chk3):
		print("Running to create control set for file: {}".format(clust_file))

	else:
		print("Error in processing, parameters not correct or empty files or files not found!!")
		print("Usage: python create_control_set.py <1.mutation_cluster_file> <2.annotation_bedfile> <3.genome_fasta_file> <4.randomization_times_integer>\n\n")
		sys.exit()


	outfile = clust_file + ".control"
	outfile_tmp = outfile + ".tmp"


	with open(outfile_tmp, 'w') as oft:
		pass

	with open(annot_file) as af:
		annot_lines = af.readlines()

	with open(clust_file,'r') as wf:
		
		for clustline in wf:	
			
			#line = "chr7	945074	945087	ADAP1	NA	-	chr7	937536	995043	893.000000	10	6.860514	5.47194007552e-17"
			
			chr_name, chr_start, chr_stop, name, val, strand = bedline_split(clustline)

			'''
			Foreach line bedfile coordinate, create random control regions of matching length 
			'''

			len_clust = int(chr_stop) - int(chr_start)

			'''
			select N random lines of annot file for each line of clust
			'''
			rand_annot_lines = []
			rand_annot_lines =  random.sample(annot_lines, rand_set)

			for bedline in rand_annot_lines:			
				#print(bedline)
				bchr_name, bchr_start, bchr_stop, bname, bval, bstrand = bedline_split(bedline)

				rand_annot_coord_start = random.randint(bchr_start, bchr_stop)

				rand_annot_coord_end = rand_annot_coord_start + len_clust

				#print("{}\t{}".format(rand_annot_coord_start, rand_annot_coord_end ))

				'''
				print to tmp bed file
				'''
				
				with open(outfile_tmp, 'a') as oft:
					oft.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(bchr_name, rand_annot_coord_start, rand_annot_coord_end, bname, bval, bstrand))
				

			oft.close()

	wf.close()
	af.close()

	'''
	Extract sequences for the control regions
	'''
	check = extract_seq.extract_kmer_seq(outfile_tmp, hg19_fa, outfile)

	if check:
		print("Control file created: {}".format(outfile))

	else:
		print("Script couldn't run properly. \n\n")

	'''
	remove tmp file
	'''

	os.remove(outfile_tmp)



def bedline_split(line):
	'''
	Returns first six fields from line in bed format
	'''
	line = line.rstrip()

	name = re.split('[\t ,]', line)[0]
	start = re.split('[\t ,]', line)[1]
	stop = re.split('[\t ,]', line)[2]
	gene = re.split('[\t ,]', line)[3]
	val  = re.split('[\t ,]', line)[4]
	strand  = re.split('[\t ,]', line)[5]

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


if __name__ == "__main__":
	main(sys.argv)