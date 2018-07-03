'''
Read maf file and gene file and overlap with bedtools
'''

import sys, os
import subprocess

def main(argv):

	'''
	get input commands
	'''
	bed_file = sys.argv[1]
	mutation_bed_file = sys.argv[2] 
	all_gene_overlap_outfile = sys.argv[3]
	gene_lib_dir = sys.argv[4]

	if len(sys.argv) < 4:
		print("Script to create gene library for mutation analysis")
		print("Usage:")
		print("create_library.py <1.bed_annotation_file> <2.Mutation_file> <3.overlap_output_file> <4.gene_lib_outdir>") 
		sys.exit()

	'''
	Example:	
	gene_file = '/projects_rg/babita/TCGA/mutation/mut_pipeline/testfiles/test_genes_of_interest.bed'
	mutation_file = '/projects_rg/babita/TCGA/mutation/larva/io_files/tcga_fred_mutations.tsv'
	overlap_outfile = '/projects_rg/babita/TCGA/mutation/mut_pipeline/testfiles/tmp_mut_all_gene_overlap'
	gene_lib_dir = '/projects_rg/babita/TCGA/mutation/mut_pipeline/testfiles/gene_lib'
	'''

	#call function
	'''
	Overlap mutation file to gene file coordinates using fjoin
	'''

	check = overlap_fjoin(mutation_file, gene_file, overlap_outfile)
	print("Done overlap_fjoin: file saved: {}".format(overlap_outfile))

	'''
	If output file from fjoin is created, for each gene create mutation library
	'''
	if(check):
		create_gene_library(bed_file, overlap_outfile, gene_lib_dir)

	else:
		print("Overlap didn't work properly, check input files {} and {}".format(mutation_file, gene_file))





def overlap_fjoin(mutation_file, gene_file, output_overlap_file):
	'''
	use fjoin to intersect file with variant (mutation) file
	important: Mutation file should be in bed format
	'''

	output_overlap_file_tmp = output_overlap_file + "tmp"


	os.system('python /soft/bio/sequence/fjoin/fjoin.py -s both -1 {} -2 {} \
				--columns1=1,2,3 --columns2=1,2,3 >{}'.format(mutation_file, gene_file, 
				output_overlap_file_tmp))

	'''
	Fjoin returns an extra first columns of overlap
	remove first column from tmp fjoin output file
	'''

	f = open(output_overlap_file_tmp, "r")
	g = open(output_overlap_file, "w")

	for line in f:
		if line.strip():
			g.write("\t".join(line.split()[1:]) + "\n")

	f.close()
	g.close()

	'''
	#remove temporary files
	'''
	os.remove(output_overlap_file_tmp)

	file_chk = check_files(output_overlap_file)
	return(file_chk) 


def create_gene_library(all_gene_id, overlap_outfile, gene_lib_dir, path_perl):
	'''
	Create gene_lib Directory to save outputs of each gene mutations
	'''
	not_found_genes = []

	if not os.path.exists(gene_lib_dir):
		os.makedirs(gene_lib_dir)

	'''
	Remove previous files, (if any) from gene_lib directory 
	'''
	#os.system('rm {}/*'.format(gene_lib_dir))


	'''
	Extract mutation for each gene

	#cnt = 0
	#num_lines = sum(1 for line in open(gene_file))
	#print("Creating gene mutation libraries\n")

	#with open(gene_file, 'r') as rf:
	#	for line in rf:
	#		gene_id = line.split()[3]

	#		cnt += 1
	#		_status_bar(gene_id, cnt, num_lines)

	#		#grep geneid from mutation file
	'''

	gene_outfile_name = all_gene_id + ".mut"
	gene_outfile_path = os.path.join(gene_lib_dir,gene_outfile_name)

	#os.system('grep -P \'\t{}\t\' {} >{}'.format(
	#	gene_id, overlap_outfile, gene_outfile_path))

	'''
	if "gene.mut" path already exist add 
	'''

	#match the gene ids in mutation file and bedfile provided
	os.system('{}/perl -nle "print if m/\t{}\t/" {} >{}'.format(
		path_perl, all_gene_id, overlap_outfile, gene_outfile_path))

	#print(gene_id)
	check = check_files(gene_outfile_path)

	if not check:
		'''
		Store Genes with no mutations found
		'''
		not_found_genes.append(all_gene_id) #print later optional
		os.remove(gene_outfile_path)
		return False

	return True

	#rf.close()


def check_files(file):
	#if file exists and is not empty return true
	if((os.path.isfile(file)) and (os.stat(file).st_size != 0) ):
		return True
	
	return False

'''
def _status_bar(gene_id, count, count_all):
	a = (float(count)/count_all)*100
	a = int(round(a,2))
	b = int(round(float(a)/2, 0)) #determines the size of bar by printing '='

	sys.stdout.write('\r')
	sys.stdout.write("{} {}/{}".format(gene_id,count,count_all) + "[%-50s] %d%%" % ('='*b, a))
	sys.stdout.flush()
'''

if __name__ == "__main__":
	main(sys.argv)
