import sys, os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import random
import pybedtools
from scipy import stats


#def main(argv):

def some_score(chr_name, chr_start, chr_stop, name, val, strand, genome_path, mutation_bed_file):
	#hg19_path = '/projects_rg/babita/ref_genome/hg19.fa'
	#bed_file = '/projects_rg/babita/TCGA/mutation/mut_pipeline/testfiles/test_genes_of_interest.bed'
	#tmp_bed_file = '/projects_rg/babita/TCGA/mutation/mut_pipeline/testfiles/temp_kmers'
	#mutation_bed_file = '/projects_rg/babita/TCGA/mutation/mut_pipeline/testfiles/test_tcga_mutation.bed'

	#with open(bed_file,'r') as wf:
	#	for line in wf:	
	#		line.rstrip()
	
	###chr_name, chr_start, chr_stop, name, val, strand = line.split()

	gene = Gene(chr_name, int(chr_start), int(chr_stop), name, val, strand)
	
	sequ = gene.sequence(genome_path)

	'''
	create temprory bedfile with gene coordinate to overlap with mutation file
	

	rand = random.random()  
	tmp_bed_file = "tmp_bed" + str(rand) + ".bed"
	with open(tmp_bed_file,'w') as bf:
		bf.write("{}\t{}\t{}\t{}\t{}\t{}".format(chr_name, chr_start, chr_stop, name, val, strand))
	bf.close()

	Get all the mutations in the gene
	send bed file and mutation file to Gene class
	
	tmp_overlap_file = mutations(tmp_bed_file, mutation_bed_file)
	print(tmp_overlap_file)
	#os.remove(tmp_bed_file)
	'''
	'''
	Calculate mutation rate per nucleotide
	R(G) = mG/N(G) #No. of mutations on G / Total no. of Gs
	steps:
		Read output overlapped-bedfile
		Count number of nucleotides found to be mutated
	'''


	base_original_1 = []
	
	with open(mutation_bed_file,'r') as tf:
		for line in tf:	
			line.rstrip()


			'''
			for TCGA and other data where mutated base is at 9th column
			'''
			if(line.rsplit()[0] == "-"):
				base_org =  line.split()[8].reverse_complement()
			else:
				base_org = line.split()[8]
			

			'''
			Only for juanlu's  SCLC data, because mutated base is at 6th column
			
			if(line.rsplit()[0] == "-"):
				base_org =  line.split()[5].reverse_complement()
			else:
				base_org = line.split()[5]
			'''
			
			base_original_1.append(base_org.upper())
			

	tf.close()

	## Conditional: remove substitutions(Mutations) with length more than S from the mutations

	S = 5    #(change, in TCGA its 5(max deletions allowed))
	
	base_original = []
	[base_original.append(x) for x in base_original_1 if len(x) <= S ]
	base_original = ''.join(base_original)

	#print(base_original_1)
	#print(base_original)

	#count mutated nucleotides
	all_mut_As = base_original.count('A')
	all_mut_Ts = base_original.count('T')
	all_mut_Gs = base_original.count('G')
	all_mut_Cs = base_original.count('C')

	'''
	#Find Rate of mutations
	Mutation found on x / Number of x in gene 
	'''
	rateA = (all_mut_As/float(gene.baseA()))
	rateT = (all_mut_Ts/float(gene.baseT()))
	rateG = (all_mut_Gs/float(gene.baseG()))
	rateC = (all_mut_Cs/float(gene.baseC()))

	'''
	Find expected value of mutation on a gene
	Si = sequqnce
	exp(Si) = Count_As*rateA + Count_Ts*rateT + Count_Gs*rateG + Count_Cs*rateC
	'''

	'''
	exp_Si = (gene.baseA()*float(rateA)) +  (gene.baseT()*float(rateT)) +  (gene.baseG()*float(rateG)) + \
			 (gene.baseC()*float(rateC))

	print("BaseA:{},BaseT:{},BaseG:{}, BaseC:{}".format(gene.baseA(), gene.baseT(), gene.baseG(), gene.baseC()))
	#print("RateA:{},RateT:{},RateG:{}, RateC:{}".format(rateA, rateT, rateG, rateC))
	#print("expected mutation: {}".format(exp_Si))

	return(exp_Si)
	'''
	#print("all_mut A:{},T:{},G:{},C:{}".format(all_mut_As,all_mut_Ts,all_mut_Gs,all_mut_Cs))
	#print("BaseA:{},BaseT:{},BaseG:{}, BaseC:{}".format(gene.baseA(), gene.baseT(), gene.baseG(), gene.baseC()))
	#print("RateA:{},RateT:{},RateG:{}, RateC:{}".format(rateA, rateT, rateG, rateC))
	
	#quit()
	return(rateA, rateT, rateG, rateC)


class Gene(object):
	'''
	A class for properties of genes, 
	1) gene length
	2) gene sequence
	3) nucleotide distribution, A/T/G/C
	4) mutations on genes 

	Attributes:
	gene coordinate line in bed format
	'''

	def __init__(self, chrmo, begin,end,name,val,strand):

		'''
		Return gene information in bed format
		'''
		#self.gene_line = gene_line
		#chrm,begin,end,name,val,strand = self.gene_line.split('\t')
		self.chrmo = chrmo
		self.begin = begin
		self.end = end
		self.name= name
		self.val = val
		self.strand = strand

	def length(self):
		'''
		Return length of the gene sequence
		'''
		if self.end > self.begin :
			return(self.end - self.begin)

		else:
			print("Not correct coordinate {}".format(self.chrmo, self.begin, self.end, self.strand))


	def sequence(self, hg19_fa_file):
		'''
		Extract fasta sequence from input genome and coordinates
		'''
		records = SeqIO.to_dict(SeqIO.parse(open(hg19_fa_file), 'fasta'))
		short_seq_records = []
		long_seq_record= []

		long_seq_record = records[self.chrmo]
		long_seq = long_seq_record.seq
		alphabet = long_seq.alphabet
		short_seq = str(long_seq)[self.begin-1:self.end]
		short_seq_record = SeqRecord(Seq(short_seq, alphabet), id=self.chrmo, description='')

		if self.strand == "-":
			self.seq = short_seq_record.seq.reverse_complement()
			self.seq = self.seq.upper() 
		else:
			self.seq = short_seq_record.seq.upper()
		
		return(self.seq)


	def baseA(self):
		'''
		Find nucleotides distribution of As/Ts/Gs/Cs (upper case string)
		'''
		return(self.seq.count('A'))

	def baseT(self):
		return(self.seq.count('T'))

	def baseG(self):
		return(self.seq.count('G'))

	def baseC(self):
		return(self.seq.count('C'))


'''
def mutations(tmp_bed, mutation_file):
	
	#Find total number of mutations on whole gene sequence
	#by overlapping it with Mutation(variant) file using 
	#Bedtools
	
	tmp_outfile = tmp_bed + ".out"
	with open(tmp_outfile, 'w') as wf:
		pass

	kmers = pybedtools.BedTool(tmp_bed)
	mutations = pybedtools.BedTool(mutation_file)

	overlap_mutations = mutations.intersect(kmers, wo = True)

	c = overlap_mutations.saveas(tmp_outfile)
	#print("tmp overlap bedfile: {}".format(c.fn))
	#return overlap mutation file
	return(c.fn)
'''



if __name__ == "__main__":
	main()
