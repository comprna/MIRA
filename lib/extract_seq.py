'''
This script extracts fasta sequences for give coordinates in bedfile format

Usage: extract_seq.py <1.bedfile> <2.fasta_file> <3.outfile>

Input:
1) bedfile coordinate
2) fasta_file (example hg19.fa)
3) outfile_name

output:
Creates Kmer sequence on 5th column of given bedfile
'''
import sys,os
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def main(argv):

	bedfile = sys.argv[1]
	fastafile = sys.argv[2]
	outfile = sys.argv[3]

	check1 = check_files(bedfile)
	check2 = check_files(fastafile)

	if len(sys.argv) == 4 and (check1 and check2):

		of = open(outfile, 'w')

		check = extract_kmer_seq(bedfile, fastafile, outfile)

		if check:
			print("sequence file created : {}".format(outfile))
		else:
			print("Error in file creation. File {} not created".format(outfile))

		of.close()

	else:
		print("Error in processing\n")
		print("Usage: extract_seq.py <1.bedfile> <2.fasta_file> <3.outfile>\n\n")
		quit()



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
			name,start,stop,gene,val,strand,rem = bedline_split(line)
			start = start
			long_seq_record = records[name]
			long_seq = long_seq_record.seq
			alphabet = long_seq.alphabet
			short_seq = str(long_seq)[start-1:stop]
			short_seq_record = SeqRecord(Seq(short_seq, alphabet))
			short_seq_record.seq.strip()

			'''
			If the sequence is in reverse strand, rev complement it
			'''
			if strand == "+":
				kmer_seq = short_seq_record.seq
			elif strand == "-":
				kmer_seq = short_seq_record.seq.reverse_complement() 

			else:
				print("Strand information not found for line:\n {}".format(line))
				#quit()

			wf.write("{}\t{}\t{}\t{}\t{}\t{}{}\n".format(name,start,stop,gene,kmer_seq,strand,rem))

	rf.close()
	wf.close()
	of.close()

	file_chk = check_files(outfile_kmers_seq)
	return(file_chk) 


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

	if len(re.split('[\t ,]', line)) > 5:
		rem  = re.split('[\t ,]', line)[6:]
		rem = "\t".join(rem)
		rem = "\t" + rem
	else:
		rem = ""

	start = int(start)
	stop = int(stop)

	return(name, start, stop, gene, val, strand,rem)


def check_files(file):
	'''
	checks if file exists and is not empty, returns true
	'''
	if((os.path.isfile(file)) and (os.stat(file).st_size != 0) ):
		return True
	
	return False




if __name__ == "__main__":
	main(sys.argv)
