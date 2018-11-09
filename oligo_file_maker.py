#!/usr/bin/env python
'''
This script reads information about qPCR primers and probes from a .csv file (i.e. qPCR_assay_oligos.csv)
and creates the various oligonucleotide files required by the qPCR Screener pipeline (i.e. oligonucleotide
definition csv file, primersearch tsv file of F/R primers, and fasta file of primer and probe.
    usage: python oligo_file_maker.py qPCR_assay_oligos.csv
    Author: Diane Eisler, Molecular Microbiology & Genomics, BCCDC Public Health Laboratory, Nov 2018
'''
import sys,string,os, time, Bio, re
from Bio import Seq, SeqIO, SeqUtils, Alphabet, SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import Bio.Data.IUPACData

def make_primersearch_tsv_file(assayName, FprimerName, FprimerSeq, RprimerSeq):
	'''Create tsv of forward and reverse primers for input into primersearch.'''
	fileName = assayName + ".tsv"
	primerName = ''
	
	if FprimerName.endswith('_F'): #ensure no '_F' suffix on forward primer name
		endpoint = FprimerName.index('_F')
		primerName = FprimerName[0:endpoint]
	else:
		primerName = FprimerName
	#open tsv file and write tab-separated line of primers
	oligo_tsv = open(fileName, 'w')
	oligo_line = primerName + '\t' + FprimerSeq + '\t' + RprimerSeq
	oligo_tsv.write(oligo_line)

def make_SeqRecord(name, sequence):
	'''Return SeqRecord object from oligo name and sequence.'''
	rec = SeqRecord(Seq(sequence, alphabet = IUPAC.ambiguous_dna), id = name, description = "")
	return rec

def make_oligo_fasta_file(assayName, seqRecord_list):
	'''Write SeqRecords for primers and probe to fasta.'''
	outFastaHandle = os.path.join(assayName + "_oligos.fasta")
	print(outFastaHandle)
	outFasta = open(outFastaHandle, 'w')
	SeqIO.write(seqRecord_list, outFasta, "fasta")
	#output_fasta.close()

def make_row(param_list):
	row = ','.join(param_list) #place comma between items in list
	return row

def make_oligo_def_csv_file(assayName, FprimerName, FprimerSeq, probeName, probeSeq, probeOrient, RprimerName, RprimerSeq):
	'''Make an oligonucleotide definition file from a list of oligo parameters and output to csv format.'''
	outputHandle = "%s_oligo_definitions.csv" % (assayName)
	oligo_def_file = open(outputHandle, 'w')
	#define the rows to print to csv
	headers_line = "Oligonucleotide_name,Nucleotide_sequence,Type(primer/probe),Strand_orientation(F/R)"
	f_row = make_row([FprimerName, FprimerSeq, "primer", "F"])
	p_row = make_row([probeName, probeSeq, "probe", probeOrient])
	r_row = make_row([RprimerName, RprimerSeq, "primer", "R"])
	rows = [headers_line, f_row, p_row, r_row] #list of rows to write to csv
	#write each row to csv and close file
	for r in rows:
		oligo_def_file.write(r)
		oligo_def_file.write('\n')
	oligo_def_file.close()

def append_suffix_to_oligoName(oligoName, suffix):
	'''Checks that forward primer name is appended with '_F', reverse primer name is appended with '_R', and probe
	name is appended with '_P' '''
	name = oligoName
	if name.endswith(suffix):
		return name
	else:
		return name + suffix
	return

def main(master_oligo_file):
	'''Reads in lines from the master csv file containing all primers and probes and creates a
	tsv, csv, and fasta file for each assay.'''
	with open(master_oligo_file, 'r') as master_file:
		master_file.readline() #skip header line
		lines = master_file.readlines() #read subsequent lines into a list
		oligo_SeqRecords = [] ##list to hold oligo SeqRecord objects

		for line in lines: #read each row and split into a list of elements separated by commas
			row_elements = line.rstrip().split(',')
			[assayName, FprimerName, FprimerSeq, probeName, probeSeq, probeOrient, RprimerName, RprimerSeq] = row_elements[:]
			#ensure oligo names have correct suffixes appended
			FprimerName = append_suffix_to_oligoName(FprimerName, '_F')
			RprimerName = append_suffix_to_oligoName(RprimerName, '_R')
			probeName = append_suffix_to_oligoName(probeName, '_P')
			#create SeqRecords from each oligo and add these to list of oligo SeqRecords
			Fprimer = make_SeqRecord(FprimerName, FprimerSeq)
			Rprimer = make_SeqRecord(RprimerName, RprimerSeq)
			Probe = make_SeqRecord(probeName, probeSeq)
			#add oligo SeqRecords to list and make fasta file from them
			oligo_SeqRecords = [Fprimer, Probe, Rprimer]
			make_oligo_fasta_file(assayName, oligo_SeqRecords)
			#create a primersearch tsv file 
			make_primersearch_tsv_file(assayName, FprimerName, FprimerSeq, RprimerSeq)
			#make oligo definition file for assay
			make_oligo_def_csv_file(assayName, FprimerName, FprimerSeq, probeName, probeSeq, probeOrient, RprimerName, RprimerSeq)

	return


if __name__ == '__main__':
	master_oligo_file = sys.argv[1]
	main(master_oligo_file)