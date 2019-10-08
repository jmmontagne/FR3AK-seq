#!/usr/bin/env python

#Algorithm to generate maximally compact primer set in FR3 region of TCRbeta chain via incorporation of universal base inosine
#Input FASTA, output list of primers

#Prompt for input FASTA file from command line

fasta_file = raw_input("Enter FASTA file name as example.fa: ") #file must be in current directory OR use full path name

print fasta_file #test for input seq

import sys
import re
import os
import numpy as np
import pandas as pd
from Bio import pairwise2

linenum = 0
tempfile = sys.argv[0] + '.tmp' #name tempfile according to fasta input file
with open(tempfile, 'w') as outfile:
	for line in open(fasta_file, 'r'):
		line = line.strip()
		if len(line) > 0:
			if line[0] == '>': #matches header line
				if linenum == 0:	
					header = re.match('^>.*', line)
					outfile.write(header.group() + '$') #Get header
					linenum +=1
				else:
					header2 = re.match('^>.*', line)
					outfile.write('\n' + header2.group() + '$')
			else:
				outfile.write(line)

#generate dictionary using tempfile	
seqs_1 = {}
with open(tempfile, 'r') as seqsfile:
	for line in seqsfile:
		line = line.split('$')
		seqs_1[line[0][1:]] = line[1:]

os.remove(tempfile) #remove tempfile

#for key, value in seqs_1.items(): #TEST
#	print key #TEST
#	print value #TEST

#print len(seqs_1) #TEST

#convert string values to lists
for key, value in seqs_1.items():
        seqs_1[key] = list(seqs_1[key][0])
	if seqs_1[key][-1] == "\n": #remove /n if present
		seqs_1[key].pop(-1)
	else:
		next
	#Test of previous, should be 111 in length
	#print len(seqs_1['L36092A|'])
	#print seqs_1['L36092A|']
	#print len(seqs_1['AF009660|'])
	#print seqs_1['AF009660|']

#print seqs_1
#Test of previous, should be 20 in length
#print len(seqs_1['L36092A|'])
#print len(seqs_1['AF009660|'])

#Put dictionary into pandas dataframe
seqs_df = pd.DataFrame.from_dict(seqs_1, orient='index')
print seqs_df

#Get most frequent nt at each position, set position with highest freq to specific nt
seq = []
freq_seq = []
headers = []
primer = [0] * 20
#print primer
remaining_seqs = len(seqs_df)
#print remaining_seqs
idx = []
two_from_last = [] #indices of last 2 positions for regex filtering
final_primers = {}
mismatch_index = {}
primer_coverage = {}
final_primer_coverage_headers = {}
primer_num = 1
primer_pos_filled = -1
nt_pos = 0
seqs_df2 = seqs_df
#print seqs_df2

while len(seqs_df) > 0:	#loop until only 1 sequence left
	while remaining_seqs > 0:
		for c in seqs_df2.columns:# get nt frequencies for each column of dataframe
                	freq = seqs_df2[c].value_counts().reset_index().rename(columns={'index': 'nt', c: 'ct'})  #put nt and count of nt into dictionary
#                       print freq
                        max_from_col = freq.loc[freq['ct'].idxmax()]
#                       print max_from_col
#                       print type(max_from_col)
                        max_freq_nt =  max_from_col.get('nt')
                        max = max_from_col.get('ct')
			freq_seq.extend(max_freq_nt) #add nt to list
#               	print freq_seq
                	seq.extend([max]) #add frequency to list
		int_seq = map(int, seq)#convert frequency to integer
#		print int_seq
		arr = np.asarray(int_seq)#put int_seq into array
#		print arr
		mask = np.zeros(arr.size, dtype=bool)
		mask[idx] = True #assign positions filled in primer to true
#		print mask
		seq_mask = np.ma.array(arr, mask=mask) #masks positions already filled in primer
#		print seq_mask
		index = np.argmax(seq_mask)
		max_freq = np.amax(seq_mask)
#		print max_freq
#               print seq
		if primer_pos_filled <15:#keep placing highest frequency nucleotide in primer until position 15 filled
			idx.extend([index])#add this index to list of filled indices
#			print idx
			third_from_last = index
#			print third_from_last
#			print primer[index]
#			print freq_seq[index]
			primer[index] = freq_seq[index]
#			print primer
			primer_pos_filled +=1
			seqs_df2 = pd.DataFrame(seqs_df2.loc[seqs_df2[index] == freq_seq[index]])#filter out sequences that have this nt at that position for next iteration
#			print seqs_df
#			print seqs_df2
			remaining_seqs = len(seqs_df2) #number of sequences left to cover in primer design
#			print remaining_seqs
#			print primer_pos_filled
			freq_seq = []
                	seq = []
			freq = {}
		elif primer_pos_filled ==15: #if len of remaining seqs = frequency of most abundant nucleotide, set to that nucleotide
			if remaining_seqs == max_freq:
                       		idx.extend([index])
#                      		print idx
                       		third_from_last = 'NoMM'
#                      		print third_from_last
                       		primer[index] = freq_seq[index]
#                      		print primer
                                seqs_df2 = pd.DataFrame(seqs_df2.loc[seqs_df2[index] == freq_seq[index]])#filter out sequences that have this nt at that position for next iteration
                        	remaining_seqs = len(seqs_df2) #number of sequences left to cover in primer design  
		     		primer_pos_filled +=1
#                      		print primer_pos_filled
                       		freq_seq = []
                       		seq = []
                       		freq = {}
			elif remaining_seqs != max_freq: #if mismatch here, remaining 3 positions also have mismatches = inosines
                        	idx.extend([index])
#                       	print idx
                        	third_from_last = index
#                       	print third_from_last
                        	primer[index] = freq_seq[index]
#	                       	print primer
                                seqs_df2 = pd.DataFrame(seqs_df2.loc[seqs_df2[index] == freq_seq[index]])#filter out sequences that have this nt at that position for next iteration
        	                remaining_seqs = len(seqs_df2) #number of sequences left to cover in primer design
	                	primer_pos_filled +=1
#                       	print primer_pos_filled
                        	freq_seq = []
                        	seq = []
                        	freq = {}
				break
		elif primer_pos_filled ==16: #if len of remaining seqs = frequency of most abundant nucleotide, set to that nucleotide
			if remaining_seqs == max_freq:
                        	idx.extend([index])
#                       	print idx
                        	third_from_last = 'NoMM'
#                       	print third_from_last
                        	primer[index] = freq_seq[index]
#                       	print primer
                                seqs_df2 = pd.DataFrame(seqs_df2.loc[seqs_df2[index] == freq_seq[index]])#filter out sequences that have this nt at that position for next iteration
                        	remaining_seqs = len(seqs_df2) #number of sequences left to cover in primer design
		        	primer_pos_filled +=1
#                       	print primer_pos_filled
                        	freq_seq = []
                        	seq = []
                        	freq = {}
			elif remaining_seqs != max_freq: #if mismatch here, remaining  position also has mismatch = inosine
		        	break
		elif primer_pos_filled ==17: #if len of remaining seqs = frequency of most abundant nucleotide, set to that nucleotide
                	if remaining_seqs == max_freq:
                        	idx.extend([index])
#                       	print idx
                        	third_from_last = 'NoMM'
#                       	print third_from_last
                        	primer[index] = freq_seq[index]
#                       	print primer
		                seqs_df2 = pd.DataFrame(seqs_df2.loc[seqs_df2[index] == freq_seq[index]])#filter out sequences that have this nt at that position for next iteration
                        	remaining_seqs = len(seqs_df2) #number of sequences left to cover in primer design
#                       	print remaining_seqs
                        	primer_pos_filled +=1
#                       	print primer_pos_filled
                        	freq_seq = []
                        	seq = []
                        	freq = {}
			elif remaining_seqs != max_freq: #if mismatch here, remaining  position also has mismatch = inosine
                        	break
		elif primer_pos_filled ==18: #if len of remaining seqs = frequency of most abundant nucleotide, set to that nucleotide
                	if remaining_seqs == max_freq:
                        	idx.extend([index])
#                       	print idx
                        	third_from_last = 'NoMM'
#                       	print third_from_last
                        	primer[index] = freq_seq[index]
#                       	print primer
                        	seqs_df2 = pd.DataFrame(seqs_df2.loc[seqs_df2[index] == freq_seq[index]])#filter out sequences that have this nt at that position for next iteration
                        	remaining_seqs = len(seqs_df2) #number of sequences left to cover in primer design
				primer_pos_filled +=1
#                       	print primer_pos_filled
                        	freq_seq = []
                        	seq = []
                        	freq = {}
			elif remaining_seqs != max_freq: #if mismatch here, remaining  position also has mismatch = inosine
                        	break

		else:
			break
#	print primer
	for i, n in enumerate(primer):
		if n==0:
			primer[i] = 'i'
			two_from_last.extend([i])
#	print two_from_last
#	print primer
	final_primers[primer_num] = primer
	mismatch_index[primer_num] = third_from_last
#	print final_primers
#	print mismatch_index

	primer_filter = list(primer)
#	print primer_filter

#	Pairwise alignment method to filter out sequences
#	convert seqs_df dataframe into dictionary
#	seqs_df = pd.DataFrame.from_dict(seqs_1, orient='index')
	seqs_df_test = seqs_df
#	print seqs_df_test
#	print pd.DataFrame.transpose(seqs_df_test)
	pw_aln = pd.DataFrame.transpose(seqs_df_test).to_dict('list')
	pw_aln2 = pd.DataFrame.transpose(seqs_df_test).to_dict('list')
#	print pw_aln
#	print type(primer_filter)
	filter_dict = {}
	
	for key, value in pw_aln.iteritems():
		for i in two_from_last:
			pw_aln[key][i] = 'i'

#	print pw_aln
#	print two_from_last
#	for key, value in pw_aln.iteritems():
#		print key, value
#	print primer

	for key, value in pw_aln.iteritems():
#		print pw_aln[key][-5:]
		if pw_aln[key][-5:] == primer_filter[-5:]:
#               	print key
#               	print value
                	alignment = pairwise2.align.globalms(primer_filter, value, 1, 0, -1000, -1000, gap_char=['-'], score_only=True)
#               	print alignment
                	filter_dict[key] = alignment #add alignment score to dictionary with header as key
		else:
			next
	for key, value in filter_dict.iteritems():
#		print key, value
		if value >=19:
			del pw_aln2[key]
		else:
			next

	seqs_df3 = pd.DataFrame.from_dict(pw_aln2, orient='index')
#	print seqs_df3
	primer_coverage[primer_num] = len(seqs_df)-len(seqs_df3)
#	print primer_coverage
	df_combo = [seqs_df, seqs_df3]
	df_concat = pd.concat(df_combo)
	df_concat_final = df_concat.drop_duplicates(keep = False)
#	print df_concat_final
	primer_coverage_headers = pd.DataFrame.transpose(df_concat_final).to_dict('list')
#	print primer_coverage_headers
	for keys, values in primer_coverage_headers.iteritems():
		headers.append(keys)
	final_primer_coverage_headers[primer_num] = headers
	primer_num +=1
	seqs_df = seqs_df3
	primer = [0]*20
	idx = []
#	print len(seqs_df)
	two_from_last = []
	primer_pos_filled = -1
	freq_seq = []
	seq = []
	pw_aln = {}
	pw_aln2 = {}
	headers = []
	seqs_df2 = seqs_df
	remaining_seqs = len(seqs_df)
#	print seqs_df2
#	print remaining_seqs
#	print seqs_df

seqs_df4 = pd.DataFrame.from_dict(final_primers, orient='index')
#print seqs_df4
seqs_df5 = seqs_df4[seqs_df4.columns[::-1]]
print "Final primer sequences:", seqs_df4
print "Number of sequences covered by each primer",primer_coverage
print "Index of mismatch position:", mismatch_index
print "Headers of primers covered:", final_primer_coverage_headers
