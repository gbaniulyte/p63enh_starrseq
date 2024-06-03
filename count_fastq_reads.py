#step 1
#create read count tables with seq and names for each sample
import sys
import pandas as pd
from datetime import datetime
import os
import gzip
os.chdir('/network/rit/lab/sammonslab/gbaniulyte/MPRA/AM_STARRSeq')
ref_path = '/network/rit/lab/sammonslab/gbaniulyte/MPRA/AM_STARRSeq/final_scripts_fastq/tables/p63enh_starrseq_info_long.csv'
fastq_path = '/network/rit/lab/sammonslab/gbaniulyte/MPRA/AM_STARRSeq/final_scripts_fastq/fastq'
out_path = '/network/rit/lab/sammonslab/gbaniulyte/MPRA/AM_STARRSeq/final_scripts_fastq/tables/read_counts'

#fastq.gz file name should be used as first imput argument
fastq = sys.argv[1]
count_names = fastq[fastq.find('q_') + 2:fastq.find('.fastq')] + '_count.csv'

df = pd.read_csv(ref_path, header=0, index_col=False)
seq_id = dict(zip(df.seq, df.id))
#create a dictionary where unique enh seq is key and ID is value, make all uppercase
seq_count = {k[:101].upper():0 for k in seq_id.keys()}
del(seq_id, df)

count_dict = seq_count.copy()
print('Empty counts dict created:', datetime.now().strftime("%H:%M:%S"), flush=True)
total_reads = 0
#read fastq.gz directly
with gzip.open(os.path.join(fastq_path, fastq), 'rt') as r:
	print('Reading:', fastq, datetime.now().strftime("%H:%M:%S"), flush=True)
	for line in r:
		total_reads += 1
		seq = next(r).rstrip() #read next line to get sequence
		try:
			count_dict[seq] += 1 #counts, if key doesn't exist adds seq:0
		except KeyError:
			pass
		next(r),
		next(r) #skip las two line so that for loop starts at the right line again
print('Writing output for individual read counts...', datetime.now().strftime("%H:%M:%S"), flush=True)
with open(os.path.join(out_path, count_names), 'w') as o:
	for key in count_dict.keys():
		o.write("%s,%s\n"%(key,count_dict[key]))
print('Total reads in file:', fastq, total_reads, flush=True)

