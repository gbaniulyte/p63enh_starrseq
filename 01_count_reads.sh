#!/bin/bash
# SBATCH -p batch
# SBATCH --cpus-per-task=4
# SBATCH --mem-per-cpu=8000
# SBATCH --mail-type=ALL
# SBATCH --time=0-23:0 # The job should take 0 days, 23 hours, 0 minutes

source /network/rit/lab/sammonslab/anaconda3/etc/profile.d/conda.sh
conda activate mpra
cd /network/rit/lab/sammonslab/gbaniulyte/MPRA/AM_STARRSeq/final_scripts_fastq

python count_fastq_reads.py p63enh_starrseq_plasmid_1.fastq.gz
python count_fastq_reads.py p63enh_starrseq_plasmid_2.fastq.gz
python count_fastq_reads.py p63enh_starrseq_MCF10A_1.fastq.gz
python count_fastq_reads.py p63enh_starrseq_MCF10A_2.fastq.gz
python count_fastq_reads.py p63enh_starrseq_MCF10Ap53KO_1.fastq.gz
python count_fastq_reads.py p63enh_starrseq_MCF10Ap53KO_2.fastq.gz
python count_fastq_reads.py p63enh_starrseq_HaCaT_1.fastq.gz
python count_fastq_reads.py p63enh_starrseq_HaCaT_2.fastq.gz
python count_fastq_reads.py p63enh_starrseq_SCC25_2.fastq.gz
python count_fastq_reads.py p63enh_starrseq_SCC25_4.fastq.gz
python count_fastq_reads.py p63enh_starrseq_MCF10ApGus.fastq.gz
python count_fastq_reads.py p63enh_starrseq_MCF10ATAp63B.fastq.gz