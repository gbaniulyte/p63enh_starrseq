#!/bin/bash
#SBATCH -p batch
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=2000
#SBATCH --mail-type=ALL
#SBATCH --time=0-01:00:00 # The job should take 0 days, 23 hours, 0 minutes

cd ../homer

# #homer motif analysis using p63RE coordinates and +/- 60 nt from the center
# #analysis using 'unchanged' enhancers as background or whole genome

# # 1. activating motifs
# ##unchanged background
findMotifsGenome.pl homer_act_MCF10Ap53KO.txt hg38 homer_act_MCF10Ap53KO_bg_unch -bg homer_unch_MCF10Ap53KO.txt -size 60 -nomotif
findMotifsGenome.pl homer_act_MCF10Ap53KO_fig5.txt hg38 homer_act_MCF10Ap53KO_fig5_bg_unch -bg homer_unch_MCF10Ap53KO_fig5.txt -size 60 -nomotif
findMotifsGenome.pl homer_act_SCC25.txt hg38 homer_act_SCC25_bg_unch -bg homer_unch_SCC25.txt -size 60 -nomotif
findMotifsGenome.pl homer_act_HaCaT.txt hg38 homer_act_HaCaT_bg_unch -bg homer_unch_HaCaT.txt -size 60 -nomotif
# ##genomic background
findMotifsGenome.pl homer_act_MCF10Ap53KO.txt hg38 homer_act_MCF10Ap53KO -size 60 -nomotif
findMotifsGenome.pl homer_act_MCF10Ap53KO_fig5.txt hg38 homer_act_MCF10Ap53KO_fig5 -size 60 -nomotif
findMotifsGenome.pl homer_act_SCC25.txt hg38 homer_act_SCC25 -size 60 -nomotif
findMotifsGenome.pl homer_act_HaCaT.txt hg38 homer_act_HaCaT -size 60 -nomotif


# # 2. repressing motifs
# ##unchanged background
findMotifsGenome.pl homer_rep_MCF10Ap53KO.txt hg38 homer_rep_MCF10Ap53KO_bg_unch -bg homer_unch_MCF10Ap53KO.txt -size 60 -nomotif
findMotifsGenome.pl homer_rep_MCF10Ap53KO_fig5.txt hg38 homer_rep_MCF10Ap53KO_fig5_bg_unch -bg homer_unch_MCF10Ap53KO_fig5.txt -size 60 -nomotif
findMotifsGenome.pl homer_rep_SCC25.txt hg38 homer_rep_SCC25_bg_unch -bg homer_unch_SCC25.txt -size 60 -nomotif
findMotifsGenome.pl homer_rep_HaCaT.txt hg38 homer_rep_HaCaT_bg_unch -bg homer_unch_HaCaT.txt -size 60 -nomotif
# ##genomic background
findMotifsGenome.pl homer_rep_MCF10Ap53KO.txt hg38 homer_rep_MCF10Ap53KO -size 60 -nomotif
findMotifsGenome.pl homer_rep_MCF10Ap53KO_fig5.txt hg38 homer_rep_MCF10Ap53KO_fig5 -size 60 -nomotif
findMotifsGenome.pl homer_rep_SCC25.txt hg38 homer_rep_SCC25 -size 60 -nomotif
findMotifsGenome.pl homer_rep_HaCaT.txt hg38 homer_rep_HaCaT -size 60 -nomotif

# # 3. all motifs in each cell type
findMotifsGenome.pl homer_all_MCF10Ap53KO_fig5.txt hg38 homer_all_MCF10Ap53KO_fig5 -size 60 -nomotif
findMotifsGenome.pl homer_all_HaCaT.txt hg38 homer_all_HaCaT -size 60 -nomotif
findMotifsGenome.pl homer_all_SCC25.txt hg38 homer_all_SCC25 -size 60 -nomotif

# 4. motif enrichment in unchanged group using genome background
findMotifsGenome.pl homer_unch_MCF10Ap53KO.txt hg38 homer_unch_MCF10Ap53KO -size 60 -nomotif
findMotifsGenome.pl homer_unch_MCF10Ap53KO_fig5.txt hg38 homer_unch_MCF10Ap53KO_fig5 -size 60 -nomotif
findMotifsGenome.pl homer_unch_SCC25.txt hg38 homer_unch_SCC25 -size 60 -nomotif
