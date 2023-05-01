#!/bin/bash
# SetBackground.sh
#Code adapted from code on Taylor lab github site

export NA=$1
cd ${NA}

# sort mutations
cat ${NA}_Muts_allSNVs.bed | sort -k1,1 -k2,2n -k3,3n >${NA}_Muts_allSNVs_sorted.bed
cat ${NA}_Muts_allDNVs.bed | sort -k1,1 -k2,2n -k3,3n >${NA}_Muts_allDNVs_sorted.bed
#cat ${NA}_Muts_nostrand_allDNVs.bed | sort -k1,1 -k2,2n -k3,3n >${NA}_Muts_nostrand_allDNVs_sorted.bed

# process transcriptional asymmetry
printf "${NA}_Muts_allSNVs_sorted.bed\n" | perl ../mutsig_orf_bins.pl
printf "${NA}_Muts_allDNVs_sorted.bed\n" | perl ../mutsig_orf_DNVs.pl
#printf "${NA}_Muts_nostrand_allDNVs_sorted.bed\n" | perl ../mutsig_orf_DNVs.pl

#tetranuc count for double mutants
#perl ../format_tetranuc_context.pl <${NA}_Muts_allDNVs_sorted.bed >${NA}_Muts_allDNVs_sorted_tetranuc.bed
#fastaFromBed -s -name -fi ../saccer3_genome.fa -bed ${NA}_Muts_allDNVs_sorted_tetranuc.bed -fo ${NA}_Muts_allDNVs_sorted_tetranuc.fa
#perl ../count_tetranucleotides.pl <${NA}_Muts_allDNVs_sorted_tetranuc.fa >${NA}_Muts_allDNVs_sorted_tetranucount.txt
 

# trinucleotide counts for mutations
printf "${NA}_Muts_allSNVs_sorted.bed\n../saccer3_trinuc_strand_background.txt" | perl ../count_trinuc_context_yeastgenome.pl

# split strands
printf "${NA}_Muts_allSNVs_sorted.bed\n" | perl ../split_strands.pl 

perl ../make_wig_mutations.pl <${NA}_Muts_allSNVs_sorted_minusstrand.bed >${NA}_Muts_allSNVs_minusstrand.wig
printf "${NA}_Muts_allSNVs_minusstrand.wig\n../initial_minus_pyrimidines.wig\n" | perl ../set_background.pl >${NA}_Muts_allSNVs_bk_minus.wig

perl ../make_wig_mutations.pl <${NA}_Muts_allSNVs_sorted_plusstrand.bed >${NA}_Muts_allSNVs_plusstrand.wig
printf "${NA}_Muts_allSNVs_plusstrand.wig\n../initial_plus_pyrimidines.wig\n" | perl ../set_background.pl >${NA}_Muts_allSNVs_bk_plus.wig

