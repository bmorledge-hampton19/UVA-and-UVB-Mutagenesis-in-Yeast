#!/bin/bash
# SetBackground.sh
#Code adapted from code on Taylor lab github site

MUT="WT_ogg1_UVA_All_mutations_comb_unique.txt"
export NA=$1
cd ${NA}

# count mutations per isolate
perl ../count_UVA_mutations_isolate.pl <${MUT}

# get fasta trinuc context for SNVs
perl ../UVA_process_trinuc_formutations.pl <${MUT}

# extract SNVs
perl ../process_UVA_mutations.pl <${MUT}

# count mut type
perl ../count_UVA_mutationclass_isolate.pl <${MUT}

# extract DNP mutations 
#perl ../UVA_process_DNPmutations.pl <${MUT}

