#!/bin/bash
# SetBackground.sh
#Code adapted from code on Taylor lab github site

MUT="UVB_New_Combined_All_v2_format_filter.txt"
export NA=$1
cd ${NA}

# count mutations per isolate
perl ../count_UVB_mutations_isolate.pl <${MUT}

# get fasta trinuc context for SNVs
perl ../UVB_process_trinuc_formutations.pl <${MUT}

# extract SNVs
perl ../process_UVB_mutations.pl <${MUT}

# extract DNP mutations 
perl ../UVB_process_DNPmutations.pl <${MUT}

# take care of mutation classes per isolate
perl ../count_trinuc_UVB_mutations_isolate.pl <${MUT}

# count dnp mutations
perl ../count_DNV_UVBmutations_isolate.pl <${MUT}
