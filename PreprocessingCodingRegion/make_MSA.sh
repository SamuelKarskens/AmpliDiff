#!/bin/sh
SEQFOLDER="$1"
GENENAME="$2"
FILE="$SEQFOLDER"/"$GENENAME"_aligned.fasta

##Check if alignment already exists in same folder as sequences, otherwise perform MSA with MAFFT and output:
##  - sequences_aligned.fasta (MSA)
##  - runtime_alignment.txt (Runtime for MSA)

if [ -f "$FILE" ]; then
	echo "$FILE exists."
else
	SECONDS=0
	echo "$FILE does not exist."
	mafft --auto --leavegappyregion "$SEQFOLDER"/"$GENENAME".fasta > "$SEQFOLDER"/"$GENENAME"_aligned.fasta
	echo "Elapsed: $(($SECONDS / 3600)) hours, $((($SECONDS / 60) % 60)) minutes, $((SECONDS % 60)) seconds" > "$SEQFOLDER"/runtime_"$GENENAME"_alignment.txt
fi