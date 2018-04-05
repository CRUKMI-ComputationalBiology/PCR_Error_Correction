#!/bin/bash

echo "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo ""
echo 'Performing Unit Tests of Build'
echo ' '
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"

if [ -e "PEC_MapReduce/pec_algorithm" ]
then
	echo "PEC tool:      has been Installed Properly"
else
	echo "PEC tool Installation appears to have FAILED"
fi
if [ -e "Fasta_Splitter_PE/Fasta_Splitter_PE" ]
then
	echo "Fasta Splitter:        has been Installed Properly"
else
	echo "C Fasta Splitter Installation appears to have FAILED"
fi
if [ -e "BWA/BWA" ]
then
        echo "BWA:             has been Installed Properly"
else
        echo "BWA Installation appears to have FAILED"
fi
