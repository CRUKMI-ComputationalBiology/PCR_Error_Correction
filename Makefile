###################################################################
#
# The default compiler is GNU mpic++.
# Run
#   make 
# to build MR-Inchworm and Fasta_Splitter 
# 

TARGETS=pec fastaSP bwampi

all: ${TARGETS}
	sh install_tests.sh

fastaSP:
	cd Fasta_Splitter_PE && $(MAKE)

bwampi:
	cd BWA && $(MAKE)

pec:
	cd PEC_MapReduce && $(MAKE) -f Makefile.mpicc

clean:
	cd PEC_MapReduce && $(MAKE) -f Makefile.mpicc clean
	cd Fasta_Splitter_PE && $(MAKE) clean 
	cd BWA && $(MAKE) clean
