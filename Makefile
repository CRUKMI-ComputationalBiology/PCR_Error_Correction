###################################################################
#
# The default compiler is GNU mpic++.
# Run
#   make 
# to build MR-Inchworm and Fasta_Splitter 
# 

TARGETS=fastaSP bwampi zlib mrmpi pec

all: ${TARGETS}
	sh install_tests.sh

fastaSP:
	cd Fasta_Splitter_PE && $(MAKE)

bwampi:
	cd BWA && $(MAKE)

zlib:
	cd LIBs && mkdir -p ZLIB
	cd LIBs/zlib-1.2.11 && ./configure
	cd LIBs/zlib-1.2.11 && make test
	cd LIBs/zlib-1.2.11 && make install prefix=../ZLIB	

mrmpi:
	cd LIBs/mrmpi-7Apr14/src/ && make mpicc

pec:
	cd PEC_MapReduce && $(MAKE) -f Makefile.mpicc

clean:
	cd LIBs/mrmpi-7Apr14/src/ && make clean-mpicc
	cd LIBs/mrmpi-7Apr14/src/ && rm libmrmpi_mpicc.a
	cd LIBs && rm -r ZLIB
	cd PEC_MapReduce && $(MAKE) -f Makefile.mpicc clean
	cd Fasta_Splitter_PE && $(MAKE) clean 
	cd BWA && $(MAKE) clean
