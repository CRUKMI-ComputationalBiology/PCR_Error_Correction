
OBJS = Fasta_Splitter_PE.o argProcessor.o string_util.o stacktrace.o Fasta_entry.o Fastq_entry.o Fasta_reader.o Fastq_reader.o
EXECUTABLE = Fasta_Splitter_PE

CXX=mpicxx
CXXFLAG=-Wall -O3 -std=c++0x -m64 -g -Wno-long-long 
COMPILE=$(CXX) $(CXXFLAG) -c

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJS) 
	$(CXX) -o $(EXECUTABLE) $(OBJS)

%.o: %.cpp
	$(COMPILE) -o $@ $<

clean:
	-rm $(OBJS) $(EXECUTABLE)

