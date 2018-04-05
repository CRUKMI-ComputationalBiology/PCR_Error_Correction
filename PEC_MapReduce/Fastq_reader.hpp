#ifndef __FASTQ_READER__
#define __FASTQ_READER__

#include <iostream>
#include <fstream>
#include <string>
#include "Fastq_entry.hpp"
#include <map>

using namespace std;

class Fastq_reader {
    
public:
    Fastq_reader(string filename); // constructor
    
    bool hasNext();  // prime for next 
    Fastq_entry getNext(); // retrieves next fasta entry, or NULL if EOF.
    
private:
    ifstream _filereader;
    bool _hasNext;
    string _lastline;
};


#endif
