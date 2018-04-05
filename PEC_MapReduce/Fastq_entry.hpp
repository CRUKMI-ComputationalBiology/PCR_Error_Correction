#ifndef __FASTQ_ENTRY__
#define __FASTQ_ENTRY__

#include <string>

using namespace std;

class Fastq_entry {

 public:

  Fastq_entry() {};
  Fastq_entry(string header, string sequence, string quality);

  string get_accession();

  string get_header();

  string get_sequence();
  
  string get_quality();
  
 private:
  string _accession;
  string _header;
  string _sequence;
  string _quality; 
 
};

#endif

