#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>

#include "LinkInfo_entry.hpp"

using namespace std;

class LinkInfo_reader {
    
public:
    LinkInfo_reader(string filename); // constructor
    
    bool hasNext();  // prime for next 
    LinkInfo_entry getNext(); // retrieves next fasta entry, or NULL if EOF.
    
private:
    ifstream _filereader;
    bool _hasNext;
    string _lastline;
};

