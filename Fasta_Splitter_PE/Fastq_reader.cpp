/*
  Copyright (C) 2017, Cancer Research UK Manchester Institute

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "Fastq_reader.hpp"
#include "string_util.hpp"
#include "stacktrace.hpp"
#include <algorithm>

//constructor
Fastq_reader::Fastq_reader (string filename) {
  
    this->_hasNext = false;
    
    if (filename == "-") {
        filename = "/dev/fd/0"; // read from stdin
    }
    this->_filereader.open(filename.c_str());
    if (! _filereader.is_open()) {
        throw(stacktrace() + "\n\nError, cannot open file " + filename );
    }
    
    getline(this->_filereader, this->_lastline);
    while ((! this->_filereader.eof()) && this->_lastline[0] != '@') {
        getline(this->_filereader, this->_lastline);
    }
    
    
}

bool Fastq_reader::hasNext() {
    bool ret;
    ret = !(this->_filereader.eof());
    return ret;
}


Fastq_entry Fastq_reader::getNext() {
   
    string quality; 
    string sequence;
    string header;
    bool ret;

        header = this->_lastline;
        
        ret = !(this->_filereader.eof());
        if (ret == true)
        {
            this->_lastline = "";
            while ((! this->_filereader.eof()) && this->_lastline[0] != '+') {
                getline(this->_filereader, this->_lastline);
                if (this->_lastline[0] != '+') {
                    sequence += this->_lastline;
                }
            }

            this->_lastline = "";
            getline(this->_filereader, this->_lastline);
            quality += this->_lastline;

            this->_lastline = "";
            getline(this->_filereader, this->_lastline);

//	    while ((! this->_filereader.eof()) && this->_lastline[0] != '@') {
//                getline(this->_filereader, this->_lastline);
//                if (this->_lastline[0] != '@') {
//                    quality += this->_lastline;
//                }
//            }

        }
    
    if (ret == true)
    {
        sequence = remove_whitespace(sequence);
        transform(sequence.begin(), sequence.end(), sequence.begin(), ::toupper);
	header   = remove_whitespace_header(header);
        Fastq_entry fe(header, sequence, quality);
        return(fe);
    } else {
        Fastq_entry fe("", "", "");
        return(fe);
    }
}



