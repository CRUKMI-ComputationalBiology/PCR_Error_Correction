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

#include "Fastq_entry.hpp"
#include "string_util.hpp"

// constructor
Fastq_entry::Fastq_entry(string header, string sequence, string quality) {
  
  // piece apart the header
  
//  if (header[0] == '@') {
//	header.erase(0,1);
//  }
  
  vector<string> toks;
  string acc;
  
  string_util::tokenize(header, toks, " \t");
  
  if (toks.size() > 0) {
	acc = toks[0];
  }
  
  this->_header = header;
  this->_accession = acc;
  this->_sequence = sequence;
  this->_quality = quality;

}

string Fastq_entry::get_accession() {
  return(this->_accession);
}

string Fastq_entry::get_header() {
  return(this->_header);
}

string Fastq_entry::get_sequence() {
  return(this->_sequence);
}

string Fastq_entry::get_quality() {
  return(this->_quality);
}

  
 
  
