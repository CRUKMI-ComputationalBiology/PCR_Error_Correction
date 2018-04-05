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

///////////////////////////////////////////////////////////////////////////////////
// This file includes the code from Trinity RNA-Seq Assembly pipeline
// by Grabherr et al. 2011 Nature Biotechnology 2011 29(7):644-652
// and is subject to their original copyright notice copied below:
////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2014, trinityrnaseq
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
//   list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
//
// * Neither the name of the {organization} nor the names of its
//   contributors may be used to endorse or promote products derived from
//   this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
///////////////////////////////////////////////////////////////////////////////////

#include "Fasta_reader.hpp"
#include "sequenceUtil.hpp"
#include "stacktrace.hpp"
#include <algorithm>

//constructor
Fasta_reader::Fasta_reader (string filename) {
  
    this->_hasNext = false;
    
    if (filename == "-") {
        filename = "/dev/fd/0"; // read from stdin
    }
    this->_filereader.open(filename.c_str());
    if (! _filereader.is_open()) {
        throw(stacktrace() + "\n\nError, cannot open file " + filename );
    }
    
    getline(this->_filereader, this->_lastline);
    while ((! this->_filereader.eof()) && this->_lastline[0] != '>') {
        getline(this->_filereader, this->_lastline);
    }
    
    
}

bool Fasta_reader::hasNext() {
    bool ret;
    
        ret = !(this->_filereader.eof());
    return ret;
}


Fasta_entry Fasta_reader::getNext() {
    
    string sequence;
    string header;
    bool ret;

        header = this->_lastline;
        
        ret = !(this->_filereader.eof());
        if (ret == true)
        {
            this->_lastline = "";
            while ((! this->_filereader.eof()) && this->_lastline[0] != '>') {
                getline(this->_filereader, this->_lastline);
                if (this->_lastline[0] != '>') {
                    sequence += this->_lastline;
                }
            }
        }
    
    if (ret == true)
    {
        sequence = remove_whitespace(sequence);
        transform(sequence.begin(), sequence.end(), sequence.begin(), ::toupper);
        Fasta_entry fe(header, sequence);
        return(fe);
    } else {
        Fasta_entry fe("", "");
        return(fe);
    }
}


Fasta_entry Fasta_reader::getNext_FQ() {

    string sequence;
    string header;
    bool ret;

        header = this->_lastline;

        ret = !(this->_filereader.eof());
        if (ret == true)
        {
            this->_lastline = "";
            while ((! this->_filereader.eof()) && this->_lastline[0] != '>') {
                getline(this->_filereader, this->_lastline);
                if (this->_lastline[0] != '>') {
                    sequence += this->_lastline;
                }
            }
        }

    if (ret == true)
    {
        sequence = remove_whitespace(sequence);
        transform(sequence.begin(), sequence.end(), sequence.begin(), ::toupper);
        Fasta_entry fe(header, sequence);
        return(fe);
    } else {
        Fasta_entry fe("", "");
        return(fe);
    }

}


map<string,string> Fasta_reader::retrieve_all_seqs_hash() {
    
    map<string,string> all_seqs_hash;
    
    while (this->hasNext()) {
        Fasta_entry f = this->getNext();
        string acc = f.get_accession();
        string seq = f.get_sequence();
        
        all_seqs_hash[acc] = seq;
    }
    
    return(all_seqs_hash);
}
