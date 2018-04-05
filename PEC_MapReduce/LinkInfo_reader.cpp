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

#include "LinkInfo_reader.hpp"
#include "sequenceUtil.hpp"
#include "string_util.hpp"
#include "stacktrace.hpp"
#include <algorithm>

//constructor
LinkInfo_reader::LinkInfo_reader (string filename) {
  
    this->_hasNext = false;
    
    if (filename == "-") {
        filename = "/dev/fd/0"; // read from stdin
    }
    this->_filereader.open(filename.c_str());

    if (! _filereader.is_open()) {
        throw(stacktrace() + "\n\nError, cannot open file " + filename );
    }
    
    // primer reader to first line 
    getline(this->_filereader, this->_lastline);
   
}

bool LinkInfo_reader::hasNext() {
    bool ret;
    ret = !(this->_filereader.eof());
    return ret;
}


LinkInfo_entry LinkInfo_reader::getNext() {
    
    int direction = 4; 
    bool ret;

    string LINE = this->_lastline;
    ret = !(this->_filereader.eof());
    if (ret == true)
    {
       vector<string> toks;
       string_util::tokenize(LINE, toks, " \t");
           
       direction = atoi(toks[2].c_str()); 
       LinkInfo_entry link(toks[0], toks[1], direction); 
       return(link);
	
       getline(this->_filereader, this->_lastline);

    } else {
       LinkInfo_entry link("", "", direction);
       return(link);
    } 

}
