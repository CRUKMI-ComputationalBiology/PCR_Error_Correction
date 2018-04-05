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

#include "LinkInfo_entry.hpp"
#include "string_util.hpp"


// constructor
LinkInfo_entry::LinkInfo_entry(string seq1, string seq2, int direction) {
  
  this->_seq1 = seq1;
  this->_seq2 = seq2;
  this->_direction = direction;
}

string LinkInfo_entry::get_seq1() {
  return(this->_seq1);
}

string LinkInfo_entry::get_seq2() {
  return(this->_seq2);
}
  
int LinkInfo_entry::get_direction() {
  return(this->_direction);
} 
  
