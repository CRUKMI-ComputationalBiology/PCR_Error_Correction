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

/////////////////////////////////////////////////////////////////////////////////
// This file includes the code derived from Trinity RNA-Seq Assembly pipeline
// by Grabherr et al. 2011 Nature Biotechnology 2011 29(7):644-652
// and is subject to their original copyright notice copied below:
//////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////
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
///////////////////////////////////////////////////////////////////////////////////////

#include "string_util.hpp"

namespace string_util {
    
    void tokenize(const string& str, vector<string>& tokens, const string& delimiters) {
        
        /* ************************************************************************************/
        /* from: http://oopweb.com/CPP/Documents/CPPHOWTO/Volume/C++Programming-HOWTO-7.html **/
        /**************************************************************************************/
        
        // Skip delimiters at beginning.
        string::size_type lastPos = str.find_first_not_of(delimiters, 0);
        // Find first "non-delimiter".
        string::size_type pos     = str.find_first_of(delimiters, lastPos);
        
        while (string::npos != pos || string::npos != lastPos) {
            // Found a token, add it to the vector.
            tokens.push_back(str.substr(lastPos, pos - lastPos));
            // Skip delimiters.  Note the "not_of"
            lastPos = str.find_first_not_of(delimiters, pos);
            // Find next "non-delimiter"
            pos = str.find_first_of(delimiters, lastPos);
        }
        
        return;
    }
    

    /* rest is my stuff */

    string join(const vector<string>& tokens, const string& delimiter) {
        
        string concat = "";
        
        for (unsigned int i = 0; i < tokens.size(); i++) {
            concat = concat + tokens[i];
            
            if (i != tokens.size()-1) {
                concat = concat + delimiter;
            }
        }

        return(concat);
    }

} // end of string_util namespace

string remove_whitespace (string s) {

  string r = "";

  for (unsigned int i = 0; i < s.length(); i++) {

        char c = s[i];

        if (c != '\t' && c != '\n' && c != ' ') {
          r += c;
        }
  }

  return(r);
}

string remove_whitespace_header (string s) {

  string r = "";

  for (unsigned int i = 0; i < s.length(); i++) {

        char c = s[i];

        if (c != '\t' && c != '\n') {
          r += c;
        }
  }

  return(r);
}

