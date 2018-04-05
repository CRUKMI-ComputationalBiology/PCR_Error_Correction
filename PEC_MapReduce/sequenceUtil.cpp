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
/////////////////////////////////////////////////////////////////////////////////////


#include "sequenceUtil.hpp"
#include <stdlib.h>
#include "stacktrace.hpp"
#include <math.h>
#include <iostream>
#include <sstream>
#include <algorithm>

static const int MAX_LINE_LENGTH = 500000;

int check_collision_base(int cbase, int base_collision) {

   int base_check;
   int Base = 0;

   switch(cbase) {

	case 0 : {

		base_check = base_collision & 7;
		if((base_check & 1) > 0) Base = Base | 4;
		if((base_check & 2) > 0) Base = Base | 2;
		if((base_check & 4) > 0) Base = Base | 1;

		return Base;	
		}

	case 1 : {

		base_check = base_collision & 56;
		base_check = base_check >> 3;
		if((base_check & 1) > 0) Base = Base | 8;
                if((base_check & 2) > 0) Base = Base | 2;
                if((base_check & 4) > 0) Base = Base | 1;	
	
		return Base;
		}	

	case 2 : {

		base_check = base_collision & 448;
		base_check = base_check >> 6;
		if((base_check & 1) > 0) Base = Base | 8;
                if((base_check & 2) > 0) Base = Base | 4;
                if((base_check & 4) > 0) Base = Base | 1;

                return Base;
		}

	case 3 : {

		base_check = base_collision & 3584;
		base_check = base_check >> 9;
                if((base_check & 1) > 0) Base = Base | 8;
                if((base_check & 2) > 0) Base = Base | 4;
                if((base_check & 4) > 0) Base = Base | 2;

                return Base;
		}

   }

}

int _consen_alt_allele(int consen, int alt) {

	switch(consen) {

	   case 0 : {
			switch(alt) {
				case 1 : return 1;
				case 2 : return 2;
				case 3 : return 3;
			}

		    }

	   case 1 : {
                        switch(alt) {
                                case 0 : return 4;
                                case 2 : return 5;
                                case 3 : return 6;
                        }

		    }

	   case 2 : {
                        switch(alt) {
                                case 0 : return 7;
                                case 1 : return 8;
                                case 3 : return 9;
                        }
		    }

	   case 3 : {
                        switch(alt) {
                                case 0 : return 10;
                                case 1 : return 11;
                                case 2 : return 12;
                        }

	   	    }
	}

}

unsigned char _base_to_Qscore [256] = {
                255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, //   0
                255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,   0,   1,   2,   3,   4,   5,   6, //  20
                  7,   8,   9,  10,  11,  12,  13,  14,  15,  16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26, //  40
                 27,  28,  29,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42, 255, 255, 255, 255, //  60
                255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, //  80
                255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 100
                255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 120
                255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 140
                255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 160
                255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 180
                255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 200
                255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 220
                255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255                      // 240
};
//               0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19

int get_Qscore(char quality) {
    int qscore = _base_to_Qscore[quality];
    if(qscore < 0)  qscore = 0;
    if(qscore > 42) qscore = 42;
    return qscore;
}

char get_QualityAscii(int qscore) {
     if(qscore < 0)  qscore = 0;
     if(qscore > 42) qscore = 42;
     int ascii_intval = qscore + 33;
     char Q = static_cast<char>(ascii_intval);
     return Q;
}

char _int_to_base [4] = {'G', 'A', 'T', 'C'};
unsigned char _base_to_int [256] = {
                255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, //   0
                255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, //  20
                255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, //  40
                255, 255, 255, 255, 255,   1, 255,   3, 255, 255, 255,   0, 255, 255, 255, 255, 255, 255, 255, 255, //  60
                255, 255, 255, 255,   2, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,   1, 255,   3, //  80
                255, 255, 255,   0, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,   2, 255, 255, 255, // 100
                255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 120
                255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 140
                255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 160
                255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 180
                255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 200
                255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 220
                255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255                      // 240
};
//               0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19

int get_int_from_base(char base) {
   return _base_to_int[base];
}

char get_base_from_int(int intval) {
  return _int_to_base[intval];
}

int count_at(string seq) {
  int at = 0;
  for(int i = 0; i < seq.size(); i++)
    if(seq[i] == 'A' || seq[i] == 'T')
      at +=  1;
  return at;
}

void count_gatc(double **gatc, string seq) {
  for(int i=0; i<seq.length(); i++) {
	if((_base_to_int[seq[i]]>=0) && (_base_to_int[seq[i]]<4)) gatc[i][ _base_to_int[seq[i]] ] += 1.0;
  }
}

void convert_cigar2bam(int32_t cigar, int32_t *output) {
	for(int i=0; i<2; i++) output[i] = 0;
        output[0] = cigar & (int32_t)15;
        cigar = cigar >> 4;
        output[1] = cigar & (int32_t)2147483647; // 2^(32-1) - 1 = 2147483647
}

int CigarOp_parse(char* s,vector<CigarOp> &v)
{
        char* p=(char*)(s);
        while(*p!=0)
            {
            char* endptr;
            CigarOp c;

            if(!isdigit(*p)) return -1;
            c.size =(int)strtol(p,&endptr,10);
            if(c.size<=0) return -1;
            p=endptr;

            c.op = p[0];
            if(!isalpha(c.op)) return -1;
            v.push_back(c);

            ++p;
            }
        return 0;
}


int smith_waterman(string seq_a, string seq_b, double mu, double delta, int *loc_a, int *loc_b, char *consensus_a, char *consensus_b) {

  int ind;
  int N_a = seq_a.length();                     // get the actual lengths of the sequences
  int N_b = seq_b.length();

  double **H = (double **) malloc(sizeof( *H ) * (N_a+1));
  if(H) { for(int i=0; i<(N_a+1); i++) H[i] = (double *) malloc(sizeof( *H[i] ) * (N_b+1)); }

  for(int i=0;i<=N_a;i++)
    for(int j=0;j<=N_b;j++) H[i][j]=0.0;

  double *temp = (double *) malloc(sizeof(double)*4);
  double **I_i = (double **) malloc(sizeof( *I_i ) * (N_a+1));
  double **I_j = (double **) malloc(sizeof( *I_j ) * (N_a+1));

  if(I_i) { for(int i=0; i<(N_a+1); i++) I_i[i] = (double *) malloc(sizeof( *I_i[i] ) * (N_b+1)); }
  if(I_j) { for(int i=0; i<(N_a+1); i++) I_j[i] = (double *) malloc(sizeof( *I_j[i] ) * (N_b+1)); }

  for(int i=1;i<=N_a;i++){
    for(int j=1;j<=N_b;j++){
      temp[0] = H[i-1][j-1]+similarity_score(seq_a[i-1],seq_b[j-1], mu);
      temp[1] = H[i-1][j]-delta;
      temp[2] = H[i][j-1]-delta;
      temp[3] = 0.;
      H[i][j] = find_array_max(temp,4, ind);
      switch(ind){
      case 0:                                  // score in (i,j) stems from a match/mismatch
        I_i[i][j] = i-1;
        I_j[i][j] = j-1;
        break;
      case 1:                                  // score in (i,j) stems from a deletion in sequence A
        I_i[i][j] = i-1;
        I_j[i][j] = j;
        break;
      case 2:                                  // score in (i,j) stems from a deletion in sequence B
        I_i[i][j] = i;
        I_j[i][j] = j-1;
        break;
      case 3:                                  // (i,j) is the beginning of a subsequence
        I_i[i][j] = i;
        I_j[i][j] = j;
        break;
      }
    }
  }

  double H_max = 0.;
  int i_max=0,j_max=0;
  for(int i=0;i<=N_a;i++){
    for(int j=0;j<=N_b;j++){
      if(H[i][j]>H_max){
        H_max = H[i][j];
        i_max = i;
        j_max = j;
      }
    }
  }

  int current_i=i_max,current_j=j_max;
  int next_i=I_i[current_i][current_j];
  int next_j=I_j[current_i][current_j];
  int tick=0;

  int *tmploc_a = (int *) malloc(sizeof(int) * (N_a+N_b+2));
  int *tmploc_b = (int *) malloc(sizeof(int) * (N_b+N_a+2));

  while(((current_i!=next_i) || (current_j!=next_j)) && (next_j!=0) && (next_i!=0)){

    if(next_i==current_i)  {
        consensus_a[tick] = '-';
        tmploc_a[tick] = -1;
    } else {
        consensus_a[tick] = seq_a[current_i-1];
        tmploc_a[tick] = current_i-1;
    }

    if(next_j==current_j)  {
        consensus_b[tick] = '-';
        tmploc_b[tick] = -1;
    } else {
        consensus_b[tick] = seq_b[current_j-1];
        tmploc_b[tick] = current_j-1;
    }

    current_i = next_i;
    current_j = next_j;
    next_i = I_i[current_i][current_j];
    next_j = I_j[current_i][current_j];
    tick++;
 }

  for(int i=tick-1;i>=0;i--) loc_a[tick-1-i] = tmploc_a[i];
  for(int j=tick-1;j>=0;j--) loc_b[tick-1-j] = tmploc_b[j];

  for(int i=0; i<(N_a+1); i++) free(H[i]);
  free(H);

  free(temp);
  for(int i=0; i<(N_a+1); i++) free(I_i[i]);
  for(int i=0; i<(N_a+1); i++) free(I_j[i]);
  free(I_i);
  free(I_j);

  free(tmploc_a);
  free(tmploc_b);

  return(tick);
}

int align_reads_assembly_start(string seq_ref, string seq_read, double mu, double delta) {

  int N_a = seq_ref.length();                     // get the actual lengths of the sequences
  int N_b = seq_read.length();

  char *consensus_a  = (char *) malloc( sizeof(char) * (N_a+N_b+2) );
  char *consensus_b  = (char *) malloc( sizeof(char) * (N_a+N_b+2) );

  int *loc_a    = (int *) malloc(sizeof(int) * (N_a+N_b+2));
  int *loc_b    = (int *) malloc(sizeof(int) * (N_b+N_a+2));

  int tick    = smith_waterman(seq_ref,  seq_read, mu, delta, loc_a, loc_b, consensus_a, consensus_b);

  int start_loc = loc_a[0] - loc_b[0];

//  int cmp1 = strcmp( seq_ref.substr(loc_a[0],1).c_str(), seq_read.substr(0,1).c_str() );
//  if(cmp1 == 0) start_loc = loc_a[0];

  free(loc_a);
  free(loc_b);
  free(consensus_a);
  free(consensus_b); 

  return(start_loc);
}

int align_query2assembly(string seq_ref, string seq_read, double mu, double delta,char *consen_a, char *consen_b, int *loc_a, int *loc_b) {

  int N_a = seq_ref.length();                      
  int N_b = seq_read.length();

  char *consensus_a  = (char *) malloc( sizeof(char) * (N_a+N_b+2) );
  char *consensus_b  = (char *) malloc( sizeof(char) * (N_a+N_b+2) );

  int tick    = smith_waterman(seq_ref,  seq_read, mu, delta, loc_a, loc_b, consensus_a, consensus_b);

  for(int i=tick-1;i>=0;i--) consen_a[tick-1-i] = consensus_a[i];
  for(int j=tick-1;j>=0;j--) consen_b[tick-1-j] = consensus_b[j];

  return(tick);
}

void align_reads_assembly(string seq_ref, string seq_read, double mu, double delta, ofstream &outFile, int DS, int *num_consensus) {

  int ind;
  int N_a = seq_ref.length();                     
  int N_b = seq_read.length();

  char *consensus_a  = (char *) malloc( sizeof(char) * (N_a+N_b+2) );
  char *consensus_b  = (char *) malloc( sizeof(char) * (N_a+N_b+2) );

  int *loc_a    = (int *) malloc(sizeof(int) * (N_a+N_b+2));
  int *loc_b    = (int *) malloc(sizeof(int) * (N_b+N_a+2));

  int tick    = smith_waterman(seq_ref,  seq_read, mu, delta, loc_a, loc_b, consensus_a, consensus_b);
 
  char *cseq_a = (char *) malloc( sizeof(char) * tick );
  char *cseq_b = (char *) malloc( sizeof(char) * tick );

  for(int i=tick-1;i>=0;i--) cseq_a[tick-1-i] = consensus_a[i];
  for(int j=tick-1;j>=0;j--) cseq_b[tick-1-j] = consensus_b[j];
  string cstring_a(cseq_a);
  string cstring_b(cseq_b); 

  free(consensus_a);
  free(consensus_b); 


  for(int i=0; i<(loc_a[0]-1); i++) outFile << " ";
  if(loc_b[0] > 0) outFile << seq_read.substr(loc_b[0]-1, 1).c_str();
  outFile <<  cstring_b.substr(0,tick).c_str();  
  for(int i=(loc_a[0]+tick+1); i<=seq_ref.length(); i++ ) outFile << " ";
  outFile  << "\t" << DS;


//  if(loc_a[0] < 5000) {

//  for(int i=0; i<(loc_a[0]-1); i++) outFile << " ";
//  outFile << seq_read.substr(0,1).c_str() << cstring_b.substr(0,tick).c_str();
//  for(int i=(loc_a[0]+tick+1); i<=seq_ref.length(); i++ ) outFile << " ";
//  outFile  << "\t" << DS << endl;

//  int loc;
//  if( (_base_to_int[ seq_read.substr(0,1).c_str()[0] ] != 255) && (loc_a[0] != -1)  ) {
//      loc = (loc_a[0]-1)*4 + _base_to_int[ seq_read.substr(0,1).c_str()[0] ];
//      num_consensus[loc]++;
//  }
//  for(int i=0; i<tick; i++) { 
//	if( (_base_to_int[cseq_b[i]] >= 0) && (_base_to_int[cseq_b[i]] <= 3) && (loc_a[i] != -1)) {
//		loc = loc_a[i]*4 + _base_to_int[cseq_b[i]];
//		num_consensus[loc]++;
//        }
//  }

//  } else {
//      outFile << " Excess ==  " << loc_a[0] << " " << loc_a[1] << endl;
//  }

 free(cseq_a);
 free(cseq_b);

 free(loc_a);
 free(loc_b);

}

int get_num_match(string seq_a, string seq_b, double mu, double delta, int range) {

  int ind;
  int N_a = seq_a.length();                      
  int N_b = seq_b.length();

  char *consensus_a = (char *) malloc( sizeof(char) * (N_a+N_b+2) );
  char *consensus_b = (char *) malloc( sizeof(char) * (N_a+N_b+2) );

  int *loc_a    = (int *) malloc(sizeof(int) * (N_a+N_b+2));
  int *loc_b    = (int *) malloc(sizeof(int) * (N_b+N_a+2));

  int tick = smith_waterman(seq_a, seq_b, mu, delta, loc_a, loc_b, consensus_a, consensus_b);

  char *cseq_a = (char *) malloc( sizeof(char) * tick );
  char *cseq_b = (char *) malloc( sizeof(char) * tick );

  for(int i=tick-1;i>=0;i--) cseq_a[tick-1-i] = consensus_a[i];
  for(int j=tick-1;j>=0;j--) cseq_b[tick-1-j] = consensus_b[j];
  string cstring_a(cseq_a);
  string cstring_b(cseq_b);

  std:size_t pos = 0;
  int ngap_b = 0;
  while(true) {
    if(pos>0) pos++;
    std::size_t found = cstring_b.find_first_of('-',pos);
    if (found!=std::string::npos) {
      pos = found;
      ngap_b++;
    }
    else {
        break;
   }
  }

  pos = 0;
  int ngap_a = 0;
  while(true) {
    if(pos>0) pos++;
    std::size_t found = cstring_a.find_first_of('-',pos);
    if (found!=std::string::npos) {
      pos = found;
      ngap_a++;
    }
    else {
        break;
   }
  }

 int match = 0;
 for(int i=0; i<tick; i++) {
        if( cstring_a.compare(i,1,cstring_b.substr(i,1)) == 0 ) match++;
 }

 free(cseq_a);
 free(cseq_b);
 free(consensus_a);
 free(consensus_b);

 free(loc_a);
 free(loc_b);

 int mis_match = tick - ngap_a - ngap_b - match;
 return(mis_match);
}


int check_strand_direction(string seq_a, string seq_b, double mu, double delta, int range, int num_mismatch) {

  int ind;
  int N_a = seq_a.length();                     
  int N_b = seq_b.length();

  char *consensus_a = (char *) malloc( sizeof(char) * (N_a+N_b+2) );
  char *consensus_b = (char *) malloc( sizeof(char) * (N_a+N_b+2) );

  int *loc_a    = (int *) malloc(sizeof(int) * (N_a+N_b+2));
  int *loc_b    = (int *) malloc(sizeof(int) * (N_b+N_a+2));

  int tick = smith_waterman(seq_a, seq_b, mu, delta, loc_a, loc_b, consensus_a, consensus_b);

  char *cseq_a = (char *) malloc( sizeof(char) * tick );
  char *cseq_b = (char *) malloc( sizeof(char) * tick );

  for(int i=tick-1;i>=0;i--) cseq_a[tick-1-i] = consensus_a[i];
  for(int j=tick-1;j>=0;j--) cseq_b[tick-1-j] = consensus_b[j];
  string cstring_a(cseq_a);
  string cstring_b(cseq_b);

  std:size_t pos = 0;
  int ngap_b = 0;
  while(true) {
    if(pos>0) pos++;
    std::size_t found = cstring_b.find_first_of('-',pos);
    if (found!=std::string::npos) {
      pos = found;
      ngap_b++;
    }
    else {
        break;
   }
  }

  pos = 0;
  int ngap_a = 0;
  while(true) {
    if(pos>0) pos++;
    std::size_t found = cstring_a.find_first_of('-',pos);
    if (found!=std::string::npos) {
      pos = found;
      ngap_a++;
    }
    else {
        break;
   }
  }

 int match = 0;
 for(int i=0; i<tick; i++) {
        if( cstring_a.compare(i,1,cstring_b.substr(i,1)) == 0 ) match++;
 }

 free(cseq_a);
 free(cseq_b);
 free(consensus_a);
 free(consensus_b);

 free(loc_a);
 free(loc_b);

 int direction = 0;
 if( ((tick - ngap_a - ngap_b - match) <= num_mismatch) && ((range-tick) <= 1) ) direction = 1;   
// if( (tick - ngap_a - ngap_b - match) <= num_mismatch ) direction = 1;
 return direction;
  
}

double similarity_score(char a,char b,double mu){

  double result;
  if(a==b){
      result=1.;
    }
  else{
      result=-mu;
    }
  return result;
}


double find_array_max(double array[],int length, int &ind){

  double max = array[0];            // start with max = first element
  ind = 0;

  for(int i = 1; i<length; i++){
      if(array[i] > max){
        max = array[i];
        ind = i;
      }
  }
  return max;                    // return highest value in array
}


// -------------------------------------------------------------------

string currAccession;


void compare_strings(int *output, int offset, int range, string seq_ref, string seq_query) {

  float error = 0.85;

  for(int i=0; i<offset; i++) {
    int cn = 0;
    for(int j=0; j<range; j++) {
        char rc = seq_ref[j+i];
        char qc = seq_query[j];
        if( _base_to_int[rc] == _base_to_int[qc] ) cn++;
    }
    if( cn >= int(error*range) ) {
        output[0] = +1;
        output[1] = i;
        break;
    }
  }

 if(output[0] == 0) {
  for(int i=0; i<offset; i++) {
    int cn = 0;
    for(int j=0; j<range; j++) {
        char rc = seq_ref[  seq_ref.length()-1-j-i];
        char qc = seq_query[seq_ref.length()-1-j];
        if( _base_to_int[rc] == _base_to_int[qc] ) cn++;
    }
    if( cn >= int(error*range) ) {
        output[0] = +2;
        output[1] = i;
        break;
    }
  }
 }

  if(output[0] == 0) {
   string rev_seq = revcomp(seq_query);
   for(int i=0; i<offset; i++) {
    int cn = 0;
    for(int j=0; j<range; j++) {
        char rc = seq_ref[j+i];
        char qc = rev_seq[j];
        if( _base_to_int[rc] == _base_to_int[qc] ) cn++;
    }
    if( cn >= int(error*range) ) {
        output[0] = -1;
        output[1] = i;
        break;
    }
   }
 }

  if(output[0] == 0) {
   string rev_seq = revcomp(seq_query);
   for(int i=0; i<offset; i++) {
    int cn = 0;
    for(int j=0; j<range; j++) {
        char rc = seq_ref[seq_ref.length()-1-j-i];
        char qc = rev_seq[rev_seq.length()-1-j];
        if( _base_to_int[rc] == _base_to_int[qc] ) cn++;
    }
    if( cn >= int(error*range) ) {
        output[0] = -2;
        output[1] = i;
        break;
    }
   }
 }

}



//void get_kmer_to_intvalArray (string kmer_forward, int *forward, string kmer_backward, int *backward, int nkmer) {
void get_kmer_to_intvalArray (string kmer_forward, int *forward, int nkmer) {
//   string rev_kmer_backward = revcomp(kmer_backward);
   for(int i=0; i<nkmer; i++) {
	char c = kmer_forward[i];
	*(forward + i) = _base_to_int[c];

//	char rc = rev_kmer_backward[i];
//	*(backward + i) = _base_to_int[rc];
   }
}

void get_kmer2intvalArray (string kmer, int *intarray, int nkmer) {
 for(int i=0; i<nkmer; i++) {
        char c = kmer[i];
        *(intarray + i) = _base_to_int[c];
 }
}
 
bool contains_non_gatc (string kmer) {

  for (unsigned int i = 0; i < kmer.size(); i++) {
	char c = kmer[i];

	if ((_base_to_int[c] > 3) || (_base_to_int[c] < 0))
	  return(false);
/*
	if (! (c == 'g' || c == 'G'
                   || c == 'a' || c == 'A'
                   || c == 't' || c == 'T'
		   || c == 'c' || c == 'C')
		) {
	  return(true);
	}
*/
  }

  return(true);
}

bool check_kmer_quality(string qkmer, int qcutoff){

  int n = 0;
  for (unsigned int i = 0; i < qkmer.size(); i++) {
      char q = qkmer[i]; 
      int qval = get_Qscore(q);
//      if(qval < qcutoff) n++;
//      if(qval < qcutoff) return(false);
  }
//  if(n > 5) return(false);

//  char q = qkmer[qkmer.size()];
//  int qval = get_Qscore(q); 
//  if(qval < qcutoff) return(false);

  return(true);
}


bool contains_non_gatc_num (string seq, int threshold) {
    int num = 0;
    for (unsigned int i = 0; i < seq.size(); i++) {
        char c = seq[i];
        if ((_base_to_int[c] > 3) || (_base_to_int[c] < 0)) num++;
	if(num > threshold) return(false);
    }

    return(true);
}

void get_error_bases(int *error_index, int *error_base, int *cseq_base, string seq_kmer, string cseq_kmer, int kmer_length) {
    for(unsigned int i = 0; i<kmer_length; i++) {
	error_index[i] = 0;
	char base = seq_kmer[i];
	char cbase = cseq_kmer[i];
	
	 if ((_base_to_int[base] > 3) || (_base_to_int[base] < 0)) throw("Rubish base composition in the kmer seq\n");
	 if ((_base_to_int[cbase] > 3) || (_base_to_int[cbase] < 0)) throw("Rubish base composition in the kmer seq\n");

	 if(_base_to_int[base] != _base_to_int[cbase]) {
		error_index[i]++;
		error_base[i] = _base_to_int[base];
		cseq_base[i]  = _base_to_int[cbase];
	 } 	
    }
}

string read_sequence_from_file (string filename) {
  
  ifstream fileReader (filename.c_str());
  if (fileReader) { 
  
    string mySequence;
    char c_line [MAX_LINE_LENGTH];
  
    fileReader.getline(c_line, MAX_LINE_LENGTH);
    while (! fileReader.eof() ) {
	if (c_line[0] != '>') {
	  string s_line (c_line);
	  
	  for (unsigned int i = 0; i < strlen(c_line); i++) {
		char c = tolower(c_line[i]);
		if (c == ' ' || c == '\t' || c == '\n') { continue;}
		mySequence += c;
	  }
	}
	
	fileReader.getline(c_line, MAX_LINE_LENGTH);
    }
  
     return (mySequence);

   } else { // couldn't open file
	throw(stacktrace() + "\n\nCould not open " + filename + "\n");
   }
}

fastaRecord readNextFastaRecord (ifstream& reader) {
  
  fastaRecord fr;

  if (reader.eof()) {
	return (fr);
  }
  
  string line;
  
  // prompt till the accession
  if (currAccession.length() == 0) {
	// must read next accession
	while (getline(reader, line)) {
	  if (line[0] == '>') {
		currAccession = line;
		break;
	  }
	}
  }
  
  fr.accession = currAccession;
  
  while (getline(reader, line)) {
	if (line[0] == '>') {
	  currAccession = line;
	  break;
	}
	else {
	  // append characters that are not whitespace
	  for (unsigned int i = 0; i < line.size(); i++) {
		char c = line[i];
                
		if (c != ' ' && c != '\t' && c != '\n') {
		  fr.sequence += c;
		}
	  }
	}
  }
   
  return(fr);
}
        
string revcomp (const string kmer) {
  
  string revstring;
  
  for (int i = kmer.size() -1; i >= 0; i--) {
	char c = kmer[i];
	char revchar;
	
	switch (c) {
	
	case 'g':
	  revchar = 'c';
	  break;
	
	case 'G':
	  revchar = 'C';
	  break;
	  
	case 'a':
	  revchar = 't';
	  break;

	case 'A':
	  revchar = 'T';
	  break;

	case 't':
	  revchar = 'a';
	  break;
	  
	case 'T':
	  revchar = 'A';
	  break;

	case 'c':
	  revchar = 'g';
	  break;
	  
	case 'C':
	  revchar = 'G';
	  break;

	default:
	  revchar = 'N';
	}
	
        
	revstring += revchar;
        
  }
  
  
  return (revstring.c_str());
}

string revquality (const string quality) {
  string rev_quality;
  for (int i = quality.size() -1; i >= 0; i--) {
        char c = quality[i];
        rev_quality += c;
  }

//  rev_quality.reserve(quality.size());
//  int i, j;
//  for(i=0, j=quality.size()-1; i<quality.size(), j>=0; i++, j--) rev_quality[i] = quality[j];

  return (rev_quality.c_str());
}

kmer_int_type_t revcomp_val(kmer_int_type_t kmer, unsigned int kmer_length)
{
	kmer_int_type_t rev_kmer = 0;
	kmer = ~kmer;
//	mpz_com(kmer.get_mpz_t(), kmer.get_mpz_t());
	for (unsigned int i = 0; i < kmer_length; i++) {

		int base = kmer & 3;
//        	kmer_int_type_t cbase;
//		kmer_int_type_t bop = 3;
//        	mpz_and(cbase.get_mpz_t(), kmer.get_mpz_t(), bop.get_mpz_t());
//        	int base = (int)mpz_get_ui(cbase.get_mpz_t());

		rev_kmer = rev_kmer << 2;
//		mpz_mul_2exp(rev_kmer.get_mpz_t(),rev_kmer.get_mpz_t(),2);

		rev_kmer+= base;
		kmer = kmer >> 2;
//		mpz_fdiv_q_2exp(kmer.get_mpz_t(),kmer.get_mpz_t(),2);

	}

	return rev_kmer;
}

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

int base_to_int_value (char nucleotide) {

  switch (nucleotide) {
	
  case 'G':
  case 'g':
	return(0);
  
  case 'A':
  case 'a':
	return(1);

  case 'T':
  case 't':
	return(2);
	
  case 'C':
  case 'c':
	return(3);
	
  default:
	return(-1);
  
  }

}

char int_to_base(int baseval) {
  
  if (baseval < 0 || baseval > 3) {
	throw (stacktrace() + "\n\nError, baseval out of range 0-3");
  }

  return(_int_to_base[baseval]);
}


bool untrusted_intersect(vector<int> untrusted_subset, vector<short> & region, int kmer_length, int read_length) {
  int start = 0;
  int end = read_length-1;

  int u;
  for(int i = 0; i < untrusted_subset.size(); i++) {
    u = untrusted_subset[i];

    if(start <= u+kmer_length-1 && u <= end) {
        start = max(start, u);
        end = min(end, u+kmer_length-1);
    } else {
        return false;
    }
 }
   
 for(short i = start; i <= end; i++)
    region.push_back(i);
 return true;
}

void untrusted_union(vector<int> untrusted_subset, vector<short> & region, int kmer_length) {
  short u;
  set<short> region_set;
  for(int i = 0; i < untrusted_subset.size(); i++) {
    u = untrusted_subset[i];

    for(short ui = u; ui < u+kmer_length; ui++)
      region_set.insert(ui);
  }

  set<short>::iterator it;
  for(it = region_set.begin(); it != region_set.end(); it++)
    region.push_back(*it);
}

vector<short> error_region_extract(vector<int> untrusted_subset, int kmer_length, int read_length) {
   vector<short> region;

if(untrusted_intersect(untrusted_subset, region, kmer_length, read_length)) {

  int right_leftkmer = untrusted_subset.front()-1;
  if(right_leftkmer >= 0) {

    vector<short> front_chop(region);
    region.clear();
    for(int i = 0; i < front_chop.size(); i++) {
      if(front_chop[i] > right_leftkmer+kmer_length-1)
        region.push_back(front_chop[i]);
    }

  } else {
    for(int i = region[0]-1; i >= 0; i--)
      region.push_back(i);
  }

  int left_rightkmer = untrusted_subset.back()+1;
  if(left_rightkmer+kmer_length-1 < read_length) {
    vector<short> back_chop(region);
    region.clear();
    for(int i = 0; i < back_chop.size(); i++) {
      if(back_chop[i] < left_rightkmer)
        region.push_back(back_chop[i]);
    }

  } else {
    for(int i = region.back()+1; i < read_length; i++)
      region.push_back(i);
  }

  return region;
 }
}



vector<short> error_region_chop(vector<int> untrusted_subset, int kmer_length, int read_length, float* prob) {
   vector<short> region;
//  if(!untrusted_intersect(untrusted_subset, region, kmer_length, read_length))
//    untrusted_union(untrusted_subset, region, kmer_length);

if(untrusted_intersect(untrusted_subset, region, kmer_length, read_length)) {

  int right_leftkmer = untrusted_subset.front()-1;
  if(right_leftkmer >= 0) {

    vector<short> front_chop(region);
    region.clear();
    for(int i = 0; i < front_chop.size(); i++) {
      if(front_chop[i] > right_leftkmer+kmer_length-1)
        region.push_back(front_chop[i]);
    }

//    for(int er = 0; er < expand_region; er++) {
//      int pre_region = region[0] - (er+1);
//      if(pre_region >= 0 && (prob[pre_region] < .99 || prob[pre_region] <= prob[pre_region+1])) {
//        vector<short>::iterator it;
//        it = region.begin();
//        region.insert(it, pre_region);
//      }
//    }

  } else {
    for(int i = region[0]-1; i >= 0; i--)
      region.push_back(i);
  }


  int left_rightkmer = untrusted_subset.back()+1;
  if(left_rightkmer+kmer_length-1 < read_length) {
    vector<short> back_chop(region);
    region.clear();
    for(int i = 0; i < back_chop.size(); i++) {
      if(back_chop[i] < left_rightkmer)
        region.push_back(back_chop[i]);
    }

  } else {
    for(int i = region.back()+1; i < read_length; i++)
      region.push_back(i);
  }

  return region;
 }
}


kmer_int_type_t get_maximum_kmer_intval(unsigned int kmer_length) {

   char c = 'C';
   kmer_int_type_t max_kmer = 0;
   for(unsigned int i = 0; i< kmer_length; i++) {
        int val = _base_to_int[c];
        max_kmer = max_kmer << 2;
        max_kmer |= val;

//	  kmer_int_type_t tval;
//	  mpz_set_ui(tval.get_mpz_t(), (unsigned long int)val );
//	  mpz_mul_2exp(max_kmer.get_mpz_t(),max_kmer.get_mpz_t(),2); 
//	  mpz_ior(max_kmer.get_mpz_t(), max_kmer.get_mpz_t(), tval.get_mpz_t());
	  	
   }
   return(max_kmer);
}


kmer_int_type_t kmer_to_intval(string kmer) {

 
  if (kmer.length() > 100) {
	return(-1);
//	throw(stacktrace() + "\n\nerror, kmer length exceeds 32");
  }

  kmer_int_type_t kmer_val = 0;
  
  for (unsigned int i = 0; i < kmer.length(); i++) {
	char c = kmer[i];
	int val = _base_to_int[c];
        /* ottmi: Don't need this here: the only possible source for non-gatc characters is
           the reads file and we already check in add_kmer() */
        /* bhaas-to-ottmi: I need to put this back in....  There are other sources where we
           don't explicitly check, and not doing so here causes problems downstream */
        /* deccles-to-bhass: In the interest of code speed / efficiency, the non-gatc check
           is better inlined here -- it's already working out the value per-char, so we get
           the check for a single if statement per character */
        if(val > 3){
 //         stringstream errstr;
 //         errstr << "\n\nerror, kmer contains nongatc: " << kmer;
 //         cerr << errstr.str();
 //         throw(stacktrace() + errstr.str());
            return(-1);
        }
   
//	 cerr << "Char " << c << "=" << val << endl;

	kmer_val = kmer_val << 2;
	kmer_val |= val;
//	mpz_mul_2exp(kmer_val.get_mpz_t(),kmer_val.get_mpz_t(),2);
//	kmer_int_type_t tval; 
//	tval = _base_to_int[c];
//	mpz_ior(kmer_val.get_mpz_t(), kmer_val.get_mpz_t(), tval.get_mpz_t());	

//	 cerr << "\tkmerval: " << kmer_val << endl;
	
  }
  
//  cerr << "Kmer: " << kmer << " = " << kmer_val << endl;
  
  return(kmer_val);
}


//uint64_t convert_mpz_uint64(kmer_int_type_t intval, unsigned int kmer_length) {

//   uint64_t kmer_uint64 = 0;
//   string kmer_seq = decode_kmer_from_intval(intval, kmer_length);
//   for (unsigned int i = 0; i < kmer_length; i++) {
//	char c = kmer_seq[i];
//        int val = _base_to_int[c];	
//	if(val > 3){
//            return(-1);
//	}

//	kmer_uint64 = kmer_uint64 << 2;
//	kmer_uint64 |= val;	
//   } 

//  return( kmer_uint64 );
//}

//string decode_kmer_from_uint64(uint64_t intval, unsigned int kmer_length) {

//  string kmer(kmer_length, ' ');
//  for (unsigned int i = 1; i <= kmer_length; i++) {
//	int base_num = intval & 3ll;
//	kmer[kmer_length-i] = _int_to_base[base_num];
//	intval = intval >> 2;	
//  }
//  return(kmer);
//
//}

string decode_kmer_from_intval(kmer_int_type_t intval, unsigned int kmer_length) {

  string kmer(kmer_length, ' ');

  kmer_int_type_t tbase = 3;
  for (unsigned int i = 1; i <= kmer_length; i++) {
	
	int base_num = intval & 3ll;
//	kmer_int_type_t tbase;	
//	kmer_int_type_t bop = 3;
//	mpz_and(tbase.get_mpz_t(), intval.get_mpz_t(), bop.get_mpz_t());
//	int base_num = (int)mpz_get_ui(tbase.get_mpz_t());	
	
	kmer[kmer_length-i] = _int_to_base[base_num];
	
	// cerr << "base: " << base << endl;
	
	intval = intval >> 2;
//	mpz_fdiv_q_2exp(intval.get_mpz_t(), intval.get_mpz_t(),2);	
  }

  return(kmer);
}

vector<kmer_int_type_t> get_Fkmer_candidates(kmer_int_type_t seed_kmer, unsigned int kmer_length) {
	kmer_int_type_t forward_prefix = (seed_kmer << (33-kmer_length)*2) >> (32-kmer_length)*2;
//	kmer_int_type_t forward_prefix;
//	mpz_mul_2exp( forward_prefix.get_mpz_t(), seed_kmer.get_mpz_t(), (101-kmer_length)*2);
//	mpz_fdiv_q_2exp( forward_prefix.get_mpz_t(), forward_prefix.get_mpz_t(), (100-kmer_length)*2);		

	vector<kmer_int_type_t> candidates;
	for (kmer_int_type_t i=0; i<4; i++) {
		candidates.push_back( forward_prefix | i );
	//	kmer_int_type_t tmp_candidate;
	//	mpz_ior(tmp_candidate.get_mpz_t(), forward_prefix.get_mpz_t(), i.get_mpz_t());
	//	candidates.push_back( tmp_candidate );
	}
	return(candidates);
}

vector<kmer_int_type_t> get_Rkmer_candidates(kmer_int_type_t seed_kmer, unsigned int kmer_length) {
	kmer_int_type_t reverse_suffix = seed_kmer >> 2;
//	kmer_int_type_t reverse_suffix;
//	mpz_fdiv_q_2exp(reverse_suffix.get_mpz_t(), seed_kmer.get_mpz_t(), 2);

	vector<kmer_int_type_t> candidates;
        for (kmer_int_type_t i=0; i<4; i++) {
		candidates.push_back( (i << (kmer_length*2-2)) | reverse_suffix );
	//	kmer_int_type_t tmp_candidate;
	//	mpz_mul_2exp(tmp_candidate.get_mpz_t(), i.get_mpz_t(), (kmer_length*2-2));
	//	mpz_ior(tmp_candidate.get_mpz_t(), tmp_candidate.get_mpz_t(), reverse_suffix.get_mpz_t());
	//	candidates.push_back( tmp_candidate );
	}
	return(candidates);	
}

float compute_entropy(kmer_int_type_t kmer, unsigned int kmer_length) {

  char counts[] = { 0, 0, 0, 0 };

  for (unsigned int i = 0; i < kmer_length; i++) {

	int c = kmer & 3;
	kmer = kmer >> 2;

//	kmer_int_type_t cmpz;
//	kmer_int_type_t bop = 3;
// 	mpz_and(cmpz.get_mpz_t(), kmer.get_mpz_t(), bop.get_mpz_t());
//	int c = (int)mpz_get_ui(cmpz.get_mpz_t());

//	mpz_fdiv_q_2exp(kmer.get_mpz_t(),kmer.get_mpz_t(),2);

	counts[c]++;
  }

  float entropy = 0;

  for (unsigned int i = 0; i < 4; i++) {

	float prob = (float)counts[i] / kmer_length;

	if (prob > 0) {
	  float val = prob * log(1/prob)/log(2.0f);
	  entropy += val;
	}
  }

  return(entropy);
}
  

float compute_entropy(string& kmer) {

  map<char,int> char_map;

  for (unsigned int i = 0; i < kmer.length(); i++) {
	
	char c = kmer[i];
	char_map[c]++;
  }

  float entropy = 0;
  
  char nucs[] = { 'G', 'A', 'T', 'C' };

  for (unsigned int i = 0; i < 4; i++) {
	
	char nuc = nucs[i];
	
	int count = char_map[nuc];
	
	float prob = (float)count / kmer.length();
	
	if (prob > 0) {
	  float val = prob * log(1/prob)/log(2.0f);
	  entropy += val;
	}
  }

  return(entropy);
}



kmer_int_type_t get_DS_kmer_val(kmer_int_type_t kmer_val, unsigned int kmer_length) {
    
    kmer_int_type_t rev_kmer = revcomp_val(kmer_val, kmer_length);
    
    if (rev_kmer > kmer_val)
        kmer_val = rev_kmer;
    
    return(kmer_val);
    
}

vector<kmer_int_type_t> sequence_string_to_kmer_int_type_vector(const string& sequence, int kmer_length) {

    vector<kmer_int_type_t> kit_vec;
    
    for (unsigned int i = 0; i <= sequence.length() - kmer_length; i++) {
        
        string kmer = sequence.substr(i, kmer_length);
                
        kmer_int_type_t k = kmer_to_intval(kmer);
        
        kit_vec.push_back(k);
    }

    return(kit_vec);
    
}

string replace_nonGATC_chars_with_A (string& str) {
    
    stringstream newStr;

    for (unsigned int i = 0; i < str.size(); i++) {
        char c = str[i];
        
        if (_base_to_int[c] > 3)
            c = 'A';
        
        newStr << c;
    }

    return(newStr.str());

}

//bool Sort_kmer_by_count_desc (const Kmer_Occurence_Pair& i, const Kmer_Occurence_Pair& j) {
//        return (i.second > j.second);
//}

