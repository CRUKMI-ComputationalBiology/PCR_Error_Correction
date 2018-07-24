/*
  Copyright (C) 2017-2018 Cancer Research UK Manchester Institute

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


#include "IRKE.hpp"
#include "KmerCounter.hpp"
#include "stacktrace.hpp"

#include "typedefs.h"
#include "mrSubroutine.h"


using namespace MAPREDUCE_NS;
using namespace std;


int _base_collision(Depth depth) {
    int Base = 0;
    if(depth.d0 == 1) Base = Base | 8;
    if(depth.d1 == 1) Base = Base | 4;
    if(depth.d2 == 1) Base = Base | 2;
    if(depth.d3 == 1) Base = Base | 1;
    return Base;
}


void fileread_RNAseq_HeadSeq_FASTA(int itask, char *fname, KeyValue *kv, void *ptr)
{
  Data *data = (Data *) ptr;
  RankID rankid;
  struct stat stbuf;
  int flag = stat(fname,&stbuf);
  if (flag < 0) {
    printf("ERROR: Could not query file size\n");
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  string filename(fname);
  std::size_t fpos = filename.find("sKmer_");
  if( fpos != std::string::npos ) {

    Fasta_reader fasta_reader(filename);
    while (true) {
        Fasta_entry fe = fasta_reader.getNext();
        string seq    = fe.get_sequence();
        string header = fe.get_header();
        if (seq == "") break;

	data->range = (int)seq.length();

        std::size_t pos = header.find_first_of(" ");
        if(pos != std::string::npos) header.replace(pos,1,"_");
	std::size_t pe_pos = header.find_first_of("_");

	HeadSeq header_seq;
	memcpy(header_seq.header, header.substr(1,pe_pos).c_str(), sizeof(char)*pe_pos);        

	if( header.compare(pe_pos+1,1,"1") == 0 ) header_seq.PE = 1;
	else	 				  header_seq.PE = 2;

	memcpy(header_seq.seq, seq.c_str(), sizeof(char)*data->range);	
	header_seq.length = data->range;

        kv->add((char *) header_seq.header, sizeof(char)*pe_pos, (char *) &header_seq, sizeof(HeadSeq) );

    }

  }
}

void fileread_Sequence_Fastq(int itask, KeyValue *kv, void *ptr)
{
  Data *data = (Data *) ptr;
  int PairEnd = data->read_side;

  string filename(data->filename);
  filename.append("/sKmer_");
  filename.append(to_string((long long int)data->lane));
  filename.append("_");
  filename.append(to_string((long long int)data->read_side));
  filename.append("_");
  filename.append(to_string((long long int)data->me));
  filename.append(".fq");

  Fastq_reader fastq_reader(filename); 
  while (true) {
        Fastq_entry fq = fastq_reader.getNext();

        string seq     = fq.get_sequence();
        string quality = fq.get_quality();
        string header  = fq.get_header();

        if (seq == "") break;
        data->range = (int)seq.length();

        if(contains_non_gatc_num(seq,5)) {

	   HeadSeq header_seq;

           std::size_t pos_1 = header.find_first_of(" ");
           std::size_t pos_2 = header.find_first_of("/");
 
           if(pos_1 != std::string::npos)      memcpy(header_seq.header, header.substr(0,(int)pos_1).c_str(), sizeof(char)*((int)pos_1));
	   else if(pos_2 != std::string::npos) memcpy(header_seq.header, header.substr(0,(int)pos_2).c_str(), sizeof(char)*((int)pos_2));
	   else {
    		printf("ERROR: Could not query headers in input fastq file, check the manual of the software!!! \n");
    		MPI_Abort(MPI_COMM_WORLD,1);
	   }

          if(     PairEnd == 1) header_seq.PE = 1;
          else if(PairEnd == 2) header_seq.PE = 2;

          if(PairEnd > 0) {
            memcpy(header_seq.seq,     seq.c_str(),     sizeof(char)*data->range);
            memcpy(header_seq.quality, quality.c_str(), sizeof(char)*data->range);
            header_seq.length = data->range;

            if(pos_1 != std::string::npos)      kv->add((char *) header_seq.header, sizeof(char)*((int)pos_1), (char *) &header_seq, sizeof(HeadSeq) );
	    else if(pos_2 != std::string::npos) kv->add((char *) header_seq.header, sizeof(char)*((int)pos_2), (char *) &header_seq, sizeof(HeadSeq) );
          }

        }

    }

}

void fileread_RNAseq_HeadSeq_FASTQ(int itask, char *fname, KeyValue *kv, void *ptr)
{
  Data *data = (Data *) ptr;
  RankID rankid;
  struct stat stbuf;
  int flag = stat(fname,&stbuf);
  if (flag < 0) {
    printf("ERROR: Could not query file size\n");
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  string filename(fname);
  std::size_t fpos = filename.find("sKmer_");
  if( fpos != std::string::npos ) {

     std::size_t fpos_1_1 = filename.find("sKmer_1_1");
     std::size_t fpos_2_1 = filename.find("sKmer_2_1");
     std::size_t fpos_3_1 = filename.find("sKmer_3_1");
     std::size_t fpos_4_1 = filename.find("sKmer_4_1");

     std::size_t fpos_1_2 = filename.find("sKmer_1_2");
     std::size_t fpos_2_2 = filename.find("sKmer_2_2");
     std::size_t fpos_3_2 = filename.find("sKmer_3_2");
     std::size_t fpos_4_2 = filename.find("sKmer_4_2");

     int PairEnd = 0;
     if( (fpos_1_1 != std::string::npos) || (fpos_2_1 != std::string::npos) || (fpos_3_1 != std::string::npos) || (fpos_4_1 != std::string::npos) ) PairEnd = 1;
     if( (fpos_1_2 != std::string::npos) || (fpos_2_2 != std::string::npos) || (fpos_3_2 != std::string::npos) || (fpos_4_2 != std::string::npos) ) PairEnd = 2;

    Fastq_reader fastq_reader(filename);
    while (true) {
	Fastq_entry fq = fastq_reader.getNext();
        string seq     = fq.get_sequence();
 	string quality = fq.get_quality();
        string header  = fq.get_header();

	if (seq == "") break;
	data->range = (int)seq.length();

	if(contains_non_gatc_num(seq,5)) {

	   std::size_t pos = header.find_first_of(" ");

	   HeadSeq header_seq;
           memcpy(header_seq.header, header.substr(0,(int)pos).c_str(), sizeof(char)*((int)pos));

          if(     PairEnd == 1) header_seq.PE = 1;
          else if(PairEnd == 2) header_seq.PE = 2;

	  if(PairEnd > 0) {
            memcpy(header_seq.seq,     seq.c_str(),     sizeof(char)*data->range);
            memcpy(header_seq.quality, quality.c_str(), sizeof(char)*data->range);
            header_seq.length = data->range;

	    kv->add((char *) header_seq.header, sizeof(char)*((int)pos), (char *) &header_seq, sizeof(HeadSeq) );
 	  }
 
	}

    }

  }   

}

void map_assign_reference(uint64_t itask, char *key, int keybytes, char *value, int valuebytes, KeyValue *kv, void *ptr)
{
  Data *data = (Data *) ptr;

  int chr = *(int *) key;
  string ref_seq;
  ref_seq.assign(value, valuebytes);

  data->reference_sequence.at(chr) = ref_seq;
}

void map_upload_reference(int itask, KeyValue *kv, void *ptr)
{
  Data *data = (Data *) ptr;

  Fasta_reader fasta_reader(data->filename);
  Fasta_entry fe;
  int loc = 0;

  while (true) {

        fe = fasta_reader.getNext();
        string seq = fe.get_sequence();
        string header = fe.get_header();
        if (seq == "") break;
        if( (loc % data->nprocs) == data->me ) {

          std::size_t pos = header.find_first_of(" ");
          int chr = 0;

          if( strcmp(header.substr(0,pos).c_str(), "1") == 0) chr = 1;
          else if( strcmp(header.substr(0,pos).c_str(), "2") == 0) chr = 2;
          else if( strcmp(header.substr(0,pos).c_str(), "3") == 0) chr = 3;
          else if( strcmp(header.substr(0,pos).c_str(), "4") == 0) chr = 4;
          else if( strcmp(header.substr(0,pos).c_str(), "5") == 0) chr = 5;
          else if( strcmp(header.substr(0,pos).c_str(), "6") == 0) chr = 6;
          else if( strcmp(header.substr(0,pos).c_str(), "7") == 0) chr = 7;
          else if( strcmp(header.substr(0,pos).c_str(), "8") == 0) chr = 8;
          else if( strcmp(header.substr(0,pos).c_str(), "9") == 0) chr = 9;
          else if( strcmp(header.substr(0,pos).c_str(),"10") == 0) chr = 10;
          else if( strcmp(header.substr(0,pos).c_str(),"11") == 0) chr = 11;
          else if( strcmp(header.substr(0,pos).c_str(),"12") == 0) chr = 12;
          else if( strcmp(header.substr(0,pos).c_str(),"13") == 0) chr = 13;
          else if( strcmp(header.substr(0,pos).c_str(),"14") == 0) chr = 14;
          else if( strcmp(header.substr(0,pos).c_str(),"15") == 0) chr = 15;
          else if( strcmp(header.substr(0,pos).c_str(),"16") == 0) chr = 16;
          else if( strcmp(header.substr(0,pos).c_str(),"17") == 0) chr = 17;
          else if( strcmp(header.substr(0,pos).c_str(),"18") == 0) chr = 18;
          else if( strcmp(header.substr(0,pos).c_str(),"19") == 0) chr = 19;
          else if( strcmp(header.substr(0,pos).c_str(),"20") == 0) chr = 20;
          else if( strcmp(header.substr(0,pos).c_str(),"21") == 0) chr = 21;
          else if( strcmp(header.substr(0,pos).c_str(),"22") == 0) chr = 22;
          else if( strcmp(header.substr(0,pos).c_str(), "X") == 0) chr = 23;
          else if( strcmp(header.substr(0,pos).c_str(), "Y") == 0) chr = 24;
          else if( strcmp(header.substr(0,pos).c_str(),"MT") == 0) chr = 25;

          if( (chr>=1) && (chr <= 25) ) {
                cerr << "Node " << data->me << " loaded chromosome " << chr << endl;
                kv->add((char *) &chr, sizeof(int), (char *) seq.substr(0,seq.length()).c_str(), sizeof(char)*seq.length() );
          }

        }

        loc++;

  }

}

void fileread_AlignmentInfo_MPI(int itask, KeyValue *kv, void *ptr)
{
  Data *data = (Data *) ptr;

  string filename(data->filename);
  filename.append("/bfiltered_");
  filename.append(to_string((long long int)data->me));
  filename.append(".out");

  struct stat stbuf;
  int flag = stat(filename.c_str(),&stbuf);
  if (flag < 0) {
    printf("ERROR: Could not query file size\n");
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  FILE* FS;
  FS = fopen(filename.c_str(), "r");

    while (!feof(FS)) {
        char temp_string[1000];
        char head[100];
        char cigar_string[500];
        char dummy[20];
        int flag, maq, lane;

        sReadInfo read;

        fgets(temp_string, 1000, FS);

        int n = sscanf(temp_string, "%s %s %d %d %s %d %d", head, dummy, &read.loc, &read.maq, cigar_string, &flag, &lane);
        if(n < 5) break;

        if( strcmp(cigar_string,"*") != 0 ) {

        read.chr = 0;
        if(      strcmp(dummy, "1") == 0) read.chr = 1;
        else if( strcmp(dummy, "2") == 0) read.chr = 2;
        else if( strcmp(dummy, "3") == 0) read.chr = 3;
        else if( strcmp(dummy, "4") == 0) read.chr = 4;
        else if( strcmp(dummy, "5") == 0) read.chr = 5;
        else if( strcmp(dummy, "6") == 0) read.chr = 6;
        else if( strcmp(dummy, "7") == 0) read.chr = 7;
        else if( strcmp(dummy, "8") == 0) read.chr = 8;
        else if( strcmp(dummy, "9") == 0) read.chr = 9;
        else if( strcmp(dummy,"10") == 0) read.chr = 10;
        else if( strcmp(dummy,"11") == 0) read.chr = 11;
        else if( strcmp(dummy,"12") == 0) read.chr = 12;
        else if( strcmp(dummy,"13") == 0) read.chr = 13;
        else if( strcmp(dummy,"14") == 0) read.chr = 14;
        else if( strcmp(dummy,"15") == 0) read.chr = 15;
        else if( strcmp(dummy,"16") == 0) read.chr = 16;
        else if( strcmp(dummy,"17") == 0) read.chr = 17;
        else if( strcmp(dummy,"18") == 0) read.chr = 18;
        else if( strcmp(dummy,"19") == 0) read.chr = 19;
        else if( strcmp(dummy,"20") == 0) read.chr = 20;
        else if( strcmp(dummy,"21") == 0) read.chr = 21;
        else if( strcmp(dummy,"22") == 0) read.chr = 22;
        else if( strcmp(dummy, "X") == 0) read.chr = 23;
        else if( strcmp(dummy, "Y") == 0) read.chr = 24;
        else if( strcmp(dummy,"MT") == 0) read.chr = 25;

      if( (read.chr > 0) && (read.chr <= 25) ) {


        string name(head);
        name = remove_whitespace(name);

        unsigned check_PE = (unsigned) flag & 64;
        if(check_PE == 64) read.pe = 1;

        check_PE = (unsigned) flag & 128;
        if(check_PE == 128) read.pe = 2;

        unsigned check_strand =  (unsigned) flag & 16;

        if(check_strand == 16) read.strand = -1;
        else                   read.strand = 1;

        vector<CigarOp_Function> cigar;
        CigarOp_Function::parse(cigar_string,cigar);

        if(cigar.size() > 0) {

          read.SOFT_S = 0;
          read.SOFT_E = 0;
          if(cigar[0].op == 'S')              read.SOFT_S += cigar[0].size;
          if(cigar[cigar.size()-1].op == 'S') read.SOFT_E += cigar[cigar.size()-1].size;

	  for(int i=1; i<cigar.size(); i++) {
		if(cigar[i].op == 'D') read.SOFT_E += cigar[i].size; 
		if(cigar[i].op == 'I') read.SOFT_E -= cigar[i].size;
	  } 

//          int nclip_hard = 0;
//          if(cigar[0].op == 'H')              nclip_hard += cigar[0].size;
//          if(cigar[cigar.size()-1].op == 'H') nclip_hard += cigar[cigar.size()-1].size;

          unsigned check_780 = flag & 780;
          unsigned check_2   = flag & 2;

          read.align_len = parse_cigar_align_length(cigar_string);

          bool isOnTarget = check_onTarget(read.chr, (read.loc - read.SOFT_S), (read.loc + read.align_len - 1 + read.SOFT_E), data->loc_vector, 1000, data->_target_loc );
          int sallow = data->range - (int)( (float)data->range * data->soft_allow);

          if( isOnTarget && (check_780 == 0) && (check_2 == 2) && ((read.pe == 1) || (read.pe == 2)) && (read.align_len >= read.align_len) )  {
            kv->add((char *) name.c_str(), sizeof(char) * name.length(), (char *) &read, sizeof(sReadInfo));
          }

        }

       }

    }

  }

  fclose(FS);

}


void fileread_AlignmentInfo(int itask, char *fname, KeyValue *kv, void *ptr)
{
  Data *data = (Data *) ptr;
  struct stat stbuf;
  int flag = stat(fname,&stbuf);
  if (flag < 0) {
    printf("ERROR: Could not query file size\n");
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  std::string filename(fname);
  std::size_t fpos = filename.find("bfiltered_");

  if( fpos != std::string::npos ) {
    FILE* FS;
    FS = fopen(filename.c_str(), "r");

    while (!feof(FS)) {
        char temp_string[1000];
        char head[100];
        char cigar_string[500];
        char dummy[20];
        int flag, maq, lane;

        sReadInfo read;

        fgets(temp_string, 1000, FS);

	int n = sscanf(temp_string, "%s %s %d %d %s %d %d", head, dummy, &read.loc, &read.maq, cigar_string, &flag, &lane);
        if(n < 5) break;

        if( strcmp(cigar_string,"*") != 0 ) {

	read.chr = 0;
        if(      strcmp(dummy, "1") == 0) read.chr = 1;
        else if( strcmp(dummy, "2") == 0) read.chr = 2;
        else if( strcmp(dummy, "3") == 0) read.chr = 3;
        else if( strcmp(dummy, "4") == 0) read.chr = 4;
        else if( strcmp(dummy, "5") == 0) read.chr = 5;
        else if( strcmp(dummy, "6") == 0) read.chr = 6;
        else if( strcmp(dummy, "7") == 0) read.chr = 7;
        else if( strcmp(dummy, "8") == 0) read.chr = 8;
        else if( strcmp(dummy, "9") == 0) read.chr = 9;
        else if( strcmp(dummy,"10") == 0) read.chr = 10;
        else if( strcmp(dummy,"11") == 0) read.chr = 11;
        else if( strcmp(dummy,"12") == 0) read.chr = 12;
        else if( strcmp(dummy,"13") == 0) read.chr = 13;
        else if( strcmp(dummy,"14") == 0) read.chr = 14;
        else if( strcmp(dummy,"15") == 0) read.chr = 15;
        else if( strcmp(dummy,"16") == 0) read.chr = 16;
        else if( strcmp(dummy,"17") == 0) read.chr = 17;
        else if( strcmp(dummy,"18") == 0) read.chr = 18;
        else if( strcmp(dummy,"19") == 0) read.chr = 19;
        else if( strcmp(dummy,"20") == 0) read.chr = 20;
        else if( strcmp(dummy,"21") == 0) read.chr = 21;
        else if( strcmp(dummy,"22") == 0) read.chr = 22;
        else if( strcmp(dummy, "X") == 0) read.chr = 23;
        else if( strcmp(dummy, "Y") == 0) read.chr = 24;
        else if( strcmp(dummy,"MT") == 0) read.chr = 25;

      if( (read.chr > 0) && (read.chr <= 25) ) {


        string name(head);
        name = remove_whitespace(name);

        unsigned check_PE = (unsigned) flag & 64;
	if(check_PE == 64) read.pe = 1;
	
	check_PE = (unsigned) flag & 128;
        if(check_PE == 128) read.pe = 2;

	unsigned check_strand =  (unsigned) flag & 16;	
	
        if(check_strand == 16) read.strand = -1;
        else 		       read.strand = 1;

        vector<CigarOp_Function> cigar;
        CigarOp_Function::parse(cigar_string,cigar);	

        if(cigar.size() > 0) {

          read.SOFT_S = 0;
	  read.SOFT_E = 0;
	  if(cigar[0].op == 'S')              read.SOFT_S += cigar[0].size;
          if(cigar[cigar.size()-1].op == 'S') read.SOFT_E += cigar[cigar.size()-1].size;

          int nclip_hard = 0;
          if(cigar[0].op == 'H')              nclip_hard += cigar[0].size;
          if(cigar[cigar.size()-1].op == 'H') nclip_hard += cigar[cigar.size()-1].size;

	  unsigned check_780 = flag & 780;
          unsigned check_2   = flag & 2;

	  read.align_len = parse_cigar_align_length(cigar_string);

          bool isOnTarget = check_onTarget(read.chr, (read.loc - read.SOFT_S), (read.loc + read.align_len - 1 + read.SOFT_E), data->loc_vector, 5000, data->_target_loc );

	  if( isOnTarget && (check_780 == 0) && ((read.pe == 1) || (read.pe == 2)) )  {	
	    kv->add((char *) name.c_str(), sizeof(char) * name.length(), (char *) &read, sizeof(sReadInfo)); 
	  }

	}

       }

    }

    }



  }

}

bool check_onTarget(int chr, int loc1, int loc2, int loc_vector, int flank, const vector<ExonLoc>& _target_loc)
{
   for(int i=0; i<loc_vector; i++) {
        if( (chr == _target_loc[i].chr) &&
            ((uint64_t)loc1 >= (_target_loc[i].loc1 - (uint64_t)flank)) && 
	    ((uint64_t)loc1 <= (_target_loc[i].loc2 + (uint64_t)flank)) && 
	    ((uint64_t)loc2 >= (_target_loc[i].loc1 - (uint64_t)flank)) && 
	    ((uint64_t)loc2 <= (_target_loc[i].loc2 + (uint64_t)flank)) 
          ) {
	      return true;
         }

   }
   return false;
}

int parse_cigar_softclipping(char *s)
{
   int num_soft = 0;
   char* p=(char*)(s);
   while(*p!=0)
   {
	 char* endptr;

         if(!isdigit(*p)) return -1;
         int c_size =strtol(p,&endptr,10);
         if(c_size<=0)    return -1;
         p=endptr;

	 int check = 0;
         if(!isalpha(*p)) return -1;
 
         if((p[0] == 'S') || (p[0] == 'H')) check++;
	 if(check>0)    num_soft += c_size; 
         ++p;
   }
   return num_soft;
}

int parse_cigar_clipping(char *s, bool start)
{
   int num_soft = 0;
   char* p=(char*)(s);
   bool check = false;

   char* endptr;
   int c_size =strtol(p,&endptr,10);
   p = endptr;
   if((p[0] == 'S') || (p[0] == 'H')) check = true;

   if(start) {
        if(check) return c_size;
	else	  return num_soft;
   } else {
     while(*p!=0)
     {
         if(!isdigit(*p)) return -1;
         c_size =strtol(p,&endptr,10);
         if(c_size<=0)    return -1;
         p=endptr;

         int check = 0;
         if(!isalpha(*p)) return -1;

         if((p[0] == 'S') || (p[0] == 'H')) check++;
         if(check>0)    num_soft = c_size;
         ++p;
     }
   }
   return num_soft;   
}

int parse_cigar_align_length(char *s)
{

   vector<CigarOp_Function> cigar;
   CigarOp_Function::parse(s,cigar);

   int length = 0;
   for(int i=0; i<cigar.size(); i++) {
      bool check = true;
      if((cigar[i].op == 'S') || (cigar[i].op == 'H') || (cigar[i].op == 'I') || (cigar[i].op == 'P')) check = false;
      if(check) length += cigar[i].size;
   }
   return length;

}

int CigarParse_Vector(char *s, char *Str, int *Len) {

   int num = 0;
   char* p=(char*)(s);
   while(*p!=0)
   {
         char* endptr;

         if(!isdigit(*p)) return -1;
         Len[num] = (int)strtol(p,&endptr,10);
         if(Len[num]<=0)    return -1;
         p=endptr;
	 Str[num] = p[0];

         ++p;
	 num++;
   }  

   return num;
}


void reduce_Link_from_Header(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr)
{

  Data *data = (Data *) ptr;
  char *value;

  int *ind_1 = (int *) malloc( sizeof(int)*nvalues );
  int *ind_2 = (int *) malloc( sizeof(int)*nvalues );

  int *loc_1 = (int *) malloc( sizeof(int)*nvalues );
  int *loc_2 = (int *) malloc( sizeof(int)*nvalues );

  int n1 = -1;
  int n2 = -1;

    string header;
    header.assign(key,keybytes);

    PairLink pair;

    uint64_t nvalues_total;
    CHECK_FOR_BLOCKS(multivalue,valuebytes,nvalues,nvalues_total)
    BEGIN_BLOCK_LOOP(multivalue,valuebytes,nvalues)

    if( nvalues == 2 ) {

      value = multivalue;
      sReadInfo *link1 = (sReadInfo *) value;

      value += valuebytes[0];
      sReadInfo *link2 = (sReadInfo *) value;

      int sloc_1 = link1->loc - link1->SOFT_S;
      int sloc_2 = link2->loc - link2->SOFT_S;

      int eloc_1 = link1->loc + link1->align_len - 1 + link1->SOFT_E;
      int eloc_2 = link2->loc + link2->align_len - 1 + link2->SOFT_E;


      if( (link1->pe == 1) && (link2->pe == 2) ) {

	   pair.chr1 = link1->chr;
           pair.chr2 = link2->chr; 

           pair.strand1 = link1->strand;
           pair.strand2 = link2->strand;

           pair.maq1 = link1->maq;
           pair.maq2 = link2->maq;

	   if(pair.strand1 > 0) {
                        pair.loc1 = sloc_1;
                        pair.loc2 = eloc_2;
	   } else {
                        pair.loc1 = eloc_1;
                        pair.loc2 = sloc_2;
	   }

           pair.ID = 0;
           if((pair.maq1 >= data->MAQ) && (pair.maq2 >= data->MAQ)) kv->add(key,keybytes,(char *) &pair, sizeof(PairLink));

      } else {

	   pair.chr1 = link2->chr;
           pair.chr2 = link1->chr;

           pair.strand1 = link2->strand;
           pair.strand2 = link1->strand;

           pair.maq1 = link2->maq;
           pair.maq2 = link1->maq;

           if(pair.strand1 > 0) {
                        pair.loc1 = sloc_1;
                        pair.loc2 = eloc_2;
           } else {
                        pair.loc1 = eloc_1;
                        pair.loc2 = sloc_2;
           }

           pair.ID = 0;
           if((pair.maq1 >= data->MAQ) && (pair.maq2 >= data->MAQ)) kv->add(key,keybytes,(char *) &pair, sizeof(PairLink));

      } 

   }

    END_BLOCK_LOOP


    free(ind_1);
    free(ind_2);
    free(loc_1);
    free(loc_2);

}

void reduce_LinkID_HeadSeq(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr)
{
   Data *data = (Data *) ptr;
   char *value;

   HeadSeq *headseq_pair = (HeadSeq *) malloc( sizeof(HeadSeq)*2 );

   uint64_t nvalues_total;
   CHECK_FOR_BLOCKS(multivalue,valuebytes,nvalues,nvalues_total)
   BEGIN_BLOCK_LOOP(multivalue,valuebytes,nvalues)

   PairLink *locInfo;
   string header;
   bool check_pair = false;

   value = multivalue;
   for(int i=0; i<nvalues; i++) {
        if(valuebytes[i] == sizeof(PairLink)) {
          locInfo = (PairLink *) value;
	
          check_pair = true;
          break;
        }
        value += valuebytes[i];
   }

   if(check_pair) {

     int nH = 0;
     value = multivalue;
     for(int i=0; i<nvalues; i++) {

       if(valuebytes[i] != sizeof(PairLink)) {

	 HeadSeq headseq = *(HeadSeq *) value;	

	 if((headseq.PE == 1) && (headseq.length == data->range)) {
		headseq.strand  = locInfo->strand1;
		headseq.maq     = locInfo->maq1;
		if(headseq.strand < 0) {
                        revcomp_sequence_carray(headseq.seq, headseq.length);
                        revcomp_quality_carray(headseq.quality, headseq.length);
                }
		headseq.type = 1;

		headseq_pair[0] = headseq;
                nH += 1;	
	 } else if((headseq.PE == 2) && (headseq.length == data->range)) {
		headseq.strand  = locInfo->strand2;
                headseq.maq     = locInfo->maq2;
                if(headseq.strand < 0) {
                        revcomp_sequence_carray(headseq.seq, headseq.length);
                        revcomp_quality_carray(headseq.quality, headseq.length);
                }
		headseq.type = 1;

                headseq_pair[1] = headseq;
                nH += 1;
	} 

      }

      value += valuebytes[i];
     }

     if(nH == 2) {
	  locInfo->ID = 0;
	  locInfo->maq1 = 0;
          locInfo->maq2 = 0;
	  locInfo->strand1 = 0;
	  locInfo->strand2 = 0;

	  if(locInfo->loc1 > locInfo->loc2) {
		int tloc = locInfo->loc1;
		locInfo->loc1 = locInfo->loc2;
		locInfo->loc2 = tloc;

		tloc = locInfo->chr1;
		locInfo->chr1 = locInfo->chr2;
		locInfo->chr2 = tloc;
	  }

          kv->add((char *) locInfo, sizeof(PairLink), (char *) headseq_pair, sizeof(HeadSeq)*2 );
     }

   }

   END_BLOCK_LOOP

   free(headseq_pair);
}


void reduce_consensus_Single_highDup(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr)
{

   PairLink *plink = (PairLink *) key;

  if((nvalues >= 3) && (plink->chr1 > 0) && (plink->chr2 == 0)) {

   Data *data = (Data *) ptr;
   char *value;
   RefConsen *refconsen;
   string seq, quality;
   int PE, Strand;

   plink->ID = nvalues;

   IRKE irke(data->kmer_length, data->MAX_RECURSION, data->MIN_SEED_ENTROPY, data->MIN_SEED_COVERAGE,
                data->min_any_entropy, data->PACMAN, data->CRAWL, data->crawl_length, data->DS);

   StripedSmithWaterman::Aligner aligner(1,4,6,1);
   StripedSmithWaterman::Filter filter;
   StripedSmithWaterman::Alignment alignment;

   int *cquality = (int *) malloc( sizeof(int)*data->range );
   for(int i=0; i<data->range; i++) cquality[i]=0;

   uint64_t nvalues_total;
   CHECK_FOR_BLOCKS(multivalue,valuebytes,nvalues,nvalues_total)
   BEGIN_BLOCK_LOOP(multivalue,valuebytes,nvalues)

    plink->maq1 = 0;

    value = multivalue;
    for(int i=0; i<nvalues; i++) {
          HeadSeq headseq = *(HeadSeq *) value;
          plink->maq1 += headseq.maq;
          value += valuebytes[i];
    }
    plink->maq1 = (int)(plink->maq1 / nvalues);

      value = multivalue;
      for(int i=0; i<nvalues; i++) {
          HeadSeq headseq = *(HeadSeq *) value;

          seq.assign( headseq.seq, data->range );
          quality.assign( headseq.quality, data->range);
          PE = headseq.PE;
          Strand = headseq.strand;
          for(int j=0; j<data->range; j++) cquality[j] += get_Qscore( quality[j] );

          int nKmer = data->range - data->kmer_length + 1;
          for(int j=0; j<nKmer; j++) {
            string seq_kmer = seq.substr(j, data->kmer_length);
            if( contains_non_gatc(seq_kmer) ) {
              kmer_int_type_t kmer_val = kmer_to_intval(seq_kmer);
              irke.add_kmer(kmer_val , 1);
            }
          }

        value += valuebytes[i];
      }

      for(int j=0; j<data->range; j++) cquality[j]=(int)(cquality[j]/nvalues);
      for(int j=0; j<data->range; j++) quality[j] = get_QualityAscii( cquality[j] );

      irke.prune_some_kmers(data->min_kmer_count, data->min_any_entropy, data->prune_error_kmers, data->min_ratio_non_error);
      irke.populate_sorted_kmers_list();

      string csequence =
            irke.compute_sequence_assemblies( data->MIN_CONNECTIVITY_RATIO, data->MIN_ASSEMBLY_LENGTH, data->MIN_ASSEMBLY_COVERAGE );

      if((int)csequence.length() == data->range) {

        HeadSeq seq_quality;

        memcpy( seq_quality.seq,   csequence.substr(0,data->range).c_str(), sizeof(char)*data->range );
        memcpy( seq_quality.quality, quality.substr(0,data->range).c_str(), sizeof(char)*data->range );
        seq_quality.strand = Strand;
        seq_quality.length = data->range;
        seq_quality.PE     = PE;

        kv->add((char *) plink, sizeof(PairLink), (char *) &seq_quality, sizeof(HeadSeq) );

      }


   END_BLOCK_LOOP

   free(cquality);

  }

}


void print_Kmer_Cov_All(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr)
{
   char *value;
   Data *data = (Data *) ptr;

   kmer_int_type_t kmer_int = *(kmer_int_type_t *) key;
   string kmer = decode_kmer_from_intval( kmer_int, data->kmer_length );
   float entp = compute_entropy( kmer_int, data->kmer_length );

   data->outFile << kmer.c_str() << "\t" << nvalues << "\t" << entp << endl;
//   data->outFile << kmer_int << "\t" << nvalues << endl; 
}

void print_Kmer_Cov(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr)
{
   char *value;
   Data *data = (Data *) ptr;

   kmer_int_type_t kmer_int = *(kmer_int_type_t *) key;
   string seq_kmer = decode_kmer_from_intval(kmer_int, data->kmer_length);
//   string seq_kmer = decode_kmer_from_uint64(kmer_int, data->kmer_length);

   int num_trusted    = 0;
   int num_untrusted  = 0;
   //double qnum_trusted   = 0.0;
   //double qnum_untrusted = 0.0;

   uint64_t nvalues_total;
   CHECK_FOR_BLOCKS(multivalue,valuebytes,nvalues,nvalues_total)
   BEGIN_BLOCK_LOOP(multivalue,valuebytes,nvalues)

   value = multivalue;
   for(int i=0; i<nvalues; i++) {
      int *value_cov = (int *) value;     

      if(value_cov[0] > 0.0) {
	num_trusted = num_trusted + 1;
	//qnum_trusted += value_cov[1];
      } else if(value_cov[0] < 0) {
        num_untrusted = num_untrusted + 1;
        //qnum_untrusted += value_cov[1];
      }

      value += valuebytes[i];
   }

   END_BLOCK_LOOP

   data->outFile << seq_kmer.c_str() << "\t" << nvalues << "\t" << num_untrusted << "\t" << num_trusted << endl; //"\t" << qnum_untrusted << "\t" << qnum_trusted << endl;
 
}


void reduce_kmers_second_resuce(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr)
{
   char *value;
   Data *data = (Data *) ptr;
   bool onError = false;

   uint64_t nvalues_total;
   CHECK_FOR_BLOCKS(multivalue,valuebytes,nvalues,nvalues_total)
   BEGIN_BLOCK_LOOP(multivalue,valuebytes,nvalues)

   for(int i=0; i<nvalues; i++) {
        if(valuebytes[i] == sizeof(kmer_int_type_t)) {
                onError = true; break;
        }
   }

   if(onError) {
     int num = 0;
     value = multivalue;
     for(int i=0; i<nvalues; i++) {
        if(valuebytes[i] != sizeof(kmer_int_type_t)) {
	   int *dummy = (int *) value;
	   num = dummy[0]; break;
	}
        value += valuebytes[i];
     }

//     kmer_int_type_t kmer_val = *(kmer_int_type_t *) key;
//     float entropy = compute_entropy(kmer_val, data->kmer_length);
     if(num >= 2) kv->add(key, keybytes, key, keybytes);
   }

   END_BLOCK_LOOP
}


void reduce_kmers_from_singleton(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr)
{
   char *value;
   Data *data = (Data *) ptr;

   int *value_cov = (int *)malloc(sizeof(int)*2);
   value_cov[0] = nvalues;
   kv->add(key,keybytes,(char *) value_cov, sizeof(int)*2 ); 
}


void reduce_error_kmers(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr)
{
   char *value;
   Data *data = (Data *) ptr;

   int num_trusted    = 0;
   int num_untrusted  = 0;

   uint64_t nvalues_total;
   CHECK_FOR_BLOCKS(multivalue,valuebytes,nvalues,nvalues_total)
   BEGIN_BLOCK_LOOP(multivalue,valuebytes,nvalues)

   value = multivalue;
   for(int i=0; i<nvalues; i++) {
//      int *value_cov = (int *) value;
//      if(value_cov[0] > 0)      num_trusted++;
//      else if(value_cov[0] < 0) num_untrusted++;

      int dummy = *(int *) value;
      if(dummy > 0)      num_trusted++;
      else if(dummy < 0) num_untrusted++;

      value += valuebytes[i];
   }

//   kmer_int_type_t kmer_val = *(kmer_int_type_t *) key;
//   float entropy = compute_entropy(kmer_val, data->kmer_length);

   if(     (num_trusted >= 2) && (num_untrusted == 1)  && (data->ftype == 0)) kv->add(key, keybytes, key, keybytes);
   else if((num_trusted == 0) && (num_untrusted  > 0)  && (data->ftype == 1)) kv->add(key, keybytes, key, keybytes);
   else if((num_trusted == 0) && (data->ftype == 2)) kv->add(key, keybytes, key, keybytes);

//   if(num_trusted == 0) kv->add(key, keybytes, key, keybytes);

   END_BLOCK_LOOP

}

void reduce_error_onAmplicons(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr)
{
   char *value;
   Data *data = (Data *) ptr;
   bool onError = false;  

   uint64_t nvalues_total;
   CHECK_FOR_BLOCKS(multivalue,valuebytes,nvalues,nvalues_total)
   BEGIN_BLOCK_LOOP(multivalue,valuebytes,nvalues)

   for(int i=0; i<nvalues; i++) {
	if(valuebytes[i] == sizeof(kmer_int_type_t)) {
		onError = true;
		break;	
	}
   }  

   if(onError) {
     value = multivalue;
     for(int i=0; i<nvalues; i++) {
        if(valuebytes[i] != sizeof(kmer_int_type_t)) kv->add(value, valuebytes[i], key, keybytes);
        value += valuebytes[i];
     } 
   }

   END_BLOCK_LOOP
}


void reduce_extract_kmers_for_error_counts(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr)
{
  PairLink *plink = (PairLink *) key;
  Data *data = (Data *) ptr;

  if((nvalues >= data->min_csize) && (plink->chr1 > 0) && (plink->chr2 > 0) && (plink->chr1 == plink->chr2)) {
   
    char *value;
    string seq1, seq2, quality1, quality2;
    int PE1, PE2, strand1, strand2;

    map<kmer_int_type_t,int> _kmer_seq;

    int nseq = 0;
    int nKmer = data->range - data->kmer_length + 1;

    uint64_t nvalues_total;
    CHECK_FOR_BLOCKS(multivalue,valuebytes,nvalues,nvalues_total)
    BEGIN_BLOCK_LOOP(multivalue,valuebytes,nvalues)

    value = multivalue;
    for(int i=0; i<nvalues; i++) {
        HeadSeq *headseq_pair = (HeadSeq *) value;
	if( (headseq_pair[0].maq >= data->MAQ)&&(headseq_pair[1].maq >= data->MAQ) ) nseq++;

        value += valuebytes[i];
    }

    if(nseq >= data->min_csize) {

      value = multivalue;
      for(int i=0; i<nvalues; i++) {
          HeadSeq *headseq_pair = (HeadSeq *) value;

	  if( (headseq_pair[0].maq >= data->MAQ)&&(headseq_pair[1].maq >= data->MAQ) )  {

	    HeadSeq headseq = headseq_pair[0];
            seq1.assign( headseq.seq, data->range );
            quality1.assign( headseq.quality, data->range );

            headseq = headseq_pair[1];
            seq2.assign( headseq.seq, data->range );
            quality2.assign( headseq.quality, data->range );

            for(int j=0; j<nKmer; j++) {
              string seq_kmer1 = seq1.substr(j, data->kmer_length);
              if( contains_non_gatc(seq_kmer1) ) {
                kmer_int_type_t kmer_val = kmer_to_intval(seq_kmer1);
                _kmer_seq[kmer_val] = 1;
              }
            }

            for(int j=0; j<nKmer; j++) {
              string seq_kmer2 = seq2.substr(j, data->kmer_length);
              if( contains_non_gatc(seq_kmer2) ) {
                kmer_int_type_t kmer_val = kmer_to_intval(seq_kmer2);
                _kmer_seq[kmer_val] = 1;
              }
            }

         }

        value += valuebytes[i];
      }

    }

    END_BLOCK_LOOP

    if(nseq >= data->min_csize) { 
        for(map<kmer_int_type_t,int>::iterator it = _kmer_seq.begin(); it != _kmer_seq.end(); ++it)
                kv->add((char *) &it->first, sizeof(kmer_int_type_t), key, keybytes);
    }

  }
}

void reduce_ekmers_from_all(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr)
{
   if(nvalues == 1) kv->add(key, keybytes, multivalue, valuebytes[0]);
}

void reduce_All_Kmers(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr)
{
  Data *data = (Data *) ptr;

  int nKmer = data->range - data->kmer_length + 1;
  string seq1, seq2;
  char *value;
  Kmer_counter_map _kmer_map;

  uint64_t nvalues_total;
  CHECK_FOR_BLOCKS(multivalue,valuebytes,nvalues,nvalues_total)
  BEGIN_BLOCK_LOOP(multivalue,valuebytes,nvalues)

  value = multivalue;
  for(int i=0; i<nvalues; i++) {
    HeadSeq *read_pair = (HeadSeq *) value;
    seq1.assign( read_pair[0].seq, read_pair[0].length );
    seq2.assign( read_pair[1].seq, read_pair[1].length );
    
    for(int j=0; j<nKmer; j++) {
       string seq_kmer1 = seq1.substr(j, data->kmer_length);
       if( contains_non_gatc(seq_kmer1) ) {
           kmer_int_type_t kmer_val = kmer_to_intval(seq_kmer1);
	   _kmer_map[kmer_val] += 1;	  
       }
    }

    for(int j=0; j<nKmer; j++) {
       string seq_kmer2 = seq2.substr(j, data->kmer_length);
       if( contains_non_gatc(seq_kmer2) ) {
           kmer_int_type_t kmer_val = kmer_to_intval(seq_kmer2);
	   _kmer_map[kmer_val] += 1;
       }
    }

    value += valuebytes[i];
  }

  END_BLOCK_LOOP 

  for(Kmer_counter_map_iterator it=_kmer_map.begin(); it!=_kmer_map.end(); ++it) {
	kv->add((char *) &it->first, sizeof(kmer_int_type_t), (char *) &it->first, sizeof(kmer_int_type_t)); 
  }  

}


void reduce_Kmer_Cov_All(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr)
{
  PairLink *plink = (PairLink *) key;
  Data *data = (Data *) ptr;

  if((nvalues >= data->min_csize) && (plink->chr1 > 0) && (plink->chr2 > 0) && (plink->chr1 == plink->chr2)) {

      char *value;
      string seq1, seq2, quality1, quality2, csequence_1, csequence_2;
      int PE1, PE2, strand1, strand2;

      Kmer_counter_map _kmer_cseq;
      Kmer_counter_map _kmer_seq;

      IRKE irke_1(data->kmer_length, data->MAX_RECURSION, data->MIN_SEED_ENTROPY, data->MIN_SEED_COVERAGE,
                data->min_any_entropy, data->PACMAN, data->CRAWL, data->crawl_length, data->DS);

      IRKE irke_2(data->kmer_length, data->MAX_RECURSION, data->MIN_SEED_ENTROPY, data->MIN_SEED_COVERAGE,
                data->min_any_entropy, data->PACMAN, data->CRAWL, data->crawl_length, data->DS);

      int nseq = 0;
      int nKmer = data->range - data->kmer_length + 1;

      uint64_t nvalues_total;
      CHECK_FOR_BLOCKS(multivalue,valuebytes,nvalues,nvalues_total)
      BEGIN_BLOCK_LOOP(multivalue,valuebytes,nvalues)

      value = multivalue;
      for(int i=0; i<nvalues; i++) {
          HeadSeq *headseq_pair = (HeadSeq *) value;

          if((headseq_pair[0].maq >= data->MAQ) && (headseq_pair[1].maq >= data->MAQ)) {

	    nseq++;

            HeadSeq headseq = headseq_pair[0];
            seq1.assign( headseq.seq, data->range );
            quality1.assign( headseq.quality, data->range );
	
            for(int j=0; j<nKmer; j++) {
              string seq_kmer1 = seq1.substr(j, data->kmer_length);
	      string qkmer1 = quality1.substr(j,data->kmer_length);
              if( contains_non_gatc(seq_kmer1) && check_kmer_quality(qkmer1, data->bQscore)) {
		   kmer_int_type_t kmer_val = kmer_to_intval(seq_kmer1);
                   if(headseq.strand > 0) irke_1.add_kmer(kmer_val , 1);
		   else		          irke_2.add_kmer(kmer_val , 1);
                   _kmer_seq[kmer_val] += 1;
              }
            }

            headseq = headseq_pair[1];
            seq2.assign( headseq.seq, data->range );
            quality2.assign( headseq.quality, data->range );

            for(int j=0; j<nKmer; j++) {
              string seq_kmer2 = seq2.substr(j, data->kmer_length);
	      string qkmer2    = quality2.substr(j,data->kmer_length);
              if( contains_non_gatc(seq_kmer2) && check_kmer_quality(qkmer2, data->bQscore) ) {
		   kmer_int_type_t kmer_val = kmer_to_intval(seq_kmer2);
		   if(headseq.strand > 0) irke_1.add_kmer(kmer_val , 1);
                   else                   irke_2.add_kmer(kmer_val , 1);
                   _kmer_seq[kmer_val] += 1;
              }
            }

         }
         value += valuebytes[i];
      }

         irke_1.prune_some_kmers(data->min_kmer_count, data->min_any_entropy, data->prune_error_kmers, data->min_ratio_non_error);
         irke_1.populate_sorted_kmers_list();
         csequence_1 = irke_1.compute_sequence_assemblies( data->MIN_CONNECTIVITY_RATIO, data->MIN_ASSEMBLY_LENGTH, data->MIN_ASSEMBLY_COVERAGE );

         irke_2.prune_some_kmers(data->min_kmer_count, data->min_any_entropy, data->prune_error_kmers, data->min_ratio_non_error);
         irke_2.populate_sorted_kmers_list();
         csequence_2 = irke_2.compute_sequence_assemblies( data->MIN_CONNECTIVITY_RATIO, data->MIN_ASSEMBLY_LENGTH, data->MIN_ASSEMBLY_COVERAGE );

      if((nseq >= data->min_csize) && ((int)csequence_1.length() == data->range) && ((int)csequence_2.length() == data->range) ) {

        value = multivalue;
        for(int i=0; i<nvalues; i++) {
          HeadSeq *headseq_pair = (HeadSeq *) value;

          if((headseq_pair[0].maq >= data->MAQ) && (headseq_pair[1].maq >= data->MAQ)) {

            HeadSeq headseq = headseq_pair[0];
            seq1.assign( headseq.seq, data->range );
            quality1.assign( headseq.quality, data->range );

            for(int j=0; j<data->range; j++) {
                    int bid = get_int_from_base( seq1[j] );
//                    if((bid>=0) && (bid<=3) && (get_Qscore(quality1[j])>=data->bQscore)) {
		      if((bid>=0) && (bid<=3)) {
	                if(headseq.strand > 0) data->ATGC1[bid]++;
                        else                   data->ATGC2[bid]++;
                    }
            }

            headseq = headseq_pair[1];
            seq2.assign( headseq.seq, data->range );
            quality2.assign( headseq.quality, data->range );

            for(int j=0; j<data->range; j++) {
                    int bid = get_int_from_base( seq2[j] );
//                    if((bid>=0) && (bid<=3) && (get_Qscore(quality2[j]) >= data->bQscore)) {
		      if((bid>=0) && (bid<=3)) {
                        if(headseq.strand > 0) data->ATGC1[bid]++;
                        else                   data->ATGC2[bid]++;
                    }
            }

         }

         value += valuebytes[i];
        }

      }

      END_BLOCK_LOOP

      if((nseq >= data->min_csize) && ((int)csequence_1.length() == data->range) && ((int)csequence_2.length() == data->range) ) {

	    //nKmer = (int)csequence_1.length() - data->kmer_length + 1;
            for(int j=0; j<nKmer; j++) {
              string cseq_kmer = csequence_1.substr(j, data->kmer_length);
              if( contains_non_gatc(cseq_kmer) ) {
                kmer_int_type_t kmer_val = kmer_to_intval(cseq_kmer);
                _kmer_cseq[kmer_val] += 1;
              }
            }

	    //nKmer = (int)csequence_2.length() - data->kmer_length + 1;
            for(int j=0; j<nKmer; j++) {
              string cseq_kmer = csequence_2.substr(j, data->kmer_length);
              if( contains_non_gatc(cseq_kmer) ) {
                kmer_int_type_t kmer_val = kmer_to_intval(cseq_kmer);
                _kmer_cseq[kmer_val] += 1;
              }
            }

           for(Kmer_counter_map::iterator it = _kmer_seq.begin(); it != _kmer_seq.end(); ++it) {

	     int dummy = 0;
	     if(_kmer_cseq.find(it->first) == _kmer_cseq.end()) dummy = -1;
	     else					        dummy =  1;
	     kv->add((char *) &it->first, sizeof(kmer_int_type_t), (char *) &dummy, sizeof(int) );

           }

       }
  }
}

void reduce_error_matrices_step2(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr)
{
  PairLink *plink = (PairLink *) key;
  Data *data = (Data *) ptr;

  string seq1, quality1, seq2, quality2;
  string cseq1, cseq2;
  int num_ekmer;
  Kmer_counter_map _ekmer_map;

  int actual_base, obs_base, qscore;
  bool onCSEQ = false;

  char *value;
  uint64_t nvalues_total;
  CHECK_FOR_BLOCKS(multivalue,valuebytes,nvalues,nvalues_total)
  BEGIN_BLOCK_LOOP(multivalue,valuebytes,nvalues)

  value = multivalue;
  for(int i=0; i<nvalues; i++) {
     if(valuebytes[i] == sizeof(HeadSeq)) {
	HeadSeq *csequence = (HeadSeq *) value;	
	  if(csequence->PE >= data->min_csize) {
	  cseq1.assign(csequence->seq,     csequence->length);
	  cseq2.assign(csequence->quality, csequence->strand);
	  num_ekmer = csequence->type;
	  onCSEQ = true;
	  break;
	}
     }
     value += valuebytes[i];
  } 

  if(onCSEQ) {

     kmer_int_type_t *ekmer = (kmer_int_type_t *) malloc( sizeof(kmer_int_type_t)*num_ekmer );
     value = multivalue;
     for(int i=0; i<nvalues; i++) {
       if((valuebytes[i] == sizeof(kmer_int_type_t)*num_ekmer) && (valuebytes[i] != sizeof(HeadSeq)) && (valuebytes[i] != sizeof(HeadSeq)*2)) {
	  memcpy( ekmer, value, sizeof(kmer_int_type_t)*num_ekmer );	
	  break;
       }
       value += valuebytes[i];
     } 
     for(int j=0; j<num_ekmer; j++) _ekmer_map[ ekmer[j] ] = 1; 
     free(ekmer);

     value = multivalue;
     for(int i=0; i<nvalues; i++) {
        if((valuebytes[i] == sizeof(HeadSeq)*2) && (valuebytes[i] != sizeof(kmer_int_type_t)*num_ekmer)) {
	   HeadSeq *headseq_pair = (HeadSeq *) value;
	   HeadSeq headseq1      = headseq_pair[0];
	   HeadSeq headseq2      = headseq_pair[1];

	if((headseq1.length == (int)cseq1.length()) && (headseq2.length == (int)cseq2.length())) {

	   seq1.assign(     headseq1.seq,     headseq1.length );
           quality1.assign( headseq1.quality, headseq1.length );

           seq2.assign(     headseq2.seq,     headseq2.length );
           quality2.assign( headseq2.quality, headseq2.length );	

           vector<int> loc1;
	   vector<int> loc2;
	   int ek1 = 0;
	   int ek2 = 0;
	   int nkmer = headseq1.length - data->kmer_length + 1;
	   for(int j=0; j<nkmer; j++) {
	     if(contains_non_gatc( seq1.substr(j,data->kmer_length).c_str()) ) {
		kmer_int_type_t kmer_val  = kmer_to_intval(  seq1.substr(j,data->kmer_length) ); 
		kmer_int_type_t ckmer_val = kmer_to_intval( cseq1.substr(j,data->kmer_length) );
		if( (_ekmer_map.find(kmer_val) != _ekmer_map.end()) && (kmer_val != ckmer_val) ) {
			if(headseq1.strand > 0) {
				ek1++;
				loc1.push_back(j);
			} else {
                        	ek2++;
                        	loc2.push_back(j);
			}
		}
	     }
	   }

           nkmer = headseq2.length - data->kmer_length + 1;
           for(int j=0; j<nkmer; j++) {
	      if(contains_non_gatc( seq2.substr(j,data->kmer_length).c_str()) ) { 
               kmer_int_type_t kmer_val  = kmer_to_intval(  seq2.substr(j,data->kmer_length) );
               kmer_int_type_t ckmer_val = kmer_to_intval( cseq2.substr(j,data->kmer_length) );
               if( _ekmer_map.find(kmer_val) != _ekmer_map.end() && (kmer_val != ckmer_val) ) {
			if(headseq2.strand > 0) {
                                ek1++;
                                loc1.push_back(j);
                        } else {
                                ek2++;
                                loc2.push_back(j);
                        }
               }
	     }
           }

          if(ek1 > 0) {
	    set<short> eloc1;
	    vector<short> chop_region;
  	    vector< vector<int> > cc_loc1;
  	    cc_loc1.push_back(vector<int>());
  	    int cc=0;
  	    cc_loc1[cc].push_back(loc1[0]);

  	    for(int j=1; j<loc1.size(); j++) {
     	      if(loc1[j-1]+data->kmer_length-1 < loc1[j]) {
        	cc++;
        	cc_loc1.push_back(vector<int>());
     	      }
     	      cc_loc1[cc].push_back(loc1[j]);
  	    }

	    for(cc = 0; cc < cc_loc1.size(); cc++) {
		chop_region.clear();
		chop_region = error_region_extract(cc_loc1[cc], data->kmer_length, headseq_pair[0].length);	
		for(int k=0; k<chop_region.size(); k++) eloc1.insert( chop_region[k] );
	    }

	    for(set<short>::iterator it=eloc1.begin(); it!=eloc1.end(); it++) {
		
		if(headseq1.strand > 0) actual_base = get_int_from_base( cseq1[*it] );
		else			actual_base = get_int_from_base( cseq2[*it] );

		if(headseq1.strand > 0) obs_base    = get_int_from_base( seq1[*it] );	
		else			obs_base    = get_int_from_base( seq2[*it] );

		if(headseq1.strand > 0) qscore = get_Qscore( quality1[*it] );
		else			qscore = get_Qscore( quality2[*it] );	

		if((actual_base != obs_base)) {
			int *error_count = (int *) malloc( sizeof(int)*3 );
	
			error_count[0] = qscore; 
			error_count[1] = actual_base; 
			error_count[2] = obs_base;
			int PE = 1;
			kv->add((char *) error_count, sizeof(int)*3, (char *) &PE, sizeof(int));

			free(error_count);
		}		
	    }  

	   }


          if(ek2 > 0) {
            set<short> eloc2;
            vector<short> chop_region;
            vector< vector<int> > cc_loc2;
            cc_loc2.push_back(vector<int>());
            int cc=0;
            cc_loc2[cc].push_back(loc2[0]);

            for(int j=1; j<loc2.size(); j++) {
              if(loc2[j-1]+data->kmer_length-1 < loc2[j]) {
                cc++;
                cc_loc2.push_back(vector<int>());
              }
              cc_loc2[cc].push_back(loc2[j]);
            }

            for(cc = 0; cc < cc_loc2.size(); cc++) {
		chop_region.clear();
		chop_region = error_region_extract(cc_loc2[cc], data->kmer_length, headseq_pair[1].length);
                for(int k=0; k<chop_region.size(); k++) eloc2.insert( chop_region[k] );
            }

            for(set<short>::iterator it=eloc2.begin(); it!=eloc2.end(); it++) {

                if(headseq1.strand < 0) actual_base = get_int_from_base( cseq1[*it] );
                else                    actual_base = get_int_from_base( cseq2[*it] );

                if(headseq1.strand < 0) obs_base    = get_int_from_base( seq1[*it] );
                else                    obs_base    = get_int_from_base( seq2[*it] );

                if(headseq1.strand < 0) qscore = get_Qscore( quality1[*it] );
                else                    qscore = get_Qscore( quality2[*it] );

                if((actual_base != obs_base)) {
                        int *error_count = (int *) malloc( sizeof(int)*3 );

			int PE = 2;
			error_count[0] = qscore;
			error_count[1] = actual_base;
                        error_count[2] = obs_base;	
			kv->add((char *) error_count, sizeof(int)*3, (char *) &PE, sizeof(int));

	                free(error_count);
                }
            }

          }

	}
        }
        value += valuebytes[i];
     }

  }

  END_BLOCK_LOOP

}


void reduce_ekmers_rescue_step2(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr)
{
  Data *data = (Data *) ptr;
  char *value;
  string seq1, quality1, seq2, quality2;

  int nKmer = data->range - data->kmer_length + 1;

  uint64_t nvalues_total;
  CHECK_FOR_BLOCKS(multivalue,valuebytes,nvalues,nvalues_total)
  BEGIN_BLOCK_LOOP(multivalue,valuebytes,nvalues)

  bool check_error = false;
  value = multivalue;
  for(int i=0; i<nvalues; i++) {
        if(valuebytes[i] == sizeof(HeadSeq)*2) {
 	   HeadSeq *headseq = (HeadSeq *) value;
	   if((headseq[0].type < 0) && (headseq[1].type < 0)) { check_error = true; break; }
        }
	value += valuebytes[i];
  }

  if(check_error) {

     string csequence_1, csequence_2;
     value = multivalue;
     for(int i=0; i<nvalues; i++) {
        if(valuebytes[i] == sizeof(HeadSeq)*2) {
           HeadSeq *headseq = (HeadSeq *) value;
           if((headseq[0].type < 0) && (headseq[1].type < 0)) { 
		csequence_1.assign( headseq[0].seq, headseq[0].length );
		csequence_2.assign( headseq[1].seq, headseq[1].length );
		break; 
	   }
        }
        value += valuebytes[i];
     }	

     KmerCounter _ckmer;
     for(int j=0; j<nKmer; j++) {
                        string seq_kmer1 = csequence_1.substr(j, data->kmer_length);
                        if( contains_non_gatc(seq_kmer1) ) {
                              kmer_int_type_t kmer_val = kmer_to_intval(seq_kmer1);
                              _ckmer.add_kmer(kmer_val,1);
                        }
                        string seq_kmer2 = csequence_2.substr(j, data->kmer_length);
                        if( contains_non_gatc(seq_kmer2) ) {
                              kmer_int_type_t kmer_val = kmer_to_intval(seq_kmer2);
                              _ckmer.add_kmer(kmer_val,1);
                        }
     }	


     vector<int> loc1;
     vector<int> loc2;
     int ek1 = 0;
     int ek2 = 0;

     value = multivalue;
     for(int i=0; i<nvalues; i++) {
          if(valuebytes[i] == sizeof(HeadSeq)*2) {
              HeadSeq *headseq = (HeadSeq *) value;

	      if((headseq[0].type > 0) && (headseq[1].type > 0) && (headseq[0].maq >= data->MAQ) && (headseq[1].maq >= data->MAQ)) { 

                seq1.assign( headseq[0].seq, data->range );
                seq2.assign( headseq[1].seq, data->range );

                for(int j=0; j<nKmer; j++) {
                        string seq_kmer1 = seq1.substr(j, data->kmer_length);
                        if( contains_non_gatc(seq_kmer1) ) {
                              kmer_int_type_t kmer_val = kmer_to_intval(seq_kmer1);
			        if(!_ckmer.kmer_exists(kmer_val)) { ek1++; loc1.push_back(j); }
                        }
                }

                for(int j=0; j<nKmer; j++) {
                        string seq_kmer2 = seq2.substr(j, data->kmer_length);
                        if( contains_non_gatc(seq_kmer2) ) {
                              kmer_int_type_t kmer_val = kmer_to_intval(seq_kmer2);
			      if(!_ckmer.kmer_exists(kmer_val)) { ek2++; loc2.push_back(j); } 
                       }
                }

                set<short> eloc1;
		set<short> eloc2;

		if(ek1 > 0) {

                    vector<short> chop_region;
                    vector< vector<int> > cc_loc1;
                    cc_loc1.push_back(vector<int>());
                    int cc=0;
                    cc_loc1[cc].push_back(loc1[0]);

                    for(int j=1; j<loc1.size(); j++) {
                     if(loc1[j-1]+data->kmer_length-1 < loc1[j]) {
                       cc++;
                       cc_loc1.push_back(vector<int>());
                     }
                     cc_loc1[cc].push_back(loc1[j]);
                    }

                    for(cc = 0; cc < cc_loc1.size(); cc++) {
                     chop_region.clear();
                     chop_region = error_region_extract(cc_loc1[cc], data->kmer_length, headseq[0].length);
                     for(int k=0; k<chop_region.size(); k++) eloc1.insert( chop_region[k] );
                    }

		}

		if(ek2 > 0) {

                    vector<short> chop_region;
                    vector< vector<int> > cc_loc2;
                    cc_loc2.push_back(vector<int>());
                    int cc=0;
                    cc_loc2[cc].push_back(loc2[0]);

                    for(int j=1; j<loc2.size(); j++) {
                       if(loc2[j-1]+data->kmer_length-1 < loc2[j]) {
                          cc++;
                          cc_loc2.push_back(vector<int>());
                       }
                       cc_loc2[cc].push_back(loc2[j]);
                    }

                    for(cc = 0; cc < cc_loc2.size(); cc++) {
                       chop_region.clear();
                       chop_region = error_region_extract(cc_loc2[cc], data->kmer_length, headseq[1].length);
                       for(int k=0; k<chop_region.size(); k++) eloc2.insert( chop_region[k] );
                    }

	        }


		if(eloc1.empty() && eloc2.empty()) {
                   PairLink plink = *(PairLink *) key;
		   plink.ID += 2*data->ftype;
                   kv->add((char *) &plink, sizeof(PairLink), value, valuebytes[i]);
		} else {
                   PairLink plink = *(PairLink *) key;
                   plink.ID += 2*data->ftype+1;
                   kv->add((char *) &plink, sizeof(PairLink), value, valuebytes[i]);
		} 

	    }

	}
	value += valuebytes[i];
    }


  } else {

    value = multivalue;
    for(int i=0; i<nvalues; i++) {
        if(valuebytes[i] == sizeof(HeadSeq)*2) kv->add(key, keybytes, value, valuebytes[i]); 
        value += valuebytes[i];
    }

  }

  END_BLOCK_LOOP

}


void reduce_rescue_by_ekmer(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr)
{
  Data *data = (Data *) ptr;
  char *value;
  string seq1, quality1, seq2, quality2;

  int nKmer = data->range - data->kmer_length + 1;
  int num_ekmer = 0;
  KmerCounter _ekmer;

  uint64_t nvalues_total;
  CHECK_FOR_BLOCKS(multivalue,valuebytes,nvalues,nvalues_total)
  BEGIN_BLOCK_LOOP(multivalue,valuebytes,nvalues)

    value = multivalue;
    for(int i=0; i<nvalues; i++) {
        if(valuebytes[i] == sizeof(kmer_int_type_t)) {
	    kmer_int_type_t kmer_val = *(kmer_int_type_t *) value;
	    num_ekmer++; _ekmer.add_kmer(kmer_val,1);
	}
	value += valuebytes[i];
    }

    if(num_ekmer > 0) {

       value = multivalue;
       for(int i=0; i<nvalues; i++) {
            if(valuebytes[i] == sizeof(HeadSeq)*2) {
              HeadSeq *headseq = (HeadSeq *) value;
	      seq1.assign( headseq[0].seq, data->range );
	      seq2.assign( headseq[1].seq, data->range );

              vector<int> loc1;
              vector<int> loc2;
              int ek1 = 0;
              int ek2 = 0;
	     
                for(int j=0; j<nKmer; j++) {
                        string seq_kmer1 = seq1.substr(j, data->kmer_length);
                        if( contains_non_gatc(seq_kmer1) ) {
                              kmer_int_type_t kmer_val = kmer_to_intval(seq_kmer1);
                                if(!_ekmer.kmer_exists(kmer_val)) { ek1++; loc1.push_back(j); }
                        }
                }

                for(int j=0; j<nKmer; j++) {
                        string seq_kmer2 = seq2.substr(j, data->kmer_length);
                        if( contains_non_gatc(seq_kmer2) ) {
                              kmer_int_type_t kmer_val = kmer_to_intval(seq_kmer2);
                              if(!_ekmer.kmer_exists(kmer_val)) { ek2++; loc2.push_back(j); }
                       }
                }

                set<short> eloc1;
                set<short> eloc2;

                if(ek1 > 0) {

                    vector<short> chop_region;
                    vector< vector<int> > cc_loc1;
                    cc_loc1.push_back(vector<int>());
                    int cc=0;
                    cc_loc1[cc].push_back(loc1[0]);

                    for(int j=1; j<loc1.size(); j++) {
                     if(loc1[j-1]+data->kmer_length-1 < loc1[j]) {
                       cc++;
                       cc_loc1.push_back(vector<int>());
                     }
                     cc_loc1[cc].push_back(loc1[j]);
                    }

                    for(cc = 0; cc < cc_loc1.size(); cc++) {
                     chop_region.clear();
                     chop_region = error_region_extract(cc_loc1[cc], data->kmer_length, headseq[0].length);
                     for(int k=0; k<chop_region.size(); k++) eloc1.insert( chop_region[k] );
                    }

                }

                if(ek2 > 0) {

                    vector<short> chop_region;
                    vector< vector<int> > cc_loc2;
                    cc_loc2.push_back(vector<int>());
                    int cc=0;
                    cc_loc2[cc].push_back(loc2[0]);

                    for(int j=1; j<loc2.size(); j++) {
                       if(loc2[j-1]+data->kmer_length-1 < loc2[j]) {
                          cc++;
                          cc_loc2.push_back(vector<int>());
                       }
                       cc_loc2[cc].push_back(loc2[j]);
                    }

                    for(cc = 0; cc < cc_loc2.size(); cc++) {
                       chop_region.clear();
                       chop_region = error_region_extract(cc_loc2[cc], data->kmer_length, headseq[1].length);
                       for(int k=0; k<chop_region.size(); k++) eloc2.insert( chop_region[k] );
                    }

                }

                if(eloc1.empty() && eloc2.empty()) {
                   kv->add(key, keybytes, value, valuebytes[i]);
                } else {
                   PairLink plink = *(PairLink *) key;
                   plink.ID += 1;
                   kv->add((char *) &plink, sizeof(PairLink), value, valuebytes[i]);
                }

            }
	    value += valuebytes[i];
       }

   } else {

    value = multivalue;
    for(int i=0; i<nvalues; i++) {
	kv->add(key, keybytes, value, valuebytes[i]);
        value += valuebytes[i];
    }

  }

  END_BLOCK_LOOP
}


void reduce_ekmers_rescue_step1(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr)
{
  Data *data = (Data *) ptr;
  char *value;
  string seq1, quality1, seq2, quality2;

  int nKmer = data->range - data->kmer_length + 1;

  int num_ekmer = 0;
  int nseq = 0;

  IRKE irke_1(data->kmer_length, data->MAX_RECURSION, data->MIN_SEED_ENTROPY, data->MIN_SEED_COVERAGE,
                data->min_any_entropy, data->PACMAN, data->CRAWL, data->crawl_length, data->DS);

  IRKE irke_2(data->kmer_length, data->MAX_RECURSION, data->MIN_SEED_ENTROPY, data->MIN_SEED_COVERAGE,
                data->min_any_entropy, data->PACMAN, data->CRAWL, data->crawl_length, data->DS);

  uint64_t nvalues_total;
  CHECK_FOR_BLOCKS(multivalue,valuebytes,nvalues,nvalues_total)
  BEGIN_BLOCK_LOOP(multivalue,valuebytes,nvalues)

    value = multivalue;
    for(int i=0; i<nvalues; i++) {
        if(valuebytes[i] == sizeof(kmer_int_type_t)) num_ekmer++;
    }


    if(num_ekmer > 0) {

      value = multivalue;
      for(int i=0; i<nvalues; i++) {
        if(valuebytes[i] == sizeof(HeadSeq)*2) {
          HeadSeq *headseq = (HeadSeq *) value;
          if((headseq[0].maq >= data->MAQ)&&(headseq[1].maq >= data->MAQ)) nseq++;
        }
        value += valuebytes[i];
      }

      if(nseq >= data->min_csize) {

          value = multivalue;
          for(int i=0; i<nvalues; i++) {
            if(valuebytes[i] == sizeof(HeadSeq)*2) {
              HeadSeq *headseq = (HeadSeq *) value;

              if(headseq[0].maq >= data->MAQ) {
                seq1.assign( headseq[0].seq, data->range );
                for(int j=0; j<nKmer; j++) {
                    string seq_kmer1 = seq1.substr(j, data->kmer_length);
                    if( contains_non_gatc(seq_kmer1) ) {
                          kmer_int_type_t kmer_val = kmer_to_intval(seq_kmer1);
                          if(headseq[0].strand > 0) irke_1.add_kmer(kmer_val, 1);
                          else                      irke_2.add_kmer(kmer_val, 1);
                    }
               }

              }

              if(headseq[1].maq >= data->MAQ) {
                seq2.assign( headseq[1].seq, data->range );
                for(int j=0; j<nKmer; j++) {
                    string seq_kmer2 = seq2.substr(j, data->kmer_length);
                    if( contains_non_gatc(seq_kmer2) ) {
                          kmer_int_type_t kmer_val = kmer_to_intval(seq_kmer2);
                          if(headseq[1].strand > 0) irke_1.add_kmer(kmer_val, 1);
                          else                      irke_2.add_kmer(kmer_val, 1);
                    }
               }
              }

            }
            value += valuebytes[i];
          }

          irke_1.prune_some_kmers(data->min_kmer_count, data->min_any_entropy, data->prune_error_kmers, data->min_ratio_non_error);
          irke_1.populate_sorted_kmers_list();
          string csequence_1 = irke_1.compute_sequence_assemblies( data->MIN_CONNECTIVITY_RATIO, data->MIN_ASSEMBLY_LENGTH, data->MIN_ASSEMBLY_COVERAGE );

          irke_2.prune_some_kmers(data->min_kmer_count, data->min_any_entropy, data->prune_error_kmers, data->min_ratio_non_error);
          irke_2.populate_sorted_kmers_list();
          string csequence_2 = irke_2.compute_sequence_assemblies( data->MIN_CONNECTIVITY_RATIO, data->MIN_ASSEMBLY_LENGTH, data->MIN_ASSEMBLY_COVERAGE );

	  if( ((int)csequence_1.length() == data->range) && ((int)csequence_2.length() == data->range) ) {

            HeadSeq seq_quality1, seq_quality2;

            memcpy( seq_quality1.seq, csequence_1.c_str(), sizeof(char)*data->range );
            memcpy( seq_quality2.seq, csequence_2.c_str(), sizeof(char)*data->range );

            memcpy( seq_quality1.quality, quality1.c_str(), sizeof(char)*data->range );
            memcpy( seq_quality2.quality, quality2.c_str(), sizeof(char)*data->range );

            seq_quality1.length = data->range;
            seq_quality2.length = data->range;

            seq_quality1.type = -1 * nvalues;
            seq_quality2.type = -1 * nvalues;

            HeadSeq *out_headseq = (HeadSeq *) malloc( sizeof(HeadSeq) * 2 );
            out_headseq[0] = seq_quality1;
            out_headseq[1] = seq_quality2;
	  
	    kv->add(key, keybytes, (char *) out_headseq, sizeof(HeadSeq)*2 );	
 
	    free(out_headseq); 

	  }

      }

    }

  END_BLOCK_LOOP

}


void reduce_error_matrices_step1(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr)
{
  PairLink *plink = (PairLink *) key;
  Data *data = (Data *) ptr;
  char *value;
  string seq1, quality1, seq2, quality2;

  vector<kmer_int_type_t> _ekmer;
  int num_ekmer = 0;
  int nseq = 0;

  IRKE irke_1(data->kmer_length, data->MAX_RECURSION, data->MIN_SEED_ENTROPY, data->MIN_SEED_COVERAGE,
                data->min_any_entropy, data->PACMAN, data->CRAWL, data->crawl_length, data->DS);

  IRKE irke_2(data->kmer_length, data->MAX_RECURSION, data->MIN_SEED_ENTROPY, data->MIN_SEED_COVERAGE,
                data->min_any_entropy, data->PACMAN, data->CRAWL, data->crawl_length, data->DS);

  uint64_t nvalues_total;
  CHECK_FOR_BLOCKS(multivalue,valuebytes,nvalues,nvalues_total)
  BEGIN_BLOCK_LOOP(multivalue,valuebytes,nvalues)

  value = multivalue;
  for(int i=0; i<nvalues; i++) {
	if(valuebytes[i] == sizeof(kmer_int_type_t)) num_ekmer++; 
  }

if(num_ekmer > 0) {

      value = multivalue;
      for(int i=0; i<nvalues; i++) {
	if(valuebytes[i] == sizeof(HeadSeq)*2) {
          HeadSeq *headseq_pair = (HeadSeq *) value;
	  if((headseq_pair[0].maq >= data->MAQ)&&(headseq_pair[1].maq >= data->MAQ)&&(headseq_pair[0].length==data->range)&&(headseq_pair[1].length==data->range)) nseq++;
        }
	value += valuebytes[i];
      }

  if(nseq >= data->min_csize) {

     num_ekmer = 0;
     value = multivalue;
     for(int i=0; i<nvalues; i++) {
        if(valuebytes[i] == sizeof(kmer_int_type_t)) {
	   kmer_int_type_t tkmer = *(kmer_int_type_t *) value;
	   _ekmer.push_back(tkmer);
	   num_ekmer++;
	} else { 
	   HeadSeq *headseq_pair = (HeadSeq *) value;
	   if((headseq_pair[0].maq >= data->MAQ)&&(headseq_pair[1].maq >= data->MAQ)&&(headseq_pair[0].length==data->range)&&(headseq_pair[1].length==data->range)) {

	       int nKmer = data->range - data->kmer_length + 1; 
	
	       HeadSeq headseq = headseq_pair[0];
               seq1.assign( headseq.seq, data->range );
               quality1.assign( headseq.quality, data->range );

               for(int j=0; j<nKmer; j++) {
                 string seq_kmer1 = seq1.substr(j, data->kmer_length);
	         string qkmer1 = quality1.substr(j,data->kmer_length);
                 if( contains_non_gatc(seq_kmer1) && check_kmer_quality(qkmer1, data->bQscore) ) {
                   kmer_int_type_t kmer_val = kmer_to_intval(seq_kmer1);
                   if(headseq.strand > 0) irke_1.add_kmer(kmer_val , 1);
		   else		          irke_2.add_kmer(kmer_val , 1);
                 }
               }

               headseq = headseq_pair[1];
               seq2.assign( headseq.seq, data->range );
               quality2.assign( headseq.quality, data->range );

               for(int j=0; j<nKmer; j++) {
                 string seq_kmer2 = seq2.substr(j, data->kmer_length);
	         string qkmer2 = quality2.substr(j, data->kmer_length);
                 if( contains_non_gatc(seq_kmer2) && check_kmer_quality(qkmer2, data->bQscore) ) {
		    kmer_int_type_t kmer_val = kmer_to_intval(seq_kmer2);
                    if(headseq.strand > 0) irke_1.add_kmer(kmer_val , 1);
                    else                   irke_2.add_kmer(kmer_val , 1);
                 }
              }

          }
	}

	value += valuebytes[i];
     }

     irke_1.prune_some_kmers(data->min_kmer_count, data->min_any_entropy, data->prune_error_kmers, data->min_ratio_non_error);
     irke_1.populate_sorted_kmers_list();
     string csequence_1 = irke_1.compute_sequence_assemblies( data->MIN_CONNECTIVITY_RATIO, data->MIN_ASSEMBLY_LENGTH, data->MIN_ASSEMBLY_COVERAGE );

     irke_2.prune_some_kmers(data->min_kmer_count, data->min_any_entropy, data->prune_error_kmers, data->min_ratio_non_error);
     irke_2.populate_sorted_kmers_list();
     string csequence_2 = irke_2.compute_sequence_assemblies( data->MIN_CONNECTIVITY_RATIO, data->MIN_ASSEMBLY_LENGTH, data->MIN_ASSEMBLY_COVERAGE );

     //int sallow = data->range - (int)( (float)data->range * data->soft_allow);

     if(((int)csequence_1.length() == data->range) && ((int)csequence_2.length() == data->range) ) {

         HeadSeq csequence;
         memcpy(csequence.seq,     csequence_1.c_str(), sizeof(char)*csequence_1.length());
         memcpy(csequence.quality, csequence_2.c_str(), sizeof(char)*csequence_2.length()); 

	 csequence.length = (int)csequence_1.length();
	 csequence.strand = (int)csequence_2.length();
	 csequence.type = num_ekmer;
	 csequence.PE = nseq;

	 kv->add(key, keybytes, (char *) &csequence, sizeof(HeadSeq));

         value = multivalue;
         for(int i=0; i<nvalues; i++) {
	    if(valuebytes[i] != sizeof(kmer_int_type_t)) kv->add(key, keybytes, value, valuebytes[i]);
	    value += valuebytes[i];
	 }

         kmer_int_type_t *error_kmers = (kmer_int_type_t *) malloc(sizeof(kmer_int_type_t)*num_ekmer);
	 for(int i=0; i<num_ekmer; i++) error_kmers[i] = _ekmer[i];
         kv->add(key, keybytes, (char *) error_kmers, sizeof(kmer_int_type_t)*num_ekmer);
         free(error_kmers);

     }

  }

}

  END_BLOCK_LOOP

}

void reduce_ntnt_count(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr)
{
  int *nt_count = (int *) key;
  Data *data =(Data *) ptr;

  char *value;
  uint64_t nvalues_total;
  CHECK_FOR_BLOCKS(multivalue,valuebytes,nvalues,nvalues_total)
  BEGIN_BLOCK_LOOP(multivalue,valuebytes,nvalues)

  value = multivalue;
  for(int i=0; i<nvalues; i++) {

      int index = *(int *) value;
      if(index == 1)      data->qs1_count[   nt_count[0] ][ nt_count[1] ][ nt_count[2] ]++;
      else if(index == 2) data->qs2_count[   nt_count[0] ][ nt_count[1] ][ nt_count[2] ]++;  
      value += valuebytes[i];
  }

  END_BLOCK_LOOP
}


void reduce_Kmer_Cov(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr)
{
  PairLink *plink = (PairLink *) key;
  Data *data = (Data *) ptr;

  if((nvalues >= data->min_csize) && (plink->chr1 > 0) && (plink->chr2 > 0) && (plink->chr1 == plink->chr2)) {

      char *value;
      string seq1, seq2, quality1, quality2;
      int PE1, PE2, strand1, strand2;

      Kmer_counter_map _kmer_cseq1, _kmer_cseq2;
      Kmer_counter_map _kmer_seq1, _kmer_seq2;

      plink->ID = nvalues;

      IRKE irke_1(data->kmer_length, data->MAX_RECURSION, data->MIN_SEED_ENTROPY, data->MIN_SEED_COVERAGE,
                data->min_any_entropy, data->PACMAN, data->CRAWL, data->crawl_length, data->DS);

      IRKE irke_2(data->kmer_length, data->MAX_RECURSION, data->MIN_SEED_ENTROPY, data->MIN_SEED_COVERAGE,
                data->min_any_entropy, data->PACMAN, data->CRAWL, data->crawl_length, data->DS);

      int nseq1 = 0;
      int nseq2 = 0;
      int nKmer = data->range - data->kmer_length + 1;

      uint64_t nvalues_total;
      CHECK_FOR_BLOCKS(multivalue,valuebytes,nvalues,nvalues_total)
      BEGIN_BLOCK_LOOP(multivalue,valuebytes,nvalues)

      value = multivalue;
      for(int i=0; i<nvalues; i++) {
          HeadSeq *headseq_pair = (HeadSeq *) value;

          HeadSeq headseq = headseq_pair[0];
          if(headseq.maq >= data->MAQ) {

            seq1.assign( headseq.seq, data->range );
            quality1.assign( headseq.quality, data->range );
            PE1 = headseq.PE;
            strand1 = headseq.strand;

            for(int j=0; j<nKmer; j++) {
              string seq_kmer1 = seq1.substr(j, data->kmer_length);
              if( contains_non_gatc(seq_kmer1) ) {
                kmer_int_type_t kmer_val = kmer_to_intval(seq_kmer1);
                irke_1.add_kmer(kmer_val , 1);
		_kmer_seq1[kmer_val] += 1;
              }
            }

            nseq1++;
          }

          headseq = headseq_pair[1];
          if(headseq.maq >= data->MAQ) {

            seq2.assign( headseq.seq, data->range );
            quality2.assign( headseq.quality, data->range );
            PE2 = headseq.PE;
            strand2 = headseq.strand;

            for(int j=0; j<nKmer; j++) {
              string seq_kmer2 = seq2.substr(j, data->kmer_length);
              if( contains_non_gatc(seq_kmer2) ) {
                kmer_int_type_t kmer_val = kmer_to_intval(seq_kmer2);
                irke_2.add_kmer(kmer_val , 1);
		_kmer_seq2[kmer_val] += 1;
              }
            }

            nseq2++;
         }

        value += valuebytes[i];
      }

      END_BLOCK_LOOP


      if((nseq1 >= data->min_csize) && (nseq2 >= data->min_csize)) {

         irke_1.prune_some_kmers(data->min_kmer_count, data->min_any_entropy, data->prune_error_kmers, data->min_ratio_non_error);
         irke_1.populate_sorted_kmers_list();
         string csequence_1 = irke_1.compute_sequence_assemblies( data->MIN_CONNECTIVITY_RATIO, data->MIN_ASSEMBLY_LENGTH, data->MIN_ASSEMBLY_COVERAGE );

         irke_2.prune_some_kmers(data->min_kmer_count, data->min_any_entropy, data->prune_error_kmers, data->min_ratio_non_error);
         irke_2.populate_sorted_kmers_list();
         string csequence_2 = irke_2.compute_sequence_assemblies( data->MIN_CONNECTIVITY_RATIO, data->MIN_ASSEMBLY_LENGTH, data->MIN_ASSEMBLY_COVERAGE );

 	 string cseq_kmer;

	 if(csequence_1.length() == data->range) {
            for(int j=0; j<nKmer; j++) {
              cseq_kmer = csequence_1.substr(j, data->kmer_length);
              if( contains_non_gatc(cseq_kmer) ) {
                kmer_int_type_t kmer_val = kmer_to_intval(cseq_kmer);
                _kmer_cseq1[kmer_val] += 1;
              }
            }
	 }	

         if(csequence_2.length() == data->range) {
            for(int j=0; j<nKmer; j++) {
              cseq_kmer = csequence_2.substr(j, data->kmer_length);
              if( contains_non_gatc(cseq_kmer) ) {
                kmer_int_type_t kmer_val = kmer_to_intval(cseq_kmer);
                _kmer_cseq2[kmer_val] += 1;
              }
            }
         }

     for(Kmer_counter_map::iterator it = _kmer_seq1.begin(); it != _kmer_seq1.end(); ++it) {
           int *buffer_value = (int *)malloc(sizeof(int)*2 );
           buffer_value[0] = 0;
           buffer_value[1] = 1;

           if(_kmer_cseq1.find(it->first) == _kmer_cseq1.end()) buffer_value[0] = -1;
           else                                                 buffer_value[0] =  1;

	   kv->add((char *) &it->first, sizeof(kmer_int_type_t), (char *) buffer_value, sizeof(int)*2 );
//	   uint64_t kmer_int = convert_mpz_uint64(it->first, data->kmer_length);
//	   kv->add((char *) &kmer_int, sizeof(uint64_t), (char *) buffer_value, sizeof(int)*2 ); 

           free(buffer_value);	
     }           


     for(Kmer_counter_map::iterator it = _kmer_seq2.begin(); it != _kmer_seq2.end(); ++it) {
           int *buffer_value = (int *)malloc(sizeof(int)*2 );
           buffer_value[0] = 0;
           buffer_value[1] = 1;

           if(_kmer_cseq2.find(it->first) == _kmer_cseq2.end()) buffer_value[0] = -1;
           else                                                 buffer_value[0] =  1;

	   kv->add((char *) &it->first, sizeof(kmer_int_type_t), (char *) buffer_value, sizeof(int)*2 ); 
//           uint64_t kmer_int = convert_mpz_uint64(it->first, data->kmer_length);
//           kv->add((char *) &kmer_int, sizeof(uint64_t), (char *) buffer_value, sizeof(int)*2 );

           free(buffer_value);
     }


  } 

 }
}


void reduce_consensus_highDup(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr)
{
   PairLink plink = *(PairLink *) key;
   Data *data = (Data *) ptr;

 if((nvalues >= data->min_csize) && (plink.chr1 > 0) && (plink.chr2 > 0)) {

   char *value;
   string seq_first, seq_second, quality_first, quality_second;

   IRKE irke_first(data->kmer_length, data->MAX_RECURSION, data->MIN_SEED_ENTROPY, data->MIN_SEED_COVERAGE,
                   data->min_any_entropy, data->PACMAN, data->CRAWL, data->crawl_length, data->DS);

   IRKE irke_second(data->kmer_length, data->MAX_RECURSION, data->MIN_SEED_ENTROPY, data->MIN_SEED_COVERAGE,
                   data->min_any_entropy, data->PACMAN, data->CRAWL, data->crawl_length, data->DS);

   StripedSmithWaterman::Aligner aligner(1,4,6,1);
   StripedSmithWaterman::Filter filter;
   StripedSmithWaterman::Alignment alignment;

   int *cquality_first  = (int *) malloc( sizeof(int)*data->range );
   int *cquality_second = (int *) malloc( sizeof(int)*data->range );
   for(int i=0; i<data->range; i++) { cquality_first[i]=0; cquality_second[i]=0; }

   int maq_first = 0;
   int maq_second = 0;

   int nseq = 0;
   int nKmer = data->range - data->kmer_length + 1;

   uint64_t nvalues_total;
   CHECK_FOR_BLOCKS(multivalue,valuebytes,nvalues,nvalues_total)
   BEGIN_BLOCK_LOOP(multivalue,valuebytes,nvalues)

    value = multivalue;
    for(int i=0; i<nvalues; i++) {
      if(valuebytes[i] == sizeof(HeadSeq)*2) {
          HeadSeq *headseq_pair = (HeadSeq *) value;
          if((headseq_pair[0].maq >= data->MAQ) && (headseq_pair[1].maq >= data->MAQ) && (headseq_pair[0].length == data->range) && (headseq_pair[1].length == data->range)) nseq++; 
      }
        value += valuebytes[i];
   }

   if(nseq >= data->min_csize) {

      value = multivalue;
      for(int i=0; i<nvalues; i++) {
          HeadSeq *headseq_pair = (HeadSeq *) value;

          HeadSeq headseq = headseq_pair[0];
          if((headseq.maq >= data->MAQ) && (headseq.length == data->range)) {

            maq_first += headseq.maq;
            seq_first.assign( headseq.seq, data->range );
            quality_first.assign( headseq.quality, data->range );

            for(int j=0; j<data->range; j++) {

		if(headseq.strand > 0) {
			if(get_Qscore(quality_first[j]) > cquality_first[j])  cquality_first[j]  = get_Qscore( quality_first[j] ); 
		} else {
			if(get_Qscore(quality_first[j]) > cquality_second[j]) cquality_second[j] = get_Qscore( quality_first[j] );
		}

	    }

            for(int j=0; j<nKmer; j++) {
              string skmer_first = seq_first.substr(j, data->kmer_length);
	      string qkmer_first = quality_first.substr(j, data->kmer_length);
              if( contains_non_gatc(skmer_first) && check_kmer_quality(qkmer_first, data->bQscore) ) {

		kmer_int_type_t kmer_val = kmer_to_intval(skmer_first);

                if(headseq.strand > 0) irke_first.add_kmer(kmer_val , 1);
		else 		       irke_second.add_kmer(kmer_val , 1);

              }
            }
          }

          headseq = headseq_pair[1];
          if((headseq.maq >= data->MAQ) && (headseq.length == data->range)) {

            maq_second += headseq.maq;

            seq_second.assign( headseq.seq, data->range );
            quality_second.assign( headseq.quality, data->range );

            for(int j=0; j<data->range; j++) {
 
                if(headseq.strand > 0) {
                        if(get_Qscore(quality_second[j]) > cquality_first[j])  cquality_first[j]  = get_Qscore( quality_second[j] );
                } else {
                        if(get_Qscore(quality_second[j]) > cquality_second[j]) cquality_second[j] = get_Qscore( quality_second[j] );
                } 

           }

            for(int j=0; j<nKmer; j++) {
              string skmer_second = seq_second.substr(j, data->kmer_length);
	      string qkmer_second = quality_second.substr(j, data->kmer_length);
              if( contains_non_gatc(skmer_second) && check_kmer_quality(qkmer_second, data->bQscore) ) {

                kmer_int_type_t kmer_val = kmer_to_intval(skmer_second);

                if(headseq.strand > 0) irke_first.add_kmer(kmer_val , 1);
                else                   irke_second.add_kmer(kmer_val , 1);

              }
            }
         }

        value += valuebytes[i];
      }

  }

  END_BLOCK_LOOP

  if(nseq >= data->min_csize) { 

      maq_first  = (int)(maq_first / nseq);
      maq_second = (int)(maq_second / nseq);

      for(int j=0; j<data->range; j++) { quality_first[j] = get_QualityAscii( cquality_first[j] ); quality_second[j] = get_QualityAscii( cquality_second[j] ); }

      irke_first.prune_some_kmers(data->min_kmer_count, data->min_any_entropy, data->prune_error_kmers, data->min_ratio_non_error);
      irke_first.populate_sorted_kmers_list();

      string csequence_first = irke_first.compute_sequence_assemblies( data->MIN_CONNECTIVITY_RATIO, data->MIN_ASSEMBLY_LENGTH, data->MIN_ASSEMBLY_COVERAGE );

      irke_second.prune_some_kmers(data->min_kmer_count, data->min_any_entropy, data->prune_error_kmers, data->min_ratio_non_error);
      irke_second.populate_sorted_kmers_list();

      string csequence_second = irke_second.compute_sequence_assemblies( data->MIN_CONNECTIVITY_RATIO, data->MIN_ASSEMBLY_LENGTH, data->MIN_ASSEMBLY_COVERAGE );

      if((int)(csequence_first.length() == data->range) && ((int)csequence_second.length() == data->range) ) {

          HeadSeq *out_headseq = (HeadSeq *) malloc( sizeof(HeadSeq) * 2 );
	  memcpy( out_headseq[0].seq, csequence_first.c_str(), sizeof(char)*data->range );
          memcpy( out_headseq[1].seq, csequence_second.c_str(), sizeof(char)*data->range );

          memcpy( out_headseq[0].quality, quality_first.c_str(), sizeof(char)*data->range );
          memcpy( out_headseq[1].quality, quality_second.c_str(), sizeof(char)*data->range );	  

          out_headseq[0].strand = 1;
          out_headseq[1].strand = -1;

          out_headseq[0].PE = 1;
          out_headseq[1].PE = 2;

          out_headseq[0].maq = maq_first;
          out_headseq[1].maq = maq_second;

          out_headseq[0].length = data->range;
          out_headseq[1].length = data->range;

          out_headseq[0].type = -1 * nvalues;
          out_headseq[1].type = -1 * nvalues;

	  kv->add((char *) &plink,sizeof(PairLink),(char *) out_headseq,sizeof(HeadSeq)*2 );
	  free(out_headseq);
      }

  }

  free(cquality_first);
  free(cquality_second);

 }

}

void reduce_ekmer_for_singleton(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr)
{
 char *value;
 bool onEkmer = false;

 uint64_t nvalues_total;
 CHECK_FOR_BLOCKS(multivalue,valuebytes,nvalues,nvalues_total)
 BEGIN_BLOCK_LOOP(multivalue,valuebytes,nvalues) 

  value = multivalue;
  for(int i=0; i<nvalues; i++) {
    if(valuebytes[i] == sizeof(kmer_int_type_t)) {
	onEkmer = true; break;
    }	
    value += valuebytes[i];
  }

  if(onEkmer) {
    int num = 0; 
    value = multivalue;
    for(int i=0; i<nvalues; i++) {
       if(valuebytes[i] != sizeof(kmer_int_type_t)) num++; 
       value += valuebytes[i];
    }

    if(num>0) kv->add(key, keybytes, key, keybytes);
  }

 END_BLOCK_LOOP

}

void reduce_ekmer_for_HighDup(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr)
{
 char *value;
 bool onEkmer = true;

 uint64_t nvalues_total;
 CHECK_FOR_BLOCKS(multivalue,valuebytes,nvalues,nvalues_total)
 BEGIN_BLOCK_LOOP(multivalue,valuebytes,nvalues)

  value = multivalue;
  for(int i=0; i<nvalues; i++) {
    if(valuebytes[i] == sizeof(kmer_int_type_t)) {
        onEkmer = false; break;
    }
    value += valuebytes[i];
  }

  if(onEkmer) {
    int num = 0;
    value = multivalue;
    for(int i=0; i<nvalues; i++) {
       if(valuebytes[i] != sizeof(kmer_int_type_t)) num++;
       value += valuebytes[i];
    }

    if(num > 0) kv->add(key, keybytes, key, keybytes);
  }

 END_BLOCK_LOOP

}

void reduce_ekmers_HiDup(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr)
{
 Data *data = (Data *) ptr;
 char *value;
 int nKmer = data->range - data->kmer_length + 1;
 string seq1, seq2;
 Kmer_counter_map _kmer_counter;

 uint64_t nvalues_total;
 CHECK_FOR_BLOCKS(multivalue,valuebytes,nvalues,nvalues_total)
 BEGIN_BLOCK_LOOP(multivalue,valuebytes,nvalues)

     value = multivalue;
     for(int i=0; i<nvalues; i++) {
      HeadSeq *read_pair = (HeadSeq *) value;
      seq1.assign( read_pair[0].seq, data->range );
      seq2.assign( read_pair[1].seq, data->range );
      for(int j=0; j<nKmer; j++) {
              string seq_kmer1 = seq1.substr(j, data->kmer_length);
              string seq_kmer2 = seq2.substr(j, data->kmer_length);
              if(contains_non_gatc(seq_kmer1)) {
                kmer_int_type_t kmer_val = kmer_to_intval(seq_kmer1);
                _kmer_counter[kmer_val]+= 1;
              }
              if(contains_non_gatc(seq_kmer2)) {
                kmer_int_type_t kmer_val = kmer_to_intval(seq_kmer2);
                _kmer_counter[kmer_val]+= 1;
              }
      }
      value += valuebytes[i];
     }

     Kmer_counter_map_iterator it;
     for (it = _kmer_counter.begin(); it != _kmer_counter.end(); it++) {
        kmer_int_type_t kmer_val = it->first;
        int dum = 1;
        kv->add((char *) &kmer_val, sizeof(kmer_int_type_t), (char *) &dum, sizeof(int));
     }

 END_BLOCK_LOOP

}

void map_kmer_from_singleton(uint64_t itask, char *key, int keybytes, char *value, int valuebytes, KeyValue *kv, void *ptr)
{
   PairLink *plink = (PairLink *) key;
   Data *data = (Data *) ptr;
   HeadSeq *head_seq = (HeadSeq *) value;
   string seq1, seq2;
   Kmer_counter_map _kmer_counter;

   seq1.assign( head_seq[0].seq, data->range );
   seq2.assign( head_seq[1].seq, data->range );

   int nKmer = data->range - data->kmer_length + 1;

   for(int j=0; j<nKmer; j++) {
              string seq_kmer1 = seq1.substr(j, data->kmer_length);
              string seq_kmer2 = seq2.substr(j, data->kmer_length);
              if(contains_non_gatc(seq_kmer1)) {
                kmer_int_type_t kmer_val = kmer_to_intval(seq_kmer1);
                _kmer_counter[kmer_val]+= 1;
              }
              if(contains_non_gatc(seq_kmer2)) {
                kmer_int_type_t kmer_val = kmer_to_intval(seq_kmer2);
                _kmer_counter[kmer_val]+= 1;
              }
   }

   Kmer_counter_map_iterator it;
   for (it = _kmer_counter.begin(); it != _kmer_counter.end(); it++) {
        kmer_int_type_t kmer_val = it->first;
        kv->add((char *) &kmer_val, sizeof(kmer_int_type_t), (char *) &kmer_val, sizeof(kmer_int_type_t));
   }

}


void reduce_extract_kmer_singleton(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr)
{
 Data *data = (Data *) ptr;
 char *value;
 //bool onSingleton = true;
 int nKmer = data->range - data->kmer_length + 1;
 string seq1, seq2;
 Kmer_counter_map _kmer_counter;

 uint64_t nvalues_total;
 CHECK_FOR_BLOCKS(multivalue,valuebytes,nvalues,nvalues_total)
 BEGIN_BLOCK_LOOP(multivalue,valuebytes,nvalues)

//  value = multivalue;
//  for(int i=0; i<nvalues; i++) {
//      HeadSeq *read_pair = (HeadSeq *) value;
//      if((read_pair[0].type < 0) && (read_pair[1].type < 0)) {
//         onSingleton = false;
//         break;
//      }
//      value += valuebytes[i];
//  }

//  if(onSingleton) {
     value = multivalue;
     for(int i=0; i<nvalues; i++) {
      HeadSeq *read_pair = (HeadSeq *) value;
      seq1.assign( read_pair[0].seq, data->range );
      seq2.assign( read_pair[1].seq, data->range );
      for(int j=0; j<nKmer; j++) {
              string seq_kmer1 = seq1.substr(j, data->kmer_length);	
	      string seq_kmer2 = seq2.substr(j, data->kmer_length);
	      if(contains_non_gatc(seq_kmer1)) {
	     	kmer_int_type_t kmer_val = kmer_to_intval(seq_kmer1);
		_kmer_counter[kmer_val]+= 1;
	      }
	      if(contains_non_gatc(seq_kmer2)) {
                kmer_int_type_t kmer_val = kmer_to_intval(seq_kmer2);
		_kmer_counter[kmer_val]+= 1;
              }
      }
      value += valuebytes[i];
     }

     Kmer_counter_map_iterator it;
     for (it = _kmer_counter.begin(); it != _kmer_counter.end(); it++) {
        kmer_int_type_t kmer_val = it->first;
	kv->add((char *) &kmer_val, sizeof(kmer_int_type_t), (char *) &it->second, sizeof(int));
     }
//  }

 END_BLOCK_LOOP

}


void reduce_sort_singleton(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr)
{
 Data *data = (Data *) ptr;
 char *value;
 bool onSingleton = true;

 uint64_t nvalues_total;
 CHECK_FOR_BLOCKS(multivalue,valuebytes,nvalues,nvalues_total)
 BEGIN_BLOCK_LOOP(multivalue,valuebytes,nvalues)

 value = multivalue;
 for(int i=0; i<nvalues; i++) {
      HeadSeq *read_pair = (HeadSeq *) value;
      if((read_pair[0].type < 0) && (read_pair[1].type < 0)) {
	 onSingleton = false;
	 break;
      }
      value += valuebytes[i];
 }

 if(onSingleton) {
     value = multivalue;
     for(int i=0; i<nvalues; i++) {
      HeadSeq *read_pair = (HeadSeq *) value;
      if((read_pair[0].type == 1) && (read_pair[1].type == 1)) kv->add(key, keybytes, value, valuebytes[i]);
      value += valuebytes[i];
     }
 }

 END_BLOCK_LOOP

}


void reduce_consensus_SeqQual_highDup(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr)
{

   PairLink *plink = (PairLink *) key;

  if((nvalues >= 3) && (plink->chr1 > 0) && (plink->chr2 > 0)) {

   Data *data = (Data *) ptr;
   char *value;
   string seq1, seq2, quality1, quality2;
   int PE1, PE2, strand1, strand2;

   IRKE irke_1(data->kmer_length, data->MAX_RECURSION, data->MIN_SEED_ENTROPY, data->MIN_SEED_COVERAGE,
                data->min_any_entropy, data->PACMAN, data->CRAWL, data->crawl_length, data->DS);

   IRKE irke_2(data->kmer_length, data->MAX_RECURSION, data->MIN_SEED_ENTROPY, data->MIN_SEED_COVERAGE,
                data->min_any_entropy, data->PACMAN, data->CRAWL, data->crawl_length, data->DS);

   IRKE irke(data->kmer_length, data->MAX_RECURSION, data->MIN_SEED_ENTROPY, data->MIN_SEED_COVERAGE,
                data->min_any_entropy, data->PACMAN, data->CRAWL, data->crawl_length, data->DS);

   StripedSmithWaterman::Aligner aligner(1,4,6,1);
   StripedSmithWaterman::Filter filter;
   StripedSmithWaterman::Alignment alignment;

   int *cquality1 = (int *) malloc( sizeof(int)*data->range );
   int *cquality2 = (int *) malloc( sizeof(int)*data->range );

   for(int i=0; i<data->range; i++) { cquality1[i]=0; cquality2[i]=0; }

   uint64_t nvalues_total;
   CHECK_FOR_BLOCKS(multivalue,valuebytes,nvalues,nvalues_total)
   BEGIN_BLOCK_LOOP(multivalue,valuebytes,nvalues)

    plink->maq1 = 0;
    plink->maq2 = 0;

    int nseq1 = 0;
    int nseq2 = 0;
    int nKmer = data->range - data->kmer_length + 1;

      value = multivalue;
      for(int i=0; i<nvalues; i++) {
          HeadSeq *headseq_pair = (HeadSeq *) value;

          HeadSeq headseq = headseq_pair[0];
	  if(headseq.maq >= data->MAQ) {

	    plink->maq1 += headseq.maq;

            seq1.assign( headseq.seq, data->range );
            quality1.assign( headseq.quality, data->range );
            for(int j=0; j<data->range; j++) cquality1[j] += get_Qscore( quality1[j] );

            for(int j=0; j<nKmer; j++) {
              string seq_kmer1 = seq1.substr(j, data->kmer_length);
              if( contains_non_gatc(seq_kmer1) ) {
                kmer_int_type_t kmer_val = kmer_to_intval(seq_kmer1);
                irke_1.add_kmer(kmer_val , 1);
                irke.add_kmer(kmer_val , 1);
              }
            }
            nseq1++;
          }


          headseq = headseq_pair[1];
	  if(headseq.maq >= data->MAQ) {

	    plink->maq2 += headseq.maq;

            seq2.assign( headseq.seq, data->range );
            quality2.assign( headseq.quality, data->range );
            for(int j=0; j<data->range; j++) cquality2[j] += get_Qscore( quality2[j] );

            for(int j=0; j<nKmer; j++) {
              string seq_kmer2 = seq2.substr(j, data->kmer_length);
              if( contains_non_gatc(seq_kmer2) ) {
                kmer_int_type_t kmer_val = kmer_to_intval(seq_kmer2);
                irke_2.add_kmer(kmer_val , 1);
	        irke.add_kmer(kmer_val , 1);
              }
            }
            nseq2++;
         }

        value += valuebytes[i];
      }

if((nseq1>0) && (nseq2>0)) {

      plink->maq1 = (int)(plink->maq1 / nseq1);
      plink->maq2 = (int)(plink->maq2 / nseq2);

      for(int j=0; j<data->range; j++) { cquality1[j]=(int)(cquality1[j]/nseq1); cquality2[j]=(int)(cquality2[j]/nseq2); }
      for(int j=0; j<data->range; j++) { quality1[j] = get_QualityAscii( cquality1[j] ); quality2[j] = get_QualityAscii( cquality2[j] ); }

      irke_1.prune_some_kmers(data->min_kmer_count, data->min_any_entropy, data->prune_error_kmers, data->min_ratio_non_error);
      irke_1.populate_sorted_kmers_list();

      string csequence_1 = irke_1.compute_sequence_assemblies( data->MIN_CONNECTIVITY_RATIO, data->MIN_ASSEMBLY_LENGTH, data->MIN_ASSEMBLY_COVERAGE );

      irke_2.prune_some_kmers(data->min_kmer_count, data->min_any_entropy, data->prune_error_kmers, data->min_ratio_non_error);
      irke_2.populate_sorted_kmers_list();

      string csequence_2 = irke_2.compute_sequence_assemblies( data->MIN_CONNECTIVITY_RATIO, data->MIN_ASSEMBLY_LENGTH, data->MIN_ASSEMBLY_COVERAGE );

      alignment.Clear();
      if(plink->loc1 <= plink->loc2) aligner.Align(csequence_2.c_str(), csequence_1.c_str(), csequence_1.length(), filter, &alignment);
      else                           aligner.Align(csequence_1.c_str(), csequence_2.c_str(), csequence_2.length(), filter, &alignment);

      int overlap = 0;
      int32_t *bam_cigar = (int32_t *) malloc( sizeof(int32_t)*3 );
      if( alignment.cigar.size() > 0) {
                convert_cigar2bam(alignment.cigar[0], bam_cigar);
                if( (bam_cigar[0] == 0) || (bam_cigar[0] == 7) || (bam_cigar[0] == 8) ) overlap += (int)bam_cigar[1];
      }
      free(bam_cigar);

      irke.prune_some_kmers(data->min_kmer_count, data->min_any_entropy, data->prune_error_kmers, data->min_ratio_non_error);
      irke.populate_sorted_kmers_list();

      string csequence =
            irke.compute_sequence_assemblies( data->MIN_CONNECTIVITY_RATIO, data->MIN_ASSEMBLY_LENGTH, data->MIN_ASSEMBLY_COVERAGE );

      int len1, len2;
      if(csequence_1.length() > data->range) len1 = data->range;
      else                                   len1 = csequence_1.length();
      if(csequence_2.length() > data->range) len2 = data->range;
      else                                   len2 = csequence_2.length();

      if( (overlap >= data->OVERLAP) && (csequence.length() >= csequence_1.length()) && (csequence.length() >= csequence_2.length()) ) {

             seq1.assign( csequence.substr(0,len1).c_str(), len1 );
             seq2.assign( csequence.substr(csequence.length() - len2,len2).c_str(), len2 );

             quality1.assign( quality1.substr(0,len1).c_str(), len1 );
             quality2.assign( quality2.substr(quality2.length() - len2,len2).c_str(), len2 );

             HeadSeq seq_quality1, seq_quality2;

             memcpy( seq_quality1.seq, seq1.c_str(), sizeof(char)*len1 );
             memcpy( seq_quality2.seq, seq2.c_str(), sizeof(char)*len2 );

             memcpy( seq_quality1.quality, quality1.c_str(), sizeof(char)*len1 );
             memcpy( seq_quality2.quality, quality2.c_str(), sizeof(char)*len2 );

             seq_quality1.strand = strand1;
             seq_quality2.strand = strand2;

	     seq_quality1.PE = PE1;
             seq_quality2.PE = PE2;

             seq_quality1.length = len1;
             seq_quality2.length = len2;

             HeadSeq *out_headseq = (HeadSeq *) malloc( sizeof(HeadSeq) * 2 );
	     out_headseq[0] = seq_quality1;
             out_headseq[1] = seq_quality2;

             kv->add((char *) plink, sizeof(PairLink), (char *) out_headseq, sizeof(HeadSeq)*2 );
	     free(out_headseq);

      } 

      if( (overlap < data->OVERLAP) && (csequence_1.length() <= data->range) && (csequence_2.length() <= data->range) ) {

          HeadSeq seq_quality1, seq_quality2;

          seq1.assign( csequence_1.substr(0,len1).c_str(), len1 );
          seq2.assign( csequence_2.substr(csequence_2.length() - len2,len2).c_str(), len2 );

          quality1.assign( quality1.substr(0,len1).c_str(), len1 );
          quality2.assign( quality2.substr(quality2.length() - len2,len2).c_str(), len2 );

          memcpy( seq_quality1.seq, csequence_1.c_str(), sizeof(char)*len1 );
          memcpy( seq_quality2.seq, csequence_2.c_str(), sizeof(char)*len2 );

          memcpy( seq_quality1.quality, quality1.c_str(), sizeof(char)*len1 );
          memcpy( seq_quality2.quality, quality2.c_str(), sizeof(char)*len2 );

          seq_quality1.strand = strand1;
          seq_quality2.strand = strand2;

	  seq_quality1.PE = PE1;
          seq_quality2.PE = PE2;

          seq_quality1.length = len1;
          seq_quality2.length = len2;


          HeadSeq *out_headseq = (HeadSeq *) malloc( sizeof(HeadSeq) * 2 );
          out_headseq[0] = seq_quality1;
          out_headseq[1] = seq_quality2;

	  plink->ID = nvalues;
          kv->add((char *) plink, sizeof(PairLink), (char *) out_headseq, sizeof(HeadSeq)*2 );
          free(out_headseq);
     }

}

   END_BLOCK_LOOP

   free(cquality1);
   free(cquality2); 

  }
}

void map_correct_singleton_step1(uint64_t itask, char *key, int keybytes, char *value, int valuebytes, KeyValue *kv, void *ptr)
{
   PairLink *plink = (PairLink *) key;
   Data *data = (Data *) ptr;
   HeadSeq *head_seq = (HeadSeq *) value;
   HeadSeq read1 = head_seq[0];
   HeadSeq read2 = head_seq[1];

   Kmer_counter_map _kmer_map;

//if((plink->chr1 == data->offside) && (plink->chr2 == data->offside)) {

   string seq;
   int nKmer = read1.length - data->kmer_length + 1; 
   seq.assign( read1.seq, read1.length );
   for(int i=0; i<nKmer; i++) {
              string seq_kmer = seq.substr(i, data->kmer_length);
              if( contains_non_gatc(seq_kmer) ) {
                  kmer_int_type_t kmer_val = kmer_to_intval(seq_kmer);
		  _kmer_map[kmer_val] += 1;
//		  kv->add((char *) &kmer_val, sizeof(kmer_int_type_t), key, keybytes);
	      }
   }

   nKmer = read2.length - data->kmer_length + 1;
   seq.assign( read2.seq, read2.length );
   for(int i=0; i<nKmer; i++) {
              string seq_kmer = seq.substr(i, data->kmer_length);
              if( contains_non_gatc(seq_kmer) ) {
                  kmer_int_type_t kmer_val = kmer_to_intval(seq_kmer);
		  _kmer_map[kmer_val] += 1;
//                  kv->add((char *) &kmer_val, sizeof(kmer_int_type_t), key, keybytes);
              }
   }

   for(Kmer_counter_map::iterator it = _kmer_map.begin(); it != _kmer_map.end(); ++it) {
        kv->add((char *) &it->first, sizeof(kmer_int_type_t), key, keybytes);
   }

//}

}

void reduce_correct_single_step2(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr)
{
   Data *data = (Data *) ptr;
   char *value;
   uint64_t nvalues_total;

   bool onError = false;
   int num = 0;

   CHECK_FOR_BLOCKS(multivalue,valuebytes,nvalues,nvalues_total)
   BEGIN_BLOCK_LOOP(multivalue,valuebytes,nvalues)    

   for(int i=0; i<nvalues; i++) {
	if(valuebytes[i] == sizeof(kmer_int_type_t)) onError = true;
	if(valuebytes[i] == sizeof(PairLink))        num++;
   }
    
   if(onError && (num >= data->nKmer)) {
     value = multivalue;
     for(int i=0; i<nvalues; i++) {
       if(valuebytes[i] == sizeof(PairLink)) kv->add(value, valuebytes[i], key, keybytes); 
       value += valuebytes[i];
     }   
   }
 
   END_BLOCK_LOOP

}

void reduce_unique_kmers(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr)
{
   Data *data = (Data *) ptr;
   char *value;
   int *buffer_value = (int *)malloc( sizeof(int)*2 );

   buffer_value[0] = 0;
   buffer_value[1] = 1;

   if(data->ftype == 1)      buffer_value[0] = -1*nvalues;
   else if(data->ftype == 2) buffer_value[0] =  1*nvalues;

   kv->add(key, keybytes, (char *) buffer_value, sizeof(int)*2 );
   free(buffer_value);
}


void reduce_ekmer_singleton(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr)
{

  PairLink *plink = (PairLink *) key;
  Data *data = (Data *) ptr;
  char *value;
  string seq1, seq2, quality1, quality2;

  int num_ekmer = 0;
  Kmer_counter_map _ekmer_map;

  StripedSmithWaterman::Aligner aligner(1,4,6,1);
  StripedSmithWaterman::Filter filter;
  StripedSmithWaterman::Alignment alignment;

  uint64_t nvalues_total;
  CHECK_FOR_BLOCKS(multivalue,valuebytes,nvalues,nvalues_total)
  BEGIN_BLOCK_LOOP(multivalue,valuebytes,nvalues)

  if(nvalues > 1) {

     value = multivalue;
     for(int i=0; i<nvalues; i++) {

       if((valuebytes[i] == sizeof(kmer_int_type_t)) && (valuebytes[i] != sizeof(HeadSeq)*2)) {
           kmer_int_type_t kmer_val = *(kmer_int_type_t *) value;
           _ekmer_map[kmer_val] += 1;
           num_ekmer++;
       }

        value += valuebytes[i];
     }

     value = multivalue;
     for(int i=0; i<nvalues; i++) {

        if((valuebytes[i] == sizeof(HeadSeq)*2) && (valuebytes[i] != sizeof(kmer_int_type_t))) {

            HeadSeq *read_pair = (HeadSeq *) value;

            seq1.assign(     read_pair[0].seq,     read_pair[0].length );
            quality1.assign( read_pair[0].quality, read_pair[0].length );

            seq2.assign(     read_pair[1].seq,     read_pair[1].length );
            quality2.assign( read_pair[1].quality, read_pair[1].length );


            vector<int> loc1;
            vector<int> loc2;
            int ek1 = 0;
            int ek2 = 0;
            int nkmer = read_pair[0].length - data->kmer_length + 1;
            for(int j=0; j<nkmer; j++) {
              if(contains_non_gatc(seq1.substr(j,data->kmer_length)) ) {
                kmer_int_type_t kmer_val  = kmer_to_intval(  seq1.substr(j,data->kmer_length) );
                if( (_ekmer_map.find(kmer_val) != _ekmer_map.end()) ) {
                        ek1++;
                        loc1.push_back(j);
                }
             }
            }

            nkmer = read_pair[1].length - data->kmer_length + 1;
            for(int j=0; j<nkmer; j++) {
             if( contains_non_gatc(seq2.substr(j,data->kmer_length)) ) {
               kmer_int_type_t kmer_val  = kmer_to_intval(  seq2.substr(j,data->kmer_length) );
               if( _ekmer_map.find(kmer_val) != _ekmer_map.end() ) {
                        ek2++;
                        loc2.push_back(j);
               }
             }
            }

	    set<kmer_int_type_t> eKmer_list;

	    if(ek1 > 0) {

                set<short> eloc1;
                vector<short> chop_region;
                vector< vector<int> > cc_loc1;
                cc_loc1.push_back(vector<int>());
                int cc=0;
                cc_loc1[cc].push_back(loc1[0]);

                for(int j=1; j<loc1.size(); j++) {
                   if(loc1[j-1]+data->kmer_length-1 < loc1[j]) {
                      cc++;
                      cc_loc1.push_back(vector<int>());
                   }
                   cc_loc1[cc].push_back(loc1[j]);
                }

                for(cc = 0; cc < cc_loc1.size(); cc++) {
                   chop_region.clear();
                   chop_region = error_region_extract(cc_loc1[cc], data->kmer_length, read_pair[0].length);

                   for(int k=0; k<chop_region.size(); k++) eloc1.insert( chop_region[k] );
                }

                for(set<short>::iterator it=eloc1.begin(); it!=eloc1.end(); it++) {
                   int obs_base = get_int_from_base( seq1[*it] );
                   int qscore   = get_Qscore( quality1[*it] );

                   if(qscore >= data->bQscore) {

			short start_loc = *it - (short)data->kmer_length + 1;
			if(start_loc < 0) start_loc = 0;
			for(short j=start_loc; j <= *it; j++) {
			   kmer_int_type_t kmer_val = kmer_to_intval(seq1.substr(j,data->kmer_length));;
			   eKmer_list.insert( kmer_val );	   
			}		

		   }

		}

	    }

	    if(ek2 > 0) {

                set<short> eloc2;
                vector<short> chop_region;
                vector< vector<int> > cc_loc2;
                cc_loc2.push_back(vector<int>());
                int cc=0;
                cc_loc2[cc].push_back(loc2[0]);

                for(int j=1; j<loc2.size(); j++) {
                   if(loc2[j-1]+data->kmer_length-1 < loc2[j]) {
                      cc++;
                      cc_loc2.push_back(vector<int>());
                   }
                   cc_loc2[cc].push_back(loc2[j]);
                }

                for(cc = 0; cc < cc_loc2.size(); cc++) {
                   chop_region.clear();
                   chop_region = error_region_extract(cc_loc2[cc], data->kmer_length, read_pair[1].length);

                  for(int k=0; k<chop_region.size(); k++) eloc2.insert( chop_region[k] );
                }

                for(set<short>::iterator it=eloc2.begin(); it!=eloc2.end(); it++) {
                   int obs_base = get_int_from_base( seq2[*it] );
                   int qscore    = get_Qscore( quality2[*it] );

                   if(qscore >= data->bQscore) {

                        short start_loc = *it - (short)data->kmer_length + 1;
                        if(start_loc < 0) start_loc = 0;
                        for(short j=start_loc; j <= *it; j++) {
                           kmer_int_type_t kmer_val = kmer_to_intval(seq2.substr(j,data->kmer_length));;
                           eKmer_list.insert( kmer_val );
                        }

		   }

	        }

	    }
	
	    for(set<kmer_int_type_t>::iterator it=eKmer_list.begin(); it!=eKmer_list.end(); it++) {
		kmer_int_type_t kmer_val = *it;
                kv->add((char *) &kmer_val, sizeof(kmer_int_type_t), (char *) &(*it), sizeof(kmer_int_type_t));	
	    }

        }
        value += valuebytes[i];

    }


  }

  END_BLOCK_LOOP

}


void reduce_correct_single_step3(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr)
{

  PairLink plink = *(PairLink *) key;
  Data *data = (Data *) ptr;
  char *value;
  string seq1, seq2, quality1, quality2;


  bool check_len = false;
  int start, len;
  if(plink.loc1 <= plink.loc2) {
        start = plink.loc1;
        len = plink.loc2 - plink.loc1 + 1;
        if( (plink.loc1+2*len) > data->reference_sequence[plink.chr1].length() ) {
                len = data->reference_sequence[plink.chr1].length() - plink.loc1 + 1;
        } else  {
                 check_len = true; len = 2*len;
        }
  } else  {
        start = plink.loc2;
        len = plink.loc1 - plink.loc2 + 1;
        if( (plink.loc2+2*len) > data->reference_sequence[plink.chr1].length() ) {
                len = data->reference_sequence[plink.chr1].length() - plink.loc2 + 1;
        } else {
                check_len = true; len = 2*len;
        }
  }


  int num_ekmer = 0;
  Kmer_counter_map _ekmer_map;

  StripedSmithWaterman::Aligner aligner(1,4,6,1);
  StripedSmithWaterman::Filter filter;
  StripedSmithWaterman::Alignment alignment;

  uint64_t nvalues_total;
  CHECK_FOR_BLOCKS(multivalue,valuebytes,nvalues,nvalues_total)
  BEGIN_BLOCK_LOOP(multivalue,valuebytes,nvalues)

  if(nvalues > 1) {

     value = multivalue;
     for(int i=0; i<nvalues; i++) {

       if((valuebytes[i] == sizeof(kmer_int_type_t)) && (valuebytes[i] != sizeof(HeadSeq)*2)) {
	   kmer_int_type_t kmer_val = *(kmer_int_type_t *) value;
	   _ekmer_map[kmer_val] += 1;
	   num_ekmer++;
       }

	value += valuebytes[i];
     } 

     value = multivalue;
     for(int i=0; i<nvalues; i++) {

        if((valuebytes[i] == sizeof(HeadSeq)*2) && (valuebytes[i] != sizeof(kmer_int_type_t))) {

            HeadSeq *read_pair = (HeadSeq *) value;

            seq1.assign(     read_pair[0].seq,     read_pair[0].length );
            quality1.assign( read_pair[0].quality, read_pair[0].length );

            seq2.assign(     read_pair[1].seq,     read_pair[1].length );
            quality2.assign( read_pair[1].quality, read_pair[1].length );

	    for(int j=0; j<read_pair[0].length; j++) {
		    int bid = get_int_from_base( seq1[j] );
                    if((bid>=0) && (bid<=3)) {
		        if(read_pair[0].strand > 0) data->sBase1[bid]++;
			else			    data->sBase2[bid]++;
		    }
	    }	

            for(int j=0; j<read_pair[1].length; j++) {
                    int bid = get_int_from_base( seq2[j] );
                    if((bid>=0) && (bid<=3)) {
                        if(read_pair[1].strand > 0) data->sBase1[bid]++;
                        else                        data->sBase2[bid]++;
		    } 
            }


            vector<int> loc1;
            vector<int> loc2;
            int ek1 = 0;
            int ek2 = 0;
            int nkmer = read_pair[0].length - data->kmer_length + 1;
            for(int j=0; j<nkmer; j++) {
	      if(contains_non_gatc(seq1.substr(j,data->kmer_length)) ) {
                kmer_int_type_t kmer_val  = kmer_to_intval(  seq1.substr(j,data->kmer_length) );
                if( (_ekmer_map.find(kmer_val) != _ekmer_map.end()) ) {
                        ek1++;
                        loc1.push_back(j);
                }
	     }
            }

            nkmer = read_pair[1].length - data->kmer_length + 1;
            for(int j=0; j<nkmer; j++) {
	     if( contains_non_gatc(seq2.substr(j,data->kmer_length)) ) { 
               kmer_int_type_t kmer_val  = kmer_to_intval(  seq2.substr(j,data->kmer_length) );
               if( _ekmer_map.find(kmer_val) != _ekmer_map.end() ) {
                        ek2++;
                        loc2.push_back(j);
               }
	     }
            }


	    string ref_sequence;
	    if(((ek1 > 0) || (ek2 > 0)) && check_len && (start > 1)) {
		ref_sequence = data->reference_sequence[plink.chr1].substr(start-1,len);
	    }

	    if(ek1 > 0) {

                set<short> eloc1;
		vector<short> chop_region;
                vector< vector<int> > cc_loc1;
                cc_loc1.push_back(vector<int>());
                int cc=0;
                cc_loc1[cc].push_back(loc1[0]);

                for(int j=1; j<loc1.size(); j++) {
                   if(loc1[j-1]+data->kmer_length-1 < loc1[j]) {
                      cc++;
                      cc_loc1.push_back(vector<int>());
                   }
                   cc_loc1[cc].push_back(loc1[j]);
                }

                for(cc = 0; cc < cc_loc1.size(); cc++) {
		   chop_region.clear();
		   chop_region = error_region_extract(cc_loc1[cc], data->kmer_length, read_pair[0].length);
                   for(int k=0; k<chop_region.size(); k++) eloc1.insert( chop_region[k] );
                }
		
		map<int,int> index1;
		int sp, sp_loc;
		if(check_len && (start > 1)) {
                  alignment.Clear();
                  aligner.Align(seq1.substr(0,read_pair[0].length).c_str(), ref_sequence.c_str(), ref_sequence.length(), filter, &alignment);
		  int sp     = (int)alignment.query_begin;
                  int sp_loc = (int)alignment.ref_begin;
		  int n = 0;
		  for(int i=0; i<sp; i++) { index1[n] = start + sp_loc - sp + i; n++; }

        	     int32_t *bam_cigar = (int32_t *) malloc( sizeof(int32_t)*2 );
                     for(size_t j = 0; j < alignment.cigar.size(); j++)  {
          		convert_cigar2bam(alignment.cigar[j], bam_cigar);

         		if( (bam_cigar[0] == 0) || (bam_cigar[0] == 3) || (bam_cigar[0] == 7) || (bam_cigar[0] == 8) ) {

                	  for(int k=0; k<(int)bam_cigar[1]; k++) {
                   		index1[n] = start + sp_loc + k;
				n++;	
			  }
                            sp     += (int)bam_cigar[1];
                            sp_loc += (int)bam_cigar[1];
                        } else if( (bam_cigar[0] == 2) || (bam_cigar[0] == 6) )  {
                	    sp_loc += (int)bam_cigar[1];
         	        } else if( bam_cigar[0] == 1) {

			   for(int k=0; k<(int)bam_cigar[1]; k++) {
                                index1[n] = start + sp_loc;
                                n++;
                           } 

			   sp += (int)bam_cigar[1];
			}

		    } 
		    free(bam_cigar);

		}

                for(set<short>::iterator it=eloc1.begin(); it!=eloc1.end(); it++) {
                   int obs_base = get_int_from_base( seq1[*it] );
                   int qscore   = get_Qscore( quality1[*it] );
		   bool check_highDup = true;


		   if(check_len && (start > 1) ) {
                     int loc = index1[*it];
                     unsigned long long point_target = (unsigned long long)plink.chr1 * (unsigned long long)1e11 + (unsigned long long)loc;
		     if(data->_loc_freq_table.count(point_target) > 0) {
			unsigned long long tloc  = data->_loc_freq_table[point_target];	
                        unsigned long long nbase = data->table_base[ tloc*4 + (unsigned long long)obs_base ];
		        if(nbase > 0) check_highDup = false;

		     }
		   }


	           if((qscore >= data->bQscore) && check_highDup && data->On_Polishing) { 
		     char rq1 = get_QualityAscii( 5 );
		     quality1.replace(*it,1,1,rq1);
		     int bind = get_int_from_base( seq1[*it] );
		     if((bind>=0) &&(bind<=3)) {
			if(read_pair[0].strand > 0) data->eBase1[bind]++;	
			else			    data->eBase2[bind]++;
		     }
	          } 


	        }

	    }


            if(ek2 > 0) {

                set<short> eloc2;
		vector<short> chop_region;
                vector< vector<int> > cc_loc2;
                cc_loc2.push_back(vector<int>());
                int cc=0;
                cc_loc2[cc].push_back(loc2[0]);

                for(int j=1; j<loc2.size(); j++) {
                   if(loc2[j-1]+data->kmer_length-1 < loc2[j]) {
                      cc++;
                      cc_loc2.push_back(vector<int>());
                   }
                   cc_loc2[cc].push_back(loc2[j]);
                }

                for(cc = 0; cc < cc_loc2.size(); cc++) {
		   chop_region.clear();
		   chop_region = error_region_extract(cc_loc2[cc], data->kmer_length, read_pair[1].length);
                   for(int k=0; k<chop_region.size(); k++) eloc2.insert( chop_region[k] );
                }

                map<int,int> index2;
                int sp, sp_loc;
                if(check_len && (start > 1)) {
                  alignment.Clear();
                  aligner.Align(seq2.substr(0,read_pair[1].length).c_str(), ref_sequence.c_str(), ref_sequence.length(), filter, &alignment);
                  int sp     = (int)alignment.query_begin;
                  int sp_loc = (int)alignment.ref_begin;
                  int n = 0;
                  for(int i=0; i<sp; i++) { index2[n] = start + sp_loc - sp + i; n++; }

                     int32_t *bam_cigar = (int32_t *) malloc( sizeof(int32_t)*2 );
                     for(size_t j = 0; j < alignment.cigar.size(); j++)  {
                        convert_cigar2bam(alignment.cigar[j], bam_cigar);

                        if( (bam_cigar[0] == 0) || (bam_cigar[0] == 3) || (bam_cigar[0] == 7) || (bam_cigar[0] == 8) ) {

                            for(int k=0; k<(int)bam_cigar[1]; k++) {
                                index2[n] = start + sp_loc + k;
                                n++;
                            }
                            sp     += (int)bam_cigar[1];
                            sp_loc += (int)bam_cigar[1];
                        } else if( (bam_cigar[0] == 2) || (bam_cigar[0] == 6) )  {
                            sp_loc += (int)bam_cigar[1];
                        } else if( bam_cigar[0] == 1) {
                           for(int k=0; k<(int)bam_cigar[1]; k++) {
                                index2[n] = start + sp_loc;
                                n++;
                           }
                           sp     += (int)bam_cigar[1];
                        }

                    }
                    free(bam_cigar);

                }


                for(set<short>::iterator it=eloc2.begin(); it!=eloc2.end(); it++) {
                   int obs_base = get_int_from_base( seq2[*it] );
                   int qscore    = get_Qscore( quality2[*it] );
                   bool check_highDup = true;

                   if(check_len && (start > 1) ) {
                     int loc = index2[*it];
                     unsigned long long point_target = (unsigned long long)plink.chr2 * (unsigned long long)1e11 + (unsigned long long)loc;
		     if(data->_loc_freq_table.count(point_target) > 0) {
		       unsigned long long tloc  = data->_loc_freq_table[point_target];
                       unsigned long long nbase = data->table_base[ tloc*4 + (unsigned long long)obs_base ];
                       if(nbase > 0) check_highDup = false;
		     }
                   }


		   if((qscore >= data->bQscore) && check_highDup && data->On_Polishing) {
                     char rq2 = get_QualityAscii( 5 );
		     quality2.replace(*it,1,1,rq2); 
		     int bind = get_int_from_base( seq2[*it] );
                     if((bind>=0) &&(bind<=3)) {
                        if(read_pair[1].strand > 0) data->eBase1[bind]++;
                        else                        data->eBase2[bind]++;
		     }
	          }


               }

            }

            memcpy(read_pair[0].quality, quality1.c_str(), sizeof(char)*read_pair[0].length);
            memcpy(read_pair[1].quality, quality2.c_str(), sizeof(char)*read_pair[1].length);

            kv->add((char *) &plink, sizeof(PairLink), (char *) read_pair, sizeof(HeadSeq)*2);
	    break;

	}

        value += valuebytes[i];
     }

   } else if(nvalues == 1) {

      if(valuebytes[0] == sizeof(HeadSeq)*2) {
	value = multivalue;
	kv->add((char *) &plink, sizeof(PairLink), value, valuebytes[0]);

        HeadSeq *read_pair = (HeadSeq *) value;
        seq1.assign(     read_pair[0].seq,     read_pair[0].length );
        seq2.assign(     read_pair[1].seq,     read_pair[1].length );

        for(int j=0; j<read_pair[0].length; j++) {
                    int bid = get_int_from_base( seq1[j] );
                    if((bid>=0) && (bid<=3)) {
                        if(read_pair[0].strand > 0) data->sBase1[bid]++;
                        else                        data->sBase2[bid]++;
                    }
        }

        for(int j=0; j<read_pair[1].length; j++) {
                    int bid = get_int_from_base( seq2[j] );
                    if((bid>=0) && (bid<=3)) {
                        if(read_pair[1].strand > 0) data->sBase1[bid]++;
                        else                        data->sBase2[bid]++;
                    }
        }


     }

  }

  END_BLOCK_LOOP

}

void map_print_like_correct(uint64_t itask, char *key, int keybytes, char *value, int valuebytes, KeyValue *kv, void *ptr)
{
   Data *data = (Data *) ptr;
   int *qscore = (int *) key;
   float *like = (float *) value;
   data->outFile << qscore[0] << "\t" << qscore[1] << "\t" << qscore[2] << "\t" << like[0] << "\t" << like[1] << endl;

//   float *like = (float *) value;
//   data->outFile << qscore << "\t" << like[0] << "\t" << like[1] << "\t" << like[2] << endl; 	   
}


void reduce_consensus_singleton(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr)
{

   Data *data = (Data *) ptr;
   char *value; 

   int qsum = -1;
   int index = -1;

   uint64_t nvalues_total;
   CHECK_FOR_BLOCKS(multivalue,valuebytes,nvalues,nvalues_total)
   BEGIN_BLOCK_LOOP(multivalue,valuebytes,nvalues)

   if(nvalues > 1) {

       value = multivalue;
       for(int i=0; i<nvalues; i++) {
	  HeadSeq *tmp_hs = (HeadSeq *) value; 

	  if((tmp_hs[0].maq >= data->MAQ) && (tmp_hs[1].maq >= data->MAQ) && (tmp_hs[0].type == 1) && (tmp_hs[1].type == 1)) {
	     int tmp_qsum = 0;
	     for(int j=0; j<tmp_hs[0].length; j++) tmp_qsum += get_Qscore(tmp_hs[1].quality[j]);  
	     for(int j=0; j<tmp_hs[1].length; j++) tmp_qsum += get_Qscore(tmp_hs[2].quality[j]); 
	     if(tmp_qsum > qsum) { index = i; qsum = tmp_qsum; } 
	  }

	  value += valuebytes[i];
       }

       if((index >= 0) && (qsum > 0)) {
	 value = multivalue;
	 for(int i=0; i<index; i++) value += valuebytes[i];
	 kv->add(key, keybytes, value, valuebytes[index]);
       }

   } else {

      value = multivalue;
      HeadSeq *tmp_hs = (HeadSeq *) value;
      if((tmp_hs[0].maq >= data->MAQ) && (tmp_hs[1].maq >= data->MAQ)&&(tmp_hs[0].type == 1) && (tmp_hs[1].type == 1)) 
		kv->add(key, keybytes, value, valuebytes[0]);

   }
   
   END_BLOCK_LOOP

}

void reduce_consensus_SeqQual_singleton(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr)
{
   PairLink *plink = (PairLink *) key;

  if((nvalues < 3) && (plink->chr1 > 0) && (plink->chr2 > 0)) {

   Data *data = (Data *) ptr;
   char *value;
   string seq1, seq2, quality1, quality2;
   int PE1, PE2, strand1, strand2;

//   plink->ID = nvalues;

   uint64_t nvalues_total;
   CHECK_FOR_BLOCKS(multivalue,valuebytes,nvalues,nvalues_total)
   BEGIN_BLOCK_LOOP(multivalue,valuebytes,nvalues)

   if(nvalues == 2) {

        int str_loc, end_loc;
        if(plink->loc1 <= plink->loc2) { str_loc = plink->loc1; end_loc = plink->loc2; }
        else                           { str_loc = plink->loc2; end_loc = plink->loc1; }

        bool check_len = false;
//        int len = end_loc - str_loc + 1;
//        if( (str_loc + len + 1) < data->reference_sequence[plink->chr1].length() ) { check_len = true; }
	if((end_loc + data->range) < data->reference_sequence[plink->chr1].length() ) { check_len = true; }

        if(check_len && (str_loc > 0) ) {

           int qsum1 = 0;

           value = multivalue;
           HeadSeq *tmp_hs = (HeadSeq *) value;
           HeadSeq tmp_hs1 = tmp_hs[0];
           HeadSeq tmp_hs2 = tmp_hs[1];
	
	   for(int i=0; i<tmp_hs1.length; i++) qsum1 += get_Qscore(tmp_hs1.quality[i]);
	   for(int i=0; i<tmp_hs2.length; i++) qsum1 += get_Qscore(tmp_hs2.quality[i]); 

	   int qsum2 = 0;
	
           value += valuebytes[0];
           tmp_hs = (HeadSeq *) value;
           tmp_hs1 = tmp_hs[0];
           tmp_hs2 = tmp_hs[1];

	   for(int i=0; i<tmp_hs1.length; i++) qsum2 += get_Qscore(tmp_hs1.quality[i]);
           for(int i=0; i<tmp_hs2.length; i++) qsum2 += get_Qscore(tmp_hs2.quality[i]);	 

           if(qsum1 >= qsum2) {
              value = multivalue;
//              HeadSeq *out_hs0 = (HeadSeq *) value;
//              plink->maq1 = out_hs0[0].maq;
//              plink->maq2 = out_hs0[1].maq;
//              kv->add((char *) plink, sizeof(PairLink), value, valuebytes[0]);
		kv->add(key, keybytes, value, valuebytes[0]);
           } else {
              value = multivalue;
              value += valuebytes[0];
//              HeadSeq *out_hs1 = (HeadSeq *) value;
//              plink->maq1 = out_hs1[0].maq;
//              plink->maq2 = out_hs1[1].maq;
//              kv->add((char *) plink, sizeof(PairLink), value, valuebytes[1]);
		kv->add(key, keybytes, value, valuebytes[1]);
           }

        }

    } 

    if(nvalues == 1) {

        value = multivalue;
//        HeadSeq *out_hs2 = (HeadSeq *) value;
//        plink->maq1 = out_hs2[0].maq;
//        plink->maq2 = out_hs2[1].maq;
//        kv->add((char *) plink, sizeof(PairLink), value, valuebytes[0]);
	kv->add(key, keybytes, value, valuebytes[0]);
    }

   END_BLOCK_LOOP

 }


}


void reduce_consensus_Single_singleton(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr)
{

   PairLink *plink = (PairLink *) key;

if((nvalues < 3) && (plink->chr1 > 0) && (plink->chr2 == 0)) {

   Data *data = (Data *) ptr;
   char *value;
   RefConsen *refconsen;
   string seq1, seq2, quality1, quality2;
   int PE1, PE2, strand1, strand2;

   plink->ID = nvalues;

   StripedSmithWaterman::Aligner aligner(1,4,6,1);
   StripedSmithWaterman::Filter filter;
   StripedSmithWaterman::Alignment alignment;

   uint64_t nvalues_total;
   CHECK_FOR_BLOCKS(multivalue,valuebytes,nvalues,nvalues_total)
   BEGIN_BLOCK_LOOP(multivalue,valuebytes,nvalues)

   if(nvalues == 2) {

        int str_loc = plink->loc1; 
	int end_loc = plink->loc1 + data->range - 1; 

        bool check_len = false;
        int len = end_loc - str_loc + 1;
        if( (end_loc - str_loc + 1) < data->reference_sequence[plink->chr1].length() ) { check_len = true; }


        if(check_len) {

           string ref_sequence = data->reference_sequence[plink->chr1].substr(str_loc-1, len);

           value = multivalue;
           HeadSeq tmp_hs = *(HeadSeq *) value;

           string seq;
           seq.assign( tmp_hs.seq, tmp_hs.length );

           alignment.Clear();
           aligner.Align(seq.c_str(), ref_sequence.c_str(), ref_sequence.length(), filter, &alignment);
           uint64_t sw_score1 = alignment.sw_score;

           value += valuebytes[0];
           tmp_hs = *(HeadSeq *) value;

           seq.assign( tmp_hs.seq, tmp_hs.length );

           alignment.Clear();
           aligner.Align(seq.c_str(), ref_sequence.c_str(), ref_sequence.length(), filter, &alignment);
           uint64_t sw_score2 = alignment.sw_score;

           if(sw_score1 >= sw_score2) {
              value = multivalue;
              HeadSeq out_hs0 = *(HeadSeq *) value;
              plink->maq1 = out_hs0.maq;
              kv->add((char *) plink, sizeof(PairLink), value, valuebytes[0]);
           } else {
              value = multivalue;
              value += valuebytes[0];
              HeadSeq out_hs1 = *(HeadSeq *) value;
              plink->maq1 = out_hs1.maq;
              kv->add((char *) plink, sizeof(PairLink), value, valuebytes[1]);
           }

       }

   } else if(nvalues == 1) {

        value = multivalue;
        HeadSeq out_hs2 = *(HeadSeq *) value;
        plink->maq1 = out_hs2.maq;
        kv->add((char *) plink, sizeof(PairLink), value, valuebytes[0]);      

   }

   END_BLOCK_LOOP

 }
}






void reduce_consensus_sequence_quality(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr)
{
   Data *data = (Data *) ptr;
   char *value;
   RefConsen *refconsen;
   string seq1, seq2, quality1, quality2;
   int PE1, PE2, strand1, strand2;

   PairLink *plink = (PairLink *) key;

   IRKE irke_1(data->kmer_length, data->MAX_RECURSION, data->MIN_SEED_ENTROPY, data->MIN_SEED_COVERAGE,
                data->min_any_entropy, data->PACMAN, data->CRAWL, data->crawl_length, data->DS);

   IRKE irke_2(data->kmer_length, data->MAX_RECURSION, data->MIN_SEED_ENTROPY, data->MIN_SEED_COVERAGE,
                data->min_any_entropy, data->PACMAN, data->CRAWL, data->crawl_length, data->DS);

   IRKE irke(data->kmer_length, data->MAX_RECURSION, data->MIN_SEED_ENTROPY, data->MIN_SEED_COVERAGE,
                data->min_any_entropy, data->PACMAN, data->CRAWL, data->crawl_length, data->DS);

   StripedSmithWaterman::Aligner aligner1(1,4,6,1);
   StripedSmithWaterman::Aligner aligner2(1,4,6,1);
   StripedSmithWaterman::Filter filter1;
   StripedSmithWaterman::Filter filter2;
   StripedSmithWaterman::Alignment alignment1;
   StripedSmithWaterman::Alignment alignment2;

   int *cquality1 = (int *) malloc( sizeof(int)*data->range );
   int *cquality2 = (int *) malloc( sizeof(int)*data->range );  

   for(int i=0; i<data->range; i++) { cquality1[i]=0; cquality2[i]=0; }

   uint64_t nvalues_total;
   CHECK_FOR_BLOCKS(multivalue,valuebytes,nvalues,nvalues_total)
   BEGIN_BLOCK_LOOP(multivalue,valuebytes,nvalues)

   if(nvalues >= 3) {

    plink->maq1 = 0;
    plink->maq2 = 0;

    value = multivalue;
    for(int i=0; i<nvalues; i++) {
          HeadSeq *headseq_pair = (HeadSeq *) value;
	  plink->maq1 += headseq_pair[0].maq;
	  plink->maq2 += headseq_pair[1].maq;

	  value += valuebytes[i];
    }
    plink->maq1 = (int)(plink->maq1 / nvalues);
    plink->maq2 = (int)(plink->maq2 / nvalues);

   }

    HeadSeq *out_headseq = (HeadSeq *) malloc( sizeof(HeadSeq) * 2 );

    if(nvalues >= 3) {

      value = multivalue;
      for(int i=0; i<nvalues; i++) {
          HeadSeq *headseq_pair = (HeadSeq *) value;

          HeadSeq headseq = headseq_pair[0];
          seq1.assign( headseq.seq, data->range );
	  quality1.assign( headseq.quality, data->range);
	  PE1 = headseq.PE;
	  strand1 = headseq.strand;
	  if(strand1 < 0) {
		string tmp = revquality( quality1 );
		quality1.assign( tmp.c_str(), data->range );
	  }

	  for(int j=0; j<data->range; j++) cquality1[j] += get_Qscore( quality1[j] );

          headseq = headseq_pair[1];
          seq2.assign( headseq.seq, data->range );
	  quality2.assign( headseq.quality, data->range);
	  PE2 = headseq.PE;
	  strand2 = headseq.strand;
	  if(strand2 < 0) {
                string tmp = revquality( quality2 );
                quality2.assign( tmp.c_str(), data->range );
          }
	  for(int j=0; j<data->range; j++) cquality2[j] += get_Qscore( quality2[j] );

	  int nKmer = data->range - data->kmer_length + 1;
          for(int j=0; j<nKmer; j++) {
            string seq_kmer1 = seq1.substr(j, data->kmer_length);
            if( contains_non_gatc(seq_kmer1) ) {
              kmer_int_type_t kmer_val = kmer_to_intval(seq_kmer1);
              if(strand1 < 0)  kmer_val = revcomp_val(kmer_val, data->kmer_length);
              irke_1.add_kmer(kmer_val , 1);
            }
          }

          for(int j=0; j<nKmer; j++) {
            string seq_kmer2 = seq2.substr(j, data->kmer_length);
            if( contains_non_gatc(seq_kmer2) ) {
              kmer_int_type_t kmer_val = kmer_to_intval(seq_kmer2);
              if(strand2 < 0)  kmer_val = revcomp_val(kmer_val, data->kmer_length);
              irke_2.add_kmer(kmer_val , 1);
            }
          }

	value += valuebytes[i];   
      }

      for(int j=0; j<data->range; j++) { cquality1[j]=(int)(cquality1[j]/nvalues); cquality2[j]=(int)(cquality2[j]/nvalues); }

      irke_1.prune_some_kmers(data->min_kmer_count, data->min_any_entropy, data->prune_error_kmers, data->min_ratio_non_error);
      irke_1.populate_sorted_kmers_list();

      string csequence_1 =
        irke_1.compute_sequence_assemblies( data->MIN_CONNECTIVITY_RATIO, data->MIN_ASSEMBLY_LENGTH, data->MIN_ASSEMBLY_COVERAGE );

      irke_2.prune_some_kmers(data->min_kmer_count, data->min_any_entropy, data->prune_error_kmers, data->min_ratio_non_error);
      irke_2.populate_sorted_kmers_list();

      string csequence_2 =
        irke_2.compute_sequence_assemblies( data->MIN_CONNECTIVITY_RATIO, data->MIN_ASSEMBLY_LENGTH, data->MIN_ASSEMBLY_COVERAGE );

      alignment1.Clear();
      if(plink->loc1 <= plink->loc2) aligner1.Align(csequence_2.c_str(), csequence_1.c_str(), csequence_1.length(), filter1, &alignment1);
      else			     aligner1.Align(csequence_1.c_str(), csequence_2.c_str(), csequence_2.length(), filter1, &alignment1);


      int overlap = 0;
      int32_t *bam_cigar = (int32_t *) malloc( sizeof(int32_t)*3 );
      if( alignment1.cigar.size() > 0) {
		convert_cigar2bam(alignment1.cigar[0], bam_cigar);
	        if( (bam_cigar[0] == 0) || (bam_cigar[0] == 7) || (bam_cigar[0] == 8) ) overlap += (int)bam_cigar[1];
      }
      free(bam_cigar);

      value = multivalue;
      for(int i=0; i<nvalues; i++) {
             if(valuebytes[i] == sizeof(HeadSeq)*2) {
                HeadSeq *headseq_pair = (HeadSeq *) value;

          	HeadSeq headseq = headseq_pair[0];
          	seq1.assign( headseq.seq, data->range );
          	quality1.assign( headseq.quality, data->range);
          	PE1 = headseq.PE;
          	strand1 = headseq.strand;

          	headseq = headseq_pair[1];
          	seq2.assign( headseq.seq, data->range );
          	quality2.assign( headseq.quality, data->range);
          	PE2 = headseq.PE;
          	strand2 = headseq.strand;

          	int nKmer = data->range - data->kmer_length + 1;
          	for(int j=0; j<nKmer; j++) {
            	string seq_kmer1 = seq1.substr(j, data->kmer_length);
            	  if( contains_non_gatc(seq_kmer1) ) {
              		kmer_int_type_t kmer_val = kmer_to_intval(seq_kmer1);
              		if(strand1 < 0)  kmer_val = revcomp_val(kmer_val, data->kmer_length);
              		irke.add_kmer(kmer_val , 1);
            	  }
          	}

          	for(int j=0; j<nKmer; j++) {
            	string seq_kmer2 = seq2.substr(j, data->kmer_length);
            	  if( contains_non_gatc(seq_kmer2) ) {
              		kmer_int_type_t kmer_val = kmer_to_intval(seq_kmer2);
              		if(strand2 < 0)  kmer_val = revcomp_val(kmer_val, data->kmer_length);
              		irke.add_kmer(kmer_val , 1);
            	  }
          	}

             }
             value += valuebytes[i];
      }

      irke.prune_some_kmers(data->min_kmer_count, data->min_any_entropy, data->prune_error_kmers, data->min_ratio_non_error);
      irke.populate_sorted_kmers_list();

      string csequence =
            irke.compute_sequence_assemblies( data->MIN_CONNECTIVITY_RATIO, data->MIN_ASSEMBLY_LENGTH, data->MIN_ASSEMBLY_COVERAGE );

      if( (overlap > data->OVERLAP) && (csequence.length() >= csequence_1.length()) && (csequence.length() >= csequence_2.length()) ) {
	
	  alignment1.Clear();
          aligner1.Align(csequence_1.c_str(), csequence.c_str(), csequence.length(), filter1, &alignment1);

	  string tmp1 = align_csequence(alignment1, csequence.c_str(), csequence_1.c_str(), csequence.length(), data->range );
	  if(plink->strand1 < 0) {
		string rev_tmp1 = revcomp( tmp1.c_str() );
	        seq1.assign(rev_tmp1.c_str(), tmp1.length());
	  } else {
		seq1.assign(tmp1.c_str(), tmp1.length());
	  }
 
	  alignment2.Clear();
          aligner2.Align(csequence_2.c_str(), csequence.c_str(), csequence.length(), filter2, &alignment2); 

	  string tmp2 = align_csequence(alignment2, csequence.c_str(), csequence_2.c_str(), csequence.length(), data->range );
	  if(plink->strand2 < 0) {
		string rev_tmp2 = revcomp( tmp2.c_str() );
          	seq2.assign(rev_tmp2.c_str(), tmp2.length());
	  } else {
		seq2.assign(tmp2.c_str(), tmp2.length());
	  }


     if( ((int)tmp1.length() == data->range) && ((int)tmp2.length() == data->range) ) {

          for(int j=0; j<data->range; j++) quality1[j] = get_QualityAscii( cquality1[j] );
          string tmp1_quality = align_quality(alignment1, csequence.c_str(), quality1.c_str(), csequence.length(), data->range );
          if(plink->strand1 < 0) {
                string rev_tmp1_quality = revquality( tmp1_quality );
                quality1.assign( rev_tmp1_quality.c_str(), tmp1_quality.length() );
          } else {
                quality1.assign( tmp1_quality.c_str(), tmp1_quality.length() );
          }


          for(int j=0; j<data->range; j++) quality2[j] = get_QualityAscii( cquality2[j] );
          string tmp2_quality = align_quality(alignment2, csequence.c_str(), quality2.c_str(), csequence.length(), data->range );
	  quality2.assign( tmp2_quality.c_str(), tmp2_quality.length() ); 
          if(plink->strand2 < 0) {
                string rev_tmp2_quality = revquality( tmp2_quality );
                quality2.assign( rev_tmp2_quality.c_str(), tmp2_quality.length() );
          } else {
                quality2.assign( tmp2_quality.c_str(), tmp2_quality.length() );
          }


	  HeadSeq seq_quality1, seq_quality2;
	
          memcpy( seq_quality1.seq, seq1.c_str(), sizeof(char)*seq1.length() );
          memcpy( seq_quality2.seq, seq2.c_str(), sizeof(char)*seq2.length() );

	  memcpy( seq_quality1.quality, quality1.c_str(), sizeof(char)*seq1.length() );
          memcpy( seq_quality2.quality, quality2.c_str(), sizeof(char)*seq2.length() );

	  seq_quality1.strand = plink->strand1;
	  seq_quality2.strand = plink->strand2;

	  seq_quality1.length = (int)seq1.length();
	  seq_quality2.length = (int)seq2.length();

	  out_headseq[0] = seq_quality1;
          out_headseq[1] = seq_quality2;

	    kv->add((char *) plink, sizeof(PairLink), (char *) out_headseq, sizeof(HeadSeq)*2 ); 

   }

       }  else {

	  HeadSeq seq_quality1, seq_quality2;

	  if(plink->strand1 < 0) {
		string ctmp1 = revcomp( csequence_1.c_str() );
		csequence_1.assign( ctmp1.c_str(), csequence_1.length() );	
	  } 
	  if(plink->strand2 < 0) {
		string ctmp2 = revcomp( csequence_2.c_str() );
		csequence_2.assign( ctmp2.c_str(), csequence_2.length() );
	  }

          memcpy( seq_quality1.seq, csequence_1.c_str(), sizeof(char)*data->range );
          memcpy( seq_quality2.seq, csequence_2.c_str(), sizeof(char)*data->range );

          for(int j=0; j<data->range; j++) quality1[j] = get_QualityAscii( cquality1[j] );
          if(plink->strand1 < 0) {
                string rev_tmp1_quality = revquality( quality1 );
                quality1.assign( rev_tmp1_quality.c_str(), data->range );
	  }

	  for(int j=0; j<data->range; j++) quality2[j] = get_QualityAscii( cquality2[j] );
          if(plink->strand2 < 0) {
                string rev_tmp2_quality = revquality( quality2 );
                quality2.assign( rev_tmp2_quality.c_str(), data->range );
          }

          memcpy( seq_quality1.quality, quality1.c_str(), sizeof(char)*data->range );
          memcpy( seq_quality2.quality, quality2.c_str(), sizeof(char)*data->range );

	  seq_quality1.strand = plink->strand1;
          seq_quality2.strand = plink->strand2;

	  seq_quality1.length = data->range;
          seq_quality2.length = data->range;	  

          out_headseq[0] = seq_quality1;
          out_headseq[1] = seq_quality2;

	  if( ((int)csequence_1.length() == data->range) && ((int)csequence_2.length() == data->range) )
                  kv->add((char *) plink, sizeof(PairLink), (char *) out_headseq, sizeof(HeadSeq)*2 );

      }
 
     } else if(nvalues == 2) {

	int str_loc, end_loc;
	if(plink->loc1 <= plink->loc2) { str_loc = plink->loc1; end_loc = plink->loc2; }
	else			       { str_loc = plink->loc2; end_loc = plink->loc1; } 

	bool check_len = false;
	int len = end_loc - str_loc + 1;
        if( (str_loc + (2*len)) > data->reference_sequence[plink->chr1].length() ) {
                len = data->reference_sequence[plink->chr1].length() - str_loc + 1;
	} else {
		len = 2 * len; check_len = true;
        }

     if(check_len) {

	string ref_sequence = data->reference_sequence[plink->chr1].substr(str_loc-1, len);  

	uint64_t sw_score1 = 0;
	uint64_t sw_score2 = 0;

	value = multivalue;
	HeadSeq *tmp_hs = (HeadSeq *) value;
	HeadSeq tmp_hs1 = tmp_hs[0];
	HeadSeq tmp_hs2 = tmp_hs[1];

	string seq11, seq12;
	seq11.assign( tmp_hs1.seq, data->range );	
	seq12.assign( tmp_hs2.seq, data->range );	

	if(tmp_hs1.strand < 0) seq11 = revcomp( seq11.c_str() );
	if(tmp_hs2.strand < 0) seq12 = revcomp( seq12.c_str() );

	alignment1.Clear();
        aligner1.Align(seq11.c_str(), ref_sequence.c_str(), ref_sequence.length(), filter1, &alignment1);	
	sw_score1 += alignment1.sw_score;

	alignment1.Clear();
        aligner1.Align(seq12.c_str(), ref_sequence.c_str(), ref_sequence.length(), filter1, &alignment1);
        sw_score1 += alignment1.sw_score;	

	value += valuebytes[0];
	tmp_hs = (HeadSeq *) value;
        tmp_hs1 = tmp_hs[0];
        tmp_hs2 = tmp_hs[1];

	seq11.assign( tmp_hs1.seq, data->range );
        seq12.assign( tmp_hs2.seq, data->range );

        if(tmp_hs1.strand < 0) seq11 = revcomp( seq11.c_str() );
        if(tmp_hs2.strand < 0) seq12 = revcomp( seq12.c_str() );

        alignment2.Clear();
        aligner2.Align(seq11.c_str(), ref_sequence.c_str(), ref_sequence.length(), filter2, &alignment2);
        sw_score2 += alignment2.sw_score;

        alignment2.Clear();
        aligner2.Align(seq12.c_str(), ref_sequence.c_str(), ref_sequence.length(), filter2, &alignment2);
        sw_score2 += alignment2.sw_score;

	if(sw_score1 >= sw_score2) {
	   value = multivalue;
	   HeadSeq *out_hs0 = (HeadSeq *) value;
	   plink->maq1 = out_hs0[0].maq;
	   plink->maq2 = out_hs0[1].maq;
	   kv->add((char *) plink, sizeof(PairLink), value, valuebytes[0]);
	} else {
	   value = multivalue;
	   value += valuebytes[0];
	   HeadSeq *out_hs1 = (HeadSeq *) value;
           plink->maq1 = out_hs1[0].maq;
           plink->maq2 = out_hs1[1].maq;
	   kv->add((char *) plink, sizeof(PairLink), value, valuebytes[1]);
	}

       }

     } else if(nvalues == 1) {
	value = multivalue;
	HeadSeq *out_hs2 = (HeadSeq *) value;
        plink->maq1 = out_hs2[0].maq;
        plink->maq2 = out_hs2[1].maq;
        if(valuebytes[0] == sizeof(HeadSeq)*2) kv->add((char *) plink, sizeof(PairLink), value, valuebytes[0]);
     }

        int nA = 0;
        int nB = 0;
        int nC = 0;
        int nD = 0;

     if(nvalues >= 3) {
        value = multivalue;
        for(int i=0; i<nvalues; i++) {
            HeadSeq *headseq_pair = (HeadSeq *) value;

            if(headseq_pair[0].header[0] == 'A') nA += 1;
            if(headseq_pair[0].header[0] == 'B') nB += 1;
            if(headseq_pair[0].header[0] == 'C') nC += 1;
            if(headseq_pair[0].header[0] == 'D') nD += 1;

           value += valuebytes[i];
        }
     }

        int check = 0;
        if(nA > 0) check += 1;
        if(nB > 0) check += 1;
        if(nC > 0) check += 1;
        if(nD > 0) check += 1;

        if(check > 1)       data->flank   += 1;
        else if(check == 1) data->offside += 1;

     free(out_headseq);
     free(cquality1);
     free(cquality2);

   END_BLOCK_LOOP

}

void find_error_base(StripedSmithWaterman::Alignment *alignment, int *error_index, int *error_base, int *ref_base, const char *seq, const char *ref_seq, int seq_length)
{
  for(int i=0; i<seq_length; i++) error_index[i] = 0;
  int32_t *bam_cigar = (int32_t *) malloc( sizeof(int32_t)*2 );

  int sp     = (int)alignment->query_begin;
  int sp_loc = (int)alignment->ref_begin;
  for(size_t j = 0; j < alignment->cigar.size(); j++)  {
       convert_cigar2bam(alignment->cigar[j], bam_cigar);

       if( (bam_cigar[0] == 0) || (bam_cigar[0] == 7) || (bam_cigar[0] == 8) ) {
		for(int k=0; k<(int)bam_cigar[1]; k++) {
			int base  = get_int_from_base( seq[sp+k] );
			int cbase = get_int_from_base( ref_seq[sp_loc+k] );
			if(base != cbase) {
				error_index[sp+k]++;
				error_base[sp+k] = base;
        			ref_base[sp+k]   = cbase;	
			}
		}
		
                sp     += (int)bam_cigar[1];
                sp_loc += (int)bam_cigar[1];
       } else if( bam_cigar[0] == 2 )  {
                sp_loc += (int)bam_cigar[1];
       } else if( bam_cigar[0] == 1) {
                sp     += (int)bam_cigar[1];
       }

  }

  free(bam_cigar);  
}


void align_seq(StripedSmithWaterman::Alignment alignment, int seq_len, string seq, string ref_sequence, int *index)
{
   for(int i=0; i<seq_len; i++) index[i] = 0;
   int32_t *bam_cigar = (int32_t *) malloc( sizeof(int32_t)*3 );

/*
   int sp = (int)alignment.query_begin;
   for(size_t j = 0; j < alignment.cigar.size(); j++)  {
	convert_cigar2bam(alignment.cigar[j], bam_cigar);

	if(bam_cigar[0] == 7) {
	  for(int k=sp; k<(int)bam_cigar[1]; k++) index[k] = 0;
	}

       if( (bam_cigar[0] == 0) || (bam_cigar[0] == 1) || (bam_cigar[0] == 7) || ( bam_cigar[0] == 8) ) sp += (int)bam_cigar[1];
   } 
   free(bam_cigar);
*/

   int sp     = (int)alignment.query_begin;
   int sp_loc = (int)alignment.ref_begin;
   for(size_t j = 0; j < alignment.cigar.size(); j++)  {
          convert_cigar2bam(alignment.cigar[j], bam_cigar);

         if( (bam_cigar[0] == 0) || (bam_cigar[0] == 3) || (bam_cigar[0] == 7) || (bam_cigar[0] == 8) ) {

	     if(bam_cigar[0] == 8) {
              for(int k=0; k<(int)bam_cigar[1]; k++) index[sp+k] = 1; 
	     }
              sp     += (int)bam_cigar[1];
              sp_loc += (int)bam_cigar[1];

         } else if( (bam_cigar[0] == 2) || (bam_cigar[0] == 6) ) sp_loc += (int)bam_cigar[1];
           else if(  bam_cigar[0] == 1 )                         sp     += (int)bam_cigar[1];

   }

   free(bam_cigar);
}

string align_csequence(StripedSmithWaterman::Alignment alignment, string ref_seq, string query_seq, int ref_len, int query_len)
{

   stringstream OUT_STR;

   if( (alignment.ref_begin - alignment.query_begin) >= 0 ) {
             if(alignment.query_begin > 0) {
                    OUT_STR << ref_seq.substr((int)(alignment.ref_begin - alignment.query_begin), (int)alignment.query_begin).c_str();
             }
   } else {
          if(alignment.ref_begin > 0)
		OUT_STR << ref_seq.substr(0,(int)(alignment.query_begin - alignment.ref_begin)).c_str();
   }
 
  int32_t *bam_cigar = (int32_t *) malloc( sizeof(int32_t)*3 );

  int sp     = (int)alignment.query_begin;
  int sp_loc = (int)alignment.ref_begin;
  for(size_t j = 0; j < alignment.cigar.size(); j++)  {
       convert_cigar2bam(alignment.cigar[j], bam_cigar);

       if( (bam_cigar[0] == 0) || (bam_cigar[0] == 7) || (bam_cigar[0] == 8) ) {
                OUT_STR << ref_seq.substr(sp_loc,(int)bam_cigar[1]).c_str();

                sp     += (int)bam_cigar[1];
                sp_loc += (int)bam_cigar[1];
       } else if( bam_cigar[0] == 2 )  {
                sp_loc += (int)bam_cigar[1];
       } else if( bam_cigar[0] == 1) {
		OUT_STR << query_seq.substr(sp,(int)bam_cigar[1]).c_str();
                sp     += (int)bam_cigar[1];
       } 

  }
  
  if(bam_cigar[0] == 4) 
          OUT_STR << ref_seq.substr(sp_loc,bam_cigar[1]);

  free(bam_cigar);

  return(OUT_STR.str().c_str());
}

string align_quality(StripedSmithWaterman::Alignment alignment, string ref_seq, string query_quality, int ref_len, int query_len)
{

   stringstream OUT_STR;

   if(alignment.query_begin > 0) 
	OUT_STR << query_quality.substr(0,alignment.query_begin).c_str();

   int32_t *bam_cigar = (int32_t *) malloc( sizeof(int32_t)*3 );

   int sp     = (int)alignment.query_begin;
   int sp_loc = (int)alignment.ref_begin;
   for(size_t j = 0; j < alignment.cigar.size(); j++)  {
       convert_cigar2bam(alignment.cigar[j], bam_cigar);

       if( (bam_cigar[0] == 0) || (bam_cigar[0] == 7) || (bam_cigar[0] == 8) ) {

                OUT_STR << query_quality.substr(sp,(int)bam_cigar[1]).c_str();

                sp     += (int)bam_cigar[1];
                sp_loc += (int)bam_cigar[1];

       } else if( bam_cigar[0] == 2 )  {
                sp_loc += (int)bam_cigar[1];
       } else if( bam_cigar[0] == 1) {

                OUT_STR << query_quality.substr(sp,(int)bam_cigar[1]).c_str();
                sp     += (int)bam_cigar[1];

       }

  }

  if(bam_cigar[0] == 4) {
          OUT_STR << query_quality.substr(sp,bam_cigar[1]);

  }

  free(bam_cigar);

  return(OUT_STR.str().c_str());
}

void map_print2samformat_singleton(uint64_t itask, char *key, int keybytes, char *value, int valuebytes, KeyValue *kv, void *ptr)
{
  Data *data = (Data *) ptr;
  PairLink *plink = (PairLink *) key;


  HeadSeq head_seq = *(HeadSeq *) value;

  bool check_len = false; 

  if( (plink->loc1+2*data->range) < data->reference_sequence[plink->chr1].length() ) check_len = true;

  if(check_len) {

    string ref_sequence = data->reference_sequence[plink->chr1].substr(plink->loc1 - 1,data->range*2);
    string seq, rev_seq, quality;

    seq.assign(head_seq.seq,        data->range);
    quality.assign(head_seq.quality,data->range);

    StripedSmithWaterman::Aligner aligner(1,4,6,1);
    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment alignment;

    alignment.Clear();
    aligner.Align(seq.substr(0,data->range).c_str(), ref_sequence.c_str(), ref_sequence.length(), filter, &alignment);


    int flag = 8; 

    if(head_seq.PE == 1) flag += 64;
    else                 flag += 128;

    if(head_seq.strand < 0) flag += 16;
    else                    flag += 32;

    stringstream head;
    head << "CRUKMI:" << plink->chr1 << ":" << plink->chr2 << ":" << plink->loc1 << ":" << plink->loc2 << ":" << plink->strand1 << ":" << plink->strand2 << ":" << plink->ID;

    int32_t TLEN = 0;

    if( (!alignment.cigar_string.empty()) && (!seq.empty()) ) {

      data->outFile << head.str().c_str() << "\t" << flag << "\t";

      if((plink->chr1 >= 1) && (plink->chr1 <= 22)) data->outFile << plink->chr1;
      else if(plink->chr1 == 23)                    data->outFile << "X";
      else if(plink->chr1 == 24)                    data->outFile << "Y";
      else if(plink->chr1 == 25)                    data->outFile << "MT";

      data->outFile << "\t" << plink->loc1 + alignment.ref_begin << "\t"
                    << plink->maq1        << "\t" << alignment.cigar_string.c_str() << "\t"
                    << "*"                << "\t" <<  0  << "\t" << 0 << "\t";

      data->outFile << seq.substr(0,data->range).c_str();

      data->outFile << "\t" << quality.substr(0,data->range).c_str() << "\t" << "PG:Z:MarkDuplicates" << "\t" << "RG:Z:AA1" << endl;

    }


  }

}





void map_print2samformat(uint64_t itask, char *key, int keybytes, char *value, int valuebytes, KeyValue *kv, void *ptr)
{
  Data *data = (Data *) ptr;
  PairLink *plink = (PairLink *) key;

 if((plink->chr1 > 0) && (plink->chr2 > 0) && (plink->chr1 == plink->chr2)) {

  HeadSeq *tmp_hs = (HeadSeq *) value;
  HeadSeq head_seq1 = tmp_hs[0];
  HeadSeq head_seq2 = tmp_hs[1];

  bool check_len = false;
  int start, len;

  if(plink->loc1 < plink->loc2) { 
	start = plink->loc1; 
	len = plink->loc2 - plink->loc1 + 1; 
	if( (plink->loc1+2*len) > data->reference_sequence[plink->chr1].length() ) {
		len = data->reference_sequence[plink->chr1].length() - plink->loc1 + 1;
	} else  {
		len = 2*len;  check_len = true;
	}
  } else  { 
	start = plink->loc2; 
	len = plink->loc1 - plink->loc2 + 1; 
	if( (plink->loc2+2*len) > data->reference_sequence[plink->chr1].length() ) {
		len = data->reference_sequence[plink->chr1].length() - plink->loc2 + 1;
        } else {
		len = 2*len; check_len = true;
	}
  } 

if(check_len && (start > 1)) { 

  string ref_sequence = data->reference_sequence[plink->chr1].substr(start-1,len); 
  string seq1, seq2, rev_seq1, rev_seq2, quality1, quality2;
 
  seq1.assign(head_seq1.seq,data->range); 
  seq2.assign(head_seq2.seq,data->range);

  quality1.assign(head_seq1.quality,data->range);
  quality2.assign(head_seq2.quality,data->range);

  StripedSmithWaterman::Aligner aligner1(1,4,6,1);
  StripedSmithWaterman::Filter filter1;
  StripedSmithWaterman::Alignment alignment1;

  StripedSmithWaterman::Aligner aligner2(1,4,6,1);
  StripedSmithWaterman::Filter filter2;
  StripedSmithWaterman::Alignment alignment2;
 
  alignment1.Clear(); 
  aligner1.Align(seq1.substr(0,data->range).c_str(), ref_sequence.c_str(), ref_sequence.length(), filter1, &alignment1);

  alignment2.Clear();
  aligner2.Align(seq2.substr(0,data->range).c_str(), ref_sequence.c_str(), ref_sequence.length(), filter2, &alignment2); 

  int flag1 = 1 + 2;
  int flag2 = 1 + 2; 

  if(head_seq1.PE == 1) flag1 += 64;
  else			flag1 += 128;

  if(head_seq2.PE == 1) flag2 += 64;
  else                  flag2 += 128;

  if(head_seq1.strand < 0) flag1 += 16;
  else	       	           flag1 += 32; 

  if(head_seq2.strand < 0) flag2 += 16;
  else			   flag2 += 32; 

  stringstream head;
  head << "PEC:" << plink->chr1 << ":" << plink->loc1 << ":" << plink->loc2 << ":" << plink->strand1 << ":" << plink->strand2 << ":" << plink->ID;

  int32_t TLEN1, TLEN2;
 
  if( alignment1.ref_begin <= alignment2.ref_begin ) { TLEN1 = alignment2.ref_end - alignment1.ref_begin + 1; TLEN2 = alignment2.ref_end - alignment1.ref_begin + 1; }
  else						     { TLEN1 = alignment1.ref_end - alignment2.ref_begin + 1; TLEN2 = alignment1.ref_end - alignment2.ref_begin + 1; } 


  int len_seq = data->range;

  if(plink->strand1 < 0) TLEN1 = -1 * TLEN1;
  else			 TLEN2 = -1 * TLEN2;

  int align_len1 = parse_cigar_align_length(const_cast<char*>(alignment1.cigar_string.c_str()));
  int align_len2 = parse_cigar_align_length(const_cast<char*>(alignment2.cigar_string.c_str()));  

if( (plink->chr1 == plink->chr2) && (plink->strand1*plink->strand2 == -1) && 
    (!alignment1.cigar_string.empty()) && (!alignment2.cigar_string.empty()) && 
    (!seq1.empty()) && (!seq2.empty())  
  ) {

    data->outFile << head.str().c_str() << "\t" << flag1 << "\t"; 

    if((plink->chr1 >= 1) && (plink->chr1 <= 22)) data->outFile << plink->chr1;       
    else if(plink->chr1 == 23)		        data->outFile << "X";
    else if(plink->chr1 == 24)		        data->outFile << "Y";
    else if(plink->chr1 == 25) 		        data->outFile << "MT";

    data->outFile << "\t" << start + alignment1.ref_begin << "\t" 
                << plink->maq1        << "\t" << alignment1.cigar_string.c_str() << "\t" 
	        << "="                << "\t" << start + alignment2.ref_begin << "\t" << TLEN1 << "\t";

         data->outFile << seq1.substr(0,len_seq).c_str();
         data->outFile << "\t" << quality1.substr(0,len_seq).c_str() << "\t" << "PG:Z:MarkDuplicates" << "\t" << "RG:Z:AA1" << endl;

    data->outFile << head.str().c_str() << "\t" << flag2 << "\t"; 

    if((plink->chr2 >= 1) && (plink->chr2 <= 22)) data->outFile << plink->chr2;       
    else if(plink->chr2 == 23)		        data->outFile << "X";
    else if(plink->chr2 == 24)			data->outFile << "Y";
    else if(plink->chr2 == 25)			data->outFile << "MT";

    data->outFile << "\t" << start + alignment2.ref_begin << "\t" 
	        << plink->maq2        << "\t" << alignment2.cigar_string.c_str() << "\t" 
	        << "="                << "\t" << start + alignment1.ref_begin << "\t" << TLEN2 << "\t";

    data->outFile << seq2.substr(0,len_seq).c_str(); 
    data->outFile << "\t" << quality2.substr(0,len_seq).c_str() << "\t" << "PG:Z:MarkDuplicates" << "\t" << "RG:Z:AA1" << endl;

 }

 }

 } else {
	kv->add(key,keybytes,value,valuebytes);
 }

}



void map_print2samformat_removeSoftClip(uint64_t itask, char *key, int keybytes, char *value, int valuebytes, KeyValue *kv, void *ptr)
{
  Data *data = (Data *) ptr;
  PairLink *plink = (PairLink *) key;

 if( (valuebytes == sizeof(HeadSeq)*2) && (plink->chr1 > 0) && (plink->chr2 > 0) && (plink->chr1 == plink->chr2) ) {

  HeadSeq *tmp_hs = (HeadSeq *) value;
  HeadSeq head_seq1 = tmp_hs[0];
  HeadSeq head_seq2 = tmp_hs[1];

  int ID = head_seq1.type;
  if(ID < 0) ID = -1 * ID;
  int strand1 = head_seq1.strand;
  int strand2 = head_seq2.strand;
  int maq1 = head_seq1.maq;
  int maq2 = head_seq2.maq;

  bool check_len = false;
  int start, len;

  if(plink->loc1 < plink->loc2) {
        len = plink->loc2 - plink->loc1 + 1;
	start = plink->loc1;
        if( (plink->loc1+2*len) > data->reference_sequence[plink->chr1].length() ) {
                len = data->reference_sequence[plink->chr1].length() - plink->loc1 + 1;
        } else  {
                len = 2*len;  check_len = true;
        }
  } else  {
        len = plink->loc1 - plink->loc2 + 1;
	start = plink->loc2;
        if( (plink->loc2+2*len) > data->reference_sequence[plink->chr1].length() ) {
                len = data->reference_sequence[plink->chr1].length() - plink->loc2 + 1;
        } else {
                len = 2*len; check_len = true;
        }
  }

if(check_len && (start > 1)) {

  string ref_sequence = data->reference_sequence[plink->chr1].substr(start-1,len);
  string seq1, seq2, rev_seq1, rev_seq2, quality1, quality2;

  int sallow = data->range - (int)( (float)data->range * data->soft_allow);
  int sallow_upper = (int)( (float)data->range * 0.02 );

  seq1.assign(head_seq1.seq,head_seq1.length);
  seq2.assign(head_seq2.seq,head_seq2.length);

  quality1.assign(head_seq1.quality,head_seq1.length);
  quality2.assign(head_seq2.quality,head_seq2.length);

  StripedSmithWaterman::Aligner aligner1(1,4,6,1);
  StripedSmithWaterman::Filter filter1;
  StripedSmithWaterman::Alignment alignment1;

  StripedSmithWaterman::Aligner aligner2(1,4,6,1);
  StripedSmithWaterman::Filter filter2;
  StripedSmithWaterman::Alignment alignment2;

  alignment1.Clear();
  aligner1.Align(seq1.substr(0,head_seq1.length).c_str(), ref_sequence.c_str(), ref_sequence.length(), filter1, &alignment1);

  alignment2.Clear();
  aligner2.Align(seq2.substr(0,head_seq2.length).c_str(), ref_sequence.c_str(), ref_sequence.length(), filter2, &alignment2);



  vector<CigarOp_Function> cigar1;
  CigarOp_Function::parse(alignment1.cigar_string.c_str() ,cigar1);

  int SOFT_S1 = 0;
  int SOFT_E1 = 0;
  if(cigar1[0].op == 'S')               SOFT_S1 += cigar1[0].size;
  if(cigar1[cigar1.size()-1].op == 'S') SOFT_E1 += cigar1[cigar1.size()-1].size;

  vector<CigarOp_Function> cigar2;
  CigarOp_Function::parse(alignment2.cigar_string.c_str() ,cigar2);

  int SOFT_S2 = 0;
  int SOFT_E2 = 0;
  if(cigar2[0].op == 'S')               SOFT_S2 += cigar2[0].size;
  if(cigar2[cigar2.size()-1].op == 'S') SOFT_E2 += cigar2[cigar2.size()-1].size;

  if(SOFT_S1 <= sallow_upper) alignment1.query_begin = 0;
  if(SOFT_E1 <= sallow_upper) alignment1.query_end = (int32_t)data->range - 1;
 
  if(SOFT_S2 <= sallow_upper) alignment2.query_begin = 0;
  if(SOFT_E2 <= sallow_upper) alignment2.query_end = (int32_t)data->range - 1;


 

  int len_seq1 = alignment1.query_end - alignment1.query_begin + 1;
  int len_seq2 = alignment2.query_end - alignment2.query_begin + 1;

  seq1.assign( seq1.substr(alignment1.query_begin,len_seq1).c_str(), len_seq1);
  seq2.assign( seq2.substr(alignment2.query_begin,len_seq2).c_str(), len_seq2); 

  quality1.assign( quality1.substr(alignment1.query_begin,len_seq1).c_str(), len_seq1);
  quality2.assign( quality2.substr(alignment2.query_begin,len_seq2).c_str(), len_seq2);


 if( (len_seq1 >= sallow) && (len_seq2 >= sallow) &&
     (!alignment1.cigar_string.empty()) && (!alignment2.cigar_string.empty()) &&
     ((int)alignment1.mismatches < (int)(len_seq1 * 0.05)) && ((int)alignment2.mismatches < (int)(len_seq2 * 0.05)) 
   ) {

  alignment1.Clear();
  aligner1.Align(seq1.substr(0,len_seq1).c_str(), ref_sequence.c_str(), ref_sequence.length(), filter1, &alignment1);

  alignment2.Clear();
  aligner2.Align(seq2.substr(0,len_seq2).c_str(), ref_sequence.c_str(), ref_sequence.length(), filter2, &alignment2);

  int flag1 = 1 + 2;
  int flag2 = 1 + 2;

  if(head_seq1.PE == 1) flag1 += 64;
  else                  flag1 += 128;

  if(head_seq2.PE == 1) flag2 += 64;
  else                  flag2 += 128;

  if(head_seq1.strand < 0) flag1 += 16;
  else                     flag1 += 32;

  if(head_seq2.strand < 0) flag2 += 16;
  else                     flag2 += 32;

  stringstream head;
  head << "PEC:" << plink->chr1 << ":" << plink->loc1 << ":" << plink->loc2 << ":" << strand1 << ":" << strand2 << ":" << ID;

  int32_t TLEN1, TLEN2;

  if( alignment1.ref_begin <= alignment2.ref_begin ) { TLEN1 = alignment2.ref_end - alignment1.ref_begin + 1; TLEN2 = alignment2.ref_end - alignment1.ref_begin + 1; }
  else                                               { TLEN1 = alignment1.ref_end - alignment2.ref_begin + 1; TLEN2 = alignment1.ref_end - alignment2.ref_begin + 1; }

  if(strand1 < 0) TLEN1 = -1 * TLEN1;
  else            TLEN2 = -1 * TLEN2;

  int align_len1 = parse_cigar_align_length(const_cast<char*>(alignment1.cigar_string.c_str()));
  int align_len2 = parse_cigar_align_length(const_cast<char*>(alignment2.cigar_string.c_str()));

  if( (plink->chr1 == plink->chr2) && (strand1*strand2 == -1) && (!seq1.empty()) && (!seq2.empty()) ) {

      stringstream cigar_out1;
      stringstream MD1;
      vector<CigarOp_Function> cigar1;
      CigarOp_Function::parse(alignment1.cigar_string.c_str(),cigar1);
      int id = 0;
      int nm1 = 0;
      for(int i=0; i<cigar1.size(); i++) {
		if(cigar1[i].op == 'M')      { MD1 << cigar1[i].size;                                                cigar_out1 << cigar1[i].size << "M"; id += cigar1[i].size;}
		else if(cigar1[i].op == 'X') { MD1 << seq1.substr(id,cigar1[i].size).c_str(); nm1 += cigar1[i].size; cigar_out1 << cigar1[i].size << "X"; id += cigar1[i].size; }
		else if(cigar1[i].op == 'I') { cigar_out1 << cigar1[i].size << "I"; id += cigar1[i].size; }
		else if(cigar1[i].op == 'D') { cigar_out1 << cigar1[i].size << "D";  }
		else if(cigar1[i].op == 'S') { cigar_out1 << cigar1[i].size << "S"; id += cigar1[i].size; }
      }

      stringstream cigar_out2;
      stringstream MD2;
      vector<CigarOp_Function> cigar2;
      CigarOp_Function::parse(alignment2.cigar_string.c_str(),cigar2);
      id = 0;
      int nm2 = 0;
      for(int i=0; i<cigar2.size(); i++) {
                if(cigar2[i].op == 'M')      { MD2 << cigar2[i].size;                                                cigar_out2 << cigar2[i].size << "M"; id += cigar2[i].size;}
                else if(cigar2[i].op == 'X') { MD2 << seq2.substr(id,cigar2[i].size).c_str(); nm2 += cigar2[i].size; cigar_out2 << cigar2[i].size << "X"; id += cigar2[i].size;}
                else if(cigar2[i].op == 'I') { cigar_out2 << cigar2[i].size << "I"; id += cigar2[i].size; }
                else if(cigar2[i].op == 'D') { cigar_out2 << cigar2[i].size << "D";  }
                else if(cigar2[i].op == 'S') { cigar_out2 << cigar2[i].size << "S"; id += cigar2[i].size; }
      }

      data->outFile << head.str().c_str() << "\t" << flag1 << "\t";

      if((plink->chr1 >= 1) && (plink->chr1 <= 22)) data->outFile << plink->chr1;
      else if(plink->chr1 == 23)                  data->outFile << "X";
      else if(plink->chr1 == 24)                  data->outFile << "Y";
      else if(plink->chr1 == 25)                  data->outFile << "MT";

      data->outFile << "\t" << start + alignment1.ref_begin << "\t"
                    << maq1 << "\t" << cigar_out1.str().c_str() << "\t"
                    << "="  << "\t" << start + alignment2.ref_begin << "\t" << TLEN1 << "\t";

      data->outFile << seq1.substr(0,len_seq1).c_str();
      data->outFile << "\t" << quality1.substr(0,len_seq1).c_str() << "\t";
      data->outFile << "MC:Z:" << cigar_out2.str().c_str() << "\t" << "MD:Z:" << MD1.str().c_str() << "\t" << "RG:Z:AA1" << "\t" << "NM:i:" << nm1 << endl;


      data->outFile << head.str().c_str() << "\t" << flag2 << "\t";

      if((plink->chr2 >= 1) && (plink->chr2 <= 22)) data->outFile << plink->chr2;
      else if(plink->chr2 == 23)                  data->outFile << "X";
      else if(plink->chr2 == 24)                  data->outFile << "Y";
      else if(plink->chr2 == 25)                  data->outFile << "MT";

      data->outFile << "\t" << start + alignment2.ref_begin << "\t"
                    << maq2 << "\t" << cigar_out2.str().c_str() << "\t"
                    << "="  << "\t" << start + alignment1.ref_begin << "\t" << TLEN2 << "\t";

      data->outFile << seq2.substr(0,len_seq2).c_str();
      data->outFile << "\t" << quality2.substr(0,len_seq2).c_str() << "\t";
      data->outFile << "MC:Z:" << cigar_out1.str().c_str() << "\t" << "MD:Z:" << MD2.str().c_str() << "\t" << "RG:Z:AA1" << "\t" << "NM:i:" << nm2 << endl;

   }

  }

 }

 } else {
        kv->add(key,keybytes,value,valuebytes);
 }

}


void reduce_base_frequency_HiDup(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr)
{
  Data *data = (Data *) ptr;
  PointLoc *ploc = (PointLoc *) key;

  unsigned long long point_target = (unsigned long long)ploc->chr * (unsigned long long)1e11 + (unsigned long long)ploc->loc;
  if(  data->_loc_freq_table.count( point_target ) > 0 ) {
     unsigned long long tloc = data->_loc_freq_table[ point_target ];

     char *value;
     uint64_t nvalues_total;
     CHECK_FOR_BLOCKS(multivalue,valuebytes,nvalues,nvalues_total)
     BEGIN_BLOCK_LOOP(multivalue,valuebytes,nvalues)

     value = multivalue;
     for(int i=0; i<nvalues; i++) {
       PointData *pData = (PointData *) value;
       int intBase = get_int_from_base( pData->base );
       if( (intBase >= 0) && (intBase <= 3) ) data->table_base[ tloc*4 + (unsigned long long)intBase ]++;
       value += valuebytes[i];
     }

     END_BLOCK_LOOP
  }
}


void reduce_summation_base_frequency(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr)
{
  TableLoc *_table_loc = (TableLoc *) ptr;

  PointLoc *ploc = (PointLoc *) key;

  unsigned long long point_target = (unsigned long long)ploc->chr * (unsigned long long)1e11 + (unsigned long long)ploc->loc;
  if(  _table_loc->_loc_freq_table.count( point_target ) > 0 ) {
     unsigned long long tloc = _table_loc->_loc_freq_table[ point_target ];

     char *value;
     uint64_t nvalues_total;
     CHECK_FOR_BLOCKS(multivalue,valuebytes,nvalues,nvalues_total)
     BEGIN_BLOCK_LOOP(multivalue,valuebytes,nvalues)

     value = multivalue;
     for(int i=0; i<nvalues; i++) {
       PointData *pData = (PointData *) value;
       int intBase = get_int_from_base( pData->base );

       if( (intBase >= 0) && (intBase <= 3) ) {
         if(pData->strand > 0) _table_loc->table_positive[ tloc*4 + (unsigned long long)intBase ]++;
         else                  _table_loc->table_negative[ tloc*4 + (unsigned long long)intBase ]++;
       }

       value += valuebytes[i];
     }

     END_BLOCK_LOOP
  }

}

void map_Base_Frequency_hGE(uint64_t itask, char *key, int keybytes, char *value, int valuebytes, KeyValue *kv, void *ptr)
{
  Data *data = (Data *) ptr;
  PairLink *plink = (PairLink *) key;

  if( (valuebytes == sizeof(HeadSeq)*2) && (plink->chr1 > 0) && (plink->chr2 > 0) && (plink->chr1 == plink->chr2) ) {

    HeadSeq *tmp_hs = (HeadSeq *) value;
    HeadSeq head_seq1 = tmp_hs[0];
    HeadSeq head_seq2 = tmp_hs[1];

    int len_cutoff = data->range  - (int)( (float)data->range * data->soft_allow );

    bool check_len = false;
    int start, len;

    if(plink->loc1 <= plink->loc2) {
        len = plink->loc2 - plink->loc1 + 1;
	start = plink->loc1;
        if( (plink->loc1+2*len) > data->reference_sequence[plink->chr1].length() ) {
                len = data->reference_sequence[plink->chr1].length() - plink->loc1 + 1;
        } else  {
                 check_len = true; len = 2*len;
        }
    } else  {
        len = plink->loc1 - plink->loc2 + 1;
	start = plink->loc2;
        if( (plink->loc2+2*len) > data->reference_sequence[plink->chr1].length() ) {
                len = data->reference_sequence[plink->chr1].length() - plink->loc2 + 1;
        } else {
                check_len = true; len = 2*len;
        }
    }

    if(check_len && (start > 1)) {

        string ref_sequence = data->reference_sequence[plink->chr1].substr(start-1,len);

        string seq1, seq2, rev_seq1, rev_seq2, quality1, quality2;

        seq1.assign(head_seq1.seq, head_seq1.length);
        seq2.assign(head_seq2.seq, head_seq2.length);

        quality1.assign(head_seq1.quality, head_seq1.length);
        quality2.assign(head_seq2.quality, head_seq2.length);

        StripedSmithWaterman::Aligner aligner(1,4,6,1);
        StripedSmithWaterman::Filter filter;
        StripedSmithWaterman::Alignment alignment1;
        StripedSmithWaterman::Alignment alignment2;

        map<unsigned long long,int> _point_base;

        PointLoc pLoc;
        PointData pData1, pData2;
        int32_t *bam_cigar = (int32_t *) malloc( sizeof(int32_t)*2 );
        int sallow = data->range - (int)( (float)data->range * data->soft_allow);
	int sallow_upper = (int)( (float)data->range * 0.02 );

        pLoc.chr = plink->chr1;
        pData1.strand = head_seq1.strand;

        alignment1.Clear();
        aligner.Align(seq1.substr(0,head_seq1.length).c_str(), ref_sequence.c_str(), ref_sequence.length(), filter, &alignment1);
        int len_seq1 = alignment1.query_end - alignment1.query_begin + 1;

        pLoc.chr = plink->chr2;
        pData2.strand = head_seq2.strand;

        alignment2.Clear();
        aligner.Align(seq2.substr(0,head_seq2.length).c_str(), ref_sequence.c_str(), ref_sequence.length(), filter, &alignment2);
        int len_seq2 = alignment2.query_end - alignment2.query_begin + 1;




        vector<CigarOp_Function> cigar1;
        CigarOp_Function::parse(alignment1.cigar_string.c_str() ,cigar1);

        int SOFT_S1 = 0;
        int SOFT_E1 = 0;
        if(cigar1[0].op == 'S')              SOFT_S1 += cigar1[0].size;
        if(cigar1[cigar1.size()-1].op == 'S') SOFT_E1 += cigar1[cigar1.size()-1].size;

        vector<CigarOp_Function> cigar2;
        CigarOp_Function::parse(alignment2.cigar_string.c_str() ,cigar2);

        int SOFT_S2 = 0;
        int SOFT_E2 = 0;
        if(cigar2[0].op == 'S')              SOFT_S2 += cigar2[0].size;
        if(cigar2[cigar2.size()-1].op == 'S') SOFT_E2 += cigar2[cigar2.size()-1].size;

        if( (len_seq1 >= sallow) && (len_seq2 >= sallow) &&
            ((int)alignment1.mismatches < (int)(len_seq1 * 0.05)) && ((int)alignment2.mismatches < (int)(len_seq2 * 0.05)) &&
            (head_seq1.maq >= data->MAQ) && (head_seq2.maq >= data->MAQ)
          ) {

        	int sp     = (int)alignment1.query_begin;
        	int sp_loc = (int)alignment1.ref_begin;

                if(SOFT_S1 <= sallow_upper) {
                   for(size_t j = 0; j < SOFT_S1; j++)  {
                        pLoc.loc = start + sp_loc - SOFT_S1 + j;
			pLoc.refbase = ref_sequence.c_str()[sp_loc-SOFT_S1+j];		

                        int qScore = get_Qscore( quality1.c_str()[j] );
                        pData1.base = seq1.c_str()[j];

                        unsigned long long point = (unsigned long long)pLoc.chr * (unsigned long long)1e11 + (unsigned long long)pLoc.loc;
                        if((qScore >= data->bQscore) && (_point_base.count(point) == 0)) {
                        	kv->add((char *) &pLoc, sizeof(PointLoc), (char *) &pData1, sizeof(PointData));
                                _point_base[point] = 1;
                        }

                   }
                }


        	for(size_t j = 0; j < alignment1.cigar.size(); j++)  {
          		convert_cigar2bam(alignment1.cigar[j], bam_cigar);

          		if( (bam_cigar[0] == 0) || (bam_cigar[0] == 3) || (bam_cigar[0] == 7) || (bam_cigar[0] == 8) ) {

                		for(int k=0; k<(int)bam_cigar[1]; k++) {
                   			pLoc.loc = start + sp_loc + k;
                   			pLoc.refbase = ref_sequence.substr(sp_loc,(int)bam_cigar[1]).c_str()[k];

                   			int qScore = get_Qscore( quality1.substr(sp,(int)bam_cigar[1]).c_str()[k] );
                   			pData1.base = seq1.substr(sp,(int)bam_cigar[1]).c_str()[k];

					unsigned long long point = (unsigned long long)pLoc.chr * (unsigned long long)1e11 + (unsigned long long)pLoc.loc;

                   			if((qScore >= data->bQscore) && (_point_base.count(point) == 0)) {
						kv->add((char *) &pLoc, sizeof(PointLoc), (char *) &pData1, sizeof(PointData));
						_point_base[point] = 1;
					}
                 		}

                 		sp     += (int)bam_cigar[1];
                 		sp_loc += (int)bam_cigar[1];

         		} else if( (bam_cigar[0] == 2) || (bam_cigar[0] == 6) )  {
                		sp_loc += (int)bam_cigar[1];
         		} else if( bam_cigar[0] == 1) {
                		sp     += (int)bam_cigar[1];
        		}

       		}


                if(SOFT_E1 <= sallow_upper) {
                   for(size_t j = 0; j < SOFT_E1; j++)  {
                        pLoc.loc = start + sp_loc + j;
                        pLoc.refbase = ref_sequence.c_str()[sp_loc+j];

                        int qScore = get_Qscore( quality1.c_str()[sp+j] );
                        pData1.base = seq1.c_str()[sp+j];

                        unsigned long long point = (unsigned long long)pLoc.chr * (unsigned long long)1e11 + (unsigned long long)pLoc.loc;
                        if((qScore >= data->bQscore) && (_point_base.count(point) == 0)) {
                                kv->add((char *) &pLoc, sizeof(PointLoc), (char *) &pData1, sizeof(PointData));
                                _point_base[point] = 1;
                        }

                   }
                }


       		sp     = (int)alignment2.query_begin;
       		sp_loc = (int)alignment2.ref_begin;

                if(SOFT_S2 <= sallow_upper) {
                   for(size_t j = 0; j < SOFT_S2; j++)  {
                        pLoc.loc = start + sp_loc - SOFT_S2 + j;
                        pLoc.refbase = ref_sequence.c_str()[sp_loc-SOFT_S2+j];

                        int qScore = get_Qscore( quality2.c_str()[j] );
                        pData2.base = seq2.c_str()[j];

                        unsigned long long point = (unsigned long long)pLoc.chr * (unsigned long long)1e11 + (unsigned long long)pLoc.loc;
                        if((qScore >= data->bQscore) && (_point_base.count(point) == 0)) {
                                kv->add((char *) &pLoc, sizeof(PointLoc), (char *) &pData2, sizeof(PointData));
                                _point_base[point] = 1;
                        }

                   }
                }


       		for(size_t j = 0; j < alignment2.cigar.size(); j++)  {
          		convert_cigar2bam(alignment2.cigar[j], bam_cigar);

         		if( (bam_cigar[0] == 0) || (bam_cigar[0] == 3) || (bam_cigar[0] == 7) || (bam_cigar[0] == 8) ) {

                		for(int k=0; k<(int)bam_cigar[1]; k++) {
                   			pLoc.loc = start + sp_loc + k;
                   			pLoc.refbase   = ref_sequence.substr(sp_loc,(int)bam_cigar[1]).c_str()[k];

                   			int qScore = get_Qscore( quality2.substr(sp,(int)bam_cigar[1]).c_str()[k] );
                   			pData2.base = seq2.substr(sp,(int)bam_cigar[1]).c_str()[k];

					unsigned long long point = (unsigned long long)pLoc.chr * (unsigned long long)1e11 + (unsigned long long)pLoc.loc;
	
                   			if((qScore >= data->bQscore) && (_point_base.count(point) == 0)) {
						 kv->add((char *) &pLoc, sizeof(PointLoc), (char *) &pData2, sizeof(PointData));
						 _point_base[point] = 1; 
					}
                		}

                		sp     += (int)bam_cigar[1];
                		sp_loc += (int)bam_cigar[1];

         		} else if( (bam_cigar[0] == 2) || (bam_cigar[0] == 6) )  {
                		sp_loc += (int)bam_cigar[1];
         		} else if( bam_cigar[0] == 1) {
                		sp     += (int)bam_cigar[1];
        		}

       		}

                if(SOFT_E2 <= sallow_upper) {
                   for(size_t j = 0; j < SOFT_E2; j++)  {
                        pLoc.loc = start + sp_loc + j;
                        pLoc.refbase = ref_sequence.c_str()[sp_loc+j];

                        int qScore = get_Qscore( quality2.c_str()[sp+j] );
                        pData2.base = seq2.c_str()[sp+j];

                        unsigned long long point = (unsigned long long)pLoc.chr * (unsigned long long)1e11 + (unsigned long long)pLoc.loc;
                        if((qScore >= data->bQscore) && (_point_base.count(point) == 0)) {
                                kv->add((char *) &pLoc, sizeof(PointLoc), (char *) &pData2, sizeof(PointData));
                                _point_base[point] = 1;
                        }

                   }
                }

	 }

	 free(bam_cigar);

   }

  }

}


void map_Base_Frequency(uint64_t itask, char *key, int keybytes, char *value, int valuebytes, KeyValue *kv, void *ptr)
{
  Data *data = (Data *) ptr;
  PairLink *plink = (PairLink *) key;

 if( (valuebytes == sizeof(HeadSeq)*2) && (plink->chr1 > 0) && (plink->chr2 > 0) && (plink->chr1 == plink->chr2) ) {

    HeadSeq *tmp_hs = (HeadSeq *) value;
    HeadSeq head_seq1 = tmp_hs[0];
    HeadSeq head_seq2 = tmp_hs[1];

    int len_cutoff = data->range  - (int)( (float)data->range * data->soft_allow );

//    plink->ID = head_seq1.type;
//    if(plink->ID < 0) plink->ID = -1 * plink->ID;
//    plink->strand1 = head_seq1.strand;
//    plink->strand2 = head_seq2.strand;
//    plink->maq1 = head_seq1.maq;
//    plink->maq2 = head_seq2.maq;

    bool check_len = false;
    int start, len;

    if(plink->loc1 <= plink->loc2) {
        len = plink->loc2 - plink->loc1 + 1;
	start = plink->loc1;
        if( (plink->loc1+2*len) > data->reference_sequence[plink->chr1].length() ) {
                len = data->reference_sequence[plink->chr1].length() - plink->loc1 + 1;
        } else  {
                 check_len = true; len = 2*len;
        }
    } else  {
        len = plink->loc1 - plink->loc2 + 1;
	start = plink->loc2;
        if( (plink->loc2+2*len) > data->reference_sequence[plink->chr1].length() ) {
                len = data->reference_sequence[plink->chr1].length() - plink->loc2 + 1;
        } else {
                check_len = true; len = 2*len;
        }
    }


   if(check_len && (start > 1)) {

        string ref_sequence = data->reference_sequence[plink->chr1].substr(start-1,len);

        string seq1, seq2, rev_seq1, rev_seq2, quality1, quality2;

        seq1.assign(head_seq1.seq, head_seq1.length);
        seq2.assign(head_seq2.seq, head_seq2.length);

        quality1.assign(head_seq1.quality, head_seq1.length);
        quality2.assign(head_seq2.quality, head_seq2.length);

        StripedSmithWaterman::Aligner aligner(1,4,6,1);
        StripedSmithWaterman::Filter filter;
        StripedSmithWaterman::Alignment alignment1;
        StripedSmithWaterman::Alignment alignment2;

        PointLoc pLoc;
        PointData pData1, pData2;
        int32_t *bam_cigar = (int32_t *) malloc( sizeof(int32_t)*2 );
	int sallow = data->range - (int)( (float)data->range * data->soft_allow);

        pLoc.chr = plink->chr1;
        pData1.strand = head_seq1.strand;

        alignment1.Clear();
        aligner.Align(seq1.substr(0,head_seq1.length).c_str(), ref_sequence.c_str(), ref_sequence.length(), filter, &alignment1);
	int len_seq1 = alignment1.query_end - alignment1.query_begin + 1;

        pLoc.chr = plink->chr2;
        pData2.strand = head_seq2.strand;

        alignment2.Clear();
        aligner.Align(seq2.substr(0,head_seq2.length).c_str(), ref_sequence.c_str(), ref_sequence.length(), filter, &alignment2);
        int len_seq2 = alignment2.query_end - alignment2.query_begin + 1;


if( (len_seq1 > sallow) && (len_seq2 > sallow) &&
    ((int)alignment1.mismatches < (int)(len_seq1 * 0.05)) && ((int)alignment2.mismatches < (int)(len_seq2 * 0.05)) &&
    (head_seq1.maq >= data->MAQ) && (head_seq2.maq >= data->MAQ)
    ) {

        int sp     = (int)alignment1.query_begin;
        int sp_loc = (int)alignment1.ref_begin;
        for(size_t j = 0; j < alignment1.cigar.size(); j++)  {
          convert_cigar2bam(alignment1.cigar[j], bam_cigar);

         if( (bam_cigar[0] == 0) || (bam_cigar[0] == 3) || (bam_cigar[0] == 7) || (bam_cigar[0] == 8) ) {

                for(int k=0; k<(int)bam_cigar[1]; k++) {
                   pLoc.loc = start + sp_loc + k;
                   pLoc.refbase = ref_sequence.substr(sp_loc,(int)bam_cigar[1]).c_str()[k];

                   int qScore = get_Qscore( quality1.substr(sp,(int)bam_cigar[1]).c_str()[k] );
                   pData1.base = seq1.substr(sp,(int)bam_cigar[1]).c_str()[k];

		   if(qScore >= data->bQscore)  kv->add((char *) &pLoc, sizeof(PointLoc), (char *) &pData1, sizeof(PointData));

                 }

                 sp     += (int)bam_cigar[1];
                 sp_loc += (int)bam_cigar[1];

         } else if( (bam_cigar[0] == 2) || (bam_cigar[0] == 6) )  {
                sp_loc += (int)bam_cigar[1];
         } else if( bam_cigar[0] == 1) {
/*
		for(int k=0; k<(int)bam_cigar[1]; k++) {
                   pLoc.loc = start + sp_loc;
                   pLoc.refbase   = ref_sequence.substr(sp_loc,(int)bam_cigar[1]).c_str()[0];

                   int qScore = get_Qscore( quality1.substr(sp,(int)bam_cigar[1]).c_str()[k] );
                   pData1.base = seq1.substr(sp,(int)bam_cigar[1]).c_str()[k];

                   if(qScore >= data->bQscore) kv->add((char *) &pLoc, sizeof(PointLoc), (char *) &pData1, sizeof(PointData));
		}
*/
                sp     += (int)bam_cigar[1];
        }

       }


       sp     = (int)alignment2.query_begin;
       sp_loc = (int)alignment2.ref_begin;
       for(size_t j = 0; j < alignment2.cigar.size(); j++)  {
          convert_cigar2bam(alignment2.cigar[j], bam_cigar);

         if( (bam_cigar[0] == 0) || (bam_cigar[0] == 3) || (bam_cigar[0] == 7) || (bam_cigar[0] == 8) ) {

                for(int k=0; k<(int)bam_cigar[1]; k++) {
                   pLoc.loc = start + sp_loc + k;
                   pLoc.refbase   = ref_sequence.substr(sp_loc,(int)bam_cigar[1]).c_str()[k];

                   int qScore = get_Qscore( quality2.substr(sp,(int)bam_cigar[1]).c_str()[k] );
                   pData2.base = seq2.substr(sp,(int)bam_cigar[1]).c_str()[k];

		   if(qScore >= data->bQscore)	kv->add((char *) &pLoc, sizeof(PointLoc), (char *) &pData2, sizeof(PointData));

                }

                sp     += (int)bam_cigar[1];
                sp_loc += (int)bam_cigar[1];

         } else if( (bam_cigar[0] == 2) || (bam_cigar[0] == 6) )  {
                sp_loc += (int)bam_cigar[1];
         } else if( bam_cigar[0] == 1) {
/*
		for(int k=0; k<(int)bam_cigar[1]; k++) {
		   pLoc.loc = start + sp_loc;
                   pLoc.refbase   = ref_sequence.substr(sp_loc,(int)bam_cigar[1]).c_str()[0];

                   int qScore = get_Qscore( quality2.substr(sp,(int)bam_cigar[1]).c_str()[k] );
                   pData2.base = seq2.substr(sp,(int)bam_cigar[1]).c_str()[k];

                   if(qScore >= data->bQscore) kv->add((char *) &pLoc, sizeof(PointLoc), (char *) &pData2, sizeof(PointData)); 
		}
*/
                sp     += (int)bam_cigar[1];
        }

       }

 }
       free(bam_cigar);

   }

 }

}


void reduce_check_DNACollision_wo_subClustering(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr)
{
   Data *data = (Data *) ptr;
   char *value;
   PairLink *plink = (PairLink *) key;

   int nPE = 0;

   int nA = 0;
   int nB = 0;
   int nC = 0;
   int nD = 0;

if(nvalues >= 3) {

   uint64_t nvalues_total;
   CHECK_FOR_BLOCKS(multivalue,valuebytes,nvalues,nvalues_total)
   BEGIN_BLOCK_LOOP(multivalue,valuebytes,nvalues)

    int ref_loc = -1;
    value = multivalue;
    for(int i=0; i<nvalues; i++) {
     if(valuebytes[i] == sizeof(RefConsen)) {
        ref_loc = i;
        break;
     }
     value += valuebytes[i];
    }

    if(ref_loc >= 0) {

      value = multivalue;
      for(int i=0; i<nvalues; i++) {
        if(valuebytes[i] == sizeof(HeadSeq)*2) { 
          HeadSeq *headseq_pair = (HeadSeq *) value;
          HeadSeq headseq = headseq_pair[0];

          if(headseq.header[0] == 'A') nA += 1;
          if(headseq.header[0] == 'B') nB += 1;
          if(headseq.header[0] == 'C') nC += 1;
          if(headseq.header[0] == 'D') nD += 1;

	  nPE += 1;
        }
        value += valuebytes[i];
      }

   } else {

      value = multivalue;
      for(int i=0; i<nvalues; i++) {
          HeadSeq *headseq_pair = (HeadSeq *) value;
          HeadSeq headseq = headseq_pair[0];

          if(headseq.header[0] == 'A') nA += 1;
          if(headseq.header[0] == 'B') nB += 1;
          if(headseq.header[0] == 'C') nC += 1;
          if(headseq.header[0] == 'D') nD += 1;

	  nPE += 1;
          value += valuebytes[i];
     }

   }

   END_BLOCK_LOOP

 }

      int check = 0;
      if(nA > 0) check += 1;
      if(nB > 0) check += 1;
      if(nC > 0) check += 1;
      if(nD > 0) check += 1;

      if(check > 1) {
        data->flank += 1;
      } else if(check == 1) {
        data->offside += 1;
      }

}


void reduce_check_DNACollision_MiniBarcode(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr)
{
   Data *data = (Data *) ptr;
   char *value;
   PairLink *plink = (PairLink *) key;

   int nA = 0;
   int nB = 0;
   int nC = 0;
   int nD = 0;

   uint64_t nvalues_total;
   CHECK_FOR_BLOCKS(multivalue,valuebytes,nvalues,nvalues_total)
   BEGIN_BLOCK_LOOP(multivalue,valuebytes,nvalues)

    value = multivalue;
    for(int i=0; i<nvalues; i++) {
      HeadSeq *headseq_pair = (HeadSeq *) value;
      HeadSeq headseq = headseq_pair[0];

      if(headseq.header[0] == 'A') nA += 1;
      if(headseq.header[0] == 'B') nB += 1;
      if(headseq.header[0] == 'C') nC += 1;
      if(headseq.header[0] == 'D') nD += 1;

      value += valuebytes[i];
    }

   END_BLOCK_LOOP

   int clen;
   if(plink->loc1 < plink->loc2) clen = plink->loc2 - plink->loc1 + 1;
   else				 clen = plink->loc1 - plink->loc2 + 1;

   data->outFile << plink->chr1 << " " << plink->chr2 << " " << plink->loc1 <<  " " << plink->loc2 << " " << plink->strand1 << " " << plink->strand2 << " " << plink->ID
		 << " " << nvalues << " " << clen << endl;

   int check = 0;
   if(nA > 0) check += 1;
   if(nB > 0) check += 1;
   if(nC > 0) check += 1;
   if(nD > 0) check += 1;

   if(check > 1) {
		  data->flank += 1;
   } else if(check == 1) {
		  data->offside += 1; 
   }

}

void revcomp_sequence_carray(char *seq, int len)
{
  string sequence;
  sequence.assign(seq, len);

  string rev_sequence = revcomp( sequence.c_str() );
  memcpy(seq, rev_sequence.c_str(),len); 
}

void revcomp_quality_carray(char *qual, int len)
{
   string quality;
   quality.assign(qual, len);

   string rev_quality = revquality( quality.c_str() ); 
   memcpy(qual, rev_quality.c_str(), len);
}


void reduce_LinkID_HeadSeq_single(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr)
{
   Data *data = (Data *) ptr;
   char *value;

   uint64_t nvalues_total;
   CHECK_FOR_BLOCKS(multivalue,valuebytes,nvalues,nvalues_total)
   BEGIN_BLOCK_LOOP(multivalue,valuebytes,nvalues)

   PairLink *locInfo;
   string header;
   bool check = false;
   bool check_single = false;

   value = multivalue;
   for(int i=0; i<nvalues; i++) {
        if(valuebytes[i] == sizeof(PairLink)) {
          locInfo = (PairLink *) value;

          check = 1;
          break;
        }
        value += valuebytes[i];
   }

   if(check && check_single) {

     HeadSeq headseq;

     int nH = 0;
     value = multivalue;
     for(int i=0; i<nvalues; i++) {

       if(valuebytes[i] != sizeof(PairLink)) {
         headseq = *(HeadSeq *) value;
       }
       value += valuebytes[i];

     }

     if(nH == 1) kv->add((char *) locInfo, sizeof(PairLink), (char *) &headseq, sizeof(HeadSeq) );

   }

   END_BLOCK_LOOP
}

