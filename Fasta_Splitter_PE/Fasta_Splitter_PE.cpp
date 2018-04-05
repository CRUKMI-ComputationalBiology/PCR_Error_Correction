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


#include <mpi.h>

#include <stdio.h>
#include <iostream>
#include <string>
#include <cstring>
#include <sstream>
#include <algorithm>
#include <map>
#include <fstream>

#include "Fasta_reader.hpp"
#include "Fastq_reader.hpp"
#include "argProcessor.hpp"

#ifndef BIG_FILE_DEFINES_H
#define BIG_FILE_DEFINES_H

#ifdef __linux
     #ifndef _LARGEFILE64_SOURCE
     #define _LARGEFILE64_SOURCE
     #endif
     #define _FILE_OFFSET_BITS 64
     #define STAT_NAME stat64
#else
     #define STAT_NAME stat
#endif

#endif //BIG_FILE_DEFINES_H

void print_split_file(std::stringstream& out_filename, string readFile, string typeFile, int rank, int numranks) {

  ofstream outFile;
  outFile.open(out_filename.str().c_str());

  if( strcmp(typeFile.substr(0,2).c_str(),"fq") == 0) {

    Fastq_reader fastq_reader(readFile);
    int line_num = 0;
    while (true) {
        Fastq_entry fq = fastq_reader.getNext();
        string seq     = fq.get_sequence();
        string quality = fq.get_quality();
        string header  = fq.get_header();
        if (seq == "") break;

        if( (line_num % numranks) == rank )
          outFile << header.c_str() << endl << seq.c_str() << endl << "+" << endl << quality.c_str() << endl;

        line_num++;
    }

  } else if( strcmp(typeFile.substr(0,2).c_str(),"fa") == 0) {

    Fasta_reader fasta_reader(readFile);
    int line_num = 0;
    while (true) {
        Fasta_entry fa = fasta_reader.getNext();
        string seq     = fa.get_sequence();
        string header  = fa.get_header();
        if (seq == "") break;

        if( (line_num % numranks) == rank )
          outFile << header.c_str() << endl << seq.c_str() << endl;

        line_num++;
    }

  } else {
        if(rank == 0) cerr << "Wrong input read file format: it should be fa or fq!!!" << endl;
        exit(EXIT_FAILURE);
  }

  outFile.close();
}



//static unsigned long long MAX_CLUSTER_SIZE = 5000;

using namespace std;

int main (int argc, char* argv[]) {

   int rank,numranks;
   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD, &numranks);
   MPI_Comm_rank(MPI_COMM_WORLD,&rank);

   string readFile;
   string sKmerDir, typeFile;
   int PE, lane_index;

  try {
       ArgProcessor in_args(argc, argv);

       if(in_args.isArgSet("-r")) {
            readFile = in_args.getStringVal("-r");
       }

       if(in_args.isArgSet("-l")) {
            lane_index = in_args.getIntVal("-l");
       }

       if(in_args.isArgSet("-o")) {
           sKmerDir = in_args.getStringVal("-o");
       }

       if(in_args.isArgSet("-t")) {
           typeFile = in_args.getStringVal("-t");
       }

       if(in_args.isArgSet("-pe")) {
            PE = in_args.getIntVal("-pe");
       }

  }

  catch(exception& e) {
        cerr << "error: " << e.what() << "\n";
        return 1;
  }


   MPI_Barrier(MPI_COMM_WORLD);

   if(!readFile.empty()) {

       stringstream out_filename;
       out_filename << sKmerDir << "/sKmer_" << lane_index;

       if(PE == 1)      out_filename << "_1_";
       else if(PE == 2) out_filename << "_2_";
       else if(PE == 0) out_filename << "_";
       out_filename << rank << "." << typeFile;

       print_split_file(out_filename, readFile, typeFile,rank, numranks);
       MPI_Barrier(MPI_COMM_WORLD);

   }


   MPI_Finalize();
   exit(0);

}


