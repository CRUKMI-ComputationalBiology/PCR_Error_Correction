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


#include <stdio.h>
#include <iostream>
#include <string>
#include <cstring>
#include <sstream>
#include <algorithm>
#include <map>
#include <fstream>
#include "argProcessor.hpp"
//#include "mutil.h"
#include <math.h>

#include<mpi.h>

using namespace std;

bool fexists(const char *filename)
{
  ifstream ifile(filename);
  if(ifile) { return true;  }
  else      { return false; }
}

void Execute_bowtie(const char * command) {
    int ret = system(command);
    if (ret != 0) {
        cout << "COMMAND: " << command << endl;
        cout << "Died with exit code " << ret << endl;
        cout << "Exiting." << endl;
        exit(-1);
    }
}

int main(int argc,char** argv)
{

  int rank,numranks, left, right;
  int buffer_send[1], buffer_receive[1];
  MPI_Request request;
  MPI_Status status;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numranks);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  string readString, bwaDIR, typeFile, refFile;
  int num_threads;
  int MAQ;
  int PE, lane_index;

  int num_args = 1;
  try {
       ArgProcessor in_args(argc, argv);

       if(in_args.isArgSet("-t")) {
            num_threads = in_args.getIntVal("-t");
            num_args += 2;
       }

       if(in_args.isArgSet("-o")) {
            bwaDIR = in_args.getStringVal("-o");
            num_args += 2;
       }

       if(in_args.isArgSet("-q")) {
            MAQ = in_args.getIntVal("-q");
            num_args += 2;
       }
	
       if(in_args.isArgSet("-f")) {
            typeFile = in_args.getStringVal("-f");
            num_args += 2;
       }

       if(in_args.isArgSet("-r")) {
            refFile = in_args.getStringVal("-r");
            num_args += 2;
       }	

       if(in_args.isArgSet("-l")) {
            lane_index = in_args.getIntVal("-l");
       }

      if(in_args.isArgSet("-pe")) {
	    PE = in_args.getIntVal("-pe");
	    num_args += 2;
      }

  }

  catch(exception& e) {
        cerr << "error: " << e.what() << "\n";
        return 1;
  }

  stringstream cmdstr;
  stringstream fstr;

  MPI_Barrier(MPI_COMM_WORLD);

if(PE == 1) {

  fstr.str("");
  fstr << bwaDIR << "/sKmer_" << lane_index << "_1_" << rank << "." << typeFile;
  if(fexists(fstr.str().c_str())) {

     cmdstr.str("");
     cmdstr << "bwa mem -t " << num_threads << " -v 0 -K 100000 -M -S " << refFile;
     cmdstr << " " << bwaDIR << "/sKmer_"<< lane_index << "_1_" << rank << "." << typeFile << " " << bwaDIR << "/sKmer_" << lane_index << "_2_" << rank << "." << typeFile;
     cmdstr << " | grep -v \"^@\" | awk '{print $1 \"\t\" $3 \"\t\" $4 \"\t\" $5 \"\t\" $6 \"\t\" $2 \"\t\" 1; }' >> " << bwaDIR << "/bfiltered_" << rank << ".out";  

     cerr << "CMD on rank " << rank << " : " << cmdstr.str() << endl;
     system(cmdstr.str().c_str());  

     MPI_Barrier(MPI_COMM_WORLD);
  } else {
	cerr << "error: splitted input read files don't exist!" << "\n";
        return 1;
  }


} else {
 
  fstr.str("");
  fstr << bwaDIR << "/sKmer_" << lane_index << "_" << rank << "." << typeFile;
  if(fexists(fstr.str().c_str())) {
 
     cmdstr.str("");
     cmdstr << "bwa mem -K 100000 -v 0 -t "  << num_threads;
     cmdstr << " -M " << refFile;
     cmdstr << " " << bwaDIR << "/sKmer_" << lane_index << "_" << rank << "." << typeFile;
     cmdstr << " | grep -v \"^@\" | awk '{print $1 \"\t\" $3 \"\t\" $4 \"\t\" $5 \"\t\" $6 \"\t\" $2 \"\t\" 1; }' >> " << bwaDIR << "/bfiltered_" << rank << ".out";

     cerr << "CMD on rank " << rank << " : " << cmdstr.str() << endl;
     system(cmdstr.str().c_str());

     MPI_Barrier(MPI_COMM_WORLD);

  } else {
        cerr << "error: splitted input read files don't exist!" << "\n";
        return 1;
  }

  fstr.str("");
  fstr << bwaDIR << "/sKmer_2_" << rank << "." << typeFile;
  if(fexists(fstr.str().c_str())) {

     cmdstr.str("");
     cmdstr << "bwa mem -K 100000 -v 0 -t "  << num_threads;
     cmdstr << " -M " << refFile;
     cmdstr << " " << bwaDIR << "/sKmer_2_" << rank << "." << typeFile;
     cmdstr << " | grep -v \"^@\" | awk '{print $1 \"\t\" $3 \"\t\" $4 \"\t\" $5 \"\t\" $6 \"\t\" $2 \"\t\" 2; }' >> " << bwaDIR << "/bfiltered_" << rank << ".out";

     cerr << "CMD on rank " << rank << " : " << cmdstr.str() << endl;
     system(cmdstr.str().c_str());

     MPI_Barrier(MPI_COMM_WORLD);

 }


  fstr.str("");
  fstr << bwaDIR << "/sKmer_3_" << rank << "." << typeFile;
  if(fexists(fstr.str().c_str())) {

     cmdstr.str("");
     cmdstr << "bwa mem -K 100000 -v 0 -t "  << num_threads;
     cmdstr << " -M " << refFile;
     cmdstr << " " << bwaDIR << "/sKmer_3_" << rank << "." << typeFile;
     cmdstr << " | grep -v \"^@\" | awk '{print $1 \"\t\" $3 \"\t\" $4 \"\t\" $5 \"\t\" $6 \"\t\" $2 \"\t\" 3; }' >> " << bwaDIR << "/bfiltered_" << rank << ".out";

     cerr << "CMD on rank " << rank << " : " << cmdstr.str() << endl;
     system(cmdstr.str().c_str());

     MPI_Barrier(MPI_COMM_WORLD);
  }

  fstr.str("");
  fstr << bwaDIR << "/sKmer_4_" << rank << "." << typeFile;
  if(fexists(fstr.str().c_str())) { 

     cmdstr.str("");
     cmdstr << "bwa mem -K 100000 -v 0 -t "  << num_threads;
     cmdstr << " -M " << refFile;
     cmdstr << " " << bwaDIR << "/sKmer_4_" << rank << "." << typeFile;
     cmdstr << " | grep -v \"^@\" | awk '{print $1 \"\t\" $3 \"\t\" $4 \"\t\" $5 \"\t\" $6 \"\t\" $2 \"\t\" 4; }' >> " << bwaDIR << "/bfiltered_" << rank << ".out";

     cerr << "CMD on rank " << rank << " : " << cmdstr.str() << endl;
     system(cmdstr.str().c_str());  

     MPI_Barrier(MPI_COMM_WORLD);

  }

}

  MPI_Finalize();
  exit(0);

}
