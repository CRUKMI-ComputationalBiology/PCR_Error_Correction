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

#include "mpi.h"
#include "zlib.h"

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "sys/stat.h"

#define __STDC_LIMIT_MACROS
#include "mapreduce.h"
#include "keyvalue.h"
#include "blockmacros.h"
#include "typedefs.h"

#include "fstream"
#include "iostream"
#include "sstream"

#include "vector"
#include "map"
#include "set"
#include "algorithm"
#include "time.h"
#include "math.h"

#include "Fasta_reader.hpp"
#include "sequenceUtil.hpp"

#include "IRKE.hpp"
#include "KmerCounter.hpp"
#include "stacktrace.hpp"
#include "argProcessor.hpp"

#include "mrSubroutine.h"

using namespace MAPREDUCE_NS;
using namespace std;

void Execute(const char * command) {
    int ret = system(command);
    if (ret != 0) {
        cout << "COMMAND: " << command << endl;
        cout << "Died with exit code " << ret << endl;
        cout << "Exiting." << endl;
        exit(-1);
    }
}

int main(int narg, char **args)
{
  MPI_Init(&narg,&args);

  Data data;

  MPI_Comm_rank(MPI_COMM_WORLD,&data.me);
  MPI_Comm_size(MPI_COMM_WORLD,&data.nprocs);

  data.kmer_length = 30;

  int page_size = 128;
  data.mismatch = 0;
  data.MAQ = 15;
  data.soft_allow = 0.1;

  data.On_Target = false;

  string SampleName, fname_read1, fname_read2, typeFile, fname_ref;

  int num_args = 1;
  try {
       ArgProcessor in_args(narg, args);

       if (in_args.isArgSet("--SoftClip")) {
            data.soft_allow = in_args.getFloatVal("--SoftClip");
            num_args += 2;
            if(data.me==0) cerr << "Allowed soft clip proportion set to: " << data.soft_allow << endl;
       }
	
       if (in_args.isArgSet("--K")) {
            data.kmer_length = in_args.getIntVal("--K");
            num_args += 2;
            if(data.me==0) cerr << "Kmer length set to: " << data.kmer_length << endl;
       }

       if(in_args.isArgSet("--PageSize")) {
            page_size = in_args.getIntVal("--PageSize");
            num_args += 2;
            if(data.me==0) cerr << "Page size for map reduce object set to: " << page_size << endl;
       }

       if(in_args.isArgSet("--MAQ")) {
            data.MAQ = in_args.getIntVal("--MAQ");
            num_args += 2;
            if(data.me==0) cerr << "Mapping Quality set to: " << data.MAQ << endl;
       } 

       if(in_args.isArgSet("--OVERLAP")) {
            data.OVERLAP = in_args.getIntVal("--OVERLAP");
            num_args += 2;
            if(data.me==0) cerr << "Seq overlapping set to: " << data.OVERLAP << endl;
       }

       if(in_args.isArgSet("--Target")) {
            data.filename = in_args.getStringVal("--Target");
            num_args += 2;
            data.On_Target = true;
            if(data.me==0) cerr << "Target Region file was given as: " << data.filename << endl;
       } 

       if(in_args.isArgSet("--fnameR1")) {
            fname_read1 = in_args.getStringVal("--fnameR1");
            num_args += 2;
            if(data.me==0) cerr << "filename for error corrected reads1 was given as: " << fname_read1 << endl;
       }

       if(in_args.isArgSet("--fnameR2")) {
            fname_read2 = in_args.getStringVal("--fnameR2");
            num_args += 2;
            if(data.me==0) cerr << "filename for error corrected reads2 was given as: " << fname_read2 << endl;
       }  

       if(in_args.isArgSet("--Ref")) {
            fname_ref = in_args.getStringVal("--Ref");
            num_args += 2;
            if(data.me==0) cerr << "filename for HG Reference was given as: " << fname_ref << endl;
       }

       if(in_args.isArgSet("--type")) {
            typeFile = in_args.getStringVal("--type");
            num_args += 2;
       }

       if(in_args.isArgSet("--SM")) {
            SampleName = in_args.getStringVal("--SM");
            num_args += 2; 
      }
      

  }

  catch(exception& e) {
        cerr << "error: " << e.what() << "\n";
        return 1;
  }

// ################################################
//

 data.MAX_RECURSION = 1;
 data.MIN_SEED_ENTROPY = 1.5;
 data.MIN_SEED_COVERAGE = 1;

 data.DS = false;

 data.PACMAN = false;
 data.CRAWL = false;
 data.crawl_length = 1;
 data.min_any_entropy = 0.4;

 data.prune_error_kmers = true;
 data.min_ratio_non_error = 0.05f;
 data.min_kmer_count = 1;

 data.MIN_CONNECTIVITY_RATIO = 0.0;
 data.MIN_ASSEMBLY_LENGTH = 151;
 data.MIN_ASSEMBLY_COVERAGE = 1;


  MapReduce *mrLink = new MapReduce(MPI_COMM_WORLD);
  mrLink->memsize = page_size;
  mrLink->verbosity = 1;
  mrLink->timer = 1;

  MapReduce *mrSeq = new MapReduce(MPI_COMM_WORLD);
  mrSeq->memsize = page_size;
  mrSeq->verbosity = 1;
  mrSeq->timer = 1;

  MapReduce *mrBuffer = new MapReduce(MPI_COMM_WORLD);
  mrBuffer->memsize = page_size;
  mrBuffer->verbosity = 1;
  mrBuffer->timer = 1;

//################################################
//

 data.min_csize = 10;
 data.range = 151;
 data.flank = 200;

// data.error_rate = 4.0;

 data.On_Alignment = true;
 data.On_Assemble = false;

// data.On_Target = true;  // OnTarget or OffTarget harvesting switch

 data.On_smallCluster = true;

 data.On_Clean = true;  // Switch for error cleaning on input reads

 bool On_Error_Clean = false;  // Error Cleaning switch
 bool On_Clustering = true;  // Clustering switch

// ##############################################

while( data.On_Target ) {

  int loc_vector = 0;

  FILE* FS;
  FS = fopen(data.filename.c_str(), "r");

  while (!feof(FS)) {
    char temp_string[300];
    fgets(temp_string, 300, FS);

    char dummy[5];
    char gene[15];
    ExonLoc _loc_info;

    int n = sscanf(temp_string, "%s %llu %llu %s", dummy, &_loc_info.loc1, &_loc_info.loc2, gene);
    if(n < 4) break;

    _loc_info.chr = 0;
    if(      strcmp(dummy, "chr1") == 0) _loc_info.chr = 1;
    else if( strcmp(dummy, "chr2") == 0) _loc_info.chr = 2;
    else if( strcmp(dummy, "chr3") == 0) _loc_info.chr = 3;
    else if( strcmp(dummy, "chr4") == 0) _loc_info.chr = 4;
    else if( strcmp(dummy, "chr5") == 0) _loc_info.chr = 5;
    else if( strcmp(dummy, "chr6") == 0) _loc_info.chr = 6;
    else if( strcmp(dummy, "chr7") == 0) _loc_info.chr = 7;
    else if( strcmp(dummy, "chr8") == 0) _loc_info.chr = 8;
    else if( strcmp(dummy, "chr9") == 0) _loc_info.chr = 9;
    else if( strcmp(dummy,"chr10") == 0) _loc_info.chr = 10;
    else if( strcmp(dummy,"chr11") == 0) _loc_info.chr = 11;
    else if( strcmp(dummy,"chr12") == 0) _loc_info.chr = 12;
    else if( strcmp(dummy,"chr13") == 0) _loc_info.chr = 13;
    else if( strcmp(dummy,"chr14") == 0) _loc_info.chr = 14;
    else if( strcmp(dummy,"chr15") == 0) _loc_info.chr = 15;
    else if( strcmp(dummy,"chr16") == 0) _loc_info.chr = 16;
    else if( strcmp(dummy,"chr17") == 0) _loc_info.chr = 17;
    else if( strcmp(dummy,"chr18") == 0) _loc_info.chr = 18;
    else if( strcmp(dummy,"chr19") == 0) _loc_info.chr = 19;
    else if( strcmp(dummy,"chr20") == 0) _loc_info.chr = 20;
    else if( strcmp(dummy,"chr21") == 0) _loc_info.chr = 21;
    else if( strcmp(dummy,"chr22") == 0) _loc_info.chr = 22;
    else if( strcmp(dummy, "chrX") == 0) _loc_info.chr = 23;
    else if( strcmp(dummy, "chrY") == 0) _loc_info.chr = 24;
    else if( strcmp(dummy, "chrM") == 0) _loc_info.chr = 25;

    if( (_loc_info.chr >= 1) && (_loc_info.chr <= 25) ) {
        string gname(gene);
        gname = remove_whitespace(gname);
        data._target_gene.push_back(gname);
        data._target_loc.push_back(_loc_info);

        loc_vector++;
    }
  }

  fclose(FS);
  break;
}

// ###############################################

 double tstart = MPI_Wtime();

 if(data.me == 0) cerr << endl << " Uploading  Reference Seq..... " << endl << endl;
 data.filename.assign(fname_ref);
 mrBuffer->map(data.nprocs,map_upload_reference,&data);

 data.reference_sequence.resize(27);
 mrBuffer->map(mrBuffer,map_assign_reference,&data);

 for(int i=0; i<25; i++) {

	int tmp_root;
	if(data.reference_sequence[i+1].empty()) tmp_root = 0;
	else	   				 tmp_root = data.me;

	int root;
	MPI_Allreduce(&tmp_root, &root, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	if(data.reference_sequence[i+1].empty()) tmp_root = 0;
	else	                                 tmp_root = data.reference_sequence[i+1].length();

	int seq_len;
	MPI_Allreduce(&tmp_root, &seq_len, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	
	char *tmp_seq = (char *) malloc( sizeof(char)*seq_len );
	if(!data.reference_sequence[i+1].empty()) memcpy( tmp_seq, data.reference_sequence[i+1].c_str(), seq_len );

	MPI_Bcast( tmp_seq, seq_len, MPI_CHAR, root, MPI_COMM_WORLD);

	if(data.reference_sequence[i+1].empty()) {
		string tmp_str;
		tmp_str.assign( tmp_seq, seq_len );
		data.reference_sequence.at(i+1) = tmp_str;
	}

	free(tmp_seq);
 }

 double tstop = MPI_Wtime();

 if(data.me == 0) cerr << endl << " Time took for uploading ref_seq : " << tstop - tstart << endl << endl;

// ###############################################

 
 if(data.me == 0) cerr << endl << " Uploading  input reads info..... " << endl << endl; 
 int nseqs;
 if(      strcmp(typeFile.substr(0,2).c_str(),"fa") == 0) { data.ftype = false; nseqs = mrSeq->map(narg-num_args,&args[num_args],0,1,0,fileread_RNAseq_HeadSeq_FASTA,&data); }
 else if( strcmp(typeFile.substr(0,2).c_str(),"fq") == 0) { data.ftype = true;  nseqs = mrSeq->map(narg-num_args,&args[num_args],0,1,0,fileread_RNAseq_HeadSeq_FASTQ,&data); }
 else {
    printf("ERROR: query file format wrong! It should be either 'fq' or 'fa' format. \n");
    MPI_Abort(MPI_COMM_WORLD,1);
 }

 int *sub_range = (int *)malloc(sizeof(int) * data.nprocs);
 MPI_Allgather(&data.range, 1, MPI_INT, sub_range, 1, MPI_INT, MPI_COMM_WORLD);

 for(int i=0; i<data.nprocs; i++) {
	if(data.range > sub_range[i]) data.range = sub_range[i];
 }

 if(data.me == 0) cerr << endl << " input seq length =  " << data.range << endl << endl;

// ############################################################################################

 if(data.me == 0) cerr << endl << " Extracting  read alignment info..... " << endl << endl;

 int nLink = mrLink->map(narg-num_args,&args[num_args],0,1,0,fileread_AlignmentInfo,&data);

 mrLink->collate(NULL);
 mrLink->reduce(reduce_Link_from_Header,&data);

 mrLink->add(mrSeq);
 mrLink->collate(NULL);

 mrLink->reduce(reduce_LinkID_HeadSeq,&data);


// ###########################################################################
//   Print consensus sequences in SAM format

 if(data.me == 0) cerr << endl << endl << "Extracting PE singleton reads....  " <<  endl << endl; 
 mrSeq = mrLink->copy();
 mrSeq->collate(NULL);
 
 mrSeq->reduce(reduce_consensus_SeqQual_singleton,&data);

 mrBuffer = mrSeq->copy();

// if(data.me == 0) cerr << endl << endl << "Extracting SE singleton reads....  " <<  endl << endl; 
// mrSeq = mrLink->copy();
// mrSeq->collate(NULL);
// mrSeq->reduce(reduce_consensus_Single_singleton,&data);

// mrBuffer->add(mrSeq);

// if(data.me == 0) cerr << endl << endl << "Extracting SE consensus....  " <<  endl << endl;
// mrSeq = mrLink->copy();
// mrSeq->collate(NULL);
// mrSeq->reduce(reduce_consensus_Single_highDup,&data);

// mrBuffer->add(mrSeq);




 if(data.me == 0) cerr << endl << endl << "Extracting PE consensus...  " <<  endl << endl;
 mrLink->collate(NULL);
 mrLink->reduce(reduce_consensus_SeqQual_highDup,&data);

 mrLink->add(mrBuffer);


 if(data.me == 0) cerr << endl << endl << "Printing result  " <<  endl;

 stringstream out_filename;
 out_filename << "cluster_" << data.me << ".fa";
 data.outFile.open(out_filename.str().c_str()); 

 if(data.me == 0) {

   for(int i=1; i<=22; i++) 
	data.outFile << "@SQ" << "\t" << "SN:" << i << "\t" << "LN:" << data.reference_sequence[i].length() << endl; 

   data.outFile << "@SQ" << "\t" << "SN:X"  << "\t" << "LN:" << data.reference_sequence[23].length() << endl;
   data.outFile << "@SQ" << "\t" << "SN:Y"  << "\t" << "LN:" << data.reference_sequence[24].length() << endl;
   data.outFile << "@SQ" << "\t" << "SN:MT" << "\t" << "LN:" << data.reference_sequence[25].length() << endl;

   data.outFile << "@RG" << "\t" << "ID:AA1" << "\t" << "SM:" << SampleName.c_str() << "\t" << "PL:ILLUMINA" << "\t" << "LB:ECPipeline" << "\t" << "PU:S2" << endl; 
   data.outFile << "@PG" << "\t" << "ID:ECPipeline" << "\t" <<  "PN:CEP" << "\t" <<  "VN:0.0.1" << "\t" << "CL:PCRError" << endl;

 }

 if(data.me == 0) cerr  << "  PE consensus reads...  "  << endl;
 mrLink->map(mrLink, map_print2samformat_removeSoftClip, &data);

 data.outFile.close();

// ###########################################################################

 MPI_Barrier(MPI_COMM_WORLD);

 delete mrSeq;
 delete mrLink;
 delete mrBuffer;

 MPI_Finalize();

 exit(0);

}

