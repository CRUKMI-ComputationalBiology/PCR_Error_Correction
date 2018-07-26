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
  data.MAQ = 10;
  data.bQscore = 20;
  data.soft_allow = 0.0;

  data.On_Target = true;

  data.On_Polishing = false;

  string SampleName, freq_filename, typeFile, fname_ref, scratch_dir;

  data.nKmer = 1;
  int num_lane = 1;
  int num_args = 1;
  try {
       ArgProcessor in_args(narg, args);

       if (in_args.isArgSet("--OnPolishing")) {
	  int dummy_onpolish = in_args.getIntVal("--OnPolishing"); 
	  num_args += 2;
	  if(dummy_onpolish == 1) {
		data.On_Polishing = false;
		if(data.me==1) cerr << "Error-Polishing mode was set to OFF " << endl;
	  } else if(dummy_onpolish == 2) {
		data.On_Polishing = true;
		if(data.me==0) cerr << "Error-Polishing mode was set to ON " << endl;
	  } else {
		if(data.me==0) cerr << "Error-Polishing mode was set to OFF by-default " << endl;
	  }
       }

       if (in_args.isArgSet("--nLane")) {
            num_lane = in_args.getIntVal("--nLane");
            num_args += 2;
            if(data.me==0) cerr << "Number of NGS lanes: " << num_lane << endl;
       }

       if (in_args.isArgSet("--nKmer")) {
            data.nKmer = in_args.getIntVal("--nKmer");
            num_args += 2;
            if(data.me==0) cerr << "eKmer number allowance: " << data.nKmer << endl;
       }

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

       if(in_args.isArgSet("--BQScore")) {
            data.bQscore = in_args.getIntVal("--BQScore");
            num_args += 2;
            if(data.me==0) cerr << "Phred Score set to: " << data.bQscore << endl;
       } 

       if(in_args.isArgSet("--Target")) {
            data.filename = in_args.getStringVal("--Target");
            num_args += 2;
            data.On_Target = true;
            if(data.me==0) cerr << "Target Region file was given as: " << data.filename << endl;
       } 

       if(in_args.isArgSet("--FreqTable")) {
            freq_filename = in_args.getStringVal("--FreqTable");
            num_args += 2;
            if(data.me==0) cerr << "filename for freq Table was given as: " << freq_filename << endl;
       }

       if(in_args.isArgSet("--Ref")) {
            fname_ref = in_args.getStringVal("--Ref");
            num_args += 2;
            if(data.me==0) cerr << "filename for HG Reference was given as: " << fname_ref << endl;
       }

       if(in_args.isArgSet("--Scratch")) {
            scratch_dir = in_args.getStringVal("--Scratch");
            num_args += 2;
            if(data.me==0) cerr << "filename for scratch dir was given as: " << scratch_dir << endl;
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

 data.OVERLAP = data.kmer_length;

 data.MAX_RECURSION = 1;
 data.MIN_SEED_ENTROPY = 1.5;
 data.MIN_SEED_COVERAGE = 2;

 data.DS = false;

 data.PACMAN = false;
 data.CRAWL = false;
 data.crawl_length = 1;
 data.min_any_entropy = 0.0;

 data.prune_error_kmers = true;
 data.min_ratio_non_error = 0.05f;
 data.min_kmer_count = 1;

 data.MIN_CONNECTIVITY_RATIO = 0.0;
 data.MIN_ASSEMBLY_LENGTH = 1; //data.kmer_length;
 data.MIN_ASSEMBLY_COVERAGE = 3;

// ##################################################
//

  data.lane = num_lane;

  MapReduce *mrLink = new MapReduce(MPI_COMM_WORLD);
  mrLink->memsize = page_size;
  mrLink->verbosity = 1;
  mrLink->timer = 1;
  mrLink->set_fpath(scratch_dir.c_str());
 
  MapReduce *mrSeq = new MapReduce(MPI_COMM_WORLD);
  mrSeq->memsize = page_size;
  mrSeq->verbosity = 1;
  mrSeq->timer = 1;
  mrSeq->set_fpath(scratch_dir.c_str());

  MapReduce *mrBuffer = new MapReduce(MPI_COMM_WORLD);
  mrBuffer->memsize = page_size;
  mrBuffer->verbosity = 1;
  mrBuffer->timer = 1;
  mrBuffer->set_fpath(scratch_dir.c_str());

  MapReduce *mrEkmer = new MapReduce(MPI_COMM_WORLD);
  mrEkmer->memsize = page_size;
  mrEkmer->verbosity = 1;
  mrEkmer->timer = 1;
  mrEkmer->set_fpath(scratch_dir.c_str());


//################################################
//

 data.min_csize = 5;
 data.flank = 200;

 data.On_Alignment = true;
 data.On_Assemble = false;

 data.On_smallCluster = true;

 data.On_Clean = true;  // Switch for error cleaning on input reads

 bool On_Error_Clean = false;  // Error Cleaning switch
 bool On_Clustering = true;  // Clustering switch

// ##############################################

data.loc_vector = 0;
while( data.On_Target ) {

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
    if(      (strcmp(dummy, "chr1") == 0) || (strcmp(dummy, "1") == 0) ) _loc_info.chr = 1;
    else if( (strcmp(dummy, "chr2") == 0) || (strcmp(dummy, "2") == 0) ) _loc_info.chr = 2;
    else if( (strcmp(dummy, "chr3") == 0) || (strcmp(dummy, "3") == 0) ) _loc_info.chr = 3;
    else if( (strcmp(dummy, "chr4") == 0) || (strcmp(dummy, "4") == 0) ) _loc_info.chr = 4;
    else if( (strcmp(dummy, "chr5") == 0) || (strcmp(dummy, "5") == 0) ) _loc_info.chr = 5;
    else if( (strcmp(dummy, "chr6") == 0) || (strcmp(dummy, "6") == 0) ) _loc_info.chr = 6;
    else if( (strcmp(dummy, "chr7") == 0) || (strcmp(dummy, "7") == 0) ) _loc_info.chr = 7;
    else if( (strcmp(dummy, "chr8") == 0) || (strcmp(dummy, "8") == 0) ) _loc_info.chr = 8;
    else if( (strcmp(dummy, "chr9") == 0) || (strcmp(dummy, "9") == 0) ) _loc_info.chr = 9;
    else if( (strcmp(dummy,"chr10") == 0) || (strcmp(dummy, "10") == 0) ) _loc_info.chr = 10;
    else if( (strcmp(dummy,"chr11") == 0) || (strcmp(dummy, "11") == 0) ) _loc_info.chr = 11;
    else if( (strcmp(dummy,"chr12") == 0) || (strcmp(dummy, "12") == 0) ) _loc_info.chr = 12;
    else if( (strcmp(dummy,"chr13") == 0) || (strcmp(dummy, "13") == 0) ) _loc_info.chr = 13;
    else if( (strcmp(dummy,"chr14") == 0) || (strcmp(dummy, "14") == 0) ) _loc_info.chr = 14;
    else if( (strcmp(dummy,"chr15") == 0) || (strcmp(dummy, "15") == 0) ) _loc_info.chr = 15;
    else if( (strcmp(dummy,"chr16") == 0) || (strcmp(dummy, "16") == 0) ) _loc_info.chr = 16;
    else if( (strcmp(dummy,"chr17") == 0) || (strcmp(dummy, "17") == 0) ) _loc_info.chr = 17;
    else if( (strcmp(dummy,"chr18") == 0) || (strcmp(dummy, "18") == 0) ) _loc_info.chr = 18;
    else if( (strcmp(dummy,"chr19") == 0) || (strcmp(dummy, "19") == 0) ) _loc_info.chr = 19;
    else if( (strcmp(dummy,"chr20") == 0) || (strcmp(dummy, "20") == 0) ) _loc_info.chr = 20;
    else if( (strcmp(dummy,"chr21") == 0) || (strcmp(dummy, "21") == 0) ) _loc_info.chr = 21;
    else if( (strcmp(dummy,"chr22") == 0) || (strcmp(dummy, "22") == 0) ) _loc_info.chr = 22;
    else if( (strcmp(dummy, "chrX") == 0) || (strcmp(dummy, "X") == 0) )     _loc_info.chr = 23;
    else if( (strcmp(dummy, "chrY") == 0) || (strcmp(dummy, "Y") == 0) )     _loc_info.chr = 24;
    else if( (strcmp(dummy, "chrM") == 0) || (strcmp(dummy, "M") == 0) )     _loc_info.chr = 25;

    if( (_loc_info.chr >= 1) && (_loc_info.chr <= 25) ) {
        string gname(gene);
        gname = remove_whitespace(gname);
        data._target_gene.push_back(gname);
        data._target_loc.push_back(_loc_info);

	data.loc_vector++;
    }
  }

  fclose(FS);
  break;
}

  data.loc_vector--;
  if(data.me == 0) cerr << endl << " The number of target exons =  \t" << data.loc_vector << endl << endl;

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

 MPI_Barrier(MPI_COMM_WORLD);
 double tstop = MPI_Wtime();

 if(data.me == 0) cerr << endl << " Time took for uploading ref_seq : " << tstop - tstart << endl << endl;

// ###############################################
//
//

 if(data.me == 0) cerr << endl << " Uploading  input reads info..... " << endl << endl; 
 int nseqs = 0;
 if(      strcmp(typeFile.substr(0,2).c_str(),"fa") == 0) { 
	data.ftype = false; 
	nseqs = mrSeq->map(narg-num_args,&args[num_args],0,1,0,fileread_RNAseq_HeadSeq_FASTA,&data); 

 } else if( strcmp(typeFile.substr(0,2).c_str(),"fq") == 0) { 
	data.ftype = true;  
//	nseqs = mrSeq->map(narg-num_args,&args[num_args],0,1,0,fileread_RNAseq_HeadSeq_FASTQ,&data);

        data.filename = scratch_dir;
        for(int i=1; i<=num_lane; i++) {
           data.lane = i;

           data.read_side = 1;
           int tmp_nseqs = mrBuffer->map(data.nprocs, fileread_Sequence_Fastq, &data);
           mrSeq->add(mrBuffer);
           nseqs += tmp_nseqs;

           data.read_side = 2;
           tmp_nseqs = mrBuffer->map(data.nprocs, fileread_Sequence_Fastq, &data);
           mrSeq->add(mrBuffer);
           nseqs += tmp_nseqs;
        }

 } else {
    printf("ERROR: query file format wrong! It should be either 'fq' or 'fa' format. \n");
    MPI_Abort(MPI_COMM_WORLD,1);
 }


 if(data.me == 0) cerr << endl << " input seq length =  " << data.range << endl << endl;

// ############################################################################################

 if(data.me == 0) cerr << endl << " Extracting  read alignment info..... " << endl << endl;

// int nLink = mrLink->map(narg-num_args,&args[num_args],0,1,0,fileread_AlignmentInfo,&data);
 int nLink = mrLink->map(data.nprocs, fileread_AlignmentInfo_MPI,&data);


 mrLink->collate(NULL);
 mrLink->reduce(reduce_Link_from_Header,&data);

 mrLink->add(mrSeq);
 mrLink->collate(NULL);

 uint64_t num_all_PEs = mrLink->reduce(reduce_LinkID_HeadSeq,&data);

// #################################################################################

 int QS_range = 43;

 data.ATGC1 = (unsigned long long *) malloc(sizeof(unsigned long long)*4);
 data.ATGC2 = (unsigned long long *) malloc(sizeof(unsigned long long)*4);

// ################################################################################

 if(data.me == 0) cerr << endl << " ############################################# " << endl << endl;
 if(data.me == 0) cerr << endl << " Extracting  Error kmers ..... " << endl << endl;

 for(int k=0; k<4; k++) { data.ATGC1[k] = 0; data.ATGC2[k] = 0; }
  
 mrEkmer = mrLink->copy();
 mrEkmer->collate(NULL);
 mrEkmer->reduce(reduce_Kmer_Cov_All,&data);

 uint64_t num_all_kmers  = mrEkmer->collate(NULL);
 data.ftype = 2;
 uint64_t num_all_ekmers = mrEkmer->reduce(reduce_error_kmers,&data);

 mrBuffer = mrLink->copy();
 mrBuffer->collate(NULL);
 mrBuffer->reduce(reduce_extract_kmers_for_error_counts, &data);

 if(data.me == 0) cerr << endl << " Adding error kmers with Read-Kmer pairs..... " << endl << endl;
 mrBuffer->add(mrEkmer);

 mrBuffer->collate(NULL);
 mrBuffer->reduce(reduce_error_onAmplicons,&data); 

 if(data.me == 0) cerr << endl << " Error counting with Read-Kmer pairs......... " << endl << endl;

 mrSeq = mrLink->copy(); 
 mrSeq->add(mrBuffer);

 mrSeq->collate(NULL);

 mrSeq->reduce(reduce_error_matrices_step1,&data);

 mrSeq ->collate(NULL);

 mrSeq->reduce(reduce_error_matrices_step2,&data);

 mrSeq->collate(NULL);

// ################################################################################
//

 unsigned long long  *atgc_all = (unsigned long long *) malloc( sizeof(unsigned long long)*4 );
  
 for(int k=0; k<4; k++) atgc_all[k] = 0;   
 MPI_Allreduce(data.ATGC1, atgc_all, 4, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
 for(int k=0; k<4; k++) data.ATGC1[k] = atgc_all[k];

 for(int k=0; k<4; k++) atgc_all[k] = 0;
 MPI_Allreduce(data.ATGC2, atgc_all, 4, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
 for(int k=0; k<4; k++) data.ATGC2[k] = atgc_all[k];

 free(atgc_all);

 data.qs1_count    = (int ***) malloc(  QS_range*sizeof(int **));
 data.qs2_count    = (int ***) malloc(  QS_range*sizeof(int **));

 for(int i=0; i<QS_range; i++) {
    data.qs1_count[i] = (int **) malloc(4*sizeof(int *)); 
    data.qs2_count[i] = (int **) malloc(4*sizeof(int *));

    for(int j=0; j<4; j++) {
	data.qs1_count[i][j] = (int *) malloc(4*sizeof(int));
	data.qs2_count[i][j] = (int *) malloc(4*sizeof(int));
    }
 }

 for(int i=0; i<QS_range; i++)
   for(int j=0; j<4; j++)
      for(int k=0; k<4; k++) {
	data.qs1_count[i][j][k] = 0;
	data.qs2_count[i][j][k] = 0;
      }

// ##############################################################################

 mrSeq->reduce(reduce_ntnt_count, &data);

 int *ncount = (int *) malloc(sizeof(int)*4);
 int *gcount = (int *) malloc(sizeof(int)*4); 

 for(int i=0; i<QS_range; i++)
   for(int j=0; j<4; j++) {

	for(int k=0; k<4; k++) { ncount[k] = data.qs1_count[i][j][k]; gcount[k] = 0; }
	MPI_Allreduce(ncount, gcount, 4, MPI_INT, MPI_SUM, MPI_COMM_WORLD); 
	for(int k=0; k<4; k++) data.qs1_count[i][j][k] = gcount[k];

        for(int k=0; k<4; k++) { ncount[k] = data.qs2_count[i][j][k]; gcount[k] = 0; }
        MPI_Allreduce(ncount, gcount, 4, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        for(int k=0; k<4; k++) data.qs2_count[i][j][k] = gcount[k];

 }

 free(ncount);
 free(gcount);

// ###########################################################################
//   Contructing consensus sequences 


  if(data.me == 0) cerr << endl << endl << "Extracting PE HighDup consensus...  " <<  endl << endl;

                              mrBuffer = mrLink->copy();
  uint64_t num_cfDNAs_input = mrBuffer->collate(NULL);
  uint64_t num_highdup      = mrBuffer->reduce(reduce_consensus_highDup,&data);


  if(data.me == 0) cerr << endl << endl << "Extracting PE singleton reads....  " <<  endl << endl;

                           mrLink->add(mrBuffer);
                           mrLink->collate(NULL);
                           mrLink->reduce(reduce_sort_singleton,&data);

                           mrLink->collate(NULL);
  uint64_t num_singleton = mrLink->reduce(reduce_consensus_singleton,&data);

// #########################################################################################
//

  if(data.me == 0) cerr << endl << endl << "Creating Lookup table for HighDup reads....  " <<  endl << endl;

  unsigned long long nloc = 0;
  for(int i=0; i<data.loc_vector; i++) {
        for(int j=data._target_loc[i].loc1; j<=data._target_loc[i].loc2; j++)  {
          unsigned long long point_target = (unsigned long long)data._target_loc[i].chr * (unsigned long long)1e11 + (unsigned long long)j;
          data._loc_freq_table[ point_target ] = nloc;
          nloc++;
      }
  }

  data.table_base = (unsigned long long *)malloc(sizeof(unsigned long long)*(nloc*4));
  for(unsigned long long i=0; i<(nloc*4); i++) data.table_base[i] = 0;

  mrSeq->map(mrBuffer, map_Base_Frequency_hGE, &data);
  mrSeq->collate(NULL);
  mrSeq->reduce(reduce_base_frequency_HiDup, &data);

  unsigned long long *tmp_point = (unsigned long long *)malloc(sizeof(unsigned long long)*(nloc*4));
  for(unsigned long long i=0; i<(nloc*4); i++) tmp_point[i] = 0; 
  MPI_Reduce(data.table_base, tmp_point, nloc*4, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  for(unsigned long long i=0; i<(nloc*4); i++) data.table_base[i] = tmp_point[i];
  free(tmp_point);


 data.eBase1 = (unsigned long long *) malloc( sizeof(unsigned long long)*4 );
 data.sBase1 = (unsigned long long *) malloc( sizeof(unsigned long long)*4 );
 data.eBase2 = (unsigned long long *) malloc( sizeof(unsigned long long)*4 );
 data.sBase2 = (unsigned long long *) malloc( sizeof(unsigned long long)*4 );
 for(int i=0; i<4; i++) { data.eBase1[i] = 0; data.eBase2[i] = 0; data.sBase1[i] = 0; data.sBase2[i] = 0; }

 if(data.me == 0) cerr << endl << endl << "Suppressing error on PE singleton reads....  " <<  endl << endl;

 mrSeq->map(mrLink,map_correct_singleton_step1,&data);

 mrSeq->add(mrEkmer);

 mrSeq->collate(NULL);

 mrSeq->reduce(reduce_correct_single_step2,&data);

 mrLink->add(mrSeq);

 mrLink->collate(NULL);

 mrLink->reduce(reduce_correct_single_step3,&data);

 free(data.table_base);

 atgc_all = (unsigned long long *) malloc( sizeof(unsigned long long)*4 );

 for(int i=0; i<4; i++) atgc_all[i] = 0; 
 MPI_Allreduce(data.eBase1, atgc_all, 4, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
 for(int i=0; i<4; i++) data.eBase1[i] = atgc_all[i];

 for(int i=0; i<4; i++) atgc_all[i] = 0;
 MPI_Allreduce(data.eBase2, atgc_all, 4, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
 for(int i=0; i<4; i++) data.eBase2[i] = atgc_all[i];

 for(int i=0; i<4; i++) atgc_all[i] = 0;
 MPI_Allreduce(data.sBase1, atgc_all, 4, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
 for(int i=0; i<4; i++) data.sBase1[i] = atgc_all[i];

 for(int i=0; i<4; i++) atgc_all[i] = 0;
 MPI_Allreduce(data.sBase2, atgc_all, 4, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
 for(int i=0; i<4; i++) data.sBase2[i] = atgc_all[i];

 free(atgc_all);

 uint64_t num_final_cfDNAs = mrLink->add(mrBuffer);


// ###############################################################################
//  Base frequency table

if( !freq_filename.empty() ) {

 
 if(data.me == 0) cerr << endl << endl << "Extracting base frequencies for target region... " << endl << endl;

 if(data.me == 0) {

    ofstream statFile;

    std::size_t loc = freq_filename.find(".txt.");
    string stat_filename = "Stats_BaseCalls.txt.";
    stat_filename.append(freq_filename,loc+5, freq_filename.length()-loc+5);    

    statFile.open(stat_filename.c_str());

    statFile << "no all PE reads = "    << num_all_PEs      << endl;
    statFile << "no. total k-mers = "   << num_all_kmers    << endl;
    statFile << "no. error k-mers = "   << num_all_ekmers   << endl;
    statFile << "no cfDNAs = "          << num_cfDNAs_input << endl;
    statFile << "no cfDNAs highdup = "  << num_highdup      << endl;
    statFile << "no cfDNAs single = "   << num_singleton    << endl;
    statFile << "no cfDNAs final = "    << num_final_cfDNAs << endl;


    int all_error_1 = 0;
    int all_error_2 = 0;

    for(int k=0; k<QS_range; k++)
        for(int i=0; i<4; i++)
                for(int j=0; j<4; j++) { all_error_1 += data.qs1_count[k][i][j]; all_error_2 += data.qs2_count[k][i][j]; }

    int sequencing_error_1 = 0;
    int sequencing_error_2 = 0;

    for(int k=0; k<data.bQscore; k++)
        for(int i=0; i<4; i++)
                for(int j=0; j<4; j++) { sequencing_error_1 += data.qs1_count[k][i][j]; sequencing_error_2 += data.qs2_count[k][i][j]; }

    int contamination_error_1 = 0;
    int contamination_error_2 = 0;

    for(int k=data.bQscore; k<QS_range; k++)
        for(int i=0; i<4; i++)
                for(int j=0; j<4; j++) { contamination_error_1 += data.qs1_count[k][i][j]; contamination_error_2 += data.qs2_count[k][i][j]; }


    int e_ga1 = 0; 
    int e_ga2 = 0;
    for(int k=data.bQscore; k<QS_range; k++) { e_ga1 += data.qs1_count[k][0][1]; e_ga2 += data.qs2_count[k][0][1]; }

    int e_gt1 = 0;
    int e_gt2 = 0;
    for(int k=data.bQscore; k<QS_range; k++) { e_gt1 += data.qs1_count[k][0][2]; e_gt2 += data.qs2_count[k][0][2]; }

    int e_gc1 = 0;
    int e_gc2 = 0;
    for(int k=data.bQscore; k<QS_range; k++) { e_gc1 += data.qs1_count[k][0][3]; e_gc2 += data.qs2_count[k][0][3]; }


    int e_ag1 = 0;
    int e_ag2 = 0;
    for(int k=data.bQscore; k<QS_range; k++) { e_ag1 += data.qs1_count[k][1][0]; e_ag2 += data.qs2_count[k][1][0]; }

    int e_at1 = 0;
    int e_at2 = 0;
    for(int k=data.bQscore; k<QS_range; k++) { e_at1 += data.qs1_count[k][1][2]; e_at2 += data.qs2_count[k][1][2]; }

    int e_ac1 = 0;
    int e_ac2 = 0;
    for(int k=data.bQscore; k<QS_range; k++) { e_ac1 += data.qs1_count[k][1][3]; e_ac2 += data.qs2_count[k][1][3]; }


    int e_tg1 = 0;
    int e_tg2 = 0;
    for(int k=data.bQscore; k<QS_range; k++) { e_tg1 += data.qs1_count[k][2][0]; e_tg2 += data.qs2_count[k][2][0]; }

    int e_ta1 = 0;
    int e_ta2 = 0;
    for(int k=data.bQscore; k<QS_range; k++) { e_ta1 += data.qs1_count[k][2][1]; e_ta2 += data.qs2_count[k][2][1]; }

    int e_tc1 = 0;
    int e_tc2 = 0;
    for(int k=data.bQscore; k<QS_range; k++) { e_tc1 += data.qs1_count[k][2][3]; e_tc2 += data.qs2_count[k][2][3]; }


    int e_cg1 = 0;
    int e_cg2 = 0;
    for(int k=data.bQscore; k<QS_range; k++) { e_cg1 += data.qs1_count[k][3][0]; e_cg2 += data.qs2_count[k][3][0]; }

    int e_ca1 = 0;
    int e_ca2 = 0;
    for(int k=data.bQscore; k<QS_range; k++) { e_ca1 += data.qs1_count[k][3][1]; e_ca2 += data.qs2_count[k][3][1]; }

    int e_ct1 = 0;
    int e_ct2 = 0;
    for(int k=data.bQscore; k<QS_range; k++) { e_ct1 += data.qs1_count[k][3][2]; e_ct2 += data.qs2_count[k][3][2]; }


    unsigned long long noBase_1 = 0; 
    unsigned long long noBase_2 = 0;
    for(int k=0; k<4; k++) { noBase_1 += data.ATGC1[k]; noBase_2 += data.ATGC2[k]; }

    statFile << endl << endl;

    double error_rate = (double)(all_error_1+all_error_2)/(double)(noBase_1 + noBase_2);
    double error_barcode = (double)(num_all_PEs * 12 * 2) * error_rate;

    statFile << "All Error Rate...... " << endl;
    statFile << "Error Rate on HiDup reads (+) " << (double)all_error_1/(double)noBase_1 << endl;
    statFile << "Error Rate on HiDup reads (-) " << (double)all_error_2/(double)noBase_2 << endl;
    statFile << "Error Rate on HiDup           " << error_rate << endl;


    statFile << endl << endl;

    error_rate = (double)(sequencing_error_1 + sequencing_error_2)/(double)(noBase_1 + noBase_2);

    statFile << "Low base quality Error Rate...... " << endl;
    statFile << "Error Rate on HiDup reads (+) " << (double)sequencing_error_1/(double)noBase_1 << endl;
    statFile << "Error Rate on HiDup reads (-) " << (double)sequencing_error_2/(double)noBase_2 << endl;
    statFile << "Error Rate on HiDup           " << error_rate << endl;


    statFile << endl << endl;

    error_rate = (double)(contamination_error_1 + contamination_error_2)/(double)(noBase_1 + noBase_2);

    statFile << "High base quality Error Rate...... " << endl;
    statFile << "Error Rate on HiDup reads (+) " << (double)contamination_error_1/(double)noBase_1 << endl;
    statFile << "Error Rate on HiDup reads (-) " << (double)contamination_error_2/(double)noBase_2 << endl;
    statFile << "Error Rate on HiDup           " << error_rate << endl;




    statFile << endl << endl;
    statFile << "G>A " << (double)e_ga1/(double)noBase_1 << "\t" << (double)e_ga2/(double)noBase_2 << endl;
    statFile << "G>T " << (double)e_gt1/(double)noBase_1 << "\t" << (double)e_gt2/(double)noBase_2 << endl;
    statFile << "G>C " << (double)e_gc1/(double)noBase_1 << "\t" << (double)e_gc2/(double)noBase_2 << endl;

    statFile << "A>G " << (double)e_ag1/(double)noBase_1 << "\t" << (double)e_ag2/(double)noBase_2 << endl;
    statFile << "A>T " << (double)e_at1/(double)noBase_1 << "\t" << (double)e_at2/(double)noBase_2 << endl;
    statFile << "A>C " << (double)e_ac1/(double)noBase_1 << "\t" << (double)e_ac2/(double)noBase_2 << endl;

    statFile << "T>G " << (double)e_tg1/(double)noBase_1 << "\t" << (double)e_tg2/(double)noBase_2 << endl;
    statFile << "T>A " << (double)e_ta1/(double)noBase_1 << "\t" << (double)e_ta2/(double)noBase_2 << endl;
    statFile << "T>C " << (double)e_tc1/(double)noBase_1 << "\t" << (double)e_tc2/(double)noBase_2 << endl;

    statFile << "C>G " << (double)e_cg1/(double)noBase_1 << "\t" << (double)e_cg2/(double)noBase_2 << endl;
    statFile << "C>A " << (double)e_ca1/(double)noBase_1 << "\t" << (double)e_ca2/(double)noBase_2 << endl;
    statFile << "C>T " << (double)e_ct1/(double)noBase_1 << "\t" << (double)e_ct2/(double)noBase_2 << endl;

    statFile << endl;
    statFile << "# Guanine  on high_dup reads (+/-):\t" << data.ATGC1[0] << "\t" << data.ATGC2[0] << endl;
    statFile << "# Adenine  on high_dup reads (+/-):\t" << data.ATGC1[1] << "\t" << data.ATGC2[1] << endl;
    statFile << "# Thymine  on high_dup reads (+/-):\t" << data.ATGC1[2] << "\t" << data.ATGC2[2] << endl;
    statFile << "# Cytosine on high_dup reads (+/-):\t" << data.ATGC1[3] << "\t" << data.ATGC2[3] << endl;

    unsigned long long GC = data.ATGC1[0] + data.ATGC1[3];
    unsigned long long tB = 0; for(int k=0; k<4; k++) tB += data.ATGC1[k];
    statFile << endl;
    statFile << "GC content on high_dup reads (+) %:\t" << (double)GC/(double)tB << endl;

    GC = data.ATGC2[0] + data.ATGC2[3];
    tB = 0; for(int k=0; k<4; k++) tB += data.ATGC2[k];
    statFile << "GC content on high_dup reads (-) %:\t" << (double)GC/(double)tB << endl;

    statFile << endl;
    statFile << "# Guanine  Error removed on singletons (+/-):\t" << data.eBase1[0] << "\t" << data.eBase2[0] << endl;
    statFile << "# Adenine  Error removed on singletons (+/-):\t" << data.eBase1[1] << "\t" << data.eBase2[1] << endl;
    statFile << "# Thymine  Error removed on singletons (+/-):\t" << data.eBase1[2] << "\t" << data.eBase2[2] << endl;
    statFile << "# Cytosine Error removed on singletons (+/-):\t" << data.eBase1[3] << "\t" << data.eBase2[3] << endl;

    statFile << endl;
    statFile << "# Guanine  bases on singletons (+/-):\t" << data.sBase1[0] << "\t" << data.sBase2[0] << endl;
    statFile << "# Adenine  bases on singletons (+/-):\t" << data.sBase1[1] << "\t" << data.sBase2[1] << endl;
    statFile << "# Thymine  bases on singletons (+/-):\t" << data.sBase1[2] << "\t" << data.sBase2[2] << endl;
    statFile << "# Cytosine bases on singletons (+/-):\t" << data.sBase1[3] << "\t" << data.sBase2[3] << endl;

    unsigned long long sB_1 = 0; unsigned long long sB_2 = 0; for(int k=0; k<4; k++) { sB_1 += data.sBase1[k]; sB_2 += data.sBase2[k]; }

    statFile << endl;
    statFile << "% error(Guanine) removed on singletons (+/-):\t"  << (double)data.eBase1[0]/(double)sB_1 << "\t" << (double)data.eBase2[0]/(double)sB_2 << endl;
    statFile << "% error(Adenine) removed on singletons (+/-):\t"  << (double)data.eBase1[1]/(double)sB_1 << "\t" << (double)data.eBase2[1]/(double)sB_2 << endl;
    statFile << "% error(Thymine) removed on singletons (+/-):\t"  << (double)data.eBase1[2]/(double)sB_1 << "\t" << (double)data.eBase2[2]/(double)sB_2 << endl;
    statFile << "% error(Cytosine) removed on singletons (+/-):\t" << (double)data.eBase1[3]/(double)sB_1 << "\t" << (double)data.eBase2[3]/(double)sB_2 << endl;

    statFile.close();

 }

 ofstream freqFile;
 if(data.me == 0) {

    freqFile.open(freq_filename.c_str());

    freqFile << "Chr" << "\t" << "Loc" << "\t" << "Depth" << "\t" << "Ref";
    freqFile << "\t" << "ref_freq_+" << "\t" << "ref_freq_-";
    freqFile << "\t" << "G+" << "\t" << "G-";
    freqFile << "\t" << "A+" << "\t" << "A-";
    freqFile << "\t" << "T+" << "\t" << "T-";
    freqFile << "\t" << "C+" << "\t" << "C-";
    freqFile << endl;

    cerr << endl << "Printing Base Frequency table" << endl << endl;
}
 
    mrBuffer->map(mrLink, map_Base_Frequency_hGE, &data);
    mrBuffer->collate(NULL);

    TableLoc _table_loc;
    unsigned long long nloc = 1;
    for(int i=0; i<data.loc_vector; i++) {
	for(int j=data._target_loc[i].loc1; j<=data._target_loc[i].loc2; j++)  {
          unsigned long long point_target = (unsigned long long)data._target_loc[i].chr * (unsigned long long)1e11 + (unsigned long long)j;
          _table_loc._loc_freq_table[ point_target ] = nloc;
          nloc++;
      }
    }	

    _table_loc.table_positive = (unsigned long long *)malloc(sizeof(unsigned long long)*(nloc*4));
    _table_loc.table_negative = (unsigned long long *)malloc(sizeof(unsigned long long)*(nloc*4));
    for(unsigned long long i=0; i<(nloc*4); i++) { _table_loc.table_positive[i] = 0; _table_loc.table_negative[i] = 0; }

    mrBuffer->reduce(reduce_summation_base_frequency, &_table_loc);

    unsigned long long *out_point_positive;
    unsigned long long *out_point_negative;

    out_point_positive = (unsigned long long *)malloc(sizeof(unsigned long long)*(nloc*4));
    out_point_negative = (unsigned long long *)malloc(sizeof(unsigned long long)*(nloc*4));

    for(unsigned long long i=0; i<(nloc*4); i++) { out_point_positive[i] = 0; out_point_negative[i] = 0; }

    if(data.me == 0) cerr << endl << endl << "Reducing point substitution frequencies " << endl << endl;

    MPI_Reduce(_table_loc.table_positive, out_point_positive, nloc*4, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(_table_loc.table_negative, out_point_negative, nloc*4, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    free(_table_loc.table_positive);
    free(_table_loc.table_negative);

    if(data.me == 0) {

    for(int i=0; i<data.loc_vector; i++) 
        for(int j=data._target_loc[i].loc1; j<=data._target_loc[i].loc2; j++)  {
           unsigned long long point_target = (unsigned long long)data._target_loc[i].chr * (unsigned long long)1e11 + (unsigned long long)j;
           unsigned long long point        = _table_loc._loc_freq_table[ point_target ];

           char ref_base = data.reference_sequence[ data._target_loc[i].chr ].at(j-1);
           int  ref_int  = get_int_from_base( ref_base );

           unsigned long long depth = 0;
           for(int k=0; k<4; k++) {
                 depth += out_point_positive[ point*4 + k ];
                 depth += out_point_negative[ point*4 + k ];
           }       

            freqFile << data._target_loc[i].chr << "\t" << j << "\t" << depth << "\t" << ref_base;
            freqFile << "\t" << out_point_positive[ point*4 + (unsigned long long)ref_int ] << "\t" << out_point_negative[ point*4 + (unsigned long long)ref_int ];

            for(int k=0; k<4; k++) {
                if(k != ref_int) freqFile << "\t" << out_point_positive[ point*4 + (unsigned long long)k ];
                else             freqFile << "\t" << 0;

                if(k != ref_int) freqFile << "\t" << out_point_negative[ point*4 + (unsigned long long)k ];
                else             freqFile << "\t" << 0;
            }

            freqFile << endl;

        }

    }


     free(out_point_positive);
     free(out_point_negative);

 if(data.me == 0) freqFile.close();


} 

 
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

 MPI_Barrier(MPI_COMM_WORLD);

// ###########################################################################

 for(int i=0; i<QS_range; i++) {
   for(int j=0; j<4; j++) {
		free(data.qs1_count[i][j]);
   		free(data.qs2_count[i][j]);
   }
 
   free(data.qs1_count[i]);
   free(data.qs2_count[i]);
 }
 free(data.qs1_count);
 free(data.qs2_count);

free(data.eBase1); free(data.eBase2);
free(data.sBase1); free(data.sBase2);
free(data.ATGC1);  free(data.ATGC2);

// ###########################################################################


 MPI_Barrier(MPI_COMM_WORLD);

 delete mrLink;
 delete mrSeq;
 delete mrEkmer;
 delete mrBuffer;

 MPI_Finalize();
 exit(0);

}

