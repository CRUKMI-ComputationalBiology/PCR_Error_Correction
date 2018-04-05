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

#include "fstream"
#include "iostream"
#include "sstream"

#include "vector"
#include "map"
#include "set"
#include "algorithm"
#include "time.h"
#include "math.h"
//#include "boost/math/special_functions/binomial.hpp"

//#include "boost/math/distributions/binomial.hpp"
//#include "boost/math/distributions/chi_squared.hpp"
//#include "boost/math/distributions/poisson.hpp"

//#include "gmp.h"
//#include "gmpxx.h"

#include "Fasta_reader.hpp"
#include "Fastq_reader.hpp"
#include "sequenceUtil.hpp"

#include "ssw_cpp.h"

using namespace MAPREDUCE_NS;
using namespace std;

//using boost::math::binomial;
//using boost::math::chi_squared;
//using boost::math::poisson;



int _base_collision(Depth depth);

void fileread_Sequence_Fastq(int itask, KeyValue *kv, void *ptr);
void fileread_AlignmentInfo_MPI(int itask, KeyValue *kv, void *ptr);


void fileread_RNAseq_HeadSeq_FASTA(int itask, char *fname, KeyValue *kv, void *ptr);
void fileread_RNAseq_HeadSeq_FASTQ(int itask, char *fname, KeyValue *kv, void *ptr);

void map_assign_reference(uint64_t itask, char *key, int keybytes, char *value, int valuebytes, KeyValue *kv, void *ptr);
void map_upload_reference(int itask, KeyValue *kv, void *ptr);

void fileread_AlignmentInfo(int itask, char *fname, KeyValue *kv, void *ptr);
bool check_onTarget(int chr, int loc1, int loc2, int loc_vector, int flank, const vector<ExonLoc>& _target_loc);
int parse_cigar_softclipping(char *s);
int parse_cigar_clipping(char *s, bool start);
int parse_cigar_align_length(char *s);

int CigarParse_Vector(char *s, char *Str, int *Len);

void reduce_Link_from_Header(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);
void reduce_LinkID_HeadSeq(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);
void reduce_LinkID_HeadSeq_single(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);
void revcomp_sequence_carray(char *seq, int len);
void revcomp_quality_carray(char *qual, int len);

void reduce_check_DNACollision_wo_subClustering(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);
void reduce_check_DNACollision_MiniBarcode(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);

void reduce_consensus_sequence_quality(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);
void reduce_consensus_SeqQual_highDup(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);
void reduce_consensus_Single_highDup(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);
void reduce_consensus_SeqQual_singleton(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);
void reduce_consensus_Single_singleton(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);


void reduce_consensus_highDup(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);
void reduce_sort_singleton(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);
void reduce_consensus_singleton(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);

void reduce_ekmers_from_all(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);
void reduce_All_Kmers(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);
void reduce_Kmer_Cov_All(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);
void reduce_Kmer_Cov(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);
void print_Kmer_Cov(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);
void print_Kmer_Cov_All(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);


void reduce_error_kmers(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);
void reduce_extract_kmers_for_error_counts(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);
void reduce_error_onAmplicons(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);
void reduce_error_matrices_step1(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);
void reduce_error_matrices_step2(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);

void reduce_rescue_by_ekmer(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);
void reduce_ekmers_rescue_step1(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);
void reduce_ekmers_rescue_step2(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);

void map_correct_singleton_step1(uint64_t itask, char *key, int keybytes, char *value, int valuebytes, KeyValue *kv, void *ptr);
void reduce_correct_single_step2(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);
void reduce_correct_single_step3(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);
void reduce_ekmer_singleton(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);

void map_print_like_correct(uint64_t itask, char *key, int keybytes, char *value, int valuebytes, KeyValue *kv, void *ptr);
void align_seq(StripedSmithWaterman::Alignment alignment, int seq_len, string seq, string ref_sequence, int *index);

void reduce_ntnt_count(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);
void find_error_base(StripedSmithWaterman::Alignment *alignment, int *error_index, int *error_base, int *ref_base, const char *seq, const char *ref_seq, int seq_length);



string align_csequence(StripedSmithWaterman::Alignment alignment, string ref_seq, string query_seq, int ref_len, int query_len);
string align_quality(StripedSmithWaterman::Alignment alignment, string ref_seq, string query_quality, int ref_len, int query_len);
void map_print2samformat_singleton(uint64_t itask, char *key, int keybytes, char *value, int valuebytes, KeyValue *kv, void *ptr);
void map_print2samformat_removeSoftClip(uint64_t itask, char *key, int keybytes, char *value, int valuebytes, KeyValue *kv, void *ptr);

void map_Base_Frequency_hGE(uint64_t itask, char *key, int keybytes, char *value, int valuebytes, KeyValue *kv, void *ptr);
void map_Base_Frequency(uint64_t itask, char *key, int keybytes, char *value, int valuebytes, KeyValue *kv, void *ptr);
void reduce_summation_base_frequency(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);
void reduce_base_frequency_HiDup(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);

void reduce_ekmers_HiDup(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);
void reduce_ekmer_for_HighDup(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);
void reduce_extract_kmer_singleton(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);
void reduce_ekmer_for_singleton(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);

void reduce_unique_kmers(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);

void map_kmer_from_singleton(uint64_t itask, char *key, int keybytes, char *value, int valuebytes, KeyValue *kv, void *ptr);
void reduce_kmers_from_singleton(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);
void reduce_kmers_second_resuce(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);


