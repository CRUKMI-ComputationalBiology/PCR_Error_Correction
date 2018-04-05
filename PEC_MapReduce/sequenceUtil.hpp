#define __STDC_LIMIT_MACROS
#include "stdint.h"

#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <vector>
#include <map>
#include <cctype>

#include <set>

#include <cstdio>
#include <cctype>
#include <cstdlib>

//#include "gmp.h"
//#include "gmpxx.h"

#ifndef __SEQUENCEUTIL__

#define __SEQUENCEUTIL__


using namespace std;

struct  CigarOp_Function {
    char op;
    int size;
    static int parse(const char* s,vector<CigarOp_Function>& v)
        {
        char* p=(char*)(s);
        while(*p!=0)
            {
            char* endptr;
            CigarOp_Function c;

            if(!isdigit(*p)) return -1;
            c.size =strtol(p,&endptr,10);
            if(c.size<=0) return -1;
            p=endptr;

            c.op = *p;
	    if(c.op == '=') c.op = 'M'; 
            if(!isalpha(c.op)) return -1;
            v.push_back(c);

            ++p;
            }
        return 0;
        }
};


struct  CigarOp {
    char op;
    int size;
/*
    static int parse(const char* s,vector<CigarOp>& v)
        {
        char* p=(char*)(s);
        while(*p!=0)
            {
            char* endptr;
            CigarOp c;

            if(!isdigit(*p)) return -1;
            c.size =strtol(p,&endptr,10);
            if(c.size<=0) return -1;
            p=endptr;

	    c.op = p[0];
            if(!isalpha(c.op)) return -1;
            v.push_back(c);

            ++p;
            }
        return 0;
        }
*/
};
int CigarOp_parse(char* s,vector<CigarOp> &v);

// misc typedefs
typedef struct { string accession; string header; string sequence; } fastaRecord;

typedef unsigned long long kmer_int_type_t;
//typedef uint64_t kmer_int_type_t;
//typedef mpz_class kmer_int_type_t;
typedef pair<kmer_int_type_t,unsigned int> Kmer_Occurence_Pair;


const static short expand_region = 1;

// function prototypes

void convert_cigar2bam(int32_t cigar, int32_t *output);

int _consen_alt_allele(int consen, int alt);
int check_collision_base(int cbase, int base_collision);

int get_num_match(string seq_a, string seq_b, double mu, double delta, int range);
int check_strand_direction(string seq_a, string seq_b, double mu, double delta, int range, int num_mismatch);
double similarity_score(char a,char b,double mu);
double find_array_max(double array[],int length, int &ind);

int align_query2assembly(string seq_ref, string seq_read, double mu, double delta,char *consen_a, char *consen_b, int *loc_a, int *loc_b);
int align_reads_assembly_start(string seq_ref, string seq_read, double mu, double delta);
void align_reads_assembly(string seq_ref, string seq_read, double mu, double delta, ofstream &outFile, int DS, int *num_consensus);
int smith_waterman(string seq_a, string seq_b, double mu, double delta, int *loc_a, int *loc_b, char *consensus_a, char *consensus_b);
char get_base_from_int(int intval);
int get_int_from_base(char base);
int count_at(string seq);
void count_gatc(double **gatc, string seq);

void compare_strings(int *output, int offset, int range, string seq_ref, string seq_query);

//void get_kmer_to_intvalArray (string kmer_forward, int *forward, string kmer_backward, int *backward, int nkmer);
void get_kmer_to_intvalArray (string kmer_forward, int *forward, int nkmer);
void get_kmer2intvalArray (string kmer, int *intarray, int nkmer);


string read_sequence_from_file (string filename);
string revcomp (const string);
string revquality (const string quality);

fastaRecord readNextFastaRecord(ifstream& reader);
bool contains_non_gatc(string kmer);
bool contains_non_gatc_num (string seq, int threshold);
string remove_whitespace(string s);

int get_Qscore(char quality);
char get_QualityAscii(int qscore);

char int_to_base(int baseval); // 0 1 2 3 => G A T C
int base_to_int_value(char nucleotide); // (GATC) = {0 1 2 3}, others = -1


void get_error_bases(int *error_index, int *error_base, int *cseq_base, string seq_kmer, string cseq_kmer, int kmer_length);

void untrusted_union(vector<int> untrusted_subset, vector<short> & region, int kmer_length);
bool untrusted_intersect(vector<int> untrusted_subset, vector<short> & region, int kmer_length, int read_length);
vector<short> error_region_chop(vector<int> untrusted_subset, int kmer_length, int read_length, float* prob);
vector<short> error_region_extract(vector<int> untrusted_subset, int kmer_length, int read_length);


kmer_int_type_t get_maximum_kmer_intval(unsigned int kmer_length);
kmer_int_type_t kmer_to_intval(string kmer); // must be less than 32 bases for 64-bit conversion
string decode_kmer_from_intval (kmer_int_type_t intval, unsigned int kmer_length);
kmer_int_type_t revcomp_val(kmer_int_type_t kmer, unsigned int kmer_length);
kmer_int_type_t get_DS_kmer_val(kmer_int_type_t kmer_val, unsigned int kmer_length);
vector<kmer_int_type_t> sequence_string_to_kmer_int_type_vector(const string& sequence, int kmer_length);

bool check_kmer_quality(string qkmer, int qcutoff);

//uint64_t convert_mpz_uint64(kmer_int_type_t intval, unsigned int kmer_length);
//string decode_kmer_from_uint64(uint64_t intval, unsigned int kmer_length);


vector<kmer_int_type_t> get_Fkmer_candidates(kmer_int_type_t seed_kmer, unsigned kmer_length);
vector<kmer_int_type_t> get_Rkmer_candidates(kmer_int_type_t seed_kmer, unsigned kmer_length);

//bool Sort_kmer_by_count_desc (const Kmer_Occurence_Pair& i, const Kmer_Occurence_Pair& j);

float compute_entropy(string& kmer);
float compute_entropy(kmer_int_type_t kmer, unsigned int kmer_length);

string replace_nonGATC_chars_with_A(string& input_seq);


#endif


