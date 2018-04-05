
#define __STDC_LIMIT_MACROS
#include "stdint.h"
#include "keyvalue.h"
#include "sequenceUtil.hpp"

using MAPREDUCE_NS::KeyValue;

#define ALLBITS UINT64_MAX
#define INT64MAX INT64_MAX
#define HIBIT UINT64_MAX-INT64_MAX

#define PI 3.14159265

typedef struct {
  map<unsigned long long,unsigned long long> _loc_freq_table;
  unsigned long long *table_positive;
  unsigned long long *table_negative;
} TableLoc;


typedef struct {
  int chr;
  int loc;
  char refbase;
} PointLoc;

typedef struct {
  int strand;
  char base;
} PointData;


typedef struct {
  int d0, d1, d2, d3;
} Depth;

typedef struct {
  int chr, loc; //ref_allele;
} ChrLoc;

typedef struct {
  int chr1, chr2, loc1, loc2, ID, strand1, strand2, maq1, maq2;
} PairLink;

typedef struct {
  int rank;
  uint64_t ID;
} RankID;

typedef struct {
  int chr, loc, strand, pe, SOFT_S, SOFT_E, maq, align_len;
} sReadInfo;

typedef struct {
  char header[100];
  int PE, strand, length, maq, type;
  char seq[200];
  char quality[200];
} HeadSeq;

struct ACGT_Count {
  int me;
  int range;
  unsigned int kmer_length;
  int acgt_count[43][4][4]; 
};

typedef struct {
  int len_consen, csize, type;
  char ref_sequence[1000];
  int num_consensus1[4100];
  int num_consensus2[4100];
  int loc_consensus[1100];
} RefConsen;

typedef struct {
  uint64_t loc1, loc2;
  int chr;
  int size;
  int *af;
} ExonLoc;

typedef struct {
    uint64_t zone,empty;
} PAD;

class ChromosomeLoc{
     int chr, loc;
   public:
     ChromosomeLoc (int a,int b) { chr = a; loc = b; } 
     void set_loc(int a, int b) { chr = a; loc = b; }
     int get_chr() { return chr; }
     int get_loc() { return loc; }
     bool operator <(const ChromosomeLoc& e) const
     {
	return (this->chr < e.chr) || ((this->chr == e.chr) && (this->loc < e.loc));
     }
     bool operator >(const ChromosomeLoc& e) const
     {
        return (this->chr > e.chr) || ((this->chr == e.chr) && (this->loc > e.loc));
     }
     bool operator ==(const ChromosomeLoc& e) const
     {
        return (this->chr == e.chr && this->loc == e.loc);
     }
};


struct LookUp {
//  map<ChromosomeLoc,Depth> Collision_Loc;
   map<uint64_t,int> Collision_Loc;
   int type, num_sclusters, num_sclusters_collision, num_sclusters_wo_collision, num_split, num_nclusters, num_singletons, nreads_sclusters, nreads_nclusters, nreads_singletons;
};

struct Data {
	uint64_t flag;	
        ofstream outFile, outFile_R1, outFile_R2;
	string filename;

        std::vector<std::string> reference_sequence; 
//	map<int,string> reference_sequence;

	std::vector<ExonLoc>     _target_loc;
	std::vector<std::string> _target_gene;
        int loc_vector;

        map<unsigned long long,unsigned long long> _loc_freq_table;
        unsigned long long *table_base;

	bool On_Alignment;
	bool On_Assemble;	
	bool On_Target;
	bool On_smallCluster;
	bool On_Clean;

	bool filter_offTarget;
	bool filter_clipping;

	bool On_Polishing;

	bool ftype;
	int  min_csize;
//	int side;
//	int depth;
//	int chromosome;

	int lane;
	int read_side;	
	int flank;
	int offside;
	int mismatch;
	int MAQ;
        int bQscore;
	int OVERLAP;
	float soft_allow;
//	float error_rate;

	unsigned int kmer_length;
	int range;
	int nKmer;

	int start;
        int end;

        int me,nprocs;
        int pshift;
        uint64_t lmask;
        int seed,nthresh;	
	PAD pad;

        bool DS;
	bool PACMAN;
	bool CRAWL;
	bool prune_error_kmers;

	unsigned int MAX_RECURSION;
	float MIN_SEED_ENTROPY;
	unsigned int MIN_SEED_COVERAGE;
	unsigned int crawl_length;
	float min_any_entropy;
	unsigned int min_kmer_count;
	float min_ratio_non_error;

        float MIN_CONNECTIVITY_RATIO;
        unsigned int MIN_ASSEMBLY_LENGTH;
        unsigned int MIN_ASSEMBLY_COVERAGE;

//	int *error_onFragment;
//	int *num_fragment;

	int ***qs1_count;
	int ***qs2_count;
	int ***read1_count;
	int ***read2_count;

	double ***qs1_prob;
	double ***qs2_prob;

	double ***read1_prob;
	double ***read2_prob;

	unsigned long long *eBase1;
	unsigned long long *eBase2;
        unsigned long long *sBase1;
	unsigned long long *sBase2;
	unsigned long long *ATGC1;
	unsigned long long *ATGC2;
	double *prior_prob;
};

typedef struct {
  int len1,len2,direct;
} LenDirect;

