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
////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////
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

#include <stdlib.h>
#include <map>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iterator>
#include <math.h>
#include <algorithm>
#include <vector>

#include "IRKE.hpp"
#include "sequenceUtil.hpp"
#include "KmerCounter.hpp"
#include "stacktrace.hpp"


IRKE::IRKE () {   // IRKE = Inchworm Recursive Kmer Extension

}

IRKE::IRKE (unsigned int kmer_length, unsigned int max_recursion, float min_seed_entropy, 
			unsigned int min_seed_coverage, float min_any_entropy, 
			bool pacman, bool crawl, unsigned int crawl_length, bool double_stranded) : kcounter(kmer_length, double_stranded) {
	
	MAX_RECURSION = max_recursion;
	MIN_SEED_ENTROPY = min_seed_entropy;
	MIN_SEED_COVERAGE = min_seed_coverage;
	MIN_ANY_ENTROPY = min_any_entropy;
	PACMAN = pacman;
	CRAWL = crawl;
	CRAWL_LENGTH = crawl_length;
	
	INCHWORM_ASSEMBLY_COUNTER = 0;
	
	DOUBLE_STRANDED_MODE = double_stranded;
	PRUNE_SINGLETON_READ_INTERVAL = 0;
	
}

bool IRKE::add_kmer (kmer_int_type_t kmer_val, unsigned int count) {
     kcounter.add_kmer(kmer_val, count);

//     cerr << kmer_val << "\t" << kcounter.get_kmer_count(kmer_val) << endl;
     return(true);
}

unsigned long IRKE::get_size() {
	return(kcounter.size());
}

void IRKE::set_prune_singleton_read_interval (unsigned long interval) {
	
	this->PRUNE_SINGLETON_READ_INTERVAL = interval;
	
}


void IRKE::prune_kmers_min_count(unsigned int min_count) {
	
	// proxy, sends message to kcounter member
	
	kcounter.prune_kmers_min_count(min_count);
}


void IRKE::prune_kmers_min_entropy(float min_entropy) {
	
	// proxy, send message to kcounter
	
	kcounter.prune_kmers_min_entropy(min_entropy);
	
	
}


bool IRKE::prune_kmer_extensions(float min_ratio_non_error) {
	
	// proxy, send message to kcounter
	
	return(kcounter.prune_kmer_extensions(min_ratio_non_error));
}


bool IRKE::prune_some_kmers(unsigned int min_count, float min_entropy, bool prune_error_kmers, float min_ratio_non_error) {

	return(kcounter.prune_some_kmers(min_count, min_entropy, prune_error_kmers, min_ratio_non_error));
}


unsigned long IRKE::get_graph_size() {
	
	// proxy call
	
	return(kcounter.size());
	
}


void IRKE::traverse_path(KmerCounter& kcounter, Kmer_Occurence_Pair seed_kmer, Kmer_visitor& visitor,
						 Kmer_visitor& place_holder, float MIN_CONNECTIVITY_RATIO, unsigned int depth) {
	
	// check to see if visited already
	if (visitor.exists(seed_kmer.first)) {
		// already visited
		
		return;
	}
	
	// check if at the end of the max recursion, and if so, don't visit it but set a placeholder.
	if (depth > MAX_RECURSION) {
		place_holder.add(seed_kmer.first);
		return;
	}
	
	visitor.add(seed_kmer.first);
    
	// try each of the forward paths from the kmer:  
	vector<Kmer_Occurence_Pair> forward_candidates = kcounter.get_forward_kmer_candidates(seed_kmer.first);
	
	for(unsigned int i = 0; i < forward_candidates.size(); i++) {
		Kmer_Occurence_Pair kmer = forward_candidates[i];
		
		if (kmer.second	&& exceeds_min_connectivity(kcounter, seed_kmer, kmer, MIN_CONNECTIVITY_RATIO)) {
			
			traverse_path(kcounter, kmer, visitor, place_holder, MIN_CONNECTIVITY_RATIO, depth + 1);
		}
	}
	
	// try each of the reverse paths from the kmer:
	vector<Kmer_Occurence_Pair> reverse_candidates = kcounter.get_reverse_kmer_candidates(seed_kmer.first);
	
	for (unsigned int i = 0; i < reverse_candidates.size(); i++) {
		Kmer_Occurence_Pair kmer = reverse_candidates[i];
		
		if (kmer.second	&& exceeds_min_connectivity(kcounter, seed_kmer, kmer, MIN_CONNECTIVITY_RATIO)) {
			
			traverse_path(kcounter, kmer, visitor, place_holder, MIN_CONNECTIVITY_RATIO, depth + 1);
		}
	}
	
	
	return;
	
}


string add_fasta_seq_line_breaks(string& sequence, int interval) {
    
    stringstream fasta_seq;
    
    int counter = 0;
    for (string::iterator it = sequence.begin(); it != sequence.end(); it++) {
        counter++;
        
        fasta_seq << *it;
        if (counter % interval == 0 && (it + 1) != sequence.end()) {
            fasta_seq << endl;
        }
    }

    return(fasta_seq.str());
}



string IRKE::compute_sequence_assemblies(float min_connectivity, unsigned int MIN_ASSEMBLY_LENGTH, unsigned int MIN_ASSEMBLY_COVERAGE) {

	string output_sequence = compute_sequence_assemblies(kcounter, min_connectivity, MIN_ASSEMBLY_LENGTH, MIN_ASSEMBLY_COVERAGE);
	return (output_sequence.c_str());
}


string IRKE::compute_sequence_assemblies(KmerCounter& kcounter, float min_connectivity,
			                 unsigned int MIN_ASSEMBLY_LENGTH, unsigned int MIN_ASSEMBLY_COVERAGE) {
	
    if (! got_sorted_kmers_flag) {
        stringstream error;
        error << stacktrace() << " Error, must populate_sorted_kmers_list() before computing sequence assemblies" << endl;
        throw(error.str());
    }

	//vector<Kmer_counter_map_iterator>& kmers = sorted_kmers; 
    	vector<Kmer_Occurence_Pair>& kmers = sorted_kmers; 
	
	unsigned long init_size = kcounter.size();
    
	unsigned int kmer_length = kcounter.get_kmer_length();

	string output_sequence("");
	
	for (unsigned int i = 0; i < kmers.size(); i++) {
		
		unsigned long kmer_counter_size = kcounter.size();
		if (kmer_counter_size > init_size) {
			
			stringstream error;
			error << stacktrace() << "Error, Kcounter size has grown from " << init_size
				  << " to " << kmer_counter_size << endl;
			throw (error.str());
		}

        	kmer_int_type_t kmer = kmers[i].first;
        	unsigned int kmer_count = kcounter.get_kmer_count(kmer);

        	//continue;
        
		if (kmer_count == 0) {
			continue;
		}

		if (kmer == revcomp_val(kmer, kmer_length)) {
			// palindromic kmer, avoid palindromes as seeds
            
            		continue;
		}
        
		if (kmer_count < MIN_SEED_COVERAGE) {
            		continue;
		}
		
		float entropy = compute_entropy(kmer, kmer_length);
		
		if (entropy < MIN_SEED_ENTROPY) {
            		continue;
		}
				
		/* Extend to the right */
		
		Kmer_visitor visitor(kmer_length, DOUBLE_STRANDED_MODE);
		Path_n_count_pair selected_path_n_pair_forward = inchworm(kcounter, 'F', kmer, visitor, min_connectivity); 
		
		visitor.clear();
		// add selected path to visitor
		
		vector<kmer_int_type_t>& forward_path = selected_path_n_pair_forward.first;

        	for (unsigned int i = 0; i < forward_path.size(); i++) {
			kmer_int_type_t kmer = forward_path[i];
			visitor.add(kmer);
		}
		
		/* Extend to the left */ 
		visitor.erase(kmer); // reset the seed
		
		Path_n_count_pair selected_path_n_pair_reverse = inchworm(kcounter, 'R', kmer, visitor, min_connectivity);
		
		unsigned int total_counts = selected_path_n_pair_forward.second + selected_path_n_pair_reverse.second + kmer_count; //kcounter.get_kmer_count(kmer); 
		
		vector<kmer_int_type_t>& reverse_path = selected_path_n_pair_reverse.first;
		
		vector<kmer_int_type_t> joined_path = _join_forward_n_reverse_paths(reverse_path, kmer, forward_path);
		
		// report sequence reconstructed from path.
		
		vector<unsigned int> assembly_base_coverage;
		string sequence = reconstruct_path_sequence(kcounter, joined_path, assembly_base_coverage);
		
		int avg_cov =  static_cast<int> ( (float)total_counts/(sequence.length()-kcounter.get_kmer_length() +1) + 0.5);
	
		
		if (sequence.length() >= MIN_ASSEMBLY_LENGTH && avg_cov >= MIN_ASSEMBLY_COVERAGE) {
			if(sequence.length() > output_sequence.length()) output_sequence.assign( sequence );
			INCHWORM_ASSEMBLY_COUNTER++;
		}
		
		// remove path
            	for (unsigned int i = 0; i < joined_path.size(); i++) {
                	kmer_int_type_t kmer = joined_path[i];
                	kcounter.clear_kmer(kmer);
            	}
		
    }
        
    	// drop sorted kmer list as part of cleanup
    	clear_sorted_kmers_list();
    
	return (output_sequence.c_str()); // end of runIRKE
}

Path_n_count_pair IRKE::inchworm (KmerCounter& kcounter, char direction, kmer_int_type_t kmer, Kmer_visitor& visitor, float min_connectivity) {
	
	// cout << "inchworm" << endl;
	
	Path_n_count_pair entire_path;
    	entire_path.second = 0; // init cumulative path coverage

	unsigned int inchworm_round = 0;
	
	unsigned long num_total_kmers = kcounter.size();
	
	Kmer_visitor eliminator(kcounter.get_kmer_length(), DOUBLE_STRANDED_MODE);
	
	while (true) {
        		
		inchworm_round++;
		eliminator.clear();
		
		if (inchworm_round > num_total_kmers) {
			throw(string ("Error, inchworm rounds have exceeded the number of possible seed kmers"));
		}
		
		visitor.erase(kmer); // seed kmer must be not visited already.
		
		Kmer_Occurence_Pair kmer_pair(kmer, kcounter.get_kmer_count(kmer));
		Path_n_count_pair best_path = inchworm_step(kcounter, direction, kmer_pair, visitor, eliminator, inchworm_round, 0, min_connectivity, MAX_RECURSION);
		

        	vector<kmer_int_type_t>& kmer_list = best_path.first;
        	unsigned int num_kmers = kmer_list.size();

		if(best_path.second > 0) {
			// append info to entire path in reverse order, so starts just after seed kmer
			
			int first_index = num_kmers - 1;
			int last_index = 0;
			if (CRAWL) {
				last_index = first_index - CRAWL_LENGTH + 1;
				if (last_index < 0) {
					last_index = 0;
				}
			}
			
			for (int i = first_index; i >= last_index; i--) {
				kmer_int_type_t kmer_extend = kmer_list[i];
				entire_path.first.push_back(kmer_extend);
				visitor.add(kmer_extend);
				//entire_path.second += kcounter.get_kmer_count(kmer_extend);
			}
			
			kmer = entire_path.first[ entire_path.first.size() -1 ];
            
            		entire_path.second += best_path.second;

		}
		else {
			// no extension possible
			break;
		}
	}
	
	return(entire_path);
}


bool compare (const Path_n_count_pair& valA, const Path_n_count_pair& valB) {
	
#ifdef _DEBUG
    if (valA.second == valB.second)
        return (valA.first > valB.first);
	else
#endif
		return(valA.second > valB.second); // reverse sort.
}



Path_n_count_pair IRKE::inchworm_step (KmerCounter& kcounter, char direction, Kmer_Occurence_Pair kmer, Kmer_visitor& visitor,
									   Kmer_visitor& eliminator, unsigned int inchworm_round, unsigned int depth, 
									   float MIN_CONNECTIVITY_RATIO, unsigned int max_recurse) {
	
	// cout << "inchworm_step" << endl;
	
	// check to see if kmer exists.  If not, return empty container
	Path_n_count_pair best_path_n_pair;
	best_path_n_pair.second = 0; // init
		
	if ( // !kmer.second || 
        	visitor.exists(kmer.first) // visited
        	|| eliminator.exists(kmer.first) // eliminated
	 ) {
        	return(best_path_n_pair);
	}
	
	visitor.add(kmer.first);
	
	if (PACMAN && depth > 0) {
		// cerr << "pacman eliminated kmer: " << kmer << endl;
		eliminator.add(kmer.first);
	}
	
	
	if (depth < max_recurse) {
		
		vector<Kmer_Occurence_Pair> kmer_candidates;
		if (direction == 'F') {
			// forward search
			kmer_candidates = kcounter.get_forward_kmer_candidates(kmer.first);
		}
		else {
			// reverse search
			kmer_candidates = kcounter.get_reverse_kmer_candidates(kmer.first);
		}
		
		bool tie = true;
		unsigned int recurse_cap = max_recurse;
		unsigned int best_path_length = 0;
		while (tie) {

            		// keep trying to break ties if ties encountered.
            		// this is done by increasing the allowed recursion depth until the tie is broken.
            		//  Recursion depth set via: recurse_cap and incremented if tie is found
            

			vector<Path_n_count_pair> paths;
			
			for (unsigned int i = 0; i < kmer_candidates.size(); i++) {
				Kmer_Occurence_Pair kmer_candidate = kmer_candidates[i];
				
				if ( kmer_candidate.second &&
				    
                    			!visitor.exists(kmer_candidate.first)  // avoid creating already visited kmers since they're unvisited below...
					&& exceeds_min_connectivity(kcounter, kmer, kmer_candidate, MIN_CONNECTIVITY_RATIO) ) {
					//cout << endl << "\ttrying " << kmer_candidate << endl;
					

                    			// recursive call here for extension
					Path_n_count_pair p = inchworm_step(kcounter, direction, kmer_candidate, visitor, eliminator, inchworm_round, depth+1, MIN_CONNECTIVITY_RATIO, recurse_cap);
					
                    		if (p.first.size() >= 1) {
                        		// only retain paths that include visited nodes.
                        		paths.push_back(p);
                    		}
					visitor.erase(kmer_candidate.first); // un-visiting
					
                	}
				
			} // end for kmer
			
			
			if (paths.size() > 1) {

                		sort(paths.begin(), paths.end(), compare);
				
				if (paths[0].second == paths[1].second
				    &&
				    paths[0].first[0] != paths[1].first[0]
                                    ) {
					
					if (paths[0].first.size() > best_path_length) {
                        		// still making progress in extending to try to break the tie.  Keep going.
                        		// note, this is the only test that keeps us in this while loop. (tie stays true)
                        			recurse_cap++;
						best_path_length = paths[0].first.size();
					}
					else {
						// cerr << "not able to delve further into the graph, though...  Stopping here." << endl;
						tie = false;
                        			best_path_n_pair = paths[0]; // pick one
					}
				}
				
				else if ((paths[0].second == paths[1].second   // same cumulative coverage values for both paths.
						  &&
						  paths[0].first[0] == paths[1].first[0] ) // same endpoint
					) {
					
					tie = false;
					best_path_n_pair = paths[0];
				}
				
				else {
					// no tie.
					tie = false;
					best_path_n_pair = paths[0];
				}
				
				
			}
			else if (paths.size() == 1) {
				tie = false;
				best_path_n_pair = paths[0];
			}
			else {
				// no extensions possible.
				tie = false;
			}
			
			
		} // end while tie
	}
	
	// add current kmer to path, as long as not the original seed kmer!
	if (depth > 0) {
		best_path_n_pair.first.push_back(kmer.first);
		best_path_n_pair.second += kmer.second;
    	}
	
	return(best_path_n_pair);
	
}


vector<kmer_int_type_t> IRKE::_join_forward_n_reverse_paths(vector<kmer_int_type_t>& reverse_path, 
         							kmer_int_type_t seed_kmer_val, 
								vector<kmer_int_type_t>& forward_path) {
	
	vector<kmer_int_type_t> joined_path;
	
	// want reverse path in reverse order
	
	for (int i = reverse_path.size()-1; i >= 0; i--) {
		joined_path.push_back( reverse_path[i] );
	}
	
	// add seed kmer
	joined_path.push_back(seed_kmer_val);
	
	// tack on the entire forward path.
	
	for (unsigned int i = 0; i < forward_path.size(); i++) {
		joined_path.push_back( forward_path[i] );
	}
	
	return(joined_path);
}


string IRKE::reconstruct_path_sequence(vector<kmer_int_type_t>& path, vector<unsigned int>& cov_counter) {
	// use kcounter member
	return(reconstruct_path_sequence(kcounter, path, cov_counter));
}


string IRKE::reconstruct_path_sequence(KmerCounter& kcounter, vector<kmer_int_type_t>& path, vector<unsigned int>& cov_counter) {
	
	if (path.size() == 0) {
		return("");
	}
	
	string seq = kcounter.get_kmer_string(path[0]);
	cov_counter.push_back( kcounter.get_kmer_count(path[0]) );
	
	for (unsigned int i = 1; i < path.size(); i++) {
		string kmer = kcounter.get_kmer_string(path[i]);
		seq += kmer.substr(kmer.length()-1, 1);
		
		cov_counter.push_back( kcounter.get_kmer_count(path[i]) );
	}
	
	return(seq);
}


bool IRKE::exceeds_min_connectivity (KmerCounter& kcounter, Kmer_Occurence_Pair kmerA, Kmer_Occurence_Pair kmerB, float min_connectivity) {
	
    if (min_connectivity < 1e5) {
        return(true); // consider test off
    }
    
	
	unsigned int kmerA_count = kmerA.second;
	if (kmerA_count == 0)
		return(false);
	unsigned int kmerB_count = kmerB.second;
	if (kmerB_count == 0)
		return(false);
	
	unsigned int minVal;
	unsigned int maxVal;
	
	if (kmerA_count < kmerB_count) {
		minVal = kmerA_count;
		maxVal = kmerB_count;
	}
	else {
		minVal = kmerB_count;
		maxVal = kmerA_count;
	}
	
	float connectivity_ratio = (float) minVal/maxVal;
	
	if (connectivity_ratio >= min_connectivity) {
		return(true);
	}
	else {
		return(false);
	}
}


bool IRKE::exceeds_min_connectivity (KmerCounter& kcounter, string kmerA, string kmerB, float min_connectivity) {

	kmer_int_type_t valA = kmer_to_intval(kmerA);
	kmer_int_type_t valB = kmer_to_intval(kmerB);

	Kmer_Occurence_Pair pairA(valA, kcounter.get_kmer_count(valA));
	Kmer_Occurence_Pair pairB(valB, kcounter.get_kmer_count(valB));

	return exceeds_min_connectivity(kcounter, pairA, pairB, min_connectivity);

}



string IRKE::thread_sequence_through_graph(string& sequence) {
	
	// describe each of the ordered kmers in the input sequence as they exist in the graph.
	
	unsigned int kmer_length = kcounter.get_kmer_length();
	
	if (sequence.length() < kmer_length) {
		cerr << "Sequence length: " << sequence.length() << " is too short to contain any kmers." << endl;
		return("");
	}
	
	stringstream s;
	
	for (unsigned int i=0; i <= sequence.length() - kmer_length; i++) {
		
		string kmer = sequence.substr(i, kmer_length);
		
		s << kcounter.describe_kmer(kmer) << endl;
	}
	
	return(s.str());
	
	
}



bool IRKE::sequence_path_exists(string& sequence, unsigned int min_coverage, float min_entropy, float min_connectivity, 
								vector<unsigned int>& coverage_counter) {
	
	unsigned int kmer_length = kcounter.get_kmer_length();
	
	if (sequence.length() < kmer_length) {
		return(false);
	}
	
	bool path_exists = true;
	
	string prev_kmer = sequence.substr(0, kmer_length);
	if (contains_non_gatc(prev_kmer) || ! kcounter.kmer_exists(prev_kmer)) {
		path_exists = false;
		coverage_counter.push_back(0);
	}
	else {
		unsigned int kmer_count = kcounter.get_kmer_count(prev_kmer);
		coverage_counter.push_back(kmer_count);

		float entropy = compute_entropy(prev_kmer);

		if (kmer_count < min_coverage || entropy < min_entropy) {
			path_exists = false;
		}
	}
	
	
	for (unsigned int i=1; i <= sequence.length() - kmer_length; i++) {
		
		string kmer = sequence.substr(i, kmer_length);
		
		if (contains_non_gatc(kmer) || ! kcounter.kmer_exists(kmer)) {
			path_exists = false;
			coverage_counter.push_back(0);
		}
		else {
			unsigned int kmer_count = kcounter.get_kmer_count(kmer);
			coverage_counter.push_back(kmer_count);

			float entropy = compute_entropy(kmer);

			if (kmer_count < min_coverage || entropy < min_entropy) {
				path_exists = false;
			}
		}

		
		if (path_exists && ! exceeds_min_connectivity(kcounter, prev_kmer, kmer, min_connectivity)) {
			path_exists = false;
		}
		
		prev_kmer = kmer;
		
	}
	
	return(path_exists);
}



void IRKE::describe_kmers() {
	
	return(describe_kmers(kcounter));
}

void IRKE::describe_kmers(KmerCounter& kcounter) {
	
	// proxy call to Kcounter method.
	
	kcounter.describe_kmers();
	
	return;
}

void IRKE::populate_sorted_kmers_list() {
    
    sorted_kmers = kcounter.get_kmers_sort_descending_counts();
    got_sorted_kmers_flag = true;

    return;
}

void IRKE::clear_sorted_kmers_list() {

    sorted_kmers.clear();
    got_sorted_kmers_flag = false;
    
    return;
}
