#ifndef __PORTHODOMO2_SCORING_H__
#define __PORTHODOMO2_SCORING_H__

#include <map>
#include <string>
#include <fstream>
#include <algorithm>
#include <iterator>
#include "Hungarian.h"
#include "GenomicNeighborhood.h"
#include "ProteinCollection.h"

/**
 *Prints the score between two genomic neighborhoods in the standard format
 */
void porthodomO2_output_score(const GenomicNeighborhood &g1, const GenomicNeighborhood &g2, double score, std::ofstream &output_file);

/**
 *Prints the chosen protein assignments to the pairings_file (treats each assignment as a pair of pairs of proteins)
 */
void porthodomO2_output_pairings(const GenomicNeighborhood &g1, const GenomicNeighborhood &g2,
                                 const std::map<std::pair<int, int>,int> &assignments, std::ofstream &pairings_file);

/**
 *Receives the porthodom assignments and a normalizing factor (length of the longest neighborhood).
 *Returns the porthodom O2 MWM score (that takes order in consideration).
 */
double porthodomO2_scoring(const std::map<std::pair<int, int>, int> &assignments, int length);

/**
 *Receives two protein_info_t vectors, a ProteinCollection, the protein stringency and a mismatch value.
 *Returns the MWM porthodom  O2 protein assignments between the two neighborhoods (vectors).
 */
std::map<std::pair<int, int>, int> porthodomO2_assignments(const std::vector<protein_info_t> &g1,
                                                           const std::vector<protein_info_t> &g2,
                                                           const ProteinCollection &clusters, double prot_stringency,
                                                           const int mismatch);
#endif
