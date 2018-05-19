#ifndef __GLOBAL_ALIGNMENT_H__
#define __GLOBAL_ALIGNMENT_H__

#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include "GenomicNeighborhood.h"
#include "ProteinCollection.h"

/**
 *Returns the dynamic programming matrix used in the global alignment algorithm (Needleman-Wunsch).
 */
std::vector<std::vector<double> > global_alignment_assignments(const std::vector<protein_info_t> &g1,
                                                               const std::vector<protein_info_t> &g2,
                                                               const ProteinCollection &clusters, double prot_stringency,
                                                               const int gap_score, const int mismatch);
/**
 *Prints the score between two genomic neighborhoods in the standard format.
 */
void global_alignment_output_score(const GenomicNeighborhood &g1, const GenomicNeighborhood &g2,
                                   double score, std::ofstream &output_file);

/**
 *Prints the aligned proteins.
 */
void global_alignment_output_pairings(const GenomicNeighborhood &g1, const GenomicNeighborhood &g2, const ProteinCollection &clusters,
                                      const std::vector<std::vector<double> > &assignments, double stringency, std::ofstream &pairings_file,
                                      const int gap_score, const int mismatch);

#endif
