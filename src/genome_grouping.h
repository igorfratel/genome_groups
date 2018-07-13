#ifndef __GENOME_GROUPING_H__
#define __GENOME_GROUPING_H__

#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include "GenomicNeighborhood.h"
#include "ProteinCollection.h"
#include "porthodom_scoring.h"
#include "porthodomO2_scoring.h"
#include "global_alignment.h"
#include "utils.h"

/**
 *Receives a vector of genomic neighborhoods,
 *a ProteinCollection, the desired genomic neighborhood clustering method, a
 *gap score value (for when an alignment method is chosen) and a mismatch value.
 *Writes the similarity between all genomic neighborhoods on the genome_sim_filename and,
 *optionally, the pairings made between their proteins on the pairings_filename.
 */
 void genome_clustering(const std::vector<GenomicNeighborhood> &neighborhoods, const ProteinCollection &clusters,
                        const std::string &method, double prot_stringency, double neigh_stringency, const std::string &genome_sim_filename,
                        const std::string &pairings_filename, const int gap_score, const int mismatch);

/**
 *Receives a genomic neighborhood filename
 *Returns a vector of genomic neighborhoods, filled with the information from the file
 */
std::vector<GenomicNeighborhood> parse_neighborhoods(const std::string &neighborhoods_filename);

#endif
