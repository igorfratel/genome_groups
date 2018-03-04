#ifndef __GENOME_GROUPING_H__
#define __GENOME_GROUPING_H__

#include <fstream>
#include <string>
#include <cstdlib>
#include "GenomicNeighborhood.h"
#include "ProteinCollection.h"
#include "porthodom_scoring.h"
#include "porthodomO2_scoring.h"

/**
 *Receives a vector of genomic neighborhoods,
 *a ProteinCollection and the desired genomic neighborhood clustering method.
 *Writes the similarity between all genomic neighborhoods on the genome_sim_filename and,
 *optionally, the pairings made between their proteins on the pairings_filename.
 */
 void genome_clustering(const std::vector<GenomicNeighborhood> &neighborhoods, const ProteinCollection &clusters,
                        const std::string &method, double prot_stringency, double neigh_stringency, const std::string &genome_sim_filename,
                        const std::string &pairings_filename);

/**
 *Receives a genomic neighborhood filename
 *Returns a vector of genomic neighborhoods, filled with the information from the file
 */
std::vector<GenomicNeighborhood> parse_neighborhoods(const std::string &neighborhoods_filename);

#endif
