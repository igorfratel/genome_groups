#ifndef __GENOME_GROUPING_H__
#define __GENOME_GROUPING_H__

#include <mutex>
#include <fstream>
#include <string>
#include <cstdlib>
#include "GenomicNeighborhood.h"
#include "ProteinCollection.h"
#include "porthodom_scoring.h"
#include "porthodomO2_scoring.h"

/*Receives a file containing all the genomic neighborhoods.
 *a ProteinCollection and the desired genomic neighborhood clustering method.
 *Writes the similarity between all genomic neighborhoods on the genome_sim_filename
 *the format "organism1 acession1 coordinates1 organism2 acession2 coordinates2 score"*/
void genome_clustering(const std::string &neighborhoods_filename, ProteinCollection &clusters,
                       const std::string &method, double stringency, const std::string &genome_sim_filename);

#endif
