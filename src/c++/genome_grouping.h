#ifndef __GENOME_GROUPING_H__
#define __GENOME_GROUPING_H__

#include <fstream>
#include <map>
#include <string>
#include <algorithm>
#include "Hungarian.h"
#include "GenomicNeighborhood.h"
#include "UndirectedEdgeWeightedGraph.h"

/*Receives a file containing all the genomic neighborhoods (as output by parse_neighborhood.py),
 *a vector of protein clusters and the desired genomic neighborhood clustering method.*/
/*INCOMPLETE*/
void genome_clustering(std::string neighborhoods_file, UndirectedEdgeWeightedGraph<std::string> &clusters,
                       std::string method, double stringency);

#endif
