#ifndef __PROTEIN_GROUPING_H__
#define __PROTEIN_GROUPING_H__

#include <string>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include "UndirectedEdgeWeightedGraph.h"

/*Receives the user's preferred protein homology/orthology detection method and runs it on a file
 *already formatted to be its input.
 *Writes the results to "homology_detection_output.txt" and returns this filename.*/
std::string homology_detection(char *format_file, char *method);

/*Receives the similarities file and stores them in a list of clusters.
 *The similarities file must be in the format "prot1 prot2 sim" in every line
 *num_prot must be the number of proteins and stringency is the minimum similarity
 *for two proteins to be considered part of the same cluster*/
std::vector<std::vector<std::string>> protein_clustering(char* prot_sim_file, int num_prot,
                                                         double stringency);

#endif
