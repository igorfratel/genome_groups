#ifndef __PROTEIN_GROUPING_H__
#define __PROTEIN_GROUPING_H__

#include <string>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include "ProteinCollection.h"

/*Receives the user's preferred protein homology/orthology detection method and runs it on a file
 *already formatted to be its input.
 *Writes the results to prot_sim_filename*/
void homology_detection(const std::string &format_file, const std::string &method, const std::string &prot_sim_filename);

/*Receives the similarities file and stores them in a ProteinCollection.
 *The similarities file must be in the format "prot1 prot2 sim" in every line
 *num_prot must be the number of proteins and stringency is the minimum similarity
 *for two proteins to be considered part of the same cluster*/
ProteinCollection protein_clustering(const std::string &prot_sim_filename, unsigned int num_prot);

#endif
