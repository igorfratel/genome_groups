#ifndef __PROTEIN_GROUPING_H__
#define __PROTEIN_GROUPING_H__

#include <string>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include "ProteinCollection.h"

/*Receives the user's preferred protein homology/orthology detection method and runs it on a file
 *already formatted to be its input, only displaying the results above the given stringency.
 *Writes the results to prot_sim_filename*/
void homology_detection(std::string format_file, std::string method, std::string prot_sim_filename, double stringency);

/*Receives the similarities file and stores them in a ProteinCollection.
 *The similarities file must be in the format "prot1 prot2 sim" in every line*/
ProteinCollection protein_clustering(std::string prot_sim_filename);

#endif
