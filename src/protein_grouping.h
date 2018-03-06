#ifndef __PROTEIN_GROUPING_H__
#define __PROTEIN_GROUPING_H__

#include <string>
#include <stdlib.h>
#include <fstream>
#include <set>
#include <vector>
#include "ProteinCollection.h"

/**
 *Receives the protein similarities file and returns the number of unique proteins in it.
 */
int total_protein_count(const std::string &prot_sim_filename);


/**
 *Receives the user's preferred protein homology/orthology detection method and runs it on a file
 *already formatted to be its input (that has the number of residues in its unique protein set equal to num_residues).
 *Writes the results to prot_sim_filename.
 */
void homology_detection(const std::string &format_file, const std::string &method, const std::string &num_residues,
                        const std::string &prot_sim_filename);

/**
 *Receives the protein similarities file and the number of unique proteins.
 *Stores the proteins and their similarity relationships in a ProteinCollection.
 *The similarities file must be in the format "prot1 prot2 sim" in every line (space separated);
 */
ProteinCollection protein_clustering(const std::string &prot_sim_filename, unsigned int num_prot);

#endif
