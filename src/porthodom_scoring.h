#ifndef __PORTHODOM_SCORING_H__
#define __PORTHODOM_SCORING_H__

#include <map>
#include <string>
#include <fstream>
#include <algorithm>
#include "Hungarian.h"
#include "GenomicNeighborhood.h"
#include "ProteinCollection.h"

/**
 *Prints the score between two genomic neighborhoods in the standard format
 */
void porthodom_output_score(GenomicNeighborhood &g1, GenomicNeighborhood &g2, double score, std::ofstream &output_file);

/**
 *Prints the chosen protein assignments to the pairings_file
 */
void porthodom_output_pairings(GenomicNeighborhood &g1, GenomicNeighborhood &g2,
                            std::map<std::pair<int, int>,int> &assignments, std::ofstream &pairings_file);

/*Receives two genomic neighborhoods a ProteinCollection and the protein stringency.
 *Returns the MWM porthodom protein assignments between the two neighborhoods
 */
std::map<std::pair<int, int>, int> porthodom_assignments(GenomicNeighborhood &g1, GenomicNeighborhood &g2,
                                                         ProteinCollection &clusters, double prot_stringency);

/*Receives the porthodom assignments and a normalizing factor (length of the longest neighborhood).
 *Returns the porthodom MWM score.
 */
double porthodom_scoring(std::map<std::pair<int, int>, int> &assignments, int length);

#endif
