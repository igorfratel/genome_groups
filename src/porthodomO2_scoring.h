#ifndef __PORTHODOMO2_SCORING_H__
#define __PORTHODOMO2_SCORING_H__

#include <map>
#include <string>
#include <algorithm>
#include <iterator>
#include "Hungarian.h"
#include "GenomicNeighborhood.h"
#include "ProteinCollection.h"

/**
 *Receives the porthodom assignments and a normalizing factor (length of the longest neighborhood).
 *Returns the porthodom O2 MWM score (that takes order in consideration).
 */
double porthodomO2_scoring(std::map<std::pair<int, int>, int> &assignments, int length);

/**
 *Receives two genomic neighborhoods, a ProteinCollection and the protein stringency.
 *Returns the MWM porthodom  O2 protein assignments between the two neighborhoods
 */
std::map<std::pair<int, int>, int> porthodomO2_assignments(GenomicNeighborhood &g1, GenomicNeighborhood &g2,
                                                           ProteinCollection &clusters, double prot_stringency);
#endif
