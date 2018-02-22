#ifndef __PORTHODOM_SCORING_H__
#define __PORTHODOM_SCORING_H__

#include <map>
#include <string>
#include <algorithm>
#include "Hungarian.h"
#include "GenomicNeighborhood.h"
#include "ProteinCollection.h"


/*Receives two genomic neighborhoods a ProteinCollection.
 *Returns the MWM porthodom protein assignments between the two neighborhoods
 */
std::map<std::pair<int, int>, int> porthodom_assignments(GenomicNeighborhood &g1, GenomicNeighborhood &g2,
                         ProteinCollection &clusters, double prot_stringency);

/*Receives the porthodom assignments and a normalizing factor (length of the longest neighborhood).
 *Returns the porthodom MWM score.
 */
double porthodom_scoring(std::map<std::pair<int, int>, int> &assignments, int length);

#endif
