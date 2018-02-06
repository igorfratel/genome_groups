#ifndef __PORTHODOMO2_SCORING_H__
#define __PORTHODOMO2_SCORING_H__

#include <map>
#include <string>
#include <algorithm>
#include <iterator>
#include "Hungarian.h"
#include "GenomicNeighborhood.h"
#include "ProteinCollection.h"

/*Receives two genomic neighborhoods and a ProteinCollection.
 *Returns the MWM_O2 porthodom score between the two neighborhoods
 *(Using the hungarian algorithm and the porthodom scoring formula).*/
double porthodomO2_scoring(GenomicNeighborhood g1, GenomicNeighborhood g2,
                             ProteinCollection &clusters, double stringency);

#endif
