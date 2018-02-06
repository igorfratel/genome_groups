#ifndef __PORTHODOM_SCORING_H__
#define __PORTHODOM_SCORING_H__

#include <map>
#include <string>
#include <algorithm>
#include "Hungarian.h"
#include "GenomicNeighborhood.h"
#include "ProteinCollection.h"


/*Receives two genomic neighborhoods and a ProteinCollection.
 *Returns the MWM porthodom score between the two neighborhoods
 *(Using the hungarian algorithm and the porthodom scoring formula).*/
double porthodom_scoring(GenomicNeighborhood g1, GenomicNeighborhood g2,
                             ProteinCollection &clusters, double stringency);


#endif
