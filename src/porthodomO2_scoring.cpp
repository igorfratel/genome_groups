#include "porthodomO2_scoring.h"

/**
 *Receives two pairs of proteins, a ProteinCollection and a threshold value.
 *The score returned is the average of the similarities between the pairs.
 */
static int clustering_value(protein_info_t prot_g1_1, protein_info_t prot_g1_2, protein_info_t prot_g2_1,
                            protein_info_t prot_g2_2, ProteinCollection &clusters, double stringency) {

  double result = (clusters.get_similarity(prot_g1_1.pid, prot_g2_1.pid) +
                   clusters.get_similarity(prot_g1_2.pid, prot_g2_2.pid))/2;
  if (result >= stringency) //Should I check if both similarities are greater and then apply the mean? !!!
    return (int)(1000000*result);
  else
    return 0;

}

/**
 *Receives two genomic neighborhoods, g1 and g2, and the ProteinCollection.
 *Fills an integer matrix where matrix[i][j] is the similarity measure between the i-th pair of proteins of g1 and
 *the j-th pair of proteins of g2.
 */
static std::vector<std::vector<int> > fill_assignment_matrix(GenomicNeighborhood &g1, GenomicNeighborhood &g2,
                                                     ProteinCollection &clusters, double stringency) {

    int i = 0;
    int j = 0;
    std::vector<std::vector<int> > matrix(g1.protein_count() - 1, std::vector<int> (g2.protein_count()));

    for(GenomicNeighborhood::iterator it = g1.begin(), it_last = --g1.end(); it != it_last; ++it) {
        j = 0;
        for(GenomicNeighborhood::iterator it2 = g2.begin(), it2_last = --g2.end(); it2 != it2_last; ++it2) {
            matrix[i][j] = clustering_value(*it, *(std::next(it)), *it2, *(std::next(it2)), clusters, stringency);
            j++;
        }
        i++;
    }
    return matrix;
}

/**
 *Receives two genomic neighborhoods, a ProteinCollection and the protein stringency.
 *Returns the MWM porthodom  O2 protein assignments between the two neighborhoods
 */
std::map<std::pair<int, int>, int> porthodomO2_assignments(GenomicNeighborhood &g1, GenomicNeighborhood &g2,
                                                           ProteinCollection &clusters, double prot_stringency) {
    //DEBUG
    /*std::cout <<"Comparing (" << g1.get_accession() << ") and "
                << "(" << g2.get_accession() << ", ):\n";*/

    std::map<std::pair<int, int>, int> assignments;
    std::vector<std::vector<int> > matrix = fill_assignment_matrix(g1, g2, clusters, prot_stringency);
    Hungarian my_hungarian (matrix, matrix.size(), matrix[0].size(), HUNGARIAN_MODE_MAXIMIZE_UTIL);

    my_hungarian.solve();
    assignments = my_hungarian.get_assignments();
    return assignments;
}

/**
 *Receives the porthodom assignments and a normalizing factor (length of the longest neighborhood).
 *Returns the porthodom MWM_O2 score (that takes order in consideration).
 */
double porthodomO2_scoring(std::map<std::pair<int, int>, int> &assignments, int length) {
    double score = 0;
    for (std::map<std::pair<int, int>,int>::iterator it = assignments.begin(); it != assignments.end(); ++it)
        score += ((double)it->second)/1000000; //Division to undo the multiplication in clustering_value()

    //apply the scoring formula
    return score/length;
}
