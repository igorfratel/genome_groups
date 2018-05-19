#include "porthodom_scoring.h"

/**
 *Prints the score between two genomic neighborhoods in the standard format
 */
void porthodom_output_score(const GenomicNeighborhood &g1, const GenomicNeighborhood &g2, double score, std::ofstream &output_file) {
    output_file << g1.get_accession() << "\t" <<
                   g1.get_first_cds() << "\t" <<
                   g1.get_last_cds() << "\t" <<
                   g2.get_accession() << "\t" <<
                   g2.get_first_cds() << "\t" <<
                   g2.get_last_cds() << "\t" <<
                   score << "\n";
}

/**
 *Prints the chosen protein assignments to the pairings_file
 */
void porthodom_output_pairings(const GenomicNeighborhood &g1, const GenomicNeighborhood &g2,
                               const std::map<std::pair<int, int>,int> &assignments, std::ofstream &pairings_file) {

    //Writes header
    pairings_file << ">" << g1.get_accession() << "\t" <<
                            g1.get_first_cds() << "\t" <<
                            g1.get_last_cds() << "\t" <<
                            g2.get_accession() << "\t" <<
                            g2.get_first_cds() << "\t" <<
                            g2.get_last_cds() << "\n";

    //Writes pairings
    for (std::map<std::pair<int, int>,int>::const_iterator it = assignments.begin(); it != assignments.end(); ++it){
        pairings_file << g1.get_pid(it->first.first) << "\t" <<
                         g2.get_pid(it->first.second) << "\t" <<
                         ((double)it->second)/1000000 << "\n";
    }
}

/**
 *Receives two proteins, a ProteinCollection, a threshold value and a mismatch value.
 *If the proteins are connected in the ProteinCollection, return an integer with 1000000x their similarity (because
 *the Hungarian class only works with integers).
 *If the proteins aren't connected or if their similarity is smaller than the threshold value, returns mismatch.
 */
static int clustering_value(protein_info_t my_prot, protein_info_t my_prot2, const ProteinCollection &clusters,
                            double stringency, const int mismatch) {
  if (!clusters.are_connected(my_prot.pid, my_prot2.pid)) return (int)(1000000*mismatch);
  double similarity = clusters.get_similarity(my_prot.pid, my_prot2.pid);
  if (similarity >= stringency) return (int)(1000000*similarity);
  return (int)(1000000*mismatch);

}

/**
 *Receives two protein_info_t vectors, g1 and g2, the ProteinCollection, the protein stringency and a mismatch value.
 *Fills an integer matrix where matrix[i][j] is the similarity measure between the i-th protein of g1 and
 *the j-th protein of g2 if they are in the same cluster (connected in the ProteinCollection) and 0 otherwise.
 */
static std::vector<std::vector<int> > fill_assignment_matrix(const std::vector<protein_info_t> &g1,
                                                             const std::vector<protein_info_t> &g2,
                                                             const ProteinCollection &clusters, double stringency,
                                                             const int mismatch) {

    std::vector<std::vector<int> > matrix(g1.size(), std::vector<int>(g2.size()));

    for(size_t i = 0; i < g1.size(); i++)
        for(size_t j = 0; j < g2.size(); j++)
            matrix[i][j] = clustering_value(g1[i], g2[j], clusters, stringency, mismatch);

    return matrix;
}

/**
 *Receives two protein_info_t vectors a ProteinCollection, the protein stringency and a mismatch value.
 *Returns the MWM porthodom protein assignments between the two neighborhoods (vectors).
 */
std::map<std::pair<int, int>, int> porthodom_assignments(const std::vector<protein_info_t> &g1,
                                                         const std::vector<protein_info_t> &g2,
                                                         const ProteinCollection &clusters, double prot_stringency,
                                                         const int mismatch) {

    std::map<std::pair<int, int>, int> assignments; //maps two protein coordinates to a score
    std::vector<std::vector<int> > matrix = fill_assignment_matrix(g1, g2, clusters, prot_stringency, mismatch);
    //The hungarian algorithm solves the MWM problem of which the porthodom algorithm is consisted.
    Hungarian my_hungarian (matrix, matrix.size(), matrix[0].size(), HUNGARIAN_MODE_MAXIMIZE_UTIL);

    my_hungarian.solve();
    return my_hungarian.get_assignments();
}

/**
 *Receives the porthodom assignments and a normalizing factor (length of the longest neighborhood).
 *Returns the porthodom MWM score.
 */
double porthodom_scoring(const std::map<std::pair<int, int>, int> &assignments, int length) {
    double score = 0;
    for (std::map<std::pair<int, int>,int>::const_iterator it = assignments.begin(); it != assignments.end(); ++it)
        score += ((double)it->second)/1000000; //Division to undo the multiplication in clustering_value()

    return score/length; //normalizing
}
