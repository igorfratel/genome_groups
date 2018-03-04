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
 *Receives two proteins, a ProteinCollection and a threshold value.
 *If the proteins are connected in the ProteinCollection, return an integer with 100x their similarity (because
 *the Hungarian class only works with integers).
 *If the proteins aren't connected or if their similarity is smaller than the threshold value, return 0.
 */
static int clustering_value(protein_info_t my_prot, protein_info_t my_prot2, const ProteinCollection &clusters,
                            double stringency) {
  //DEBUG
  //std::cout << "clustering value between " << my_prot.pid << " and " << my_prot2.pid << " " << similarity << "\n";
  double similarity = clusters.get_similarity(my_prot.pid, my_prot2.pid);
  if (similarity >= stringency)
    return (int)(1000000*similarity);
  else
    return 0;

}

/**
 *Receives two protein_info_t vectors, g1 and g2, and the ProteinCollection.
 *Fills an integer matrix where matrix[i][j] is the similarity measure between the i-th protein of g1 and
 *the j-th protein of g2 are in the same cluster (connected in the ProteinCollection) and 0 otherwise.
 */
static std::vector<std::vector<int> > fill_assignment_matrix(const std::vector<protein_info_t> &g1,
                                                            const std::vector<protein_info_t> &g2,
                                                            const ProteinCollection &clusters, double stringency) {

    std::vector<std::vector<int> > matrix(g1.size(), std::vector<int>(g2.size()));

    for(size_t i = 0; i < g1.size(); i++) {
        for(size_t j = 0; j < g2.size(); j++) {
            //DEBUG
            /*std::cout <<"matrix: " << i << " " << j << " " << it->pid << " " << it2->pid << "\n";*/
            matrix[i][j] = clustering_value(g1[i], g2[j], clusters, stringency);
            //DEBUG
            //std::cout <<"matrix: " << i << " " << j << " " << it->pid << " " << it2->pid << " score = " << matrix[i][j]<<"\n";
        }
    }
    return matrix;
}

/*Receives two protein_info_t vectors and a ProteinCollection.
 *Returns the MWM porthodom protein assignments between the two neighborhoods (vectors)
 */
std::map<std::pair<int, int>, int> porthodom_assignments(const std::vector<protein_info_t> &g1,
                             const std::vector<protein_info_t> &g2,
                             const ProteinCollection &clusters, double prot_stringency) {

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

/*Receives the porthodom assignments and a normalizing factor (length of the longest neighborhood).
 *Returns the porthodom MWM score.
 */

double porthodom_scoring(const std::map<std::pair<int, int>, int> &assignments, int length) {
    double score = 0;
    for (std::map<std::pair<int, int>,int>::const_iterator it = assignments.begin(); it != assignments.end(); ++it)
        score += ((double)it->second)/1000000; //Division to undo the multiplication in clustering_value()

    //apply the scoring formula
    return score/length;

}
