#include "porthodomO2_scoring.h"

/**
 *Prints the score between two genomic neighborhoods in the standard format
 */
void porthodomO2_output_score(const GenomicNeighborhood &g1, const GenomicNeighborhood &g2, double score, std::ofstream &output_file) {
    output_file << g1.get_accession() << "\t" <<
                   g1.get_first_cds() << "\t" <<
                   g1.get_last_cds() << "\t" <<
                   g2.get_accession() << "\t" <<
                   g2.get_first_cds() << "\t" <<
                   g2.get_last_cds() << "\t" <<
                   score << "\n";
}

/**
 *Prints the chosen protein assignments to the pairings_file (treats each assignment as a pair of pairs of proteins)
 */
void porthodomO2_output_pairings(const GenomicNeighborhood &g1, const GenomicNeighborhood &g2,
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
                         g1.get_pid(it->first.first + 1) << "\t" <<
                         g2.get_pid(it->first.second) << "\t" <<
                         g2.get_pid(it->first.second + 1) << "\t" <<
                         ((double)it->second)/1000000 << "\n";
    }
}

/**
 *Receives two pairs of proteins, a ProteinCollection and a threshold value.
 *The score returned is the average of the similarities between the pairs.
 */
static int clustering_value(protein_info_t prot_g1_1, protein_info_t prot_g1_2, protein_info_t prot_g2_1,
                            protein_info_t prot_g2_2, const ProteinCollection &clusters, double stringency) {

  double result = (clusters.get_similarity(prot_g1_1.pid, prot_g2_1.pid) +
                   clusters.get_similarity(prot_g1_2.pid, prot_g2_2.pid))/2;
  if (result >= stringency) //Should I check if both similarities are greater and then apply the mean? !!!
    return (int)(1000000*result);
  else
    return 0;

}

/**
 *Receives two protein_info_t vectors, g1 and g2, and the ProteinCollection.
 *Fills an integer matrix where matrix[i][j] is the similarity measure between the i-th pair of proteins of g1 and
 *the j-th pair of proteins of g2.
 */
static std::vector<std::vector<int> > fill_assignment_matrix(const std::vector<protein_info_t> &g1,
                                                     const std::vector<protein_info_t> &g2,
                                                     const ProteinCollection &clusters, double stringency) {

    std::vector<std::vector<int> > matrix(g1.size(), std::vector<int>(g2.size()));

    for(size_t i = 0; i < g1.size() - 1; i++) {
        for(size_t j = 0; j < g2.size() - 1; j++) {
            //DEBUG
            /*std::cout <<"matrix: " << i << " " << j << " " << it->pid << " " << it2->pid << "\n";*/
            matrix[i][j] = clustering_value(g1[i], g1[i+1], g2[j], g2[j+1], clusters, stringency);
            //DEBUG
            //std::cout <<"matrix: " << i << " " << j << " " << it->pid << " " << it2->pid << " score = " << matrix[i][j]<<"\n";
        }
    }
    return matrix;
}

/**
 *Receives two protein_info_t vectors, a ProteinCollection and the protein stringency.
 *Returns the MWM porthodom  O2 protein assignments between the two neighborhoods (vectors)
 */
std::map<std::pair<int, int>, int> porthodomO2_assignments(const std::vector<protein_info_t> &g1,
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

/**
 *Receives the porthodom assignments and a normalizing factor (length of the longest neighborhood).
 *Returns the porthodom MWM_O2 score (that takes order in consideration).
 */
double porthodomO2_scoring(const std::map<std::pair<int, int>, int> &assignments, int length) {
    double score = 0;
    for (std::map<std::pair<int, int>,int>::const_iterator it = assignments.begin(); it != assignments.end(); ++it)
        score += ((double)it->second)/1000000; //Division to undo the multiplication in clustering_value()

    //apply the scoring formula
    return score/length;
}
