#include "porthodom_scoring.h"

/**
 *Prints the score between two genomic neighborhoods in the standard format
 */
void porthodom_output_score(GenomicNeighborhood &g1, GenomicNeighborhood &g2, double score, std::ofstream &output_file) {
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
void porthodom_output_pairings(GenomicNeighborhood &g1, GenomicNeighborhood &g2,
                            std::map<std::pair<int, int>,int> &assignments, std::ofstream &pairings_file) {

    //Writes header
    pairings_file << ">" << g1.get_accession() << "\t" <<
                            g1.get_first_cds() << "\t" <<
                            g1.get_last_cds() << "\t" <<
                            g2.get_accession() << "\t" <<
                            g2.get_first_cds() << "\t" <<
                            g2.get_last_cds() << "\n";

    //Writes pairings
    for (std::map<std::pair<int, int>,int>::iterator it = assignments.begin(); it != assignments.end(); ++it){
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
static int clustering_value(protein_info_t my_prot, protein_info_t my_prot2, ProteinCollection &clusters,
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
 *Receives two genomic neighborhoods, g1 and g2, and the ProteinCollection.
 *Fills an integer matrix where matrix[i][j] is the similarity measure between the i-th protein of g1 and
 *the j-th protein of g2 are in the same cluster (connected in the ProteinCollection) and 0 otherwise.
 */
static std::vector<std::vector<int> > fill_assignment_matrix(GenomicNeighborhood &g1, GenomicNeighborhood &g2,
                                                     ProteinCollection &clusters, double stringency) {

    int i = 0;
    int j = 0;
    std::vector<std::vector<int> > matrix(g1.protein_count(), std::vector<int>(g2.protein_count()));

    for(GenomicNeighborhood::iterator it = g1.begin(); it != g1.end(); ++it) {
        j = 0;
        for(GenomicNeighborhood::iterator it2 = g2.begin(); it2 != g2.end(); ++it2) {

            //DEBUG
            /*std::cout <<"matrix: " << i << " " << j << " " << it->pid << " " << it2->pid << "\n";*/
            matrix[i][j] = clustering_value(*it, *it2, clusters, stringency);
            //DEBUG
            //std::cout <<"matrix: " << i << " " << j << " " << it->pid << " " << it2->pid << " score = " << matrix[i][j]<<"\n";
            j++;
        }
        i++;
    }
    return matrix;
}

/*Receives two genomic neighborhoods and a ProteinCollection.
 *Returns the MWM porthodom protein assignments between the two neighborhoods
 */
std::map<std::pair<int, int>, int> porthodom_assignments(GenomicNeighborhood &g1, GenomicNeighborhood &g2,
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

/*Receives the porthodom assignments and a normalizing factor (length of the longest neighborhood).
 *Returns the porthodom MWM score.
 */

double porthodom_scoring(std::map<std::pair<int, int>, int> &assignments, int length) {
    double score = 0;
    for (std::map<std::pair<int, int>,int>::iterator it = assignments.begin(); it != assignments.end(); ++it)
        score += ((double)it->second)/1000000; //Division to undo the multiplication in clustering_value()

    //apply the scoring formula
    return score/length;

}
