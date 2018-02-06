#include "porthodom_scoring.h"

static int clustering_value(protein_info_t my_prot, protein_info_t my_prot2, ProteinCollection &clusters,
                            double stringency) {
  /*Receives two proteins, a ProteinCollection and a threshold value.
   *If the proteins are connected in the ProteinCollection, return an integer with 100x their similarity (because
    the Hungarian class only works with integers).
   *If the proteins aren't connected or if their similarity is smaller than the threshold value, return 0.*/
  double similarity = clusters.get_similarity(my_prot.pid, my_prot2.pid);
  if (similarity >= stringency)
    return (int)(100*similarity);
  else
    return 0;

}

static std::vector<std::vector<int> > fill_assignment_matrix(GenomicNeighborhood g1, GenomicNeighborhood g2,
                                                     ProteinCollection &clusters, double stringency) {
    /*Receives two genomic neighborhoods, g1 and g2, and the ProteinCollection.
     *Fills an integer matrix where matrix[i][j] is the similarity measure between the i-th protein of g1 and
     *the j-th protein of g2 are in the same cluster (connected in the ProteinCollection) and 0 otherwise.*/

    int i = 0;
    int j = 0;
    std::vector<std::vector<int> > matrix;

    matrix.resize(g1.protein_count());
    for (unsigned int k = 0; k < matrix.size(); k++)
        matrix[k].resize(g2.protein_count());

    for(GenomicNeighborhood::iterator it = g1.begin(); it != g1.end(); ++it) {
        j = 0;
        for(GenomicNeighborhood::iterator it2 = g2.begin(); it2 != g2.end(); ++it2) {
            matrix[i][j] = clustering_value(*it, *it2, clusters, stringency);
	        //DEBUG
            //std::cout <<"matrix: " << i << " " << j << " " << it->pid << " " << it2->pid << " score = " << matrix[i][j]<<"\n";
            j++;
        }
        i++;
    }
    return matrix;
}

double porthodom_scoring(GenomicNeighborhood g1, GenomicNeighborhood g2,
                             ProteinCollection &clusters, double stringency) {
    /*Receives two genomic neighborhoods and a ProteinCollection.
     *Returns the MWM porthodom score between the two neighborhoods
     *(Using the hungarian algorithm and the porthodom scoring formula).*/

    double sum_temp = 0;
    std::map<std::pair<int, int>, int> assignments;
    std::vector<std::vector<int> > matrix = fill_assignment_matrix(g1, g2, clusters, stringency);
    Hungarian my_hungarian (matrix, matrix.size(), matrix[0].size(), HUNGARIAN_MODE_MAXIMIZE_UTIL);

    my_hungarian.solve();
    assignments = my_hungarian.get_assignments();
    //DEBUG
    /*std::cout <<"Assignments between (" << g1.get_accession() << ", " << g1.get_organism() << ") and "
    << "(" << g2.get_accession() << ", " << g2.get_organism() << "):\n";
    my_hungarian.print_assignment();
    my_hungarian.print_cost();
    fprintf(stderr, "\n");*/
    //PORTHODOM similarity measure calculation, adapted for genome neighborhoods
    //The similarity between to neighborhoods is given by the normalized sum of the
    //similarities between each pair of "assigned" proteins (as given by the Hungarian algorithm)
    //TO compute:
    //    1. add the assignment values
    //    2. divide by the number of proteins in the "largest" genome
    for (std::map<std::pair<int, int>,int>::iterator it = assignments.begin(); it != assignments.end(); ++it)
        sum_temp += ((double)it->second)/100; //Division to undo the multiplication in clustering_value()
        //DEBUG
        /*std::cout << "SUM_TEMP: " << sum_temp << "\n";*/
    return sum_temp/std::max(g1.protein_count(), g2.protein_count());
}
