#include "global_alignment.h"

/**
 *Receives two proteins, a ProteinCollection, a threshold value and the mismatch value.
 *If the proteins are connected in the ProteinCollection, returns their similarity
 *If the proteins aren't connected or if their similarity is smaller than the threshold value, return mismatch.
 */
static double clustering_value(protein_info_t my_prot, protein_info_t my_prot2, const ProteinCollection &clusters,
                               double stringency, const int mismatch) {

  if (!clusters.are_connected(my_prot.pid, my_prot2.pid)) return (double)mismatch;
  double similarity = clusters.get_similarity(my_prot.pid, my_prot2.pid);
  if (similarity >= stringency) return similarity;
  return (double)mismatch;
}

/**
 *Receives two protein_info_t vectors, a ProteinCollection, the protein stringency, gap score and mismatch value.
 *Computes the global similarity (Needleman-Wunsch algorithm) between the sequences (vectors) using dynamic programming.
 *Returns the dynamic programming matrix.
 */
static std::vector<std::vector<double> > compute_sim (const std::vector<protein_info_t> &seq_1,
                                                      const std::vector<protein_info_t> &seq_2,
                                                      const ProteinCollection &clusters, double stringency,
                                                      const int gap_score, const int mismatch) {

    double match;
    double options[3];
    std::vector<std::vector<double> > dyn_matrix(seq_1.size() + 1, std::vector<double>(seq_2.size() + 1));

    for (size_t i = 0; i <= seq_1.size(); i++)
        dyn_matrix[i][0] = (double)i * gap_score;
    for (size_t i = 0; i <= seq_2.size(); i++)
        dyn_matrix[0][i] = (double)i * gap_score;
    for (size_t i = 1; i <= seq_1.size(); i++) {
        for (size_t j = 1; j <= seq_2.size(); j++) {
            match = clustering_value(seq_1[i-1], seq_2[j-1], clusters, stringency, mismatch);
            options[0] = dyn_matrix[i-1][j] + (double)gap_score;
            options[1] = dyn_matrix[i-1][j-1] + match;
            options[2] = dyn_matrix[i][j-1] + (double)gap_score;
            dyn_matrix[i][j] = *std::max_element(std::begin(options), std::end(options));
        }
    }
    return dyn_matrix;
}

/**
 *Recursive.
 *Receives the dynamic programming matrix,
 *two integers (the sizes of the first and second sequences),
 *an int len = 0, two string vectors where the aligned sequences will be
 *stored, the two sequences to be aligned (protein_info_t vectors), a ProteinCollection, the protein stringency,
 *the gap score and a mismatch value.
 *Recovers the upmost optimal alignment between the two sequences and stores it in align_1 and align_2.
 */
static void align(const std::vector<std::vector<double> > &dyn_matrix, int i, int j, int len, std::vector<std::string> &align_1,
                  std::vector<std::string> &align_2,  const std::vector<protein_info_t> &seq_1,
                  const std::vector<protein_info_t> &seq_2, const ProteinCollection &clusters, double stringency,
                  const int gap_score, const int mismatch) {

    if (i == 0 && j == 0)
        len = 0;

    else if (i > 0 && dyn_matrix[i][j] == dyn_matrix[i-1][j] + gap_score) {
        align(dyn_matrix, i - 1, j, len, align_1, align_2, seq_1, seq_2, clusters, stringency, gap_score, mismatch);
        len++;
        align_1.push_back(seq_1[i-1].pid);
        align_2.push_back("-");
    }
    else if (i > 0 && j > 0 && dyn_matrix[i][j] == dyn_matrix[i-1][j-1] +
             clustering_value(seq_1[i-1], seq_2[j-1], clusters, stringency, mismatch)) {
        align(dyn_matrix, i - 1, j - 1, len, align_1, align_2, seq_1, seq_2, clusters, stringency, gap_score, mismatch);
        len++;
        align_1.push_back(seq_1[i-1].pid);
        align_2.push_back(seq_2[j-1].pid);
    }
    else {
        align(dyn_matrix, i, j - 1, len, align_1, align_2, seq_1, seq_2, clusters, stringency, gap_score, mismatch);
        len++;
        align_1.push_back("-");
        align_2.push_back(seq_2[j-1].pid);
    }
}

/**
 *Returns the dynamic programming matrix used in the global alignment algorithm (Needleman-Wunsch).
 */
std::vector<std::vector<double> > global_alignment_assignments(const std::vector<protein_info_t> &g1,
                                                               const std::vector<protein_info_t> &g2,
                                                               const ProteinCollection &clusters, double prot_stringency,
                                                               const int gap_score, const int mismatch) {
    return compute_sim(g1, g2, clusters, prot_stringency, gap_score, mismatch);
}


/**
 *Prints the score between two genomic neighborhoods in the standard format.
 */
void global_alignment_output_score(const GenomicNeighborhood &g1, const GenomicNeighborhood &g2, double score,
                                   std::ofstream &output_file) {

    output_file << g1.get_accession() << "\t" <<
                   g1.get_first_cds() << "\t" <<
                   g1.get_last_cds() << "\t" <<
                   g2.get_accession() << "\t" <<
                   g2.get_first_cds() << "\t" <<
                   g2.get_last_cds() << "\t" <<
                   score << "\n";
}

/**
 *Prints the aligned proteins.
 */
void global_alignment_output_pairings(const GenomicNeighborhood &g1, const GenomicNeighborhood &g2,
                                      const ProteinCollection &clusters, const std::vector<std::vector<double> > &assignments,
                                      double stringency, std::ofstream &pairings_file, const int gap_score, const int mismatch) {

    std::vector<std::string> align_1, align_2;
    align(assignments, g1.protein_count(), g2.protein_count(), 0, align_1, align_2, g1.get_protein_vector(),
          g2.get_protein_vector(), clusters, stringency, gap_score, mismatch);

    //Writes header
    pairings_file << ">" << g1.get_accession() << "\t" <<
                            g1.get_first_cds() << "\t" <<
                            g1.get_last_cds() << "\t" <<
                            g2.get_accession() << "\t" <<
                            g2.get_first_cds() << "\t" <<
                            g2.get_last_cds() << "\n";

    //Writes pairings
    for (size_t i = 0; i < align_1.size(); i++)
        pairings_file << align_1[i] << "\t" << align_2[i] << "\n";
}
