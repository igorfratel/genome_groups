#include <vector>
#include <thread>

#include "genome_grouping.h"

// Dummy filename for the pairings option
const std::string DUMMY_FILENAME = "&";

/**
 * Receives a string and delimiters.
 * Splits the string in a vector according to delimiters
 */
static std::vector<std::string> split(const std::string& in, const std::string& delim) {
  std::string::size_type start = in.find_first_not_of(delim), end = 0;

  std::vector<std::string> out;
  while(start != in.npos) {
    end = in.find_first_of(delim, start);
    if (end == in.npos) {
      out.push_back(in.substr(start));
      break;
    } else {
      out.push_back(in.substr(start, end-start));
    }
    start = in.find_first_not_of(delim, end);
  }
  return out;
}

/**
 * Receives a genomic neighborhood filename.
 * Returns a vector of genomic neighborhoods, filled with the information from the file.
 */
std::vector<GenomicNeighborhood> parse_neighborhoods(const std::string &neighborhoods_filename) {
  std::ifstream file;
  std::string line;
  std::vector<std::string> split_line;
  std::vector<GenomicNeighborhood> neighborhoods;

  file.open(neighborhoods_filename.c_str());
  if (file.fail()) {
    std::cerr << "ERROR: trouble opening the neighborhoods file\n";
    exit(1);
  }

  int organism_index = -1;
  while(std::getline(file, line)) {
    split_line = split(line, " \t");

    if (split_line[0] == "ORGANISM") { // beginning of organism
      organism_index++;
      for (unsigned int accession_index = 0; accession_index < split_line.size(); accession_index++){ // find index of accession code
        if (split_line[accession_index] == "accession") {
          accession_index += 3;
          neighborhoods.emplace_back(GenomicNeighborhood(split_line[accession_index]));
          break;
        }
      }
    }
    else if (split_line[0] == "." && split_line[1] != "cds") // protein
      neighborhoods[organism_index].add_protein(split_line[7], split_line[4], split_line[1]); // locus pid cds

    else if (split_line[0] == "-->") { // seed protein
      neighborhoods[organism_index].add_protein(split_line[7], split_line[4], split_line[1]); // locus pid cds
      neighborhoods[organism_index].add_seed(split_line[7], split_line[4], split_line[1]); // locus pid cds
    }
  }
  file.close();
  return neighborhoods;
}
/**
 *Prints the score between two genomic neighborhoods in the standard format
 */
static void output_score(GenomicNeighborhood &g1, GenomicNeighborhood &g2, double score, std::ofstream &output_file) {
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
static void output_pairings(GenomicNeighborhood &g1, GenomicNeighborhood &g2,
    std::map<std::pair<int, int>,int> &assignments, std::ofstream &pairings_file) {

  // Writes header
  pairings_file << ">" << g1.get_accession() << "\t" <<
    g1.get_first_cds() << "\t" <<
    g1.get_last_cds() << "\t" <<
    g2.get_accession() << "\t" <<
    g2.get_first_cds() << "\t" <<
    g2.get_last_cds() << "\n";

  // Writes pairings
  for (std::map<std::pair<int, int>,int>::iterator it = assignments.begin(); it != assignments.end(); ++it){
    pairings_file << g1.get_pid(it->first.first) << "\t" <<
      g2.get_pid(it->first.second) << "\t" <<
      ((double)it->second)/1000000 << "\n";
  }
}

/**
 *Prints the chosen protein assignments to the pairings_file (treats each assignment as a pair of pairs of proteins)
 */
static void output_pairingsO2(GenomicNeighborhood &g1, GenomicNeighborhood &g2,
    std::map<std::pair<int, int>,int> &assignments, std::ofstream &pairings_file) {

  // Writes header
  pairings_file << ">" << g1.get_accession() << "\t" <<
    g1.get_first_cds() << "\t" <<
    g1.get_last_cds() << "\t" <<
    g2.get_accession() << "\t" <<
    g2.get_first_cds() << "\t" <<
    g2.get_last_cds() << "\n";

  // Writes pairings
  for (std::map<std::pair<int, int>,int>::iterator it = assignments.begin(); it != assignments.end(); ++it){
    pairings_file << g1.get_pid(it->first.first) << "\t" <<
      g1.get_pid(it->first.first + 1) << "\t" <<
      g2.get_pid(it->first.second) << "\t" <<
      g2.get_pid(it->first.second + 1) << "\t" <<
      ((double)it->second)/1000000 << "\n";
  }
}

/**
 *Receives a vector of genomic neighborhoods and returns the number of unique proteins in them.
 */
int total_protein_count(std::vector<GenomicNeighborhood> &neighborhoods) {
  std::set<std::string> protein_set;
  for (unsigned int i = 0; i < neighborhoods.size(); i++) {
    for (GenomicNeighborhood::iterator it = neighborhoods[i].begin(); it != neighborhoods[i].end(); ++it) {
      protein_set.emplace(it->pid);
    }
  }
  return protein_set.size();
}

/**
 *Receives a vector of genomic neighborhoods,
 *a ProteinCollection and the desired genomic neighborhood clustering method.
 *Writes the similarity between all genomic neighborhoods on the genome_sim_filename and,
 *optionally, the pairings made between their proteins on the pairings_filename.
 */
void genome_clustering(std::vector<GenomicNeighborhood> &neighborhoods, ProteinCollection &clusters,
    const std::string &method, double prot_stringency, double neigh_stringency, size_t num_threads, const std::string &genome_sim_filename,
    const std::string &pairings_filename) {

  std::ofstream output_file;
  if (genome_sim_filename == "-")
    output_file.basic_ios<char>::rdbuf(std::cout.rdbuf());
  else
    output_file = std::ofstream(genome_sim_filename.c_str());

  std::ofstream pairings_file;
  if (pairings_filename != DUMMY_FILENAME)
    pairings_file = std::ofstream(pairings_filename.c_str());

  if (method == "porthodom") {
    using Assignment = std::map<std::pair<int,int>, int>;

    size_t num_sequences = 1600/* neighborhoods.size() */;

    std::vector<std::vector<Assignment> > assignments(
        num_sequences, std::vector<Assignment>(num_sequences));
    std::vector<std::vector<double> > score(
        num_sequences, std::vector<double>(num_sequences));

    size_t m = 0;
    size_t n = 0;

    // #pragma omp parallel for collapse(2) private(m, n)
    // for (m = 0; m < num_sequences#<{(| neighborhoods.size() |)}>#; m++) {
    //   for (n = 0#<{(| m + 1 |)}>#; n < num_sequences#<{(| neighborhoods.size() |)}>#; n++) {
    //     assignments[m][n] = porthodom_assignments(
    //         neighborhoods[m], neighborhoods[n], clusters, prot_stringency);
    //
    //     score[m][n] = porthodom_scoring(assignments[m][n],
    //         std::max(neighborhoods[m].protein_count(),
    //           neighborhoods[n].protein_count()));
    //   }
    // }

    auto max = num_threads;

    std::vector<std::thread > threads;
    for (size_t k = 0; k < max; k++) {
      threads.emplace_back([k, &num_sequences, &neighborhoods, &max,
                               &clusters, &prot_stringency,
                               &score, &assignments](){
        auto num_jobs = (num_sequences*(num_sequences-1))/2;
        auto left_limit = (num_jobs/max)*k;
        auto right_limit = (num_jobs/max)*(k+1);

        for (size_t l = left_limit; l < right_limit; l++) {
          auto m = l / num_sequences;
          auto n = l % num_sequences;

          assignments[m][n] = porthodom_assignments(
              neighborhoods[m], neighborhoods[n], clusters, prot_stringency);

          score[m][n] = porthodom_scoring(assignments[m][n],
              std::max(neighborhoods[m].protein_count(),
                neighborhoods[n].protein_count()));
        }

        // auto left_limit = (num_sequences/max)*k;
        // auto right_limit = (num_sequences/max)*(k+1);
        //
        // for (size_t m = left_limit; m < right_limit; m++) {
        //   for (size_t n = m + 1; n < num_sequences; n++) {
        //     assignments[m][n] = porthodom_assignments(
        //         neighborhoods[m], neighborhoods[n], clusters, prot_stringency);
        //
        //     score[m][n] = porthodom_scoring(assignments[m][n],
        //         std::max(neighborhoods[m].protein_count(),
        //           neighborhoods[n].protein_count()));
        //   }
        // }
      });
    }

    for (auto& thread : threads)
      thread.join();

    for (size_t m = 0; m < num_sequences; m++) {
      for (size_t n = m + 1; n < num_sequences; n++) {
        if (score[m][n] >= neigh_stringency) {
          output_score(neighborhoods[m], neighborhoods[n],
              score[m][n], output_file);
        }

        if (pairings_filename != DUMMY_FILENAME) {
          output_pairings(neighborhoods[m], neighborhoods[n],
            assignments[m][n], pairings_file);
        }
      }
    }
  }

  else if (method == "porthodomO2") {
    std::map<std::pair<int,int>, int> assignments;
    double score;
    for (unsigned int m = 0; m < neighborhoods.size(); m++) {

      if (neighborhoods[m].protein_count() == 1) continue; // Ignores neighborhoods with less than 2 proteins

      for (unsigned int n = m + 1; n < neighborhoods.size(); n++) {

        if (neighborhoods[n].protein_count() == 1) continue; // Ignores neighborhoods with less than 2 proteins

        // Edges chosen by the algorithm
        assignments = porthodomO2_assignments(neighborhoods[m], neighborhoods[n], clusters, prot_stringency);

        // apply the scoring formula
        score = porthodomO2_scoring(assignments,
            std::max(neighborhoods[m].protein_count(), neighborhoods[n].protein_count()) - 1);

        if (score < neigh_stringency) continue; // ignore scores below stringency
        // Writes scores to output_file
        output_score(neighborhoods[m], neighborhoods[n], score, output_file);

        if (pairings_filename == "&") continue; // Dummy filename indicating this option was not chosen
        // Writes pairing to pairings_file
        output_pairingsO2(neighborhoods[m], neighborhoods[n], assignments, pairings_file);
      }
    }
  }
}
