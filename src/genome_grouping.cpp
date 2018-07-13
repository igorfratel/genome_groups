#include "genome_grouping.h"

/**
 *Receives a genomic neighborhood filename.
 *Returns a vector of genomic neighborhoods, filled with the information from the file.
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

         if (split_line[0] == "ORGANISM") { //beginning of organism
             organism_index++;
             for (unsigned int accession_index = 0; accession_index < split_line.size(); accession_index++){ //find index of accession code
                 if (split_line[accession_index] == "accession") {
                     accession_index += 3;
                     neighborhoods.emplace_back(GenomicNeighborhood(split_line[accession_index]));
                     break;
                }
             }
         }
 		else if (split_line[0] == "." && split_line[1] != "cds") //protein
             neighborhoods[organism_index].add_protein(split_line[7], split_line[4], split_line[1]); //locus pid cds

 		else if (split_line[0] == "-->") { //seed protein
             neighborhoods[organism_index].add_protein(split_line[7], split_line[4], split_line[1]); //locus pid cds
             neighborhoods[organism_index].add_seed(split_line[7], split_line[4], split_line[1]); //locus pid cds
         }
     }
     file.close();
     return neighborhoods;
 }

 /**
  *Receives a vector of genomic neighborhoods,
  *a ProteinCollection, the desired genomic neighborhood clustering method, a
  *gap score value (for when an alignment method is chosen) and a mismatch value.
  *Writes the similarity between all genomic neighborhoods on the genome_sim_filename and,
  *optionally, the pairings made between their proteins on the pairings_filename.
  */
void genome_clustering(const std::vector<GenomicNeighborhood> &neighborhoods, const ProteinCollection &clusters,
                       const std::string &method, double prot_stringency, double neigh_stringency,
                       const std::string &genome_sim_filename, const std::string &pairings_filename, const int gap_score,
                       const int mismatch) {

    std::ofstream output_file;
    if(genome_sim_filename == "-")
        output_file.basic_ios<char>::rdbuf(std::cout.rdbuf());
    else
        output_file = std::ofstream(genome_sim_filename.c_str());

    std::ofstream pairings_file;
    if(pairings_filename != "&") //Dummy filename indicating this option was not chosen
        pairings_file = std::ofstream(pairings_filename.c_str());

    if (method == "porthodom") {
        std::map<std::pair<int,int>, int> assignments;
        double score;
        for(size_t m = 0; m < neighborhoods.size(); m++) {
            for (size_t n = m + 1; n < neighborhoods.size(); n++) {
                //Edges chosen by the algorithm
                assignments = porthodom_assignments(neighborhoods[m].get_protein_vector(),
                                                    neighborhoods[n].get_protein_vector(), clusters, prot_stringency,
                                                    mismatch);

                //apply the scoring formula
                score = porthodom_scoring(assignments,
                                          std::max(neighborhoods[m].protein_count(), neighborhoods[n].protein_count()));

                if (score < neigh_stringency) continue; //ignore scores below stringency
                //Writes scores to output_file
                porthodom_output_score(neighborhoods[m], neighborhoods[n], score, output_file);

                if (pairings_filename == "&") continue; //Dummy filename indicating this option was not chosen
                //Writes pairing to pairings_file
                porthodom_output_pairings(neighborhoods[m], neighborhoods[n], assignments, pairings_file);
            }
        }
    }

    else if (method == "porthodomO2") {
        std::map<std::pair<int,int>, int> assignments;
        double score;
        for (unsigned int m = 0; m < neighborhoods.size(); m++) {

            if(neighborhoods[m].protein_count() == 1) continue; //Ignores neighborhoods with less than 2 proteins

            for (unsigned int n = m + 1; n < neighborhoods.size(); n++) {

                if(neighborhoods[n].protein_count() == 1) continue; //Ignores neighborhoods with less than 2 proteins

                //Edges chosen by the algorithm
                assignments = porthodomO2_assignments(neighborhoods[m].get_protein_vector(),
                                                      neighborhoods[n].get_protein_vector(), clusters, prot_stringency,
                                                      mismatch);

                //apply the scoring formula
                score = porthodomO2_scoring(assignments,
                                          std::max(neighborhoods[m].protein_count(), neighborhoods[n].protein_count()) - 1);

                if (score < neigh_stringency) continue; //ignore scores below stringency
                //Writes scores to output_file
                porthodomO2_output_score(neighborhoods[m], neighborhoods[n], score, output_file);

                if (pairings_filename == "&") continue; //Dummy filename indicating this option was not chosen
                //Writes pairing to pairings_file
                porthodomO2_output_pairings(neighborhoods[m], neighborhoods[n], assignments, pairings_file);
            }
        }
    }

    else if (method == "global_alignment") {
        std::vector<std::vector<double> > assignments;
        double score;
        for (unsigned int m = 0; m < neighborhoods.size(); m++) {
            for (unsigned int n = m + 1; n < neighborhoods.size(); n++) {

                //dynamic programming matrix outputed by the alignment algorithm
                assignments = global_alignment_assignments(neighborhoods[m].get_protein_vector(),
                                                           neighborhoods[n].get_protein_vector(), clusters, prot_stringency,
                                                           gap_score, mismatch);

                //Retrieves the alignment score
                score = assignments[neighborhoods[m].protein_count()][neighborhoods[n].protein_count()];

                if (score < neigh_stringency) continue; //ignore scores below stringency
                //Writes scores to output_file
                global_alignment_output_score(neighborhoods[m], neighborhoods[n], score, output_file);

                if (pairings_filename == "&") continue; //Dummy filename indicating this option was not chosen
                //Retrieves alignment and writes pairings to pairings_file
                global_alignment_output_pairings(neighborhoods[m], neighborhoods[n], clusters, assignments,
                                                 prot_stringency, pairings_file, gap_score, mismatch);
            }
        }
    }
}
