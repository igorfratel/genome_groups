#include "genome_grouping.h"

/**
 *Receives a string and delimiters.
 *Splits the string in a vector according to delimiters
 */
static std::vector<std::string> split(const std::string& in, const std::string& delim) {
   std::string::size_type start = in.find_first_not_of(delim), end = 0;

   std::vector<std::string> out;
   while(start != in.npos) {
      end = in.find_first_of(delim, start);
      if(end == in.npos) {
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
 *Receives a genomic neighborhood filename.
 *Returns a vector of genomic neighborhoods, filled with the information from the file.
 */
 static std::vector<GenomicNeighborhood> parse_neighborhoods(const std::string &neighborhoods_filename) {
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
 *Receives a file containing all the genomic neighborhoods,
 *a ProteinCollection and the desired genomic neighborhood clustering method.
 *Writes the similarity between all genomic neighborhoods on the genome_sim_filename in
 *the format "acession1 coordinates1 acession2 coordinates2 score".
 */
void genome_clustering(const std::string &neighborhoods_filename, ProteinCollection &clusters,
                       const std::string &method, double prot_stringency, double neigh_stringency, const std::string &genome_sim_filename,
                       const std::string &pairings_filename) {

    std::ofstream output_file;
    if(genome_sim_filename == "-")
        output_file.basic_ios<char>::rdbuf(std::cout.rdbuf());
    else
        output_file = std::ofstream(genome_sim_filename.c_str());

    std::ofstream pairings_file;
    if(pairings_filename != "&") //Dummy filename indicating this option was not chosen
        pairings_file = std::ofstream(pairings_filename.c_str());

    std::vector<GenomicNeighborhood> neighborhoods = parse_neighborhoods(neighborhoods_filename);

    if (method == "porthodom") {
        std::map<std::pair<int,int>, int> assignments;
        double sum_temp;
        double score;
        for(unsigned int m = 0; m < neighborhoods.size(); m++) {
            for (unsigned int n = m + 1; n < neighborhoods.size(); n++) {
                sum_temp = 0;
                //Edges chosen by the algorithm
                assignments = porthodom_scoring(neighborhoods[m], neighborhoods[n], clusters, prot_stringency);

                for (std::map<std::pair<int, int>,int>::iterator it = assignments.begin(); it != assignments.end(); ++it)
                    sum_temp += ((double)it->second)/1000000; //Division to undo the multiplication in clustering_value()

                //apply the scoring formula
                score = sum_temp/std::max(neighborhoods[m].protein_count(), neighborhoods[n].protein_count());

                if (score < neigh_stringency) continue; //ignore scores below stringency

                //Writes scores to output_file
                output_file << neighborhoods[m].get_accession() << "\t" <<
                        neighborhoods[m].get_first_cds() << "\t" <<
                        neighborhoods[m].get_last_cds() << "\t" <<
                        neighborhoods[n].get_accession() << "\t" <<
                        neighborhoods[n].get_first_cds() << "\t" <<
                        neighborhoods[n].get_last_cds() << "\t" <<
                        score << "\n";

                if (pairings_filename == "&") continue; //Dummy filename indicating this option was not chosen
                //Writes pairing header to pairings_file
                pairings_file << ">" << neighborhoods[m].get_accession() << "\t" <<
                        neighborhoods[m].get_first_cds() << "\t" <<
                        neighborhoods[m].get_last_cds() << "\t" <<
                        neighborhoods[n].get_accession() << "\t" <<
                        neighborhoods[n].get_first_cds() << "\t" <<
                        neighborhoods[n].get_last_cds() << "\n";

                //Writes pairings to pairings_file
                for (std::map<std::pair<int, int>,int>::iterator it = assignments.begin(); it != assignments.end(); ++it){
                    pairings_file << neighborhoods[m].get_pid(it->first.first) << "\t" <<
                                     neighborhoods[n].get_pid(it->first.second) << "\t" <<
                                     ((double)it->second)/1000000 << "\n";
                }



            }
        }
    }

    else if (method == "porthodomO2") {
        for(unsigned int m = 0; m < neighborhoods.size(); m++) {
            for (unsigned int n = m + 1; n < neighborhoods.size(); n++) {
                output_file << neighborhoods[m].get_accession() << "\t" <<
                        neighborhoods[m].get_first_cds() << "\t" <<
                        neighborhoods[m].get_last_cds() << "\t" <<
                        neighborhoods[n].get_accession() << "\t" <<
                        neighborhoods[n].get_first_cds() << "\t" <<
                        neighborhoods[n].get_last_cds() << "\t" <<
                        porthodomO2_scoring(neighborhoods[m], neighborhoods[n], clusters, prot_stringency) << "\n";
            }
        }

    }
}
