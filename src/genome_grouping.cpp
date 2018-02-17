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
                     neighborhoods.emplace_back(GenomicNeighborhood(split_line[accession_index], split_line[1] + " " +
                                                split_line[2]));
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
 *the format "organism1 acession1 coordinates1 organism2 acession2 coordinates2 score".
 */
void genome_clustering(const std::string &neighborhoods_filename, ProteinCollection &clusters,
                       const std::string &method, double stringency, const std::string &genome_sim_filename) {

    std::ofstream file(genome_sim_filename.c_str());
    std::cout.precision(2); //configure output limiting to 2 decimal points
    std::vector<GenomicNeighborhood> neighborhoods = parse_neighborhoods(neighborhoods_filename);

    //DEBUG print the genomic neighborhoods we are considering
    /*for (int i = 0; i < neighborhoods.size(); i++) {
        std::cout << neighborhoods[i].get_accession() << " " << neighborhoods[i].get_organism() <<
                     " " << neighborhoods[i].protein_count() << "\n";
        for (GenomicNeighborhood::iterator it = neighborhoods[i].begin(); it != neighborhoods[i].end(); ++it) {
            std::cout << it->pid << " " << it->locus << " " << it->cds << "\n";
        }
    }
    std::cout << "\n"*/;
    if (method == "porthodom") {
        std::mutex mtx;
        double tmp;
        unsigned int n;
        for(unsigned int m = 0; m < neighborhoods.size(); m++) {
            #pragma omp parallel for private(n, tmp)
            for (n = m + 1; n < neighborhoods.size(); n++) {
                tmp = porthodom_scoring(neighborhoods[m], neighborhoods[n], clusters, stringency);
                mtx.lock();
                file << neighborhoods[m].get_organism() << " " <<
                        neighborhoods[m].get_accession() << " " <<
                        neighborhoods[m].get_cds_string() << " " <<
                        neighborhoods[n].get_organism() << " " <<
                        neighborhoods[n].get_accession() << " " <<
                        neighborhoods[n].get_cds_string() << " " <<
                        tmp << "\n";
                mtx.unlock();
            }
        }
    }

    else if (method == "porthodomO2") {
        for(unsigned int m = 0; m < neighborhoods.size(); m++) {
            for (unsigned int n = m + 1; n < neighborhoods.size(); n++) {
                file << neighborhoods[m].get_organism() << " " <<
                        neighborhoods[m].get_accession() << " " <<
                        neighborhoods[m].get_cds_string() << " " <<
                        neighborhoods[n].get_organism() << " " <<
                        neighborhoods[n].get_accession() << " " <<
                        neighborhoods[n].get_cds_string() << " " <<
                        porthodomO2_scoring(neighborhoods[m], neighborhoods[n], clusters, stringency) << "\n";
            }
        }

    }
}
