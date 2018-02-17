#include "protein_grouping.h"

/**
 *Receives the user's preferred protein homology/orthology detection method and runs it on a file
 *already formatted to be its input.
 *Writes the results to prot_sim_filename
 */
void homology_detection(const std::string &format_file, const std::string &method, const std::string &prot_sim_filename) {

    if (method == "nc") {
        std::string command = "NC_standalone -f " + format_file + " -o " + prot_sim_filename;
        system(command.c_str()); //I'm using the default threshold and other defaut parameters(!!!)
    }
}


/**
* Receives the protein similarities file and stores them in a ProteinCollection.
* The similarities file must be in the format "prot1 prot2 sim" in every line
* num_prot must be the number of proteins and stringency is the minimum similarity
* for two proteins to be considered part of the same cluster
*/
ProteinCollection protein_clustering(const std::string &prot_sim_filename, unsigned int num_prot) {

    ProteinCollection my_proteins (num_prot);
    std::ifstream file;
    std::string prot1;
    std::string prot2;
    std::string similarity;

    file.open(prot_sim_filename.c_str());
    if (file.fail()) {
        std::cerr << "ERROR: trouble opening the protein similarities file\n";
        exit(1);
    }

    while(std::getline(file, prot1, ' ')) {
        std::getline(file, prot2, ' ');
        std::getline(file, similarity);
        my_proteins.add_protein(prot1);
        my_proteins.add_protein(prot2);

        //DEBUG
        //std::cout << prot1 << " " << prot2 << " similarity: " << similarity_aux << "\n";
        my_proteins.connect_proteins(prot1, prot2, std::atof(similarity.c_str()));
    }
    file.close();
    return my_proteins;
}
