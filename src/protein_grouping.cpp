#include "protein_grouping.h"

/**
 *Receives the protein similarities file and returns the number of unique proteins in it.
 */
int total_protein_count(const std::string &prot_sim_filename) {
    std::set<std::string> protein_set;
    std::ifstream file;
    std::string prot;
    std::string test_string;

    file.open(prot_sim_filename.c_str());

    if (file.fail()) {
        std::cerr << "ERROR: trouble opening the protein similarities file\n";
        exit(1);
    }

    std::getline(file, test_string);
    if(test_string.find_first_of(' ') != std::string::npos) {
        std::vector<std::string> splitted = split(test_string, " ");
        protein_set.emplace(splitted[0]);
        protein_set.emplace(splitted[1]);
        while(std::getline(file, prot, ' ')) {
            protein_set.emplace(prot);
            std::getline(file, prot, ' ');
            protein_set.emplace(prot);
            std::getline(file, prot); //similarity is ignored
        }
    }

    else {
        std::vector<std::string> splitted = split(test_string, "\t");
        protein_set.emplace(splitted[0]);
        protein_set.emplace(splitted[1]);
        while(std::getline(file, prot, '\t')) {
            protein_set.emplace(prot);
            std::getline(file, prot, '\t');
            protein_set.emplace(prot);
            std::getline(file, prot); //similarity is ignored
        }
    }
    
    file.close();
    return protein_set.size();
}

/**
 *Receives the user's preferred protein homology/orthology detection method and runs it on a file
 *already formatted to be its input (that has the number of residues in its unique protein set equal to num_residues).
 *Writes the results to prot_sim_filename.
 */
void homology_detection(const std::string &format_file, const std::string &method, const std::string &num_residues,
                        const std::string &prot_sim_filename) {

    if (method == "nc") {
        std::string command = "NC_standalone -f " + format_file + " --num_residues " + num_residues + " -o " + prot_sim_filename;
        system(command.c_str()); //I'm using the default threshold and other defaut parameters(!!!)
    }
}

/**
 *Receives the protein similarities file and the number of unique proteins.
 *Stores the proteins and their similarity relationships in a ProteinCollection.
 *The similarities file must be in the format "prot1 prot2 sim" in every line (space separated);
 */
ProteinCollection protein_clustering(const std::string &prot_sim_filename, unsigned int num_prot) {

    ProteinCollection my_proteins (num_prot);
    std::ifstream file;
    std::string prot1;
    std::string prot2;
    std::string similarity;
    std::string test_string;

    file.open(prot_sim_filename.c_str());
    if (file.fail()) {
        std::cerr << "ERROR: trouble opening the protein similarities file\n";
        exit(1);
    }

    std::getline(file, test_string);
    if(test_string.find_first_of(' ') != std::string::npos) {
        std::vector<std::string> splitted = split(test_string, " ");
        prot1 = splitted[0];
        prot2 = splitted[1];
        similarity = splitted[2];
        my_proteins.add_protein(prot1);
        my_proteins.add_protein(prot2);
        my_proteins.connect_proteins(prot1, prot2, std::stod(similarity));
        while(std::getline(file, prot1, ' ')) {
            std::getline(file, prot2, ' ');
            std::getline(file, similarity);
            my_proteins.add_protein(prot1);
            my_proteins.add_protein(prot2);
            my_proteins.connect_proteins(prot1, prot2, std::stod(similarity));
        }
    }

    else {
        std::vector<std::string> splitted = split(test_string, "\t");
        prot1 = splitted[0];
        prot2 = splitted[1];
        similarity = splitted[2];
        my_proteins.add_protein(prot1);
        my_proteins.add_protein(prot2);
        my_proteins.connect_proteins(prot1, prot2, std::stod(similarity));
        while(std::getline(file, prot1, '\t')) {
            std::getline(file, prot2, '\t');
            std::getline(file, similarity);
            my_proteins.add_protein(prot1);
            my_proteins.add_protein(prot2);
            my_proteins.connect_proteins(prot1, prot2, std::stod(similarity));
        }
    }
    file.close();
    return my_proteins;
}
