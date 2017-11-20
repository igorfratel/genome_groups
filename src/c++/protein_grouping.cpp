#include "protein_grouping.h"

std::string homology_detection(char *format_file, char *method) {
    /*Receives the user's preferred protein homology/orthology detection method and runs it on a file
    /*already formatted to be its input.
    /*Writes the results to "homology_detection_output.txt" and returns this filename.*/
    std::string method_aux = method;

    if (method_aux == "nc") {
        std::string format_file_aux = format_file;
        std::string command = "NC_standalone -f " + format_file_aux + " -o homology_detection_output.txt";
        system(command.c_str()); //I'm using the default threshold and other defaut parameters(!!!)
    }

    return "homology_detection_output.txt";
}


std::vector<std::vector<std::string>> protein_clustering(char* prot_sim_file, int num_prot,
                                                         double stringency) {
    /*Receives the protein similarities file and stores them in a list of clusters.
    /*The similarities file must be in the format "prot1 prot2 sim" in every line
    /*num_prot must be the number of proteins and stringency is the minimum similarity
    /*for two proteins to be considered part of the same cluster*/
    UndirectedEdgeWeightedGraph my_graph (num_prot);
    std::ifstream file;
    std::string prot1;
    std::string prot2;
    std::string weight;
    std::vector<std::vector<std::string>> components;
    double weight_aux;

    file.open(prot_sim_file);
    while(std::getline(file, prot1, ' ')) {
        std::getline(file, prot2, ' ');
        std::getline(file, weight);
        my_graph.add_node(prot1);
        my_graph.add_node(prot2);
        weight_aux = ::atof(weight.c_str());
        my_graph.add_edge(prot1, prot2, weight_aux);
    }
    file.close();

    components = my_graph.connected_components(stringency);
    return components;
}
