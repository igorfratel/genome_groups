#include <vector>
#include <string.h>
#include <cstring>
#include "protein_grouping.h"
#include "genome_grouping.h"
#include "ProteinCollection.h"
#include "cxxopts.hpp"


int main(int argc, char *argv[]) {
	/*Main program, coordinates all the modes of execution calling the apropriate functions*/

	/*Dealing with command line arguments*/
	cxxopts::Options options("neighborhood_comparer", "genomic neighborhood comparison software");
	options.add_options()
		("h, help", "Prints help message")
		("e,execution_mode", "full or partial execution mode (default: full)", cxxopts::value<std::string>()->default_value("full"))
		("n,neighborhoods_filename", "File containing the genomic neighborhoods", cxxopts::value<std::string>())
		("s,prot_sim_filename", "File containing pairs of proteins and their similarities", cxxopts::value<std::string>())
		("l,normalize_prot_sim", "Indicates that the protein similarities file should be normalized (Used in the partial execution mode)")
		("f,formatted_prot_filename", "File already formatted as the input for the homology detection method", cxxopts::value<std::string>())
		("p,protein_comparing", "Method for comparing proteins (default: nc)", cxxopts::value<std::string>()->default_value("nc"))
		("m,num_prot", "Number of unique proteins", cxxopts::value<int>())
		("t,prot_stringency", "Minimum similarity required to treat two proteins as a related pair", cxxopts::value<double>()->default_value("0.0"))
		("r,neigh_stringency", "Minimum threshold to display the similarity between two neighborhoods", cxxopts::value<double>()->default_value("0.0"))
		("g,neigh_comparing","Method for comparing genomic neighborhoods (default: porthodom method)", cxxopts::value<std::string>()->default_value("porthodom"))
		("o,output", "Where the neighborhood similarities should be written", cxxopts::value<std::string>()->default_value("-"))
		("a,pairings_filename", "Where the chosen pairings between proteins in the neighborhoods should be written", cxxopts::value<std::string>()->default_value("&"))
		;

	auto result = options.parse(argc, argv);

	if (result.count("help")){
		/*Shows help message*/
		std::cout << "HOW TO RUN:\n"
				<< "  - 1: ./make\n"
				<<"  - 2: ./python parse_neighborhood.py <file with neighborhoods>\n"
				<<"  - 3: ./neighborhood_comparer <full or partial> <args according to chosen mode>\n\n"

				<<"ARGUMENTS FOR EACH EXECUTION MODE:\n"
				<<"  full --> default execution\n"
		    	<<"    -e --execution_mode full\n"
			    <<"    -n --neighborhoods_filename\n"
			    <<"    -s --prot_sim_filename\n"
			    <<"    -f --formatted_prot_filename\n"
			    <<"    -p --protein_comparing\n"
			    <<"    -m --num_prot\n"
			    <<"    -t --prot_stringency\n"
				<<"    -r --neigh_stringency\n"
			    <<"    -g --neigh_comparing\n"
			    <<"    -o --output\n"
				<<"    -a --pairings_filename\n"

				<<"  partial --> Already has the similarities between the proteins.\n"
			    <<"    -e --execution_mode partial\n"
			    <<"    -n --neighborhoods_filename\n"
			    <<"    -s --prot_sim_filename\n"
			    <<"    -m --num_prot\n"
			    <<"    -t --prot_stringency\n"
				<<"    -r --neigh_stringency\n"
			    <<"    -g --neigh_comparing\n"
			    <<"    -o --output\n"
				<<"-a --pairings_filename\n";

		return 0;
	}

	/*Actual program execution*/
	std::string execution_mode = result["execution_mode"].as<std::string>();
	std::string neighborhoods_filename = result["neighborhoods_filename"].as<std::string>();
	std::string prot_sim_filename = result["prot_sim_filename"].as<std::string>();
	int num_prot = result["num_prot"].as<int>();
	double prot_stringency = result["prot_stringency"].as<double>();
	double neigh_stringency = result["neigh_stringency"].as<double>();

	std::string neigh_comparing = result["neigh_comparing"].as<std::string>();
	std::string output = result["output"].as<std::string>() ;
	std::string pairings_filename = result["pairings_filename"].as<std::string>();

	ProteinCollection prot_clusters;

	if (execution_mode == "full") {
		//default execution
		std::string formatted_prot_filename = result["formatted_prot_filename"].as<std::string>();
		std::string protein_comparing = result["protein_comparing"].as<std::string>();

		std::cout << "Applying homology detection method...\n";
		homology_detection(formatted_prot_filename, protein_comparing, prot_sim_filename);

		std::cout << "\nClustering proteins...\n";
		prot_clusters = protein_clustering(prot_sim_filename, num_prot);

		std::cout << "\nClustering genomic neighborhoods...\n";
		genome_clustering(neighborhoods_filename, prot_clusters, neigh_comparing, prot_stringency, neigh_stringency, output, pairings_filename);

		std::cout << "\nDone!";
	}

	else if (execution_mode == "partial") {
		//Already has the similarities between the proteins.
		int normalize_prot_sim = result.count("normalize_prot_sim");

		std::cout << "\nClustering proteins...\n";
		prot_clusters = protein_clustering(prot_sim_filename, num_prot);

		if (normalize_prot_sim)
			prot_clusters.normalize();

		std::cout << "\nClustering genomic neighborhoods...\n";
		genome_clustering(neighborhoods_filename, prot_clusters, neigh_comparing, prot_stringency, neigh_stringency, output, pairings_filename);

		std::cout << "\nDone!\n";

	}
}
