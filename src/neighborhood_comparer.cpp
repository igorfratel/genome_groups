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
		("d,num_residues", "Total number of residues in unique protein set of formatted_prot_filename", cxxopts::value<std::string>()->default_value("537"))
		("t,prot_stringency", "Minimum similarity required to treat two proteins as a related pair", cxxopts::value<double>()->default_value("0.0"))
		("r,neigh_stringency", "Minimum threshold to display the similarity between two neighborhoods", cxxopts::value<double>()->default_value("0.0"))
		("m,mismatch", "Similarity value used when protein pair similarity is below prot_stringency", cxxopts::value<int>()->default_value("-1"))
		("g,neigh_comparing","Method for comparing genomic neighborhoods (default: porthodom method)", cxxopts::value<std::string>()->default_value("porthodom"))
		("G,gap_score", "Value added to similarity when a gap is chosen. Used when the neigh_comparing option is an alignment algorithm", cxxopts::value<int>()->default_value("0"))
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
				<<"    -d --num_residues\n"
			    	<<"    -t --prot_stringency\n"
				<<"    -r --neigh_stringency\n"
				<<"    -m --mismatch\n"
			    	<<"    -g --neigh_comparing\n"
				<<"    -G --gap_score\n"
			    	<<"    -o --output\n"
				<<"    -a --pairings_filename\n"
				<<"  partial --> Already has the similarities between the proteins.\n"
			    	<<"    -e --execution_mode partial\n"
			    	<<"    -n --neighborhoods_filename\n"
			    	<<"    -s --prot_sim_filename\n"
				<<"    -l --normalize_prot_sim\n"
			    	<<"    -t --prot_stringency\n"
				<<"    -r --neigh_stringency\n"
				<<"    -m --mismatch\n"
			    	<<"    -g --neigh_comparing\n"
				<<"    -G --gap_score\n"
			    	<<"    -o --output\n"
				<<"    -a --pairings_filename\n";
		return 0;
	}

	/*Actual program execution*/
	std::string execution_mode = result["execution_mode"].as<std::string>();
	std::string neighborhoods_filename = result["neighborhoods_filename"].as<std::string>();
	std::string prot_sim_filename = result["prot_sim_filename"].as<std::string>();
	double prot_stringency = result["prot_stringency"].as<double>();
	double neigh_stringency = result["neigh_stringency"].as<double>();

	std::string neigh_comparing = result["neigh_comparing"].as<std::string>();
	std::string output = result["output"].as<std::string>() ;
	std::string pairings_filename = result["pairings_filename"].as<std::string>();

	int gap_score = result["gap_score"].as<int>();
	int mismatch = result["mismatch"].as<int>();

	ProteinCollection prot_clusters;
	int num_prot;

	if (execution_mode == "full") {
		//default execution
		std::string formatted_prot_filename = result["formatted_prot_filename"].as<std::string>();
		std::string protein_comparing = result["protein_comparing"].as<std::string>();
		std::string num_residues = result["num_residues"].as<std::string>();

		std::cout << "Applying homology detection method...\n";
		homology_detection(formatted_prot_filename, protein_comparing, num_residues, prot_sim_filename);

		num_prot = total_protein_count(prot_sim_filename);

		std::cout << "\nClustering proteins...\n";
		prot_clusters = protein_clustering(prot_sim_filename, num_prot);

		std::vector<GenomicNeighborhood> neighborhoods = parse_neighborhoods(neighborhoods_filename);

		std::cout << "\nClustering genomic neighborhoods...\n";
		genome_clustering(neighborhoods, prot_clusters, neigh_comparing, prot_stringency, neigh_stringency, output,
			 	          pairings_filename, gap_score, mismatch);

		std::cout << "\nDone!";
	}

	else if (execution_mode == "partial") {
		//Already has the similarities between the proteins.
		int normalize_prot_sim = result.count("normalize_prot_sim");

		num_prot = total_protein_count(prot_sim_filename);

		std::cout << "\nClustering proteins...\n";
		prot_clusters = protein_clustering(prot_sim_filename, num_prot);

		if (normalize_prot_sim)
			prot_clusters.normalize();

		std::vector<GenomicNeighborhood> neighborhoods = parse_neighborhoods(neighborhoods_filename);

		std::cout << "\nClustering genomic neighborhoods...\n";
		genome_clustering(neighborhoods, prot_clusters, neigh_comparing, prot_stringency, neigh_stringency, output,
			              pairings_filename, gap_score, mismatch);

		std::cout << "\nDone!\n";

	}
}
