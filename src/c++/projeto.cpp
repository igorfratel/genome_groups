#include <vector>
#include <string.h>
#include <cstring>
#include "protein_grouping.h"
#include "genome_grouping.h"


int main(int argc, char *argv[]) {
	/*Main program, coordinates all the modes of execution calling the apropriate functions*/

	std::string neighborhoods_file; //Contains the genomic neighborhoods
	std::string genome_clustering_method;
	std::string prot_sim_filename;
	std::string formatted_prot_file; //File already formatted for the protein homology method
	std::string protein_homology_method;
	UndirectedEdgeWeightedGraph<std::string> prot_clusters;
	int num_prot; //Total number of unique proteins
	double stringency; //Baseline protein similarity score to consider two proteins as grouped

	switch (atoi(argv[1])) {

		case 0:
			//default execution
			//neighborhoods_file, prot_sim_filename, formatted_prot_file, protein_homology_method, num_prot, stringency,
			//genome_clustering_method

			neighborhoods_file = argv[2];
			prot_sim_filename = argv[3];
			formatted_prot_file = argv[4];
			protein_homology_method = argv[5];
			num_prot = atoi(argv[6]);
			stringency = atof(argv[7]);
			genome_clustering_method = argv[8];

			std::cout << "Applying homology detection method...\n";
			homology_detection(formatted_prot_file, protein_homology_method, prot_sim_filename);

			std::cout << "\nClustering proteins...\n";
			prot_clusters = protein_clustering(prot_sim_filename, num_prot);

			std::cout << "\nClustering genomic neighborhoods...\n";
			genome_clustering(neighborhoods_file, prot_clusters, genome_clustering_method, stringency);
			break;

		case 1:
			//Already has the similarities between the proteins. Needs to cluster them and the
	  		//genomic neighborhoods.
			//neighborhoods_file, prot_sim_filename, num_prot, stringency, genome_clustering_method

			neighborhoods_file = argv[2];
			prot_sim_filename = argv[3];
			num_prot = atoi(argv[4]);
			stringency = atof(argv[5]);
			genome_clustering_method = argv[6];

			std::cout << "\nClustering proteins...\n";
			prot_clusters = protein_clustering(prot_sim_filename, num_prot);

			std::cout << "\nClustering genomic neighborhoods...\n";
			genome_clustering(neighborhoods_file, prot_clusters, genome_clustering_method, stringency);
			break;


	}
}
