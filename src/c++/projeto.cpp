#include <vector>
#include <string.h>
#include <cstring>
#include "protein_grouping.h"
#include "genome_grouping.h"


int main(int argc, char *argv[]) {
	/*Main program, coordinates all the modes of execution calling the apropriate functions*/

	char *neighborhoods_file;
	char *genome_clustering_method;
	char *prot_sim_file;
	std::string prot_sim_file_aux;
	std::vector<std::vector<std::string>> prot_clusters;
	int num_prot;
	double stringency;

	switch (atoi(argv[1])) {
		case 0:
		//Already has clustered proteins. Just needs to cluster the genomic neighborhoods
	  //neighborhoods_file, prot_clusters_file, genome_clustering_method

			break;

		case 1:
		//Already has the similarities between the proteins. Needs to cluster them and the
	  //genomic neighborhoods.
		//neighborhoods_file, prot_sim_file, num_prot, stringency, genome_clustering_method
			neighborhoods_file = argv[2];
			prot_sim_file = argv[3];
			num_prot = atoi(argv[4]);
			stringency = atof(argv[5]);
			genome_clustering_method = argv[6];

			std::cout << "Clustering proteins...\n";
			prot_clusters = protein_clustering(prot_sim_file, num_prot, stringency);

			std::cout << "Protein clusters: \n";
			for (unsigned int i = 0; i < prot_clusters.size(); i++) {
				std::cout << "[";
				for (unsigned int j = 0; j < prot_clusters[i].size(); j++) {
					if (j == prot_clusters[i].size() - 1) std::cout << prot_clusters[i][j];
					else std::cout << prot_clusters[i][j] << ", ";
				}
				std::cout << "]\n";
			}

			std::cout << "Clustering genomic neighborhoods...\n";
			genome_clustering(neighborhoods_file, prot_clusters, genome_clustering_method);
			break;

		case 2:
		  //default execution
			//neighborhoods_file, formatted_prot_file, protein_homology_method, num_prot, stringency,
			//genome_clustering_method
			neighborhoods_file = argv[2];
			char *formatted_prot_file = argv[3];
			char *protein_homology_method = argv[4];
			num_prot = atoi(argv[5]);
			stringency = atof(argv[6]);
			genome_clustering_method = argv[7];

			std::cout << "Applying homology detection method...\n";
			prot_sim_file_aux = homology_detection(formatted_prot_file, protein_homology_method);
			prot_sim_file = new char[prot_sim_file_aux.length() + 1];
			std::strcpy(prot_sim_file, prot_sim_file_aux.c_str());

			std::cout << "Clustering proteins...\n";
			prot_clusters = protein_clustering(prot_sim_file, num_prot, stringency);
			std::cout << "Protein clusters: \n";
			for (unsigned int i = 0; i < prot_clusters.size(); i++) {
				std::cout << "[";
				for (unsigned int j = 0; j < prot_clusters[i].size(); j++) {
					if (j == prot_clusters[i].size() - 1) std::cout << prot_clusters[i][j];
					else std::cout << prot_clusters[i][j] << ", ";
				}
				std::cout << "]\n";
			}

			std::cout << "Clustering genomic neighborhoods...\n";
			genome_clustering(neighborhoods_file, prot_clusters, genome_clustering_method);
			break;
	}
}
