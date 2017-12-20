#ifndef __PROTEIN_COLLECTION_H__
#define __PROTEIN_COLLECTION_H__

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <stack>

/*Undirected edge-weighted graph implementation*/
class ProteinCollection {

	typedef struct {
		int index;
		int visited;
	} node_info_t;

	std::map<std::string, node_info_t> nodes; //nodes are proteins
	std::vector<std::vector<double> > adj; //Adjacency matrix

	public:

		ProteinCollection();

		/*Creates object with known number of proteins to be added*/
		ProteinCollection(int n_nodes);

		/*Adds protein with string identifier*/
		void add_protein(std::string node);

		/*Adds connection between two existing proteins with given similarity*/
		void connect_proteins(std::string node1, std::string node2, double weight);

		/*Adds two proteins and connects them with given similarity*/
		void add_connected_proteins(std::string node1, std::string node2, double weight);

		/*Returns true if given proteins are directly connected and false otherwise*/
		bool are_connected(std::string node1, std::string node2);

		/*Returns similarity value between two existing and directly connected proteins.
		 *If not connected, returns 0.0*/
		double get_similarity(std::string node1, std::string node2);

		/*Returns number of proteins in object*/
		int get_number_proteins();

		/*Returns vector of connected components where each position is a vector of nodes in the same
		 *component*/
		std::vector<std::vector<std::string> > connected_components(double weight);

	private:

		/*connected_components auxiliar function*/
		void DFS_vector_fill(typename std::map<std::string, node_info_t>::iterator n,
							 std::vector<std::vector<std::string> > &components, double weight);
};

#endif
