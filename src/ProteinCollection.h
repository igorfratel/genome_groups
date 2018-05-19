#ifndef __PROTEIN_COLLECTION_H__
#define __PROTEIN_COLLECTION_H__

#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <stack>

/*Undirected edge-weighted graph with adjacency lists implementation*/
class ProteinCollection {

	struct node_info_t {
		size_t index;
		bool visited;
	};

	std::unordered_map<std::string, node_info_t> nodes; //nodes are proteins
	std::vector<std::unordered_map<int, double>> adj; //Adjacency list

	public:

		/**
		 *Does not allocate anything.
		 *Used only to declare a ProteinCollection that will be assigned to a properly constructed one.
		 */
		ProteinCollection() = default;

		/**
	 	 *Creates object with known number of proteins to be added.
		 */
		ProteinCollection(size_t n_nodes);

		/**
		 *Adds protein with string identifier.
		 */
		void add_protein(const std::string& node);

		/**
		 *Adds connection between two existing proteins with given similarity.
		 */
		void connect_proteins(const std::string& node1, const std::string& node2, double weight);

		/*
		 *Returns true if given proteins are directly connected and false otherwise.
		 */
		bool are_connected(const std::string& node1, const std::string& node2) const;

		/**
		 *Returns similarity value between two existing and directly connected proteins.
		 *If not connected, returns 0.0.
		 */
		double get_similarity(const std::string& node1, const std::string& node2) const;

		/**
		 *Normalizes all the similarity scores, that is, divides them by the highest score in the object.
		 */
		void normalize();

		/**
		 *Returns vector of connected components where each position is a vector of nodes in the same
		 *component.
		 */
		std::vector<std::vector<std::string> > connected_components(double weight);

	private:

		/**
	 	 *Runs DFS from a starting node "n" and adds a vector to components with all the nodes in the same
		 *connected component as n.
		 */
		void DFS_vector_fill(typename std::unordered_map<std::string, node_info_t>::iterator n,
							 std::vector<std::vector<std::string> > &components, double weight);
};

#endif
