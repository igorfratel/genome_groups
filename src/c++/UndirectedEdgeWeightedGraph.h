#ifndef __UNDIRECTED_EDGE_WEIGHTED_GRAPH_H__
#define __UNDIRECTED_EDGE_WEIGHTED_GRAPH_H__

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <stack>

class UndirectedEdgeWeightedGraph {

	typedef struct {
		int index;
		int visited;
	} node;

	std::map<std::string, node> nodes;
	std::vector<std::vector<double>> adj;

	public:

		UndirectedEdgeWeightedGraph();

		/*Creates object with known number of nodes to be added*/
		UndirectedEdgeWeightedGraph(int n_nodes);

		/*Adds node with string identifier*/
		void add_node(std::string n);

		/*Adds edge connecting two existing nodes with given weight*/
		void add_edge(std::string n1, std::string n2, double weight);

		/*Returns 1 if given nodes are directly connected and 0 otherwise*/
		int is_connected(std::string n1, std::string n2);

		/*Returns weight of the edge connecting two existing and directly connected nodes*/
		double get_edge_weight(std::string n1, std::string n2);

		/*Returns number of nodes in graph*/
		int get_number_nodes();

		/*Returns vector of connected components where each position is a vector of nodes in the same
		 *component*/
		std::vector<std::vector<std::string>> connected_components(double weight);

	private:

		/*connected_components auxiliar function*/
		void DFS_vector_fill(std::map<std::string, node>::iterator n,
												 std::vector<std::vector<std::string>> &components, double weight);
};

#endif
