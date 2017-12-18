#ifndef __UNDIRECTED_EDGE_WEIGHTED_GRAPH_H__
#define __UNDIRECTED_EDGE_WEIGHTED_GRAPH_H__

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <stack>

template <typename T> //Value represented by each vertex

class UndirectedEdgeWeightedGraph {

	typedef struct {
		int index;
		int visited;
	} node_info_t;

	std::map<T, node_info_t> nodes;
	std::vector<std::vector<double> > adj; //Adjacency matrix

	public:

		UndirectedEdgeWeightedGraph();

		/*Creates object with known number of nodes to be added*/
		UndirectedEdgeWeightedGraph(int n_nodes);

		/*Adds node with string identifier*/
		void add_node(T node);

		/*Adds edge connecting two existing nodes with given weight*/
		void add_edge(T node1, T node2, double weight);

		/*Adds two nodes and connects them with given weight*/
		void add_connected_nodes(T node1, T node2, double weight);

		/*Returns true if given nodes are directly connected and false otherwise*/
		bool are_connected(T node1, T node2);

		/*Returns weight of the edge connecting two existing and directly connected nodes.
		 *If not connected, returns 0.0*/
		double get_edge_weight(T node1, T node2);

		/*Returns number of nodes in graph*/
		int get_number_nodes();

		/*Returns vector of connected components where each position is a vector of nodes in the same
		 *component*/
		std::vector<std::vector<T> > connected_components(double weight);

	private:

		/*connected_components auxiliar function*/
		void DFS_vector_fill(typename std::map<T, node_info_t>::iterator n,
							 std::vector<std::vector<T> > &components, double weight);
};

#endif
