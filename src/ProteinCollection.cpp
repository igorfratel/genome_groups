/*Undirected edge-weighted graph implementation*/

#include "ProteinCollection.h"

ProteinCollection::ProteinCollection() {}

ProteinCollection::ProteinCollection(int n_nodes) {
	/*Creates object with known number of nodes to be added*/

	adj.resize(n_nodes);
	for(unsigned int i = 0; i < adj.size(); i++)
		adj[i].resize(n_nodes);
}

void ProteinCollection::add_protein(std::string node) {
	/*Adds node with string identifier*/

	node_info_t my_node_info;
	my_node_info.index = nodes.size();
	my_node_info.visited = 0;
	nodes.insert(std::make_pair(node, my_node_info));

	if(adj.size() < nodes.size()){
		adj.resize(adj.size() + 1);
		adj[adj.size() - 1].resize(adj.size());
	}
}

void ProteinCollection::connect_proteins(std::string node1, std::string node2, double weight) {
	/*Adds edge connecting two existing nodes with given weight*/

	int x = (nodes.at(node1)).index;
	int y = (nodes.at(node2)).index;
	adj[x][y] = weight;
	adj[y][x] = weight;
}

void ProteinCollection::add_connected_proteins(std::string node1, std::string node2, double weight) {
	/*Adds two nodes and connects them with given weight*/
	node_info_t my_node_info1, my_node_info2;

	my_node_info1.index = nodes.size();
	my_node_info1.visited = 0;
	nodes.insert(std::make_pair(node1, my_node_info1));

	my_node_info2.index = nodes.size();
	my_node_info2.visited = 0;
	nodes.insert(std::make_pair(node2, my_node_info2));

	if(adj.size() < nodes.size()){
		adj.resize(adj.size() + 2);
		adj[adj.size() - 1].resize(adj.size());
		adj[adj.size() - 2].resize(adj.size());
	}

	adj[my_node_info1.index][my_node_info2.index] = weight;
	adj[my_node_info2.index][my_node_info1.index] = weight;
}

bool ProteinCollection::are_connected(std::string node1, std::string node2) {
	/*Returns 1 if given nodes are directly connected and 0 otherwise*/
	node_info_t aux_x;
	node_info_t aux_y;
	try { //In case some node does not exist
		aux_x = (nodes.at(node1));
		aux_y = (nodes.at(node2));
	}
	catch (...) {
		return false;
	}
	int x = aux_x.index;
	int y = aux_y.index;
	if (adj[x][y] > 0) return true;
	return false;
}

double ProteinCollection::get_similarity(std::string node1, std::string node2) {
	/*Returns weight of the edge connecting two existing and directly connected nodes.
	 *If not connected, returns 0.0*/
	if (!are_connected(node1, node2)) return 0.0;
	int x = (nodes.at(node1)).index;
	int y = (nodes.at(node2)).index;
	return adj[x][y];
}

int ProteinCollection::get_number_proteins() {
	return nodes.size();
}

std::vector<std::vector<std::string> > ProteinCollection::connected_components(double weight) {
	/*Returns vector of connected components where each position is a vector of nodes in the
	/*same component*/

	std::vector<std::vector<std::string> > components;
	for(typename std::map<std::string, node_info_t>::iterator it = nodes.begin(); it != nodes.end(); ++it)
		(it->second).visited = 0;

	for(typename std::map<std::string, node_info_t>::iterator it = nodes.begin(); it != nodes.end(); ++it)
		if (!(it->second).visited)
			//For every unvisited node, runs DFS to get all connected nodes
			DFS_vector_fill(it, components, weight);
	return components;
}

void ProteinCollection::DFS_vector_fill(typename std::map<std::string, node_info_t>::iterator n,
											   	  std::vector<std::vector<std::string> > &components,
												  double weight) {
	/*Runs DFS from a starting node "n" and adds a vector to components with all the nodes in the same
	/*connected component as n*/

	std::stack<typename std::map<std::string, node_info_t>::iterator> my_stack;
	std::vector<std::string> aux;
	my_stack.push(n);
	while (!my_stack.empty()) {
		n = my_stack.top();
		my_stack.pop();
		if (!(n->second).visited) {
			aux.push_back(n->first);
			(n->second).visited = 1;
		}
		for(typename std::map<std::string, node_info_t>::iterator it = nodes.begin(); it != nodes.end(); ++it)
			if (adj[(n->second).index][(it->second).index] >= weight && !(it->second).visited)
				my_stack.push(it);
	}
	components.push_back(aux);
}
