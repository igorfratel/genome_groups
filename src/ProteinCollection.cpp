/*Undirected edge-weighted graph implementation*/

#include "ProteinCollection.h"
#include <iostream>

/**
 * Creates object with known number of nodes to be added
 */
ProteinCollection::ProteinCollection(size_t n_nodes)
	: nodes(n_nodes), adj(n_nodes) {
}

/**
 * Adds node with string identifier
 */
void ProteinCollection::add_protein(const std::string& node) {
	auto new_id = nodes.size();
	nodes.emplace(node, node_info_t { new_id, false });
	adj.emplace_back();
}

/**
 * Adds edge connecting two existing nodes with given weight
 */
void ProteinCollection::connect_proteins(const std::string& node1, const std::string& node2, double weight) {
    auto x = nodes.at(node1).index;
 	auto y = nodes.at(node2).index;

	//DEBUG
 	/*std::cerr << x << std::endl;
 	std::cerr << y << std::endl;
 	std::cerr << adj.size() << std::endl;*/
	if (x >= y)
  		adj[x].emplace(y, weight);
	else
  		adj[y].emplace(x, weight);
}

/**
 * Adds two nodes and connects them with given weight
 */
void ProteinCollection::add_connected_proteins(const std::string& node1, const std::string& node2, double weight) {
	add_protein(node1);
	add_protein(node2);
	connect_proteins(node1, node2, weight);
}

/**
 * @returns True if given nodes are directly connected and false otherwise.
 */
bool ProteinCollection::are_connected(const std::string& node1, const std::string& node2) {

	try{
    	auto x = nodes.at(node1).index;
  		auto y = nodes.at(node2).index;
		return x >= y ? adj[x].count(y) : adj[y].count(x);
	}
	catch(...){
		return false;
	}
}

/**
 * @returns Weight of the edge connecting two existing and directly
 * connected nodes or 0.0 otherwise.
 */
double ProteinCollection::get_similarity(const std::string& node1,
	 									 const std::string& node2) {

	try{
		auto x = nodes.at(node1).index;
		auto y = nodes.at(node2).index;
		try{
			return x >= y ? adj[x].at(y) : adj[y].at(x);
		}
		catch(...){
			return 0.0;
		}
	}
	catch(...){
		return 0.0;
	}
}

/**
* @returns Vector of connected components where each position is a vector
* of nodes in the same component.
*/
std::vector<std::vector<std::string>> ProteinCollection::connected_components(double weight) {

	std::vector<std::vector<std::string> > components;
	for(typename std::unordered_map<std::string, node_info_t>::iterator it = nodes.begin(); it != nodes.end(); ++it)
		(it->second).visited = false;

 	for(typename std::unordered_map<std::string, node_info_t>::iterator it = nodes.begin(); it != nodes.end(); ++it)
		if (!(it->second).visited)
	  	//For every unvisited node, runs DFS to get all connected nodes
	  	DFS_vector_fill(it, components, weight);
  	return components;
}

/**
* Runs DFS from a starting node "n" and adds a vector to components with all the nodes in the same
* connected component as n.
*/
void ProteinCollection::DFS_vector_fill(typename std::unordered_map<std::string, node_info_t>::iterator n,
										std::vector<std::vector<std::string> > &components,
	  									double weight) {

	std::stack<typename std::unordered_map<std::string, node_info_t>::iterator> my_stack;
	std::vector<std::string> aux;
	my_stack.push(n);
	while (!my_stack.empty()) {
		n = my_stack.top();
		my_stack.pop();
		if (!(n->second).visited) {
			aux.push_back(n->first);
			(n->second).visited = true;
		}
		for(typename std::unordered_map<std::string, node_info_t>::iterator it = nodes.begin(); it != nodes.end(); ++it){
			if (!(it->second).visited) {

				try {
					if (adj[(n->second).index].at((it->second).index) >= weight)
					my_stack.push(it);
				}

				catch(...){}
		/*for(size_t i = 0; i < adj[(n->second).index].size(); i++){
			if(adj[(n->second).index][i].first == (it->second).index && adj[(n->second).index][i].second >= weight) {
			my_stack.push(it);
			break;
			}*/
			}
		}
	}
	components.push_back(aux);
}
