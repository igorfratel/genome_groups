#include "UndirectedEdgeWeightedGraph.h"

UndirectedEdgeWeightedGraph::UndirectedEdgeWeightedGraph() {}

UndirectedEdgeWeightedGraph::UndirectedEdgeWeightedGraph(int n_nodes) {
	/*Creates object with known number of nodes to be added*/
	adj.resize(n_nodes);
	for(unsigned int i = 0; i < adj.size(); i++)
		adj[i].resize(n_nodes);
}

void UndirectedEdgeWeightedGraph::add_node(std::string n) {
	/*Adds node with string identifier*/
	node my_node;
	my_node.index = nodes.size();
	my_node.visited = 0;
	nodes.insert(std::make_pair(n, my_node));
	adj.resize(adj.size() + 1);
	adj[adj.size() - 1].resize(adj.size());
}

void UndirectedEdgeWeightedGraph::add_edge(std::string n1, std::string n2, double weight) {
	/*Adds edge connecting two existing nodes with given weight*/
	int x = (nodes.at(n1)).index;
	int y = (nodes.at(n2)).index;
	adj[x][y] = weight;
	adj[y][x] = weight;
}

int UndirectedEdgeWeightedGraph::is_connected(std::string n1, std::string n2) {
	/*Returns 1 if given nodes are directly connected and 0 otherwise*/
	int x = (nodes.at(n1)).index;
	int y = (nodes.at(n2)).index;
	if (adj[x][y] > 0) return 1;
	return 0;
}

double UndirectedEdgeWeightedGraph::get_edge_weight(std::string n1, std::string n2) {
	/*Returns weight of the edge connecting two existing and directly connected nodes*/
	int x = (nodes.at(n1)).index;
	int y = (nodes.at(n2)).index;
	return adj[x][y];
}

int UndirectedEdgeWeightedGraph::get_number_nodes() {
	return nodes.size();
}

void UndirectedEdgeWeightedGraph::DFS_vector_fill(std::map<std::string, node>::iterator n,
																									std::vector<std::vector<std::string>> &components,
																								  double weight) {
	/*Runs DFS from a starting node "n" and adds a vector to components with all the nodes in the same
	/*connected component as n*/
	std::stack<std::map<std::string, node>::iterator> my_stack;
	std::vector<std::string> aux;
	my_stack.push(n);
	while (!my_stack.empty()) {
		n = my_stack.top();
		my_stack.pop();
		if (!(n->second).visited) {
			aux.push_back(n->first);
			(n->second).visited = 1;
		}
		for(std::map<std::string, node>::iterator it = nodes.begin(); it != nodes.end(); ++it)
			if (adj[(n->second).index][(it->second).index] >= weight && !(it->second).visited)
				my_stack.push(it);
	}
	components.push_back(aux);
}

std::vector<std::vector<std::string>> UndirectedEdgeWeightedGraph::connected_components(double weight) {
	/*Returns vector of connected components where each position is a vector of nodes in the
	/*same component*/
	std::vector<std::vector<std::string>> components;
	for(std::map<std::string, node>::iterator it = nodes.begin(); it != nodes.end(); ++it)
		(it->second).visited = 0;

	for(std::map<std::string, node>::iterator it = nodes.begin(); it != nodes.end(); ++it)
		if (!(it->second).visited)
			//For every unvisited node, runs DFS to get all connected nodes
			DFS_vector_fill(it, components, weight);
	return components;
}
