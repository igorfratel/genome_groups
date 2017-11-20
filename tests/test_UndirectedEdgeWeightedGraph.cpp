#include "UndirectedEdgeWeightedGraph.h"

int main() {
	UndirectedEdgeWeightedGraph my_graph;
	UndirectedEdgeWeightedGraph my_graph2 (5);
	std::vector<std::vector<std::string>> components;

	my_graph.add_node("primeiro");
	my_graph.add_node("segundo");
	my_graph.add_node("terceiro");
	my_graph.add_node("quarto");
	my_graph.add_node("quinto");
	my_graph.add_edge("primeiro", "terceiro", 0.05);
	my_graph.add_edge("quarto", "quinto", 1);
	std::cout << my_graph.is_connected("primeiro", "segundo") << "\n";
	std::cout << my_graph.is_connected("primeiro", "terceiro") << "\n";
	std::cout << my_graph.is_connected("quarto", "quinto") << " " << my_graph.get_edge_weight("quarto", "quinto") << "\n";

	std::cout << "\n";

	my_graph2.add_node("primeiro");
	my_graph2.add_node("segundo");
	my_graph2.add_node("terceiro");
	my_graph2.add_node("quarto");
	my_graph2.add_node("quinto");
	my_graph.add_edge("primeiro", "terceiro", 0.05);
	my_graph.add_edge("quarto", "quinto", 1);
	std::cout << my_graph.is_connected("primeiro", "segundo") << "\n";
	std::cout << my_graph.is_connected("primeiro", "terceiro") << " " << my_graph.get_edge_weight("primeiro", "terceiro") << "\n";
	std::cout << my_graph.is_connected("quarto", "quinto") << "\n";


	components = my_graph.connected_components(0.05);

	std::cout << "\n";

	for (std::vector<std::vector<std::string>>::iterator it = components.begin(); it != components.end(); ++it) {
		for (std::vector<std::string>::iterator it2 = it->begin(); it2 != it->end(); ++it2)
			std::cout << *it2 << " ";
		std::cout << "\n";
	}
			
}