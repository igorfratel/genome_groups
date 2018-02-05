#include "../src/ProteinCollection.h"

int main() {
	ProteinCollection my_graph;
	ProteinCollection my_graph2 (5);
	std::vector<std::vector<std::string>> components;

	my_graph.add_protein("primeiro");
	my_graph.add_protein("segundo");
	my_graph.add_protein("terceiro");
	my_graph.add_protein("quarto");
	my_graph.add_protein("quinto");
	my_graph.connect_proteins("primeiro", "terceiro", 0.05);
	my_graph.connect_proteins("quarto", "quinto", 1);
	std::cout << my_graph.are_connected("primeiro", "segundo") << "\n";
	std::cout << my_graph.are_connected("primeiro", "terceiro") << "\n";
	std::cout << my_graph.are_connected("quarto", "quinto") << " " << my_graph.get_similarity("quarto", "quinto") << "\n";

	std::cout << "\n";

	my_graph2.add_protein("primeiro");
	my_graph2.add_protein("segundo");
	my_graph2.add_protein("terceiro");
	my_graph2.add_protein("quarto");
	my_graph2.add_protein("quinto");
	my_graph2.connect_proteins("primeiro", "terceiro", 0.05);
	my_graph2.connect_proteins("quarto", "quinto", 1);
	std::cout << my_graph2.are_connected("primeiro", "segundo") << "\n";
	std::cout << my_graph2.are_connected("primeiro", "terceiro") << " " << my_graph2.get_similarity("primeiro", "terceiro") << "\n";
	std::cout << my_graph2.are_connected("quarto", "quinto") << "\n";


	components = my_graph2.connected_components(0.05);

	std::cout << "\n";

	for (std::vector<std::vector<std::string>>::iterator it = components.begin(); it != components.end(); ++it) {
		for (std::vector<std::string>::iterator it2 = it->begin(); it2 != it->end(); ++it2)
			std::cout << *it2 << " ";
		std::cout << "\n";
	}

}
