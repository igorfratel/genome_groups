#include "../src/ProteinCollection.h"
#include <iostream>

int main() {
	ProteinCollection my_graph;
    for(int i = 0; i < 1000000; i++){
        my_graph.add_protein("protein"+std::to_string(i));
    }
    my_graph.connect_proteins("protein1", "protein500");
	std::cout << my_graph.get_similarity("protein1", "protein500") << "\n";
}
