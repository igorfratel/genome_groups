#include "../src/ProteinCollection.h"

int main() {
	ProteinCollection my_graph;
    for(int i = 0; i < 10000000; i++){
        my_graph.add_protein("protein"+std::to_string(i));
    }
    my_graph.connect_proteins("protein1", "protein500", 0.05);
	std::cout << my_graph.get_similarity("protein1", "protein500") << "\n";
}
