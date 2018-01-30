#include "../src/ProteinCollection.h"
#include <iostream>
#include <cstdlib>

int main() {
	ProteinCollection my_graph(10000);
    for(int i = 0; i < 10000; i++){
        my_graph.add_protein("protein"+std::to_string(i));
    }

    int x = 0;
    int y = 0;
    for(int i = 0; i < 100000000; i++){
        x = rand()%100000000;
        y = rand()%100000000;
        my_graph.connect_proteins("protein"+std::to_string(x), "protein"+std::to_string(y), 0.05);
        my_graph.get_similarity("protein"+std::to_string(x), "protein"+std::to_string(y));
    }
}
