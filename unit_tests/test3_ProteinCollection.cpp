#include "../src/ProteinCollection.h"
#include <iostream>
#include <cstdlib>

int main() {
	ProteinCollection my_graph;
    for(int i = 0; i < 10000; i++){
        my_graph.add_protein("protein"+std::to_string(i));
    }

    int x = 0;
    int y = 0;
    for(int i = 0; i < 177463380; i++){
        x = rand()%10000;
        y = rand()%10000;
        my_graph.connect_proteins("protein"+std::to_string(x), "protein"+std::to_string(y));
        my_graph.get_similarity("protein"+std::to_string(x), "protein"+std::to_string(y));
    }
}
