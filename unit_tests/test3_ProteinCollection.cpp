#include "../src/ProteinCollection.h"
#include <iostream>
#include <cstdlib>

int main() {
	std::string input;
    std::cin >> input;
    int test_size = std::stoi(input);
	ProteinCollection my_graph(10000);
    for(int i = 0; i < 10000; i++){
        my_graph.add_protein("protein"+std::to_string(i));
    }

    int x = 0;
    int y = 0;
    for(int i = 0; i < test_size; i++){
        x = rand()%10000;
        y = rand()%10000;
        my_graph.connect_proteins("protein"+std::to_string(x), "protein"+std::to_string(y), 0.05);
        my_graph.get_similarity("protein"+std::to_string(x), "protein"+std::to_string(y));
    }
}
