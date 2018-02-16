#include "../src/ProteinCollection.h"

int main() {
  std::string input;
  std::cin >> input;
  unsigned int test_size = std::stoul(input);

  ProteinCollection my_graph(test_size);
  for (unsigned int i = 0; i < test_size; i++)
    my_graph.add_protein("protein" + std::to_string(i));

  my_graph.connect_proteins("protein1", "protein2", 0.05);
}
