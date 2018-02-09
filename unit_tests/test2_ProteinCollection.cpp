#include "../src/ProteinCollection.h"

int main() {
  const size_t test_size = 1000000;

  ProteinCollection my_graph(test_size);
  for (unsigned int i = 0; i < test_size; i++)
    my_graph.add_protein("protein" + std::to_string(i));

  my_graph.connect_proteins("protein1", "protein500", 0.05);

  std::cout << my_graph.get_similarity("protein1", "protein500") << "\n";
}
