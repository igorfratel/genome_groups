#include "GenomicNeighborhood.h"

int main() {
	protein_info_t my_prot;
	GenomicNeighborhood my_neighborhood ("AABW01000001.1", "Rickettsia sibirica");
	my_neighborhood.set_seed("rsib_orf887", "EAA26079.1", "815461..817893");
	my_neighborhood.add_protein("rsib_orf877", "EAA26069.1", "803087..804568");
	my_neighborhood.add_protein("rsib_orf878", "EAA26070.1", "804722..806695");
	my_neighborhood.add_protein("rsib_orf879", "EAA26071.1", "806695..807009");
	my_prot = my_neighborhood.get_seed();
	std::cout << "printing seed: \n";
	std::cout << my_prot.locus << " ";
	std::cout << my_prot.pid << " ";
	std::cout << my_prot.cds << "\n";

	std::cout << "printing proteins: \n";
	for(GenomicNeighborhood::iterator it = my_neighborhood.begin(); it != my_neighborhood.end(); ++it) {
		my_prot = *it;
		std::cout << my_prot.locus << " ";
		std::cout << my_prot.pid << " ";
		std::cout << my_prot.cds << "\n";
	}
}