#include "GenomicNeighborhood.h"

int main() {
	protein_info_t my_prot;
	std::vector<protein_info_t> seeds;
	GenomicNeighborhood my_neighborhood ("AABW01000001.1", "Rickettsia sibirica");
	my_neighborhood.add_seed("rsib_orf887", "EAA26079.1", "815461..817893");
	my_neighborhood.add_protein("rsib_orf877", "EAA26069.1", "803087..804568");
	my_neighborhood.add_protein("rsib_orf878", "EAA26070.1", "804722..806695");
	my_neighborhood.add_protein("rsib_orf879", "EAA26071.1", "806695..807009");
	seeds = my_neighborhood.get_seeds();
	std::cout << "printing seeds: \n";
	for(int i = 0; i < seeds.size(); i++){
		std::cout << seeds[i].locus << " ";
		std::cout << seeds[i].pid << " ";
		std::cout << seeds[i].cds_begin << "..";
		std::cout << seeds[i].cds_end << "\n";
	}

	std::cout << "printing proteins: \n";
	for(GenomicNeighborhood::iterator it = my_neighborhood.begin(); it != my_neighborhood.end(); ++it) {
		my_prot = *it;
		std::cout << my_prot.locus << " ";
		std::cout << my_prot.pid << " ";
		std::cout << my_prot.cds_begin << "..";
		std::cout << my_prot.cds_end << "\n";
	}

	std::cout << "printing cds using integers: \n";
	std::vector<int> cds_vector = my_neighborhood.get_cds();
	std::cout << cds_vector[0] << ".." << cds_vector[1] << "\n";

	std::cout << "printing cds using string: \n";
	std::cout << my_neighborhood.get_cds_string() << "\n";
}
