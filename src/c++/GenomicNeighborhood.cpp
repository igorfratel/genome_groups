#include "GenomicNeighborhood.h"

typedef std::vector<protein_info_t>::iterator iterator;

GenomicNeighborhood::GenomicNeighborhood (std::string accession, std::string organism) {
	/*Receives the name and acession code of the genomic neighborhood. Creates a new object.*/
	this->accession = accession;
	this->organism = organism;
}

void GenomicNeighborhood::set_seed(std::string locus, std::string pid, std::string cds) {
	/*Receives the locus, pid and cds of the anchor/seed protein and adds that info to the object*/
	protein_info_t my_prot;
	my_prot.locus = locus;
	my_prot.pid = pid;
	my_prot.cds = cds;
	seed = my_prot;
}

void GenomicNeighborhood::add_protein(std::string locus, std::string pid, std::string cds) {
	/*Receives the locus, pid and cds of a protein and adds the protein to the object*/
	protein_info_t my_prot;
	my_prot.locus = locus;
	my_prot.pid = pid;
	my_prot.cds = cds;
	proteins.push_back(my_prot);
}

std::string GenomicNeighborhood::get_accession() {return accession;}

std::string GenomicNeighborhood::get_organism() {return organism;}

protein_info_t GenomicNeighborhood::get_seed() {return seed;}

int GenomicNeighborhood::protein_count() {return proteins.size();}

iterator GenomicNeighborhood::begin() {return proteins.begin();}

iterator GenomicNeighborhood::end() {return proteins.end();}
