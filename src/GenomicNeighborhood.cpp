#include "GenomicNeighborhood.h"

GenomicNeighborhood::GenomicNeighborhood (const std::string &accession_code) {
	this->accession = accession_code;
}

void GenomicNeighborhood::add_seed(const std::string &locus, const std::string &pid, const std::string &cds) {
	/*Receives the locus, pid and cds of an anchor/seed protein and adds that info to the object*/
	protein_info_t my_prot;
	my_prot.locus = locus;
	my_prot.pid = pid;
	std::vector<int> cds_vector = parse_cds(cds);
	my_prot.cds_begin = cds_vector[0];
	my_prot.cds_end = cds_vector[1];
	seeds.push_back(my_prot);
}

void GenomicNeighborhood::add_protein(const std::string &locus, const std::string &pid, const std::string &cds) {
	/*Receives the locus, pid and cds of a protein and adds the protein to the object*/
	protein_info_t my_prot;
	my_prot.locus = locus;
	my_prot.pid = pid;
	std::vector<int> cds_vector = parse_cds(cds);
	my_prot.cds_begin = cds_vector[0];
	my_prot.cds_end = cds_vector[1];
	proteins.push_back(my_prot);
}

std::vector<int> GenomicNeighborhood::parse_cds(std::string cds) {
	/*Receives the cds string in a format like "534..345" and splits it in the two composing numbers*/
	/*Returns a vector with two positions.*/
	std::vector<int> cds_array;
	int temp;
	std::replace(cds.begin(), cds.end(), '.', ' ');  // replace '.' by ' '
	std::stringstream ss(cds);
    while (ss >> temp)
    cds_array.push_back(temp);
    return cds_array;
}

int GenomicNeighborhood::get_first_cds() const {
	/*Returns the first coordinate of the genomic neighborhood*/
	return proteins[0].cds_begin;
}

int GenomicNeighborhood::get_last_cds() const{
	/*Returns the last coordinate of the genomic neighborhood*/
	return proteins.back().cds_end;
}

/*Receives an index and returns the corresponding protein in the neighborhood sequence*/
std::string GenomicNeighborhood::get_pid(int index) const{
	if ((unsigned int)index >= proteins.size()) return ".";
	return proteins[index].pid;
}

/*returns a vector containing the proteins in order*/
std::vector<protein_info_t> GenomicNeighborhood::get_protein_vector() const{
	return proteins;
}

std::string GenomicNeighborhood::get_accession() const{return accession;}

std::vector<protein_info_t> GenomicNeighborhood::get_seeds() const{return seeds;}

int GenomicNeighborhood::protein_count() const{return proteins.size();}
