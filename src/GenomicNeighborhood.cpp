#include "GenomicNeighborhood.h"

GenomicNeighborhood::GenomicNeighborhood (const std::string &accession_code) {
	this->accession = accession_code;
}

/**
 *Receives the locus, pid and cds of an anchor/seed protein and adds that info to the object
 */
void GenomicNeighborhood::add_seed(const std::string &locus, const std::string &pid, const std::string &cds) {
	protein_info_t my_prot;
	my_prot.locus = locus;
	my_prot.pid = pid;
	std::vector<int> cds_vector = parse_cds(cds);
	my_prot.cds_begin = cds_vector[0];
	my_prot.cds_end = cds_vector[1];
	seeds.push_back(my_prot);
}

/**
 *Receives the locus, pid and cds of a protein and adds the protein to the object
 */
void GenomicNeighborhood::add_protein(const std::string &locus, const std::string &pid, const std::string &cds) {
	protein_info_t my_prot;
	my_prot.locus = locus;
	my_prot.pid = pid;
	std::vector<int> cds_vector = parse_cds(cds);
	my_prot.cds_begin = cds_vector[0];
	my_prot.cds_end = cds_vector[1];
	proteins.push_back(my_prot);
}

/**
 *Receives the cds string in a format like "534..345" and splits it in the two composing numbers
 *Returns a vector with two positions.
 */
std::vector<int> GenomicNeighborhood::parse_cds(std::string cds) {
	std::vector<int> cds_array;
	int temp;
	std::replace(cds.begin(), cds.end(), '.', ' ');  // replace '.' by ' '
	std::stringstream ss(cds);
    while (ss >> temp)
    cds_array.push_back(temp);
    return cds_array;
}

/**
 *Returns the first coordinate of the genomic neighborhood
 */
int GenomicNeighborhood::get_first_cds() const {
	return proteins[0].cds_begin;
}

/**
 *Returns the last coordinate of the genomic neighborhood
 */
int GenomicNeighborhood::get_last_cds() const{
	return proteins.back().cds_end;
}

/**
 *Receives an index and returns the corresponding protein in the neighborhood sequence
 */
std::string GenomicNeighborhood::get_pid(int index) const{
	if ((unsigned int)index >= proteins.size()) return ".";
	return proteins[index].pid;
}

/**
 *Returns a vector containing the proteins in order
 */
std::vector<protein_info_t> GenomicNeighborhood::get_protein_vector() const{
	return proteins;
}

/**
 *Returns a vector containing the proteins in reverse order
 */
std::vector<protein_info_t> GenomicNeighborhood::get_reverse_protein_vector() const{
    std::vector<protein_info_t> reverse_proteins = proteins;
    std::reverse(reverse_proteins.begin(), reverse_proteins.end());
	return reverse_proteins;
}


std::string GenomicNeighborhood::get_accession() const{return accession;}

std::vector<protein_info_t> GenomicNeighborhood::get_seeds() const{return seeds;}

int GenomicNeighborhood::protein_count() const{return proteins.size();}
