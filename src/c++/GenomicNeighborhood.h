#ifndef __GENOMIC_NEIGHBORHOOD_H__ 
#define __GENOMIC_NEIGHBORHOOD_H__ 

#include <iostream>
#include <string>
#include <vector>

typedef struct {
	std::string locus;
	std::string pid;
	std::string cds;
} protein_info_t;


class GenomicNeighborhood {

	private:
		std::string accession;
		std::string organism;
		protein_info_t seed;
		std::vector<protein_info_t> proteins;


	public:
		typedef std::vector<protein_info_t>::iterator iterator;
		
		/*Receives the name and acession code of the genomic neighborhood. Creates a new object.*/
		GenomicNeighborhood (std::string, std::string);
		
		/*Receives the locus, pid and cds of the anchor/seed protein and adds that info to the object*/
		void set_seed(std::string, std::string, std::string);
		
		/*Receives the locus, pid and cds of a protein and adds the protein to the object*/
		void add_protein(std::string, std::string, std::string);
		
		/*Returns genomic neighborhood accession code*/
		std::string get_accession();
		
		/*Returns name of the genomic neighborhood's organism*/
		std::string get_organism();
		
		/*Returns anchor/seed protein*/
		protein_info_t get_seed();
		
		/*Returns number of proteins in the genomic neighborhood*/
		int protein_count();
		
		iterator begin();
		
		iterator end();
};

#endif