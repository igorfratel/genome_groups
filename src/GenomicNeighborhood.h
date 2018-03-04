#ifndef __GENOMIC_NEIGHBORHOOD_H__
#define __GENOMIC_NEIGHBORHOOD_H__

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>

typedef struct {
	std::string locus;
	std::string pid;
	int cds_begin;
	int cds_end;
} protein_info_t;


class GenomicNeighborhood {

	private:
		std::string accession;
		std::vector<protein_info_t> seeds;
		std::vector<protein_info_t> proteins;

		/*Receives the cds string in a format like "534..345" and splits it in the two composing numbers*/
		/*Returns a vector with two positions.*/
		std::vector<int> parse_cds(const std::string cds);


	public:

		GenomicNeighborhood (const std::string &acession_code);

		/*Receives the locus, pid and cds of an anchor/seed protein and adds that info to the object*/
		void add_seed(const std::string &locus, const std::string &pid, const std::string &cds);

		/*Receives the locus, pid and cds of a protein and adds the protein to the object*/
		void add_protein(const std::string &locus, const std::string &pid, const std::string &cds);

		/*Returns genomic neighborhood accession code*/
		std::string get_accession() const;

		/*Returns anchor/seed proteins*/
		std::vector<protein_info_t> get_seeds() const;

		/*Receives an index and returns the corresponding protein in the neighborhood sequence*/
		std::string get_pid(int index) const;

		/*returns a vector containing the proteins in order*/
		std::vector<protein_info_t> get_protein_vector() const;

		/*Returns number of proteins in the genomic neighborhood*/
		int protein_count() const;

		/*Returns the first coordinate of the genomic neighborhood*/
		int get_first_cds() const;

		/*Returns the last coordinate of the genomic neighborhood*/
		int get_last_cds() const;
};

#endif
