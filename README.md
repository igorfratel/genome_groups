# Genome Comparer

Software that calculates a similarity score between genomic neighborhoods using the similarities between their proteins as
an additional information.  

## How to run:  

- 1: Go to /src folder
- 2: ./make  
- 3: ./neighborhood_comparer <full or partial> <args according to chosen mode>  


## Command line options
```
full --> default execution  
    -e --execution_mode full    
    -n --neighborhoods_filename  "File containing the genomic neighborhoods"   
    -s --prot_sim_filename  "File containing pairs of proteins and their similarities"  
    -f --formatted_prot_filename  "File already formatted as the input for the homology detection method (eg: nc method receives a blast file containing the pid's and the bitscore separated by spaces)"  
    -p --protein_comparing  "Method for comparing proteins (default: nc)"  
    -d --num_residues  "Total number of residues in unique protein set of formatted_prot_filename. Used in the nc method. Defaults to 537, which is the default for the NC_Standalone program."  
    -t --prot_stringency  "Minimum similarity required to treat two proteins as a related pair (default 0)"  
    -r --neigh_stringency  "Minimum threshold to display the similarity between two neighborhoods (default 0)"  
    -m --mismatch  "Similarity value used when protein pair similarity is below prot_stringency (default -1)"
    -g --neigh_comparing  "Method for comparing genomic neighborhoods (default: porthodom method)"  
    -G --gap_score  "Value added to similarity when a gap is chosen. Used when the neigh_comparing option is an alignment algorithm (default 0)"  
    -o --output  "Where the neighborhood similarities should be written (outputs to stdout if the filename is - or if not used)"  
    -a --pairings_filename  "Where the chosen pairings between proteins in the neighborhoods should be written (if not specified, does not generate a pairings file)"  


partial --> Already has the similarities between the proteins.  
    -e --execution_mode partial    
    -n --neighborhoods_filename    
    -s --prot_sim_filename    
    -l --normalize_prot_sim" "Indicates that the protein similarities file should be normalized"  
    -t --prot_stringency  
    -r --neigh_stringency  
    -m --mismatch
    -g --neigh_comparing  
    -G --gap_score  
    -o --output  
    -a --pairings_filename  

Help option: -h --help  
```
## Additional information
* protein scoring methods: nc or raw score in the partial mode (just use a file with the desired scores and normalize if needed).

* neighborhood scoring methods: porthodom, porthodomO2, global_alignment.

* prot_sim_filename format: "prot1 prot2 score" <-- tabs or whitespaces (nc uses whitespaces)

* Output format: "accession1    cds_begin1    cds_end1    accession2    cds_begin2    cds_end2    score" <-- tabs  
    Note: If the neigh_comparing option is set to global_alignment, this file will have an additional last column consisting
    of a "+" if the highest scoring alignment is the regular one or a "-" if the highest scoring alignment was reached
    inverting the order of the second neighborhood.

* Pairings format:  
  header: ">accession1    cds_begin1    cds_end1    accession2    cds_begin2    cds_end2" <-- tabs  
  pairings: "prot1 prot2 sim" <-- tabs

* To add a new neighborhood scoring method: include the filename containing the scoring method functions in genome_grouping.h,
add a new "if else" clause at the genome_clustering function in the genome_grouping.cpp file comparing the genomic neighborhoods using the new scoring function. Add new files to Makefile.
