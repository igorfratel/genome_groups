HOW TO RUN:  

- 1: ./make  
- 2: ./python parse_neighborhood.py "file with neighborhoods"  
- 3: ./neighborhood_comparator <full or partial> <args according to chosen mode>  


ARGUMENTS FOR EACH EXECUTION MODE:  

full --> default execution  
    -e --execution_mode full  
    -n --neighborhoods_filename  
    -s --prot_sim_filename  
    -f --formatted_prot_filename  
    -p --protein_comparing  
    -m --num_prot  
    -t --stringency  
    -g --genome_comparing  
    -o --output  

partial --> Already has the similarities between the proteins.  
    -e --execution_mode partial  
    -n --neighborhoods_filename  
    -s --prot_sim_filename  
    -m --num_prot  
    -t --stringency  
    -g --genome_comparing  
    -o --output  
