COMO EXECUTAR:

- 1: ./make
- 2: ./python parse_neighborhood.py "arquivo com vizinhanças gênicas"
- 3: ./projeto <0 ou 1 ou 2> <args referentes ao modo escolhido>

Modo 0 --> //default execution
            //neighborhoods_file, prot_sim_filename, formatted_prot_file, protein_homology_method, num_prot, stringency,
            //genome_clustering_method, genome_sim_filename

Modo 1 --> 	//Already has the similarities between the proteins. Needs to cluster them and the neighborhoods.
            //neighborhoods_file, prot_sim_filename, num_prot, stringency, genome_clustering_method,
            //genome_sim_filename
