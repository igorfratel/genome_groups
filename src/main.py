import sys #sys.argv[], sys.exit()
import subprocess
import argparse 
import os.path #os.path.isfile()
def main():
    #User should be able to specify:
    # 1 - the preferred method for protein homology detection.
    # 2 - the name of a format file, to be used as input to the chosen protein homology detection method
    # 4 - the preferred method for genome neighborhoods clustering.
    # 5 - the name of the genomic neighborhoods file.
    # 6 - the name of the database file.
    # 7 - the types of visualization, neighborhoods(1), graph(2), sets(3).
    # 8 - the name of a similarities file, if the user has already used a protein homology detection method on his dataset
    # 9 - the name of a clusters file, if the user has already used a protein clustering method on his dataset 

    #Definition of all the arguments using the argparse library
    parser = argparse.ArgumentParser(description = "Receives a number of genomic neighborhoods, uses a protein "
                                                   "clustering method to determine similarity information on these "
                                                   "neighborhoods and clusters them using a genome clustering method. "
                                                   "Displays results in the specified format.")

    parser.add_argument('-g', '--genome_clustering', help = "Method for grouping genomic neighborhoods."
                        "(default: simple method)", default = "simple")
    parser.add_argument('-p', '--protein_homology', help = "Method for homology/othology detection on proteins."
                        "(default: Neighborhood Correlation)", default = "nc")
    parser.add_argument('-f', '--formatted_prot_file', help = "File already formatted as the input for the homology/orthology "
                        "detection method")
    parser.add_argument('-v', '--visualization', help = "How information should be displayed: neighborhoods[1]," 
                        "graph[2] and/or sets[3] (default: neighborhoods[1])", type = int, choices = [1, 2, 3],
                         nargs = "*", default = 1)
    parser.add_argument('-n', '--neighborhoods_file', help = "File containing the genomic neighborhoods with the locus "
                        "tags and pids of their proteins.", required = True)
    parser.add_argument('-d', '--database_file', help = "File containing the pid of each protein and their fasta "
                        "sequences")
    parser.add_argument('-s', '--prot_sim_file', help = "File containing pairs of proteins and their similarities "
                        "calculated using some homology/orthology detection method")
    parser.add_argument('-c', '--prot_clusters_file', help = "File containing the already clustered proteins with their locus "
                        "tags and pids")
    parser.add_argument('-m', '--num_prot', help = "Number of unique proteins", type = int, required = True)
    parser.add_argument('-t', '--stringency', help = "Minimum similarity required to treat two proteins as a related pair", type = float, default = 0.0)
    

    #Parsing the arguments
    args = parser.parse_args()
    #Prints this if arguments are valid
    print("Executing with following arguments: \n" + 
          "genome clustering method: " + args.genome_clustering + "\n"
          "protein homology method: " + args.protein_homology + "\n" 
          "formattted proteins file: " + str(args.formatted_prot_file) + "\n"
          "visualization options: " + str(args.visualization) + "\n"
          "neighborhoods file: " + str(args.neighborhoods_file) + "\n"
          "database file: " + str(args.database_file) + "\n"
          "protein clusters file: " + str(args.prot_clusters_file)
          "protein similarities file: " + str(args.prot_sim_file) + "\n"
          "number of proteins: " + str(args.num_prot) + "\n"
          "Stringency for protein clustering: " + str(args.stringency) + "\n"
          )


    #Check validity of arguments (Only checking valid filenames!!!)
    if (not os.path.isfile(args.neighborhoods_file))  :
        sys.exit("ERROR: Genomic neighborhoods file does not exist")
    if (args.database_file and not os.path.isfile(args.database_file)):
        sys.exit("ERROR: Database file does not exist")
    if (args.prot_clusters_file and not os.path.isfile(args.prot_clusters_file))  :
    if (args.prot_sim_file and not os.path.isfile(args.prot_sim_file)):
        sys.exit("ERROR: Protein similarities file does not exist")
    if (args.formatted_prot_file and not os.path.isfile(args.formatted_prot_file)):
        sys.exit("ERROR: Formatted proteins file does not exist")


    if (args.prot_clusters_file):
      print("Não implementado")
      #Already has clustered proteins. Just need to cluster the genomic neighborhoods
      #genome_clustering_output = gengroup.genome_clustering(args.neighborhoods_file, args.prot_clusters_file, args.genome_clustering)


    elif (args.prot_sim_file):
      #Already has the similarities between the proteins. Need to cluster them and the genomic neighborhoods.
      #protein_clustering_output = protgroup.protein_clustering(args.prot_sim_file)
      #genome_clustering_output = gengroup.genome_clustering(args.neighborhoods_file, protein_clustering_output, args.genome_clustering)
      p = subprocess.run(["./projeto", "1", "args.neighborhoods_file", "args.prot_sim_file", "args.num_prot", "args.stringency", "args.genome_clustering"], stdout=subprocess.PIPE);
    elif (args.formatted_prot_file):
      #default execution
      #homology_detection_output = protgroup.homology_detection(args.formatted_prot_file, args.protein_homology)
      #protein_clustering_output = protgroup.protein_clustering(homology_detection_output)
      #genome_clustering_output = gengroup.genome_clustering(args.neighborhoods_file, protein_clustering_output, args.genome_clustering)
      subprocess.Popen(["./projeto", "2", "args.neighborhoods_file", "args.formatted_prot_file", "protein_homology",
                      "args.num_prot", "args.stringency", "args.genome_clustering"], stdout=subprocess.PIPE);

    else:
      sys.exit("ERROR: Invalid combination of parameters")


    #vis.visualize(genome_clustering_output, args.visualization)


if __name__ == "__main__":
    main()