import subprocess
import networkx as nx
def homology_detection(format_file, method):
    #Receives the user's preferred protein homology/orthology detection method and runs it on a file already formatted to be its input.
    #Writes the results to "homology_detection_output.txt" and returns this filename.
    if (method == 'nc'):
        #Will the file contain the fasta sequences or the Blast bit-score? Should I run Blast too?
        subprocess.run(["NC_standalone", "-f " + format_file, " -o homology_detection_output.txt"]) #I'm using the default threshold (!!!)
    #What about the other methods? (!!!)

    return ("homology_detection_output.txt")

def protein_clustering(similarities_file):
    #Receives the similarities file and stores them in a dictionary of pairs of proteins
    sim_dict = {}
    with open(similarities_file, 'r') as f:
        for line in f:
            (key1, key2, val) = line.split()
            sim_dict[key1, key2] = float(val)
    return sim_dict

    #Binary clustering using connected components
    #clusters = nx.Graph()
    #with open(similarities_file, 'r') as f:
    #    for line in f:
    #        line = line.split()
    #        clusters.add_edge(line[0], line[1])
    #components = nx.connected_components(clusters) 
    #with open("protein_clustering_output.txt", 'w') as out:
    #    for component in components:
    #        out.write(str(component) + "\n")
    #return ("protein_clustering_output")