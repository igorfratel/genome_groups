import sys

#Takes an input file, creates a dictionary and parses a blast file to delete the proteins not in the dictionary
def parse_neighbourhood(filename, blast_file):
	proteins = []
	with open(filename) as f:
		for line in f:
			line = line.split()
			if ((line[0] == "." and line[1] != "cds") or line[0] == "-->"): #protein
				pid = line[4]
				proteins.append(pid)

	with open(blast_file) as f:
		for line in f:
			line = line.split()
			if (line[0] in proteins and line[1] in proteins):
				print(line[0]+" "+line[1]+" "+line[2])


def main():
	#argv[1] <--- neighborhoods file name
	#argv[1] <--- blast file name
	parse_neighbourhood(sys.argv[1], sys.argv[2])

main();
