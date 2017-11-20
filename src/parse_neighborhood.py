import sys

#Takes an input file and formats it to be better parsed by the c++ part of the project
def parse_neighbourhood(filename):
	with open(filename) as f:
		for line in f:
			line = line.split()
			if (line[0] == "ORGANISM"): #beginning of organism
				index = line.index("accession")
				print(line[1] + "\n" + line[2] + "\n" + line[index + 3]) #genus species accession

			elif (line[0] == "." and line[1] != "cds"): #protein
				print(line[7] + "\n" + line[4] + "\n" + line[1]) #locus pid cds

			elif (line[0] == "-->"): #seed protein
				print("seed\n" + line[7] + "\n" + line[4] + "\n" + line[1]) #"seed" locus pid cds 

			elif (line[0][0] == "-"): #end of organism
				print("end")

def main():
	#argv[1] <--- file name
	parse_neighbourhood(sys.argv[1])

main();