import sys

#Takes an input neighborhoods file and prints the pid of all proteins
def parse_neighbourhood(filename):
	with open(filename) as f:
		for line in f:
			line = line.split()
			if ((line[0] == "." and line[1] != "cds") or line[0] == "-->"): #protein
				pid = line[4]
				if(pid != "."):
					print(pid, end=' ')


def main():
	#argv[1] <--- file name
	parse_neighbourhood(sys.argv[1])

main();
