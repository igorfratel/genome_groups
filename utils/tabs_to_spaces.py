import sys
def tabs_to_spaces(file):
	with open(file, 'r') as f:
		for line in f:
			line = line.split()
			print(line[0] + " " + line[1] + " " + line[2] + " " + line[3] + " " + line[4] + " " + line[5] + " " + line[6])


def main():
	tabs_to_spaces(sys.argv[1])

main()
