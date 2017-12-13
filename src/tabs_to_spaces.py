import sys
def tabs_to_spaces(file):
	with open(file, 'r') as f:
		for line in f:
			line = line.split()
			print(line[0] + " " + line[1] + " " + line[2])


def main():
	tabs_to_spaces(sys.argv[1])

main()
