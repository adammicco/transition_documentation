import argparse

def reverse_compliment(seq):
	base_map = {
		'a' : 't',
		't' : 'a',
		'c' : 'g',
		'g' : 'c',
		'A' : 'T',
		'T' : 'A',
		'C' : 'G',
		'G' : 'C'
	}
	ret = []
	for base in reversed(seq):
		try:
			ret.append(base_map[base])
		except KeyError:
			ret.append(base)
	return ''.join(ret)


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('-c', '--columns', action='store_true', help='set this flag if the input file has multiple colums that should be considered separate sequences')
	parser.add_argument('input_file')

	args = parser.parse_args()
	columns = args.columns
	input_file = args.input_file

	with open(input_file, 'r') as f:
		if columns:
			lines = [list(map(str.strip, line.strip().split())) for line in f]
		else:
			lines = [line.strip() for line in f]

	if columns:
		rev_comp_line = ['\t'.join([reverse_compliment(col) for col in line]) for line in lines]
	else:
		rev_comp_line = [reverse_compliment(line) for line in lines]

	for line in rev_comp_line:
		print(line)
