from argparse import ArgumentParser

LEGACY_COVERAGE_LINE_KEY = "mean depth:"
COV_TABLE_LABEL = "coverage report (autosomes):"

def extract_covs(log_file):
	with open(log_file, 'r') as f:
		while True:
			line = f.readline()
			if LEGACY_COVERAGE_LINE_KEY in line:
				break
		legacy_line = line.strip().split()
		output = {
			"sample_id" : legacy_line[0],
			"legacy_snps" : legacy_line[6],
			"legacy_mean_depth" : legacy_line[4]
			}
		while True:
			line = f.readline()
			if COV_TABLE_LABEL in line:
				break
		headers = f.readline().strip().replace('Mean Depth', 'mean_depth').split()
		cov_lines = []
		while True:
			line = f.readline()
			if "coverage" not in line:
				break
			cov_lines.append(line.strip().split())
	for idx, header in enumerate(headers):
		if idx == 0:
			continue
		for cov_line in cov_lines:
			output.update({f'{cov_line[1]}_{header}'.lower() : cov_line[idx+1]})
	return output

if __name__ == '__main__':
	parser = ArgumentParser()
	parser.add_argument('--header', action="store_true")
	parser.add_argument('log_file', type=str, nargs='+')
	args = parser.parse_args()

	output_list = []
	for log_file in args.log_file:
		output_list.append(extract_covs(log_file=log_file))

	if args.header:
			print('\t'.join(output_list[0].keys()))
	for output_table in output_list:
		print('\t'.join(output_table.values()))
