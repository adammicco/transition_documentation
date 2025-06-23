import argparse
from os import rename

def read_driver(driver_filename):
	driver_dict = {}
	with open(driver_filename, 'r') as driver:
		for line in driver:
			fields = line.split()
			pulldownID = fields[0]
			geneticID = fields[1]
			sex = fields[2]
			groupID = fields[3]

			if pulldownID in driver_dict.keys():
				if driver_dict[pulldownID] != (geneticID, sex, groupID):
					raise KeyError("Pulldown ID {} repeated in driver file with differing attributes!")
			else:
				driver_dict.update({pulldownID : (geneticID, sex, groupID)})
	return driver_dict

def update_ind_file(ind_filename, driver_dict):
	with open(ind_filename, 'r') as old_ind:
		with open("{}.tmp".format(ind_filename), 'w') as new_ind:
			for line in old_ind:
				fields = line.split()
				pulldownID = fields[0]
				sex = fields[1]
				
				if pulldownID in driver_dict.keys():
					new_ind.write('\t'.join(driver_dict[pulldownID]) + '\n')
				else:
					new_ind.write('\t'.join([pulldownID, sex, "Ignore"]) + '\n')
	rename(ind_filename, '{}.pdOrig'.format(ind_filename))
	rename('{}.tmp'.format(ind_filename), ind_filename)
	pass

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="set up repulled ind files for merge, ignoring unneeded samples")

	parser.add_argument('-d', '--driver', help="driver file of samples to keep in [pulldown ID]\t[new genetic ID]\ts[ex]\t[group ID] format", required=True)
	parser.add_argument("ind_files", help="ind files to modify", nargs="+")

	args = parser.parse_args()

	driver_filename = args.driver
	ind_files = args.ind_files

	driver_dict = read_driver(driver_filename)

	for ind_file in ind_files:
		update_ind_file(ind_file, driver_dict)