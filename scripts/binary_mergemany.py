from genericpath import exists
from glob import glob
from hashlib import md5
from multiprocessing import Pool
from os import mkdir
from shutil import rmtree
import subprocess
import argparse
from os.path import realpath

def read_input_file(filepath):
	with open(filepath, 'r') as f:
		return list(map(str.strip, f.readlines()))

def generate_instances(stems):
	ret = [list(x) for x in zip(stems[0::2], stems[1::2])]
	if len(stems) % 2 != 0:
		ret[-1].append(stems[-1])
	return ret

def run_merge(input_stems, suffix, exe):
	cmd = [exe]
	cmd.extend(input_stems)
	input_stems_str = "\t".join(input_stems)
	output_md5 = md5(input_stems_str.encode('utf-8')).hexdigest()
	output_stem = "{}_{}".format(suffix, output_md5)
	cmd.append(output_stem)
	with open(suffix + '.md5_key', 'a') as f:
		f.write('\t'.join([output_md5, input_stems_str]) + '\n')
	oe = open(output_stem + ".oe", 'w')
	subprocess.run(cmd, shell=False, stdout=oe, stderr=oe)
	oe.close()
	pass

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="run pairwise (or n-wise) genotype merge iterations")
	parser.add_argument('input_stems', help="file with stems to merge")
	# parser.add_argument('-n', '--instances_per_merge', help="number of stems to include in each merge", default=2)
	parser.add_argument('-e', '--executable', help="pointer to mergemany", default="/home/np29/o2bin/mergemany") # TODO: make a switch to handle packed/eigenstrat
	parser.add_argument('-p', "--max_processes", help="Maximun concurrent processes to run", default=20)
	parser.add_argument('-o', "--overwrite", help="overwrite merge directories", action="store_true")

	args = parser.parse_args()

	input_stems = realpath(args.input_stems)
	# instances = int(args.instances_per_merge)
	mergemany_path = realpath(args.executable)
	max_processes = int(args.max_processes)
	overwrite = args.overwrite

	stems = read_input_file(input_stems)
	
	suffix = 0
	while len(stems) > 1:
		instances = generate_instances(stems)
		output = "merge{}".format(suffix)
		output_w_dir = "{}/{}".format(output, output)
		if exists(output) and overwrite:
			rmtree(output)
		mkdir(output)
		# new_stems = []
		p = Pool(processes=min(len(stems), max_processes))
		# for inputs in instances:
		p.starmap(run_merge, [(x, output_w_dir, mergemany_path) for x in instances])
			# new_stems.append(output_w_dir)
		p.close()
		trash_dirs = glob("./trashdir*")
		for trash_dir in trash_dirs:
			try:
				rmtree(trash_dir)
			except:
				print("Error while deleting file : ", trash_dir)
		suffix += 1
		stems = [x[:-4] for x in glob((output_w_dir + '*.ind'))]
	print("Done! Final merge in: {}".format(output))