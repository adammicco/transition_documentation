import argparse
import re
from os import mkdir, chdir
from os.path import exists, abspath, basename
import subprocess

# TODO: Update to real filepaths
SNP_SET_FILENAMES = {
	"HO" : "/n/groups/reich/array_snp_lists/nextGenHO/v0/nextGenHO_v0.snp",
	"1240k_classic" : "/n/groups/reich/matt/pipeline/static/snp_split/1240K.snp",
	"2M" : "/n/groups/reich/array_snp_lists/nextGenCapture/cov-0.1__maf-0.01/v3/nextGenCapture_cov0.1_maf0.01_v3.snp",
	"3M" : "/n/groups/reich/matt/pipeline/static/twist_plus.v6.snp",
	"Master" : "/n/groups/reich/matt/pipeline/static/twist_plus.v6.snp"
}
# TODO: replace par_template, sbatch_template with real, absolute filepaths
PAR_DEFAULT = "/n/groups/reich/matt/pipeline/static/snp_split/par_template"
SBATCH_TEMPLATE = "/n/groups/reich/matt/pipeline/static/snp_split/sbatch_template"
CONVERTF_DEFAULT = '/home/np29/o2bin/convertf'
SUBMIT_COMMAND = 'sbatch'

def build_par_file(master_stem, snp_file_path, driver_file_path, output_name, par_template=PAR_DEFAULT):
	"""
	Constructs a par file to feed to convertf.

	RETURNS:
		par_file_path: Path pointing to newly created par file
	ACCEPTS:
		master_stem: Path to master Eigenstrat file.
		snp_file_path: Path to snp_file to use for convertf
		driver_file_path: Path to driver file to use in convertf.
		output_name: Name to use for output Eigenstrat files. Should be in SNPSET/name_SNPSET format.
		par_template: Path to par file template to use. Defaults to boilerplate par file specified in PAR_DEFAULT
	"""
	with open(par_template, 'r') as template:
		par_file_content = template.read()
	
	par_file_content = par_file_content.replace("MASTER_STEM", master_stem)
	par_file_content = par_file_content.replace("SNP_FILE_PATH", snp_file_path)
	par_file_content = par_file_content.replace("DRIVER_PATH", driver_file_path)
	par_file_content = par_file_content.replace("OUTPUT_NAME", output_name)
	
	par_file_path = output_name + "_par"
	with open(par_file_path, 'w') as par_file:
		par_file.write(par_file_content)
	return abspath(par_file_path)

def build_sbatch_script(output_name, par_file, sbatch_template=SBATCH_TEMPLATE, convertf_execuable=CONVERTF_DEFAULT):
	"""
	Constructs an sbatch script or similar scheduler submission script using a template
	
	RETURNS:
		sbatch_path: Path to the newly created submission script file
	ACCEPTS:
		output_name: Name to use for output submission script files. Should be in SNPSET/name_SNPSET format.
		par_file: Path to the parfile to run convertf using.
		sbatch_template: Path to submission script template to use. Defaults to boilerplate submission script specified in SBATCH_TEMPLATE
		convertf_execuable: Path to the convertf executable file to use. Defaults to the path in CONVERTF_DEFAULT.
		"""
	with open(sbatch_template, 'r') as template:
		sbatch_file_content = template.read()
	# TODO: Scale resource requirements to size of snpset. Verify with Matt that this is worth it due to Cokie's changes.
	sbatch_file_content = sbatch_file_content.replace("OUTPUT_NAME", output_name.split("/")[-1])
	sbatch_file_content = sbatch_file_content.replace("CONVERTF_EXE", convertf_execuable)
	sbatch_file_content = sbatch_file_content.replace("PAR_FILE", par_file)
	
	sbatch_path = output_name + '.sh'
	with open(sbatch_path, 'w') as sbatch_file:
		sbatch_file.write(sbatch_file_content)
	return abspath(sbatch_path)

def get_file_len(target_file, header=False):
	"""
	Given a file path, return the line count. Can account for a header.

	RETURNS:
		line_count: Number of non-header lines in the target file
	ACCEPTS:
		target_file: Path to the file.
		header: Set to True if the file has a header that should not be counted in the line count. Default is False.
	"""
	with open(target_file) as file:
		for i, l in enumerate(file):
			pass
	if header == True:
		line_count = i
	else:
		line_count = i + 1
	return line_count

def exclude_by_data_type(data_types, anno_data):
	"""
	Finds samples to exclude based on a passed in data type(s).

	RETURNS:
		exclude_set: Set of Genetic IDs that should be dropped from the output Eigenstrat files.
	ACCEPTS:
		data_types: List of data sources that should be removed from the output datasets.
		anno_data: Dict of lists representing anno file data. Keys are column headers and list indicies correspond to anno file row numbers.
	"""
	# Iterate through each data type to exclude
	exclude_set = set()
	for data_type in data_types:
		# Add the Genetic ID to the exclude set if the excluded data type is found in the Data source column in the anno file
		exclude_set.update([anno_data['Genetic ID'][i] for i, v in enumerate(anno_data['Data types in bam']) if v.lower().find(data_type.lower()) != -1])
	return exclude_set

def exclude_by_assessment(assessment_threshold, anno_data):
	"""
	Finds samples to exclude based on anno file QC assessments.

	RETURNS:
		exclude_set: Set of Genetic IDs that should be dropped from the output Eigenstrat files.
	ACCEPTS:
		assessment_threshold: Set to "QUESTIONABLE" to exclude only QUESTIONABLE_CRITICAL samples. Set to "PASS" to exclude QUESTIONABLE AND QUESTIONABLE_CRITICAL samples. Any other value will cause an error.
		anno_data: Dict of lists representing anno file data. Keys are column headers and list indicies correspond to anno file row numbers.
	"""
	exclude_set = set()
	assessment_threshold = assessment_threshold.upper()
	if assessment_threshold == 'QUESTIONABLE':
		exclude_assessment = 'QUESTIONABLE_CRITICAL'
	elif assessment_threshold == 'PASS':
		exclude_assessment = 'QUESTIONABLE'
	else:
		raise ValueError("Invalid assessment threshold! Must be either 'questionable' or 'pass'.")
	exclude_set.update([anno_data['Genetic ID'][i] for i, v in enumerate(anno_data['ASSESSMENT']) if v.upper().find(exclude_assessment) != -1])
	return exclude_set

def exclude_by_udg(udg_treatments, anno_data):	
	"""
	Finds samples to exclude based on anno file UDG treatments.

	RETURNS:
		exclude_set: Set of Genetic IDs that should be dropped from the output Eigenstrat files.
	ACCEPTS:
		udg_treatments: List of UDG treatments to exclude from the output datasets. Valid UDG treatments: 'half', 'plus', 'minus', 'mixed'. 'mixed' will exclude any sample that is a merge of more than one UDG treatment.
		anno_data: Dict of lists representing anno file data. Keys are column headers and list indicies correspond to anno file row numbers.
	"""
	exclude_set = set()
	# Iterate through each udg_treatment to exclude
	for udg_treatment in udg_treatments:
		if udg_treatment == 'mixed':
			# Exclude mixed udg treatment samples
			exclude_set.update([anno_data['Genetic ID'][i] for i,v in enumerate(anno_data['Library type (minus=no.damage.correction, half=damage.retained.at.last.position, plus=damage.fully.corrected, ds=double.stranded.library.preparation, ss=single.stranded.library.preparation)']) if len(set(filter(None, v.replace("ss.", "").replace("ds.", "").replace("..", "").split(',')))) > 1])
		else:
			exclude_set.update([anno_data['Genetic ID'][i] for i, v in enumerate(anno_data['Library type (minus=no.damage.correction, half=damage.retained.at.last.position, plus=damage.fully.corrected, ds=double.stranded.library.preparation, ss=single.stranded.library.preparation)']) if v.lower().find(udg_treatment.lower()) != -1])
	return exclude_set

def exclude_unpublished(anno_data):
	"""
	Exludes unpublished samples. Unpublished samples are determined by scanning the Publication column in the anno file for values that contain either the string "prepub" or "unpub.
	
	RETURNS:
		exclude_set: Set of Genetic IDs that should be dropped from the output Eigenstrat files.
	ACCEPTS:
		anno_data: Dict of lists representing anno file data. Keys are column headers and list indicies correspond to anno file row numbers.
	"""
	exclude_set = set()
	# TODO figure out if there is potentially a better way to ascertain publication status.
	# This is just a Python implementation of how I've been doing it when constructing the HO anno file...
	exclude_set.update([anno_data['Genetic ID'][i] for i, v in enumerate(anno_data['Publication abbreviation or plan']) if v.lower().find("unpub") != -1 or v.lower().find("prepub") != -1])
	return exclude_set

def exclude_by_custom_list(exclusion_list_path):
	"""
	Excludes samples from the output datasets based on a custom list of samples.

	RETURNS:
		exclude_set: Set of Genetic IDs that should be dropped from the output Eigenstrat files.
	ACCEPTS:
		exclusion_list_path: Path to a new line-delimited list of Genetic IDs to be dropped from the output datasets. These should match those in the master ind file.
	"""
	exclude_set = set()
	with open(exclusion_list_path, 'r') as exclusion_file:
		exclude_set.update(exclusion_file.read().split('\n'))
	return exclude_set

def include_by_custom_list(inclusion_list_path, exclude_set):
	"""
	Includes samples in the output datasets based on a custom list of samples. This fucntion is used to restore specific samples that would be otherwise filtered out with other filtering functions.
	
	RETURNS:
		exclude_set: Set of Genetic IDs that should be dropped from the output Eigenstrat files.
	ACCEPTS:
		inclusion_list_path: Path to a new line-delimited list of Genetic IDs to be restored to the output datasets if they would be otherwise dropped. These should match those in the master ind file.
		exclude_set: Existing set of Genetic IDs to be dropped from the output Eigenstrat files.
	"""
	with open(inclusion_list_path, 'r') as inclusion_file:
		for ind in inclusion_file.readlines():
			try:
				exclude_set.remove(ind.strip())
			except:
				pass
	return exclude_set

def build_driver_file(master_ind_filepath, exclusion_set, driver_filepath):
	"""
	Constructs a .driver file to be specified in convertf par files given the input master dataset and a set of samples to be excluded from the output datasets.

	RETURNS:
		driver_filepath: Path to the newly created driver file
	ACCEPTS:
		master_ind_filepath: Path to the master ind file from which samples will be read in order to ensure ordinal integrity with output driver file.
		exclusion_set: Set of Genetic IDs in the master ind file that should be marked as 'Ignore' in the driver file.
		driver_filepath: Path, including filename and extension, where the driver file should be saved.
	"""
	with open(master_ind_filepath, 'r') as ind_template:
		with open(driver_filepath, 'w') as driver:
			for line in ind_template:
				line = list(filter(None,re.split("\t| ", line.strip(), )))
				ID = line[0]
				if ID in exclusion_set:
					driver.write("{}\t{}\tIgnore\n".format(line[0], line[1]))
				else:
					driver.write("{}\t{}\t{}\n".format(line[0], line[1], line[2]))
	return abspath(driver_filepath)

def read_anno_file(anno_file_path):
	"""
	Reads in a tab-separated anno file and parses this data into a dict of lists where the keys represent column names and the lists contain column data in the same order as the input anno file. Indicies in these lists therefore correspond to anno file (as well as ind file) rows.

	RETURNS:
		anno_file_data: Dict of lists representing anno file data. Keys are column headers and list indicies correspond to anno file row numbers.
	ACCEPTS:
		anno_file_path: Path to the anno file to read in.
	"""
	with open(anno_file_path, "r") as f:
		header_line = f.readline()
		headers = header_line.split('\t')
		anno_file_data = {header : [] for header in headers}
		for line in f:
			fields = line.split('\t')
			for header in headers:
				anno_file_data[header].append(fields[headers.index(header)])
	return anno_file_data

	# Override mode does not respect filters and generates a best representative sampleset  regardless of any specified filters. WARN USER IN STDOUT.
	# filter_union mode excludes a union of filtered samples and those that are filtered out of a best representative sample. WARN USER IN STDOUT
	# Next_best mode will try to find the next best sample for each Master ID if the otherwise best representative is in the exclude_set. If no next best sample is defined, the highest coverage sample is still accepted.
	# Next_best_strict mode functions similarly to Next_best mode. It will attempt to find the next best version of each Master ID but will NOT accept the highest coverage version if no next best, non excluded, version is found.

def get_best_representatives(anno_data, preferred_data_source="", filtered_set = set(), strict = False):
	"""
	Determine 'best representative' samples for each unique Master ID in the input data.
	The best representative sample is defined as the sample with highest SNP count for a given Master ID.
	This function will find all Genetic IDs for each unique Master ID and loop through them, adding all but the best representative verison to the exclude_set.
	
	Under normal conditions this will drop all but the sample with the most SNPs hit.
	However, a set of pre-filtered samples may be provded to the function. These will be skipped over while determining the best representative sample.

	Additonally, if a preferred_data_source is set, the function will try to find the higest coverage sample for the Master ID with the specified data source.
	If one cannot be found (or it is otherwise pre-filtered out), the highest coverage sample will be accepted and all others dropped for the Master ID.

	If the strict flag is set, however, and no sample matching the preferred_data_source that is not filtered out is found, the all samples for the Master ID will be dropped.

	RETURNS:
		exclude_set: Existing set of Genetic IDs to be dropped from the output Eigenstrat files.
	ACCEPTS:
		anno_data: Dict of lists representing anno file data. Keys are column headers and list indicies correspond to anno file row numbers.
		preferred_data_source: Data source, corresponding to those in the Data Source column in the anno file to prefer while making best representative determinations. Defaults to "", meaning no preferred data source.
		filtered_set: Set of Genetic IDs that should be skipped over while making best representative determinations. Defaults to an empty set, meaning no samples are pre-filtered.
		strict: Boolean indicating if best representative determination should be run in strict mode. Defaults to False.

	"""
	exclude_set = set()

	# This was too slow
	# master_id_set = set(anno_data['Master ID'])

	# Loop through unique master IDs after building a reverse index of unique Master IDs to occurance indicies in the anno/ind files
	master_id_dict = {}
	for ix, master_id in enumerate(anno_data['Master ID']):
		if master_id in master_id_dict:
			master_id_dict[master_id].append(ix)
		else:
			master_id_dict.update({master_id : [ix]})
	for master_id in master_id_dict.keys():
		# ...for each master ID, construct a dictionary of samples, keyed on the anno file index, containing SNP count and data source information'
		match_ixs = {i : (anno_data['SNPs hit on autosomal targets (Computed using easystats on 3.2M snpset)'][i], anno_data['Data source'][i].lower()) for i in master_id_dict[master_id]}
		
		# This was way too slow...
		# match_ixs = {i : (anno_data['SNPs hit on autosomal targets'][i], anno_data['Data source'][i].lower()) for i, v in enumerate(anno_data['Master ID']) if v == master_id}
		
		# Get keys sorted by SNP count
		keys_by_cov = [k for k, v in sorted(match_ixs.items(), key = lambda x: float(x[1][0]), reverse=True)]
		if preferred_data_source != "":
			for i, k in enumerate(keys_by_cov):
				# Attempt to find the higest coverage sample with our preferred data source and drop it from the exclusion set...
				if match_ixs[k][1] == preferred_data_source.lower() and anno_data['Genetic ID'][k] not in filtered_set:
					keys_by_cov.remove(k)
					exclude_set.update([anno_data['Genetic ID'][id] for id in keys_by_cov])
					break
			# ...if we reach the end of keys_by_cov without finding a sample with our preferred data source just accept the highest coverage one and drop the rest:
				elif i == len(keys_by_cov)-1:
					if strict:
						exclude_set.update([anno_data['Genetic ID'][id] for id in keys_by_cov])
					else:
						exclude_set.update([anno_data['Genetic ID'][id] for id in keys_by_cov[1:]])
		# If there's no preferred data source, just accept the highest coverage sample and drop the rest
		else:
			for i, k in enumerate(keys_by_cov):
				if anno_data['Genetic ID'][k] not in filtered_set:
					try:
						exclude_set.update([anno_data['Genetic ID'][id] for id in keys_by_cov[i+1:]])
						break
					except IndexError:
						pass
				elif i == len(keys_by_cov)-1:
					if strict:
						exclude_set.update([anno_data['Genetic ID'][id] for id in keys_by_cov])
					else:
						exclude_set.update([anno_data['Genetic ID'][id] for id in keys_by_cov[1:]])
	return exclude_set

if __name__ == "__main__":
	Parser = argparse.ArgumentParser(description="Generate usable release dataset(s) based on master Eigenstrat file.")

	# Specify anno file for this release. This will drive filtering fuctionality.
	Parser.add_argument('--anno_file', help="Specify path to a tsv formatted anno file.", type=str)
	Parser.add_argument('--master', help="Specify path to master Eigentrat files, pulled down on largest SNP union.", required=True, type=str)
	Parser.add_argument('--label', help="Label to use on output files.", type=str, required=True)
	# Parser.add_argument('--HO', help="Generate genotype files on HO SNPs", action='store_false')
	# Parser.add_argument('--1240k_classic', help="Generate genotype files on 1240k classic SNPs", action='store_false')
	# Parser.add_argument('--2M', help="Generate genotype files on the union of 1240k-enhanced + Twist-enhanced + Yfull SNPs.", action='store_false')
	# Parser.add_argument('--3M', help="Generate genotype files on the union of 1240k-enhanced + Twist-enhanced + BigYoruba + Yfull SNPs", action='store_false')
	Parser.add_argument('--snp_sets', nargs='+', default=['Master'], choices=['HO', '1240k_classic', '2M', '3M', 'Master'])
	Parser.add_argument('--custom_list', nargs='+', help="Generate a genotype files based on a custom .snp file. Use this option to specify the path to this file.", type=str, default=[])
	Parser.add_argument('--foreground', help="Set this flag to run convertf in the foreground insted of submitting jobs to SLURM to run in parallel. NOT RECOMMENDED.", action="store_true")
	Parser.add_argument('--convertf_exe', help="Specify custom convertf executable. Default is stored in CONVERTF_DEFAULT.", type=str, default=CONVERTF_DEFAULT)
	Parser.add_argument('--sbatch_template', help="Specify a custom submission script. This could be used to adapt this script to work on other schedulers. Default is stored in SBATCH_TEMPLATE", type=str, default=SBATCH_TEMPLATE)
	Parser.add_argument('--working_dir', help="Specify the working directory. Defaults to the currect directory.", default=".")

	# Filter arguments
	# TODO make a list of anno file data types to use a choices for --data_exclude
	Parser.add_argument('--data_exclude', help="Exclude samples based on thier data source. Requires a valid anno file.", nargs = '+', default=[])
	Parser.add_argument('--assessment_threshold', help="Threshold at which to exclude samples due to QC assessment. Defaults to include all samples, 'questionable' will exclude QUESTIONABLE_CRITICAL samples and allow all others, 'pass' will exlude both QESTIONABLE and QUESTIONABLE_CRITICAL samples. Requires a valid anno file.", choices=['questionable', 'pass'], type=str, default="")
	Parser.add_argument('--UDG_exclude', help="Exlude samples based on UDG treatment. Merged samples including libraries with the specified treatment will be excluded. Defaults to include all samples. Choices are 'minus', 'half', 'plus', and 'mixed'. 'mixed' filters out all individuals with libraires of mixed UDG treatments. Requires a valid anno file.", nargs='+', choices=['minus', 'half', 'plus', 'mixed'], default=[])
	Parser.add_argument('--published_only', help="Only include samples that have been published. Publication status is inferred based on the publication field in the anno file. Defaults to include all samples. Requires a valid anno file.", action='store_true')

	# Best representative sampleset arguments
	Parser.add_argument('--best_representative', help="Creates a list of best representative samples by taking highest coverage of any master ID.", action='store_true')
	Parser.add_argument('--preferred_data_source', help="Specify a preferred data source to use when generating a best representative sample set. If this option is not specified, the higest coverage sample will allways be used. If this option is specified, a sample with the specified data source will be preferred over other data sources, regardless of coverage. This option has no effect if the --best_representative flag is not set.", default="")
	Parser.add_argument('--best_representative_mode', help="Determines if best representative sample determination overrides anno filters or if it respects these filters, either applied before or after determination.", choices=['override', 'filter_union', 'next_best', 'next_best_strict'], default='override')

	# Custom sample exclusion/inclusion arguments
	Parser.add_argument('--custom_inclusion', help="Specify a text file with new-line delimited samples to include in the output. Overrides all filters.", type=str, default="")
	Parser.add_argument('--custom_exclusion', help ="Specify a text file with new-line delimited samples to exclude from the output. Overrides all filters.", type=str, default="")

	args = Parser.parse_args()

	# Set working directory
	# TODO Move this to just before driver file generation so we dont break relative filepaths in args
	# OR just resolve the filepaths before... but this would be a pain in the ass.
	chdir(args.working_dir)

	# Ensure that the preferred_data_source flag does not conflict with the data_exclude flag
	if args.preferred_data_source.lower() in list(map(str.lower, args.data_exclude)):
		raise argparse.ArgumentError("--data_exclude and --preferred_data_source flags conflict!")

	# Make sure the input geno, ind files exist
	if args.master[-4:] in ['.ind', '.snp']:
		master_stem = args.master[:-4]
	elif args.master[-5:] == '.geno':
		master_stem = args.master[:-5]
	else:
		master_stem = args.master
	master_geno = abspath(master_stem + '.geno')
	master_ind = abspath(master_stem + '.ind')
	if not(exists(master_geno)) or not(exists(master_ind)):
		raise FileNotFoundError("{} does not exist!".format(args.master))

	# Check if an anno file has been explicitly specified. If not, infer one based on the 'master' argument.
	if args.anno_file == None:
		anno_path = master_stem + '.anno'
	else:
		anno_path = args.anno_file		

	exclusion_set = set()

	# Run custom exclusion if a custom exclusion list is provided unless we're in best_representative/override mode.
	if args.custom_exclusion != "" and (args.best_representative == False or args.best_representative_mode != 'override'):
		list_filepath = abspath(args.custom_exclusion)
		if exists(list_filepath):
			exclusion_set.update(exclude_by_custom_list(list_filepath))
		else:
			raise FileNotFoundError("Custom exclusion list filepath invalid!")

	# Disable sample based filters if we don't have an anno file.
	if exists(anno_path):
		# Make sure it's the right anno file based on line count
		if get_file_len(anno_path, header=True) != get_file_len(master_ind, header=False):
			raise IndexError("Length of input ind file does not match that of the specified anno file.")
		anno_data = read_anno_file(anno_path)

		# Exclusion filters
		# These will not be applied in best representative override mode
		if args.best_representative == False or args.best_representative_mode != 'override':
			
			# Data type filter
			if len(args.data_exclude) > 0:
				exclusion_set.update(exclude_by_data_type(data_types=args.data_exclude, anno_data=anno_data))
			
			# Assessment threshold filter
			if args.assessment_threshold != "":
				exclusion_set.update(exclude_by_assessment(assessment_threshold=args.assessment_threshold, anno_data=anno_data))
			
			# UDG treatment filter
			if len(args.UDG_exclude) > 0:
				exclusion_set.update(exclude_by_udg(udg_treatments=args.UDG_exclude, anno_data=anno_data))
			
			# Exclude unpublished samples
			if args.published_only == True:
				exclusion_set.update(exclude_unpublished(anno_data=anno_data))
		
		# Define best representative sample set if --best_representative flag is set
		if args.best_representative:
			
			# If the override mode is set, exclude should be empty. If we're in filter_union mode, we'll take the union of the filtered samples and the non-best representative samples
			if args.best_representative_mode == 'override' or args.best_representative_mode == 'filter_union':
				exclusion_set.update(get_best_representatives(anno_data=anno_data, preferred_data_source=args.preferred_data_source))
			
			# If the next_best mode is set, also pass in the pre-filtered samples.
			elif args.best_representative_mode == 'next_best':
				exclusion_set = get_best_representatives(anno_data=anno_data, preferred_data_source=args.preferred_data_source, filtered_set=exclusion_set)

			# If the next_best_strict mode is set, pass in the pre-filtered samples and set the strict argument to True.
			elif args.best_representative_mode == 'nest_best_strict':
				exclusion_set = get_best_representatives(anno_data=anno_data, preferred_data_source=args.preferred_data_source, filtered_set=exclusion_set, strict=True)
	else:
		print('Anno file not found. Filters and/or best representative mode will be NOT be applied.')
	
	# Restore samples if a custom inclusion list has been specified. This will override all exclusion filters as well as best representative determinations.
	if args.custom_inclusion != "":
		list_filepath = abspath(args.custom_inclusion)
		if exists(list_filepath):
			exclusion_set = include_by_custom_list(list_filepath, exclusion_set)
		else:
			raise FileNotFoundError("Custom inclusion list filepath invalid!")

	# Build the driver file based on the exclusion set
	driver_file = build_driver_file(master_ind, exclusion_set, '{}.driver'.format(args.label))
	
	# Determine full list of snp sets to output to, including custom ones.
	snp_sets = args.snp_sets.extend(args.custom_list)
	for snp_set in args.snp_sets:
		# Pull standard snp set filepaths from the global...
		try:
			snp_file_path = SNP_SET_FILENAMES[snp_set]
		# ...and ensure custom snp set files exist and derive thier names.
		except:
			if exists(snp_set):
				snp_file_path = abspath(snp_set)
				snp_set = basename(snp_set)[:-4]
			else:
				raise FileNotFoundError("{} custom SNP set does not exist!".format(snp_set))
		# Make a directory for each output dataset by snp set
		try:
			mkdir(snp_set)
		except:
			pass
		output_name = "{}/{}_{}".format(snp_set, args.label, snp_set)

		par_file = build_par_file(master_stem=master_stem, snp_file_path=snp_file_path, driver_file_path=driver_file, output_name=output_name, par_template=PAR_DEFAULT)
		
		# Build command to run submission script or convertf directly if the forground flag is set and run.
		if args.foreground:
			bash_command = "{} -p {}".format(args.convertf_exe, par_file)
		else:
			submit = build_sbatch_script(output_name=output_name, convertf_execuable=args.convertf_exe, par_file=par_file, sbatch_template=args.sbatch_template)
			bash_command = "{} {}".format(SUBMIT_COMMAND, submit)
		subprocess.run(bash_command, shell=True)
