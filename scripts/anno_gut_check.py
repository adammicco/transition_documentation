import argparse
import numpy as np
import re
import jellyfish
from numpy.lib.arraysetops import unique


IGNORE_PATTERNS = [
	"lc",
	"contam",
	"published",
	"in.preparation",
	"o[0-9]+",
	"o[A-Z][a-z]+",
	"^o$"
	"brother.+",
	"sister.+",
	"mother.+",
	"father.+",
	"son.+",
	"daughter.+",
	"sibling.+",
	"parent.+",
	"dup.+",
	"duplicate.+",
	"twin.+",
	"halfbro.+",
	"halfbrother.+",
	"halfsis.+",
	"halfsister.",
	"halfsib.+",
	"halfsiblng.+",
	"in.preparation"
]

Parser = argparse.ArgumentParser(description="Run a series of checks on an anno file to ensure quality before release.")

Parser.add_argument('--anno_file', help="Specify path to a tsv formatted anno file.", required=True)

args = Parser.parse_args()

def read_anno_file(anno_file_path):
	with open(anno_file_path, "r") as f:
		header_line = f.readline()
		headers = header_line.split('\t')
		anno_file_data = {header : [] for header in headers}
		for line in f:
			fields = line.split('\t')
			for header in headers:
				anno_file_data[header].append(fields[headers.index(header)])
	return anno_file_data

def detect_cov_outliers(anno_data_set, extreme = False):
	version_ID = anno_data_set['Version ID']
	SNP_count = anno_data_set['SNPs hit on autosomal targets']
	autosomal_cov = anno_data_set['Coverage on autosomal targets']

	missing_SNP_count = set([index for index, value in enumerate(SNP_count) if (value == "" or value == ".." or value == "#N/A")])
	missing_autosomal_cov = set([index for index, value in enumerate(autosomal_cov) if (value == "" or value == ".." or value == "#N/A")])

	bad_indices = missing_SNP_count.union(missing_autosomal_cov)

	cleaned_version_ID = [value for index, value in enumerate(version_ID) if index not in bad_indices]
	cleaned_SNP_count = [int(value) for index, value in enumerate(SNP_count) if index not in bad_indices]
	cleaned_autosomal_cov = [float(value) for index, value in enumerate(autosomal_cov) if index not in bad_indices]

	data_array = np.array(list(zip(cleaned_SNP_count, cleaned_autosomal_cov)))

	covariance_matrix = np.cov(data_array, rowvar=False)
	if is_pos_def(covariance_matrix):
		inv_covariance_matrix = np.linalg.inv(covariance_matrix)
		if is_pos_def(inv_covariance_matrix):
			vars_mean = []
			for i in range(data_array.shape[0]):
				vars_mean.append(list(data_array.mean(axis=0)))
			diff = data_array - vars_mean
			md = []
			for i in range(len(diff)):
				md.append(np.sqrt(diff[i].dot(inv_covariance_matrix).dot(diff[i])))

			print("Covariance Matrix:\n {}\n".format(covariance_matrix))
			print("Inverse of Covariance Matrix:\n {}\n".format(inv_covariance_matrix))
			#print("Variables Mean Vector:\n {}\n".format(vars_mean))
			#print("Variables - Variables Mean Vector:\n {}\n".format(diff))
			print("Mahalanobis Distance:\n {}\n".format(md))
		else:
			print("Error: Inverse of Covariance Matrix is not positive definite!")
			pass
	else:
		print("Error: Covariance Matrix is not positive definite!")
		pass
	
	std = np.std(md)
	k = 3. * std if extreme else 2. * std
	m = np.mean(md)
	up_t = m + k
	low_t = m - k
	outliers = []
	for i in range(len(md)):
		if (md[i] >= up_t) or (md[i] <= low_t):
			outliers.append(i)  # index of the outlier
	outliers = {value : (cleaned_SNP_count[index], cleaned_autosomal_cov[index]) for index, value in enumerate(cleaned_version_ID) if index in outliers}
	return outliers

def is_pos_def(A):
	if np.allclose(A, A.T):
		try:
			np.linalg.cholesky(A)
			return True
		except np.linalg.LinAlgError:
			return False
	else:
		return False

def group_ID_set_match(set1, set2, similarity_threshold = 0.85, exact_match = False):
	# Jaro-Winkler similarity to determine pairwise similarity scores of groupID elements regardless of order
	for element1 in set1:
		for element2 in set2:
			if jellyfish.jaro_winkler_similarity(element1, element2) < similarity_threshold:
				return False
	return True

def find_base_ids_with_common_elements(query_tokenized_id, tokenized_ID_list):
	query_tokenized_id_set = set(query_tokenized_id)
	tokenized_group_ID_sets = map(set, tokenized_ID_list)
	matching_ID_ixs = []
	for ix, id in enumerate(tokenized_group_ID_sets):
		# This if statement should run if we essentially evaluate as true the Jaro-Winkler-ized version of `set1 == set2`.
		if len(query_tokenized_id_set) == len(id) and group_ID_set_match(query_tokenized_id_set, id, 0.85):
			matching_ID_ixs.append(ix)
	return matching_ID_ixs

def find_group_ID_typos(anno_data_set, ignore_patterns = IGNORE_PATTERNS):
	# Create parallel lists of version ID, group ID, tokenized group ID, base group ID, and base group ID occruance count
	version_ID = anno_data_set['Version ID']
	group_ID = anno_data_set['Group ID']
	
	# Implement a negative lookaround to exclude IGNORE_PATTERNS in the following filter
	group_ID_tags = re.compile("(?! {})".format("|".join(ignore_patterns)))

	tokenized_group_ID = [list(filter(group_ID_tags.match, re.split('_', ID_str.replace(".SG", "").replace(".DG", "")))) for ID_str in group_ID]
	
	base_group_ID = ["_".join(ID_tokens) for ID_tokens in tokenized_group_ID]
	base_group_ID_dict = {ID : [ix for ix, v in enumerate(base_group_ID) if v==ID] for ID in set(base_group_ID)}

	# TODO: refactor to reference base_group_IDs. May need to swap key and value in base_group_ID_dict to support lookup of base IDs by index...

	# Find instances where components of group IDs are swapped. We should generally accept the one with the higher higher entry count

	# putative_matches is a dictionary with keys representing indicies of "source" group IDs and values that are lists of indicies of matches
	putative_matches = {}
	for ix, id in enumerate(tokenized_group_ID):
		matching_ID_ixs = find_base_ids_with_common_elements(id, tokenized_group_ID)
		if len(matching_ID_ixs) > 0:
			putative_matches.update({ix : matching_ID_ixs})

	# Assemble results based on putative match indicies giving priority to the "base" ID with higher occurance count
	for source_id in putative_matches:
		if len(putative_matches[source_id]) > 0:
			pass

	pass

data = read_anno_file(args.anno_file)
outliers = detect_cov_outliers(data, extreme=True)
find_group_ID_typos(data)
print(outliers)