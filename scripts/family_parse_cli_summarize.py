import argparse
import networkx as nx
import re

import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt

RELATIONSHIP_DEGREE_MAP = {
	"dup" : 0,
	"dup/1d" : 0.5,
	"1d" : 1,
		"mother-son" : 1,
		"mother-daughter" : 1,
		"mother-child" : 1,
		"father-son" : 1,
		"father-daughter" : 1,
		"father-child" : 1,
		"parent-son" : 1,
		"parent-daughter" : 1,
		"parent-child" : 1,
		"son-mother" : 1,
		"daughter-mother" : 1,
		"child-mother" : 1,
		"son-father" : 1,
		"daughter-father" : 1,
		"child-father" : 1,
		"son-parent" : 1,
		"daughter-parent" : 1,
		"child-parent" : 1,
		"brothers" : 1,
		"brother" : 1,
		"brother-brother" : 1,
		"brother-sister" : 1,
		"sisters" : 1,
		"sister" : 1,
		"sister-sister" : 1,
		"sister-brother" : 1,
		"siblings" : 1,
		"sibling" : 1,
		"sibs" : 1,
		"sib" : 1,
		"sibling-sibling" : 1,
	"1/2d" : 1.5,
	"2d" : 2,
		"grandparent-grandchild" : 2,
		"grandparent-grandson" : 2,
		"grandparent-granddaughter" : 2,
		"grandfather-grandchild" : 2,
		"grandfather-grandson" : 2,
		"grandfather-granddaughter" : 2,
		"grandmother-grandchild" : 2,
		"grandmother-grandson" : 2,
		"grandmother-granddaughter" : 2,
		"grandchild-grandparent" : 2,
		"grandson-grandparent" : 2,
		"granddaughter-grandparent" : 2,
		"grandchild-grandfather" : 2,
		"grandson-grandfather" : 2,
		"granddaughter-grandfather" : 2,
		"grandchild-grandmother" : 2,
		"grandson-grandmother" : 2,
		"granddaughter-grandmother" : 2,
		"parsib-nibling" : 2,
		"parsib-nephew" : 2,
		"parsib-niece" : 2,
		"uncle-nibling" : 2,
		"uncle-nephew" : 2,
		"uncle-niece" : 2,
		"aunt-nibling" : 2,
		"aunt-nephew" : 2,
		"aunt-niece" : 2,
		"nibling-parsib" : 2,
		"nephew-parsib" : 2,
		"niece-parsib" : 2,
		"nibling-uncle" : 2,
		"nephew-uncle" : 2,
		"niece-uncle" : 2,
		"nibling-aunt" : 2,
		"nephew-aunt" : 2,
		"niece-aunt" : 2,
	"2/3d" : 2.5,
	"3d" : 3,
	"3d+" : 3.5,
	"4d" : 4,
	"5d" : 5,
	"6d" : 6
}

parser = argparse.ArgumentParser(description = "Builds a directed graph from pairwise relatedness information. This graph can be used to generate and plot subgraphs representing putative family trees for individuals in the dataset.")
parser.version = '1.0'
parser.add_argument("--csv", action='store', required=True)
parser.add_argument("--suppress_output", "-s", action='store_true')
parser.add_argument("--verbose_output", "-v", action='store_true')
parser.add_argument("--ind_summary", "-i", action='store_true')
parser.add_argument("--rel_summary", "-r", action='store_true')
parser.add_argument("--plot", "-p", action='store_true')

args = parser.parse_args()

csv_path = args.csv
SUPPRESS = args.suppress_output
PLOT = args.plot
VERBOSE = args.verbose_output
if not(VERBOSE):
	if args.ind_summary:
		SUMMARY_LEVEL = 'ind'
	elif args.rel_summary:
		SUMMARY_LEVEL = 'rel'
	else:
		SUMMARY_LEVEL = ''
		VERBOSE = True

# Reads in pairwise relationships in ind1	ind2	loc1	loc2	country1	country2	relationship format
# Yields: 7-tuple with above fields
def read_pairwise_relationships(csv_path):
	pair_set = set()
	with open(csv_path, 'r', encoding="utf8") as pairs:
		for line in pairs:
			feilds = line.split("\t")
			ind1 = feilds[0].strip()
			ind2 = feilds[1].strip()
			loc1 = feilds[2].strip()
			loc2 = feilds[3].strip()
			country1 = feilds[4].strip()
			country2 = feilds[5].strip()
			relationship = feilds[6].strip()
			order_unknown = bool(feilds[7].strip())
			if frozenset([ind1, ind2]) in pair_set:
				continue
			pair_set.add(frozenset([ind1, ind2]))
			yield (ind1, ind2, loc1, loc2, country1, country2, relationship, order_unknown)

# Constructs directed graph representing all family relationships in the loaded dataset.
# Returns:
#	G: DiGraph object representing family relationships in a directed graph
#	ind_set: Set of unique individuals represented in this graph
# 	dup_set: Set of individuals detected as duplicates. These are represented as alises on nodes in G.
# TODO It will be a bit of a PITA to handle arbitrary relative clusters I think but that could be a future enhancement. In the meantime, clusters should be decomposed into pairwise relationships.
def construct_graph(csv_path, locality_match_tolerance = 1):
	# Define graph
	G = nx.DiGraph()
	# Keep track of the unique individuals and duplicates
	ind_set = set()
	dup_set = set()
	# Get pairwise generator object and loop through reach relationship
	pairwise_relationships = read_pairwise_relationships(csv_path)
	for rel_line in pairwise_relationships:
		# If this is a duplicate relationship...
		if rel_line[6] == "dup":
			# Try to contract existing nodes
			try:
				G.nodes[rel_line[0]]
			except:
				for x,y in G.nodes(data=True):
					try:
						if rel_line[0] in y['aliases']:
							rel_line[0] = x
							break
					except:
						pass
			try:
				G.nodes[rel_line[1]]
			except:
				for x,y in G.nodes(data=True):
					try:
						if rel_line[1] in y['aliases']:
							rel_line[1] = x
							break
					except:
						pass
			try:
				nx.contracted_nodes(G, rel_line[0], rel_line[1])
				# If the nodes exists and we're able to contract them, either add an alias to the existing list in node attributes...
				try:
					G.nodes[rel_line[0]]['aliases'].append(rel_line[1])
					G.nodes[rel_line[1]]['aliases'].append(rel_line[0])
					try:
						G.nodes[rel_line[0]]['aliases'].extend(G.nodes[rel_line[1]]['aliases'])
						G.nodes[rel_line[1]]['aliases'].extend(G.nodes[rel_line[0]]['aliases'])
					except:
						pass
				# ...or define the attribue and add the alias to it.
				except:
					G.nodes[rel_line[0]].update({'aliases' : [rel_line[1]]})
					G.nodes[rel_line[1]].update({'aliases' : [rel_line[0]]})
					G.nodes[rel_line[0]]['aliases'].extend(G.nodes[rel_line[1]]['aliases'])
					G.nodes[rel_line[1]]['aliases'].extend(G.nodes[rel_line[0]]['aliases'])

				G = nx.contracted_nodes(G, rel_line[0], rel_line[1])
			# If we don't have nodes to contract, add the node for the ind we want to 'keep' and add the alias
			except:
				G.add_node(rel_line[0], aliases = [rel_line[1]])
			# Add the duplicates' IDs to the dup_set
			dup_set.add(rel_line[0])
			dup_set.add(rel_line[1])
		# If this is not a duplicate relationship...
		else:
			# Check if the locations match with respect to tolerance option
			loc0_element_set = set(re.split(", |\s", rel_line[2]))
			loc1_element_set = set(re.split(", |\s", rel_line[3]))
			shared_loc_elements = loc0_element_set.intersection(loc1_element_set)
			if len(shared_loc_elements) >= locality_match_tolerance:
				location = rel_line[2]
			# Call it 'cross-site' if the localities don't match
			else:
				location = "cross-site"
			# For non-duplicate relationships, add an edge to the graph connecting the inds and storing locality, country, and relationship info in the attributes.
			# Loop through nodes...
			source_node_0 = []
			source_node_1 = []
			for x,y in G.nodes(data=True):
					try:
						# ...to find nodes with an alias that matches our ind and append to source_node
						if rel_line[0] in y['aliases']:
							source_node_0.append(x)
						if rel_line[1] in y['aliases']:
							source_node_1.append(x)
					except:
						pass
			if len(source_node_0) == 0:
				rel0 = rel_line[0]
			else:
				rel0 = source_node_0[0]
			if len(source_node_1) == 0:
				rel1 = rel_line[1]
			else:
				rel1 = source_node_1[0]
			G.add_edge(rel0, rel1, relationship = rel_line[6], location = location, country = rel_line[4], order_unknown = rel_line[7])
			if "-" in rel_line[6]:
				relationship_peices = rel_line[6].split('-')
				relationship_peices.reverse()
				invert = '-'.join(relationship_peices)
			else:
				invert = rel_line[6]
			G.add_edge(rel1, rel0, relationship = invert, location = location, country = rel_line[4], order_unknown = rel_line[7])
				
			# Add the inds' IDs to the ind_set
			ind_set.add(rel_line[0])
			ind_set.add(rel_line[1])
	return G, ind_set, dup_set

# Traverse a family subgraph and return a generator that yeilds all involved edges using nx.bfs_edges() or nx.dfs_edges(). Handles referencing the source node by both it's own ID and aliases.
# Accepts:
#	G: The full dataset's graph
#	source_ind: An individual to use as the starting point of the family subgraph
#	search_type: Either "bfs" to do a breadth-first search or "dfs" to do a depth-first search (Default: "bfs")
# Returns:
#	family: A generator object that yeilds all edges in the family subgraph
#	nodes: A set of all unique nodes in the family subtree
#	master_node: The node used to build out the subgraph. This will euqal the source_ind unless we're referencing a node by an alias - as is the case with duplicates.
def get_family_generator(G, source_ind, search_type = "bfs"):
	if source_ind in G.nodes():
		master_node = source_ind
	else:
		for x,y in G.nodes(data=True):
			try:
				if source_ind in y['aliases']:
					master_node = x
					break
			except:
				pass
	if search_type == "bfs":
		family = nx.edge_bfs(G, source = master_node)
	elif search_type == "dfs":
		family = nx.edge_dfs(G, source = master_node)
	else:
		raise Exception('Invalid search type! The supported search types are "bfs" and "dfs".')
	
	return family, master_node

def get_family_members(family, master_node):
	nodes = [master_node] + [v for u,v in family]
	return set(nodes)

# Checks a node for aliases saved in the 'aliases' attribute and returns a list of the node's ID plus these aliases.
# If there are no aliases, a list with only the node's ID is returned.
# Accepts:
#	G: The dataset's full graph
#	node_ID: The node ID to check
# Returns:
#	IDs: A list containing the node's ID plus any aliases
def get_node_aliases(G, node_ID):
	try:
		IDs = [node_ID] + G.nodes[node_ID]['aliases']
	except:
		IDs = [node_ID]
	return IDs

# Setup dictionary to create family string summarizing relationships by individual, ordered by relationship degree.
# Accepts:
#   running_dict: Dictionary containing relationship information. Grows as the family subgraph is traversed.
#   addressed_IDs: List to keep track of IDs addressed. Grows as the family subgraph is traversed.
#   rel: Tuple defining pairwise relationship from family subgraph.
def summarize_rels_by_ind(rel, running_dict = {}, addressed_IDs = []):
	rel_edge = G.edges[rel]
	if rel_edge['order_unknown']:
		rel_order = " (order unknown)"
	else:
		rel_order = ""
	rel_relationship = rel_edge['relationship'] + rel_order

	rel0_IDs = get_node_aliases(G, rel[0])
	rel1_IDs = get_node_aliases(G, rel[1])
	rel0_str = "/".join(rel0_IDs)
	rel1_str = "/".join(rel1_IDs)
	addressed_IDs = addressed_IDs + rel0_IDs + rel1_IDs

	if rel0_str in running_dict.keys():
		if rel_relationship in running_dict[rel0_str].keys():
			running_dict[rel0_str][rel_relationship].append(rel1_str)
		else:
			running_dict[rel0_str].update({rel_relationship : [rel1_str]})
	else:
		running_dict.update({rel0_str : {rel_relationship : [rel1_str]}})

	return running_dict, addressed_IDs

# Setup dictionary to create family string summarizing relationships by degree, ordered by relationship degree.
# Accepts:
#   rel: Tuple defining pariwise relationship from family subgraph
#   running_dict: Dictionary containing relationship infromation. Grows as the family subgraph is traversed.
#   addressed_IDs: List to keep track of IDs addressed. Grows as the family subgraph is traversed.
def summarize_rels_by_degree(rel, running_dict = {}, addressed_IDs = []):
	rel_edge = G.edges[rel]
	if rel_edge['order_unknown']:
		rel_order = " (order unknown)"
	else:
		rel_order = ""
	rel_relationship = rel_edge['relationship'] + rel_order

	rel0_IDs = get_node_aliases(G, rel[0])
	rel1_IDs = get_node_aliases(G, rel[1])
	rel0_str = "/".join(rel0_IDs)
	rel1_str = "/".join(rel1_IDs)
	addressed_IDs = addressed_IDs + rel0_IDs + rel1_IDs

	if rel_relationship in running_dict.keys():
		running_dict[rel_relationship].append("{}-{}".format(rel0_str, rel1_str))
	else:
		running_dict.update({rel_relationship : ["{}-{}".format(rel0_str, rel1_str)]})
	
	return running_dict, addressed_IDs

def define_families(G, ind_set):
	# Keep track of the localities used in family names so we can iterate family name letters
	locality_list = []
	# Initialize done_list to track inds already addressed by way of being in an alredy traversed family subgraph -- We'll return this as a set at the end of this function
	done_list = []
	family_member_dict = {}
	for ind in ind_set:
		family, master_node = get_family_generator(G, ind, search_type = "bfs")
		if master_node not in done_list:
			nodes = get_family_members(family, master_node)
			num_members = len(nodes)
			# Initialize list to store relationship strings
			rel_string_list = list()
			# Initialize set to track if we end up with a cross-site family.
			location_set = set()
			# 'Reset' family generator
			family, master_node_check = get_family_generator(G, master_node)
			if master_node_check != master_node:
				raise Exception('Ran family generator twice using alias {} and recieved incorrect source node! Expected {}, Got {}.'.format(ind, master_node, master_node_check))
			addressed_IDs = []
			family_dict = {}
			reversed_rels = []
			for rel in family:
				reversed_rels.append(rel[::-1])
				if rel not in reversed_rels:
					if not(VERBOSE):
						if SUMMARY_LEVEL == 'ind':
							family_dict, addressed_IDs = summarize_rels_by_ind(rel, family_dict, addressed_IDs)
						elif SUMMARY_LEVEL == 'rel':
							family_dict, addressed_IDs = summarize_rels_by_degree(rel, family_dict, addressed_IDs)
					if VERBOSE:
					# Get relationship and determine appropriate substring syntax
						rel0_IDs = get_node_aliases(G, rel[0])
						rel1_IDs = get_node_aliases(G, rel[1])
						addressed_IDs = addressed_IDs + rel0_IDs + rel1_IDs
						rel_relationship = G.edges[rel]['relationship']
						if G.edges[rel]['order_unknown']:
							order = " (order unknown)"
						else:
							order = ""
						if re.match(r'\d', rel_relationship):
							rel_string = "{}-{} have a {} relationship".format("/".join(rel0_IDs), "/".join(rel1_IDs), rel_relationship)
						else:
							rel_string = "{}-{} are {}{}".format("/".join(rel0_IDs), "/".join(rel1_IDs), rel_relationship, order)
						rel_string_list.append(rel_string)
					rel_loc = G.edges[rel]['location']
					rel_country = G.edges[rel]['country']
					location_set.add(rel_loc)
			if not(VERBOSE):
				if SUMMARY_LEVEL == 'ind':
					for member in family_dict:
						rel_string_list.append("{}: {}".format(member, "; ".join(["{}{} {}".format(x, " rel" if bool(re.search(r'\d', x)) else "", ", ".join(family_dict[member][x])) for x in sorted(family_dict[member].keys(), key=lambda rel_deg: RELATIONSHIP_DEGREE_MAP[rel_deg.split(" ")[0]])])))
				elif SUMMARY_LEVEL == 'rel':
					for degree in family_dict.keys():
						rel_string_list.append("{}: {}".format(degree, ", ".join(family_dict[degree])))

			if len(location_set) > 1 or "cross-site" in location_set:
				locality = "{} cross-site".format(rel_country)
			else:
				try:
					locality = location_set.pop()
				except:
					raise(ValueError("empty location_set for family of ind: {}".format(master_node)))
			locality_list.append(locality)
			locality_occurrence_num = locality_list.count(locality)
			if locality_occurrence_num <= 26:
				family_letter = chr(ord('@') + locality_occurrence_num)
			else:
				family_letter ="{}{}".format(chr(ord('@')+(locality_occurrence_num//26)), chr(ord('A')+(locality_occurrence_num%26)))
			family_name = "{} Family {}".format(locality, family_letter)
			# family_member_dict.update({family_name : })
			if not(VERBOSE) and SUMMARY_LEVEL == 'rel':
				# print(rel_string_list)
				rel_string_list = sorted(rel_string_list, key=lambda x: RELATIONSHIP_DEGREE_MAP[x[:(x.find(':'))].split(" ")[0]])
			if not(VERBOSE):
				delimiter = ' | '
			else:
				delimiter = ', '
			family_string = "{} ({} members) ({})".format(family_name, num_members, delimiter.join(rel_string_list))

			if PLOT:
				#TODO improve plotting (https://towardsdatascience.com/customizing-networkx-graphs-f80b4e69bedf)
				family_subgraph = G.subgraph(nodes)
				family_subgraph_edges = family_subgraph.edges()
				weights = [1/RELATIONSHIP_DEGREE_MAP[family_subgraph[u][v]['relationship']]*3 for u,v in family_subgraph_edges]
				nx.draw(family_subgraph, with_labels = True, width = weights)
				plt.savefig("{} Family {}.pdf".format(locality, family_letter))
				print("{} Family {}.pdf saved.".format(locality, family_letter))

			for ID in set(addressed_IDs):
				done_list.append(ID)
				# if -s, --suppress_output set, do not print strings
				if not(SUPPRESS):
					print("{}\t{}".format(ID, family_string))
	return set(done_list)

G, ind_set, dup_set = construct_graph(csv_path)
define_families(G, ind_set)
