# Scripts
## family_parse_cli_summarize.py
This script is what I had been using to generate family group calls from Iñigo's pair(ish)-wise kinship calls. It uses [networkx](https://networkx.org) (`nx`), a graph/network analysis python package, to handle parsing pairwise calls into edges/nodes of a graph. Each disjoin subgraph is treated as a family and is traversed to generate the family definition string that is saved into the anno file.

This will likely become less directly relevant as adna2 becomes more robust and the anno file moves towards being a dynamically generated view rather than a discrete document. However, the code implemented here could be repurposed to drive a graph-based representation of kinship relationships in adna2, perhaps via a package like [Apache AGE](https://age.apache.org), and ever integrating directly with kinship calling software in order to take some of the burden off of Iñigo.

### Implementation
The script works by reading in pairwise relatiponships from a tsv file in following format in order to build out the edges of the graph:

master ID 1 | master ID 2 | locality 1 | locality 2 | country 1 | country 2 | relationship | order_unknown (optional)
--- | --- | --- | --- | --- | --- | --- | --- 
I12345 | I54321 | Massachusetts, Cambridge, Central | Massachusetts, Cambridge, Inman | USA | USA | parent-child | 1
I51423 | I15324 | Massachusetts, Somerville, Davis | Massachusetts, Boston, Fenway | USA | USA | 2d

If the individuals are duplicates and nodes already exist for them, we will contract the nodes with `nx.contracted_nodes()` and record the IDs into an attribute on the node called `"aliases"`. If the nodes do not exist in the graph, a new node will be created with the first master ID and we will save off the other in the `"aliases"` attribute.

For non-duplicate relationships, the logic is simpler. We just create an edge -- this implicitly creates the nodes is they do not already exist -- and save as edge attributesd the locality, country, relationship, and order_unknown boolean. We also save off each ID into the `ind_set` for use in defining families.

We then iterate through each individual in the `ind_set` and, depending on the search type (we use a breadth-first search by default in the `define_families()` function), we create an `nx.edge_bfs()` or `nx.edge_dfs()` generator that is used to iterate through the edges of the disjoint subgraph that represents the family to which that individual belongs and we then can construct the family string for that family, assigning a name either based on the individuals' shared locality or country, in the case of cross-site families. Matching/comparison of locality names is done by breaking each locatity string into components (delimited by ', ') and making sure there are not more than `locality_match_tolerance` # of mismatching components, allowing for locality strings with grave # designations to not trigger the 'cross-site' moniker. For example, the first line in the example table above would not be a cross-site family, but the second line would be.

We also log the individuals in all already processed families in the `done_list` so that we do not repeat them as we keep iterating.

There is logic for printing relationships in verbose forms as well as summarisiung them by relationship type or by individual. I have also included very redimentray plotting logic but the results are very ugly right now and could stand to be substantially improved upon before this is actually useful.

### Shortcomings and Future Improvements

Right now, the script does not retain its graph database beteen runs or have tha ability to accept a gml graph database file as input to build off of or report from. It only takes a tsv of all pairwise relationships and constructs the database in memory with each run. We could resolve this by taking advantage of Posgres graph plugins and simply operating on a graph database in adna2 instead of an in-memory database. 

Because of the way that the graph is constructed each time and traversed, the family letters and relationship 'directions' (e.g. mother-son vs. son-mother) are not retained between runs. This can mean that, while the constituents of the families and their relationships are always the same, the family names and strings may chnage between runs. I have been manually cleaning this up as part of the kinship call process but it is cretainly not ideal. Adding support of tracking family names between runs could be a big improvement and should ideally be part of an integration with adna2.

Inputs must be pairwise. However, Iñigo often sends calls that reference relationships between clusters of individuals such as...

rel | inds
--- | --- 
siblings | IndA, IndB, IndC

...In this case I have been exploding these into pairwise relationships manually. Making the script able to handle these would be nice but it would be best to integrate it directly with the kinship caller.

This script does not currently have any support for assisting with pedigree construction. However, if this is of interest, the edge and node attribues could be very useful for storing the data required to fully define pedigree relatiponships. This could also be made more robust by integration with a graph database exttension in adna2.

Graphing functionality is curently very rudimentary and produces ugly results. This is something that could be improved upon and potentailly even expanded to assist with pedigree construction and/or programatic figure creation.

## PCA_clustering.R
This is a basic starting point for a project I started working on a long time ago to try to add rigor to the PCA-based pre-release curation/QA process.

Assuming that the process has not changed, we are currently relying on manual checks to determine if individuals cluster as expected in PCA space, both with respect to other individuals in their same population/group ID and to the background data. This is useful for QA, outlier detection, and subdivision of group IDs were applicable.

This script attempts to assist in detecting outliers and instances where subpopulations may be defined. It does not address comparison with the background PCA data.

### Implementation
To address the first problem of defining subpopulations within a specific group ID, I first attempt to calculate the ideal number of clusters that best fit the data in `calc_gap_stat_opt_clust` by using a gap statistic. My naïve approach here does not tune itself to better suit different populations that map more/less densely into PCA space.

Once the ideal number of clusters is determined, we then run k-means to assign individuals to a cluster in `run_kmeans`.

From there, `detect_outliers` handles calcuating euclidean distance between cluster centers and each individual and assigning a boolean outlier designation to each at various z-score thresholds (2, 1, 0.75, 0.5). These are then returned in a dataframe which then used for plotting clusters and outliers.

### Shortcomings and Future Improvements

I never got particularly far in this project but I think that it could still be a valuable improvement for the QA/curation process.

The most basic improvements could center around tuning the gap stat logic or implementing other methods (elbow method? hierarchical clustering?) for determining the ideal number of clusters to map individuals to.

Outlier detection is also done pretty naïvely here. During manual curation, it is normal to adjust outlier tolerace based on the expected spread of individuals within PCA space as some populations tend to map more densely and others less so. This script simply relies on analyzing individuals' PC_1 and PC_2 coordinates once the PCA mapping has already been performed externally and does not consider the PCA background data. This poses a problem in outlier detection across different levels of expected spread. I've thought about this a bit but have not had a chance to implement anything to attempt to solve it. I do not fully understand the complexities here, but I would expect that the differing spread across populations could possibly be accounted for by running clustering more PCs in higher dimentional space in order to normalize for spread across different populations/clusters when doing outlier detection.

Both for improved clustering/normalization and to streamline the pipeline, it could also be beneficial to integrate this with the PCA transform and plotting logic.

## rev_complement.py
This is a super simple script that simply takes in sequences and produces thier reverse complements. Self-explainatory and probably inferior to another tool. I don't remember why I wrote this!

## anno_gut_check.py
This is part of an unfinished project to implement various quality checks that would run on the anno file before a release. With the anno file likely moving to being a view extracted from adna2, this is likely not something that we would want to continue developing as a stand alone utility. However, a lot of the logic here could be reused and run against adna2 with some modification.

### Implementation
`detect_cov_outliers` takes in SNP count and coverage data (column names will need to be adjusted) in order to detect samples for which these are not correlated as expected as this is a huge red flag of some sort of sample scramble or data processing issue. It ignores any samples where there is missing or invalid data for SNP count, builds a covariance matrix and checks if it is positive definite with `is_pos_def()`. From there, it calculates the Mahalanobis distance, which takes into account the covariance between these two variables, to measure how far each data point deviates from the expected relationship between SNP count and coverage. Outliers are flagged based on how many standard deviations they fall away from the mean, with thresholds being adjustable through the extreme parameter, althrough this may need adjustment as we have different capture methods and pulldown parameters now.

`find_group_ID_typos()` is a work in progress for detecting typos in group IDs as these were (and likely still are) a huge pain point. It first toeknizes group IDs from the anno file, exluding patterns defined in `IGNORE_PATTERNS`. These are suffixes and infixes that I found to be common in the anno file (circa 2021ish) and should be excluded when comparing group IDs to eachother. It then looks for similarities in the tokenized group IDs to detect instances where components of group IDs may have been inadvertently swapped or misformatted and suggest correction to the version of each group ID that appears most frequently, based on the assumption that this is the one most likely to be correct.

There are two helper function that drive this:

1. `group_ID_set_match()` uses the Jaro-Winkler similarity between sets of tokens from two different group IDs to determine their similarity regardless of order. If this similarity is above `similarity_threshold`, the two tokenized group IDs are flagged as being possible matches, allowing for slight variations (swapped components) to be captured.

2. `find_base_ids_with_common_elements()` uses `group_ID_set_match()` to compare ach tokenized group ID to the others in the dataset, returning the indices of those that are determined to be above the specified Jaro-Winkler `similarity_threshold`.

### Shortcomings and Future Improvements
This is an unfinished script and will not work out of the box, especially given changes in anno file formatting, and methodology impacting SNP counts and coverage depth. That said, this could be useful as a juming off point for implementing validation checks in adna2 -- though this might be better implemented via a tool like [Great Expectations](https://docs.greatexpectations.io/docs/0.18/oss/get_started/get_started_with_gx_and_sql/) or similar tools.

In terms of methodology, the outliter detection could probably use some work and may even be able to be shifted to a simpler method.

Locality matching/validation could also be improved and even integrated with a mapping API in the future (see below).

Additonal checks could be easily performed on filepaths pointing to bams, pulldowns/logs, etc. in the anno file/adna2 to flag issues with dead pointers.

## lat_long_verify.py
In much the same vein as the above script, this script is a very basic starting point to vaildate latitude and longitude based on locality data.

### Implementation
This script was meant to use [geopy](https://geopy.readthedocs.io/en/stable/) to handle hitting a geocoding API, [nominatim](https://nominatim.org/release-docs/develop/). This uses Open Street Map data to attempt to geocode based on a query made up of the tokenized locality string components and the country.

### Shortcomings and Future Improvements

This did not go very well in my testing as many of the locality names do not appear in Open Street Map data. This will likely continue to be a barrier as many locality strings reference archaeological site and/or histroical features. Additionally, many names are non-standard Romanizations/Anglicanizations of names in other languages and this makes querying against a datbase mostly focused on modern points of interest difficult.

This may still be an idea worth pursuing in the future, but it will prove challenging.

## scratch_prepare.sh
This is a basic script that handles creating a folder on /n/sratch3 to be used for running the demultiplexing and analysis pipeline and copying over the fastq files given the flowcell and lane.

This is not really relevant anymore given the sequencing platform and data delivery changes.

## c_eig_reorder
This is a program that I wrote out of desperation and I think it is thankfully not really needed anymore but I will document it here in case it does come in handy at some point.

This program reorders .geno and .ind files (eigenstrat) to match a user-defined sample order specified in a "driver" file (a list of individual IDs). The output consists of reordered .geno and .ind files, along with a copy of the original .snp file. 

### Dependencies
* C++17

### Implementation
`read_geno_2_transpose` handles reading in the .geno file into its transpose, with individuals being on the rows and snps on the columns axis. `read_ind_file` then handles reading the .ind file.

The geno transpose and the ind vectors are then destructively reordered in place in `reorder_destructive` based on the order defined in the driver file. This saves a lot of overhead by mutating the vecotrs in place.

The now reordered geno vector is now writted back out to disk in its untransposed form, putting snps back on the rows and individuals on the columns, in `write_untransposed_geno` while `write_ind_file` handles writing the reordered ind file. The snp file is copied with no updates.

## c_snp_reorder
This program operates very similarly to c_eig_reorder but it handles reordering the snps in an eigenstrat geno file and its accompanying snp file.

### Dependencies
* C++17

### Implementation
`read_geno` handles reading the eigenstrat geno file into a vector while `read_snp_file` does this for the snp file. In this case, no transposing is needed.

The geno transpose and the snp vectors are then destructively reordered in place in `reorder_destructive` based on the order defined in the driver file. This saves a lot of overhead by mutating the vecotrs in place.

The now reordered geno and snp files are then written to disk with `write_geno` and `write_snp_file`, respectively, while the ind file is copied with out chnages.