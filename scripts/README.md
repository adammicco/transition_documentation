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