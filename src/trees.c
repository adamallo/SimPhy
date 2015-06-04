/**
 *
 * \file trees.c
 * Source of the library of tree-manipulation functions (and tree types definitions).
 *******************************************************************************/

/**
 * \defgroup trees_priv Trees managing module -> PRIVATE FUNCTIONS
 *
 * Strict functions within trees managing module. For duplossim I´ve unlocked some
 * of them.
 *******************************************************************************/
///@{

#include "trees.h"

// **** Prototipes of private functions **** //

/** \name Recursive functions **/ 
///@{

// ** Node deletion ** //

/**
 * Frees a group of s_node with tree structure.
 *
 * It frees all the allocated memory of each s_node of a group of s_nodes,
 * including s_node::children.
 *
 * \param node
 *  s_node pointer (root of the group of s_nodes in the first call).
 *
 *******************************************************************************/
static void FreeSNodes (s_node * node);

/**
 * Frees a group of l_node with tree structure.
 *
 * It frees all the allocated memory of each l_node of a group of l_nodes,
 * including l_node::children and l_node::g_nodes (this could be avoid using the
 * free_gnodes tag).
 *
 * \param node
 *  l_node pointer (root of the group of l_nodes in the first call).
 * \param free_gnodes
 *  Logical flag. If =0 the g_node pointers (\ref l_node::g_nodes) are not freed.
 *
 *******************************************************************************/
static void FreeLNodes (l_node * node, int free_gnodes);

/**
 * Frees a group of g_node with tree structure.
 *
 * It frees all the allocated memory of each g_node of a group of g_node with
 * tree structure. It does not free the root memory if complete==0.
 *
 * \param node
 *  g_node pointer (root of the group of g_node in the first call).
 * \param complete
 *  Logical flag. If ==0 the root node is not freed.
 *******************************************************************************/
static void FreeGNodes (g_node * node, int complete);

// ** Node reorganitation ** //

/**
 * Reorders a group of s_node with tree structure in a post-order.
 *
 * It puts a consecutive index in each s_node following a post-order.
 *
 * \param node
 *  s_node (root in the first call of the function).
 * \param index
 *  Pointer where updating the index in each function call. Its value should be
 *  0 in the first call.
 *******************************************************************************/
static void PostReorderSNodes (s_node * node, int * index);

/**
 * Reorders a group of l_node with tree structure in a post-order.
 *
 * It puts a consecutive index in each l_node following a post-order.
 *
 * \param node
 *  l_node (root in the first call of the function).
 * \param index
 *  Pointer where updating the index in each function call. Its value should be
 *  0 in the first call.
 *******************************************************************************/
static void PostReorderLNodes (l_node * node, int * index);

/**
 * Reorders a group of g_node with tree structure in a post-order.
 *
 * It puts a consecutive index in each g_node following a post-order.
 *
 * \param node
 *  g_node (root in the first call of the function).
 * \param index
 *  Pointer where updating the index in each function call. Its value should be
 *  0 in the first call.
 *******************************************************************************/
static void PostReorderGNodes (g_node * node, int * index); 

/**
 * Reorders a group of s_node with tree structure in a pre-order.
 *
 * It puts a consecutive index in each s_node following a pre-order.
 *
 * \param node
 *  s_node (root in the first call of the function).
 * \param index
 *  Pointer where updating the index in each function call. Its value should be
 *  0 in the first call.
 *******************************************************************************/
static void PreReorderSNodes (s_node * node, int * index);

/**
 * Reorders a group of l_node with tree structure in a pre-order.
 *
 * It puts a consecutive index in each l_node following a pre-order.
 *
 * \param node
 *  l_node (root in the first call of the function).
 * \param index
 *  Pointer where updating the index in each function call. Its value should be
 *  0 in the first call.
 *******************************************************************************/
static void PreReorderLNodes (l_node * node, int * index);

/**
 * Reorders a group of g_node with tree structure in a pre-order.
 *
 * It puts a consecutive index in each g_node following a pre-order.
 *
 * \param node
 *  g_node (root in the first call of the function).
 * \param index
 *  Pointer where updating the index in each function call. Its value should be
 *  0 in the first call.
 *******************************************************************************/
static void PreReorderGNodes (g_node * node, int * index);

// ** Node copy ** //

/**
 * Copies a group of s_nodes (with tree structure) in an array of s_nodes following
 * a post-order.
 *
 * It clones a tree/subtree (group of \ref s_node "s_nodes" with tree structure) in an array
 * of \ref s_node "s_nodes". This is an important part of \ref CollapseSTree public function.
 *
 * \param output
 *  Pre-allocated memory of s_nodes (with enough memory for the number of
 *  \ref s_node "s_nodes" present in the input).
 * \param input
 *  Node of the original group of nodes, to be copied into the output (root in
 *  the first call).
 *
 * \attention It requires that the group of s_nodes have their index in a post-order.
 *******************************************************************************/
static void PostCollapseSNodes(s_node *output,s_node *input);

/**
 * Copies a group of s_nodes (with tree structure) in an array of l_nodes
 *
 * \param output
 *  Pre-allocated array of l_nodes (with enough memory for the number of
 *  \ref s_node "s_nodes" present in the input).
 * \param input
 *  Node of the original group of nodes, to be copied into the output (root in
 *  the first call).
 * \param l_id
 *  Pointer to set the l_nodes ids.
 *******************************************************************************/
static void CopyStoLNodes(l_node *output,s_node *input,int *l_id);

/**
 * Copies a group of l_nodes (with tree structure) in an array of l_nodes following
 * a post-order.
 *
 * It clones a tree/subtree (group of \ref l_node "l_nodes" with tree structure) in an array
 * of \ref l_node "l_nodes". This is an important part of \ref CollapseLTree public function.
 * The l_node::g_nodes are also copied if there is allocated memory for them.
 *
 * \param output
 *  Pre-allocated memory of l_nodes (with enough memory for the number of
 *  \ref l_node "l_nodes" present in the input.
 * \param input
 *  Node of the original group of nodes, to be copied into the output (root in
 *  the first call).
 * \param n_gleaves
 *  Number of leaves of the g_tree pointed by these \ref l_node "l_nodes". It is used to copy
 *  l_node::g_nodes.
 * \param retain_lateral
 *  Logical flag. If retain_lateral = 1, the \ref l_node::lat_node will not be set as NULL, retained an outdated value (DANGEROUS).
 * \param probs
 *  Logical flag. If probs = 1,l_tree::i_probs, l_tree::o_probs and ltree::i_combprobs will be copied with the rest of the tree.
 *
 * \attention It requires that the group of l_nodes have their index in a post-order.
 *******************************************************************************/
static void PostCollapseLNodes(l_node *output,l_node *input,int n_gleaves, int retain_lateral, int probs);

/**
 * Copies a group of g_nodes (with tree structure) in an array of g_nodes following
 * a post-order.
 *
 * It clones a tree/subtree (group of \ref g_node "g_nodes" with tree structure) in an array
 * of \ref g_node "g_nodes". This is an important part of \ref CollapseGTree public function.
 *
 * \param output
 *  Pre-allocated memory of g_nodes (with enough memory for the number of
 *  \ref g_node "g_nodes" present in the input).
 * \param input
 *  Node of the original group of nodes, to be copied into the output (root in
 *  the first call).
 *
 * \attention It requires that the group of g_nodes have their index in a post-order.
 *******************************************************************************/
static void PostCollapseGNodes(g_node *output,g_node *input);

/**
 * Copies a group of s_nodes (with tree structure) in an array of s_nodes following
 * a pre-order.
 *
 * It clones a tree/subtree (group of \ref s_node "s_nodes" with tree structure) in an array
 * of \ref s_node "s_nodes". This is an important part of \ref CollapseSTree public function.
 *
 * \param output
 *  Pre-allocated memory of s_nodes (with enough memory for the number of
 *  \ref s_node "s_nodes" present in the input).
 * \param input
 *  Node of the original group of nodes, to be copied into the output (root in
 *  the first call).
 *
 * \attention It requires that the group of s_nodes have their index in a pre-order.
 *******************************************************************************/
static void PreCollapseSNodes(s_node *output,s_node *input);

/**
 * Copies a group of l_nodes (with tree structure) in an array of l_nodes following
 * a pre-order.
 *
 * It clones a tree/subtree (group of \ref l_node "l_nodes" with tree structure) in an array
 * of \ref l_node "l_nodes". This is an important part of \ref CollapseLTree public function.
 * The l_node::g_nodes are also copied if there is allocated memory for them.
 *
 * \param output
 *  Pre-allocated memory of l_nodes (with enough memory for the number of
 *  \ref l_node "l_nodes" present in the input.
 * \param input
 *  Node of the original group of nodes, to be copied into the output (root in
 *  the first call).
 * \param n_gleaves
 *  Number of leaves of the g_tree pointed by these \ref l_node "l_nodes". It is used to copy
 *  l_node::g_nodes.
 * \param retain_lateral
 *  Logical flag. If retain_lateral = 1, the \ref l_node::lat_node will not be set as NULL, retained an outdated value (DANGEROUS).
 * \param probs
 *  Logical flag. If probs = 1,l_tree::i_probs, l_tree::o_probs and ltree::i_combprobs will be copied with the rest of the tree.
 *
 * \attention It requires that the group of l_nodes have their index in a pre-order.
 *******************************************************************************/
static void PreCollapseLNodes(l_node *output,l_node *input,int n_gleaves, int retain_lateral, int probs);

/**
 * Copies a group of g_nodes (with tree structure) in an array of g_nodes following
 * a pre-order.
 *
 * It clones a tree/subtree (group of \ref g_node "g_nodes" with tree structure) in an array
 * of \ref g_node "g_nodes". This is an important part of \ref CollapseGTree public function.
 *
 * \param output
 *  Pre-allocated memory of g_nodes (with enough memory for the number of
 *  \ref g_node "g_nodes" present in the input).
 * \param input
 *  Node of the original group of nodes, to be copied into the output (root in
 *  the first call).
 *
 * \attention It requires that the group of g_nodes have their index in a pre-order.
 *******************************************************************************/
static void PreCollapseGNodes(g_node *output,g_node *input);

// ** Node features obtention ** //

/**
 * Gets the number of duplications of a bunch of l_nodes following a post-order recursion.
 *
 * \param node
 *  l_node to analyze.
 * \return
 *  Number of duplications.
 *******************************************************************************/
static int Count_duplications_lnodes(l_node *node);


/**
 * Gets the number of transfers of a bunch of l_nodes following a post-order recursion.
 *
 * \param node
 *  l_node to analyze.
 * \return
 *  Number of transfers.
 *******************************************************************************/
static int Count_transfers_lnodes(l_node *node);

/**
 * Gets the number of gene conversions of a bunch of l_nodes following a post-order recursion.
 *
 * \param node
 *  l_node to analyze.
 * \return
 *  Number of gene conversions.
 *******************************************************************************/
static int Count_gc_lnodes(l_node *node);

/**
 * Gets the number of losses of a bunch of l_nodes following a post-order recursion.
 *
 * \param node
 *  l_node to analyze.
 * \return
 *  Number of losses.
 *******************************************************************************/
static int Count_losses_lnodes(l_node *node);

/**
 * Gets the n_gen of the tree in coalescent units (n_gen/Ne).
 *
 * \param node
 *  s_node to analyze.
 * \param g_Ne
 *  Global effective population size.
 * \return Tree n_gen in coalescent units.
 *******************************************************************************/
static long double Measure_s_node_cu_height(s_node *node, int g_Ne);

/**
 * Gets the n_gen of the tree in generations.
 *
 * \param node
 *  g_node to analyze.
 * \return Tree n_gen in number of generations.
 *******************************************************************************/
static long double Measure_s_node_gl_height(s_node *node);

/**
 * Gets the length of the tree in coalescent units (n_gen/Ne).
 *
 * \param node
 *  g_node to analyze.
 * \param g_Ne
 *  Global effective population size.
 * \return Tree length in coalescent units.
 *******************************************************************************/
static long double Measure_s_node_cu_length(s_node *node, int g_Ne);

/**
 * Gets the length of the tree in expected number of substitutions per site.
 *
 * \param node
 *  g_node to analyze.
 * \return Tree length in units of expected number of substitutions per site.
 *******************************************************************************/
static long double Measure_s_node_gl_length(s_node *node);

/**
 * Checks the ultrametricity of the tree in time units.
 *
 * \param node
 *  l_node to analyze.
 * \return Time if the tree is ultrametric, -1 if it is not ultrametric.
 *******************************************************************************/
static double CheckUltrametricityLNodes(l_node *node);

/**
 * Checks the ultrametricity of the tree in time units.
 *
 * \param node
 *  s_node to analyze.
 * \return Time if the tree is ultrametric, -1 if it is not ultrametric.
 *******************************************************************************/
static double CheckUltrametricitySNodes(s_node *node);

/**
 * Gets the n_gen of the tree in coalescent units (n_gen/Ne).
 *
 * \param node
 *  g_node to analyze.
 * \param g_Ne
 *  Global effective population size.
 * \return Tree n_gen in coalescent units.
 *******************************************************************************/
static long double Measure_g_node_cu_height(g_node *node, int g_Ne);

/**
 * Gets the n_gen of the tree in expected number of substitutions per site.
 *
 * \param node
 *  g_node to analyze.
 * \return Tree n_gen in units of expected number of substitutions per site.
 *******************************************************************************/
static long double Measure_g_node_bl_height(g_node *node);

/**
 * Gets the lenght of the tree in coalescent units (n_gen/Ne).
 *
 * \param node
 *  g_node to analyze.
 * \param g_Ne
 *  Global effective population size.
 * \return Tree length in coalescent units.
 *******************************************************************************/
static long double Measure_g_node_cu_length(g_node *node, int g_Ne);

/**
 * Gets the lenght of the tree in expected number of substitutions per site.
 *
 * \param node
 *  g_node to analyze.
 * \return Tree length in units of expected number of substitutions per site.
 *******************************************************************************/
static long double Measure_g_node_bl_length(g_node *node);


// ** Node modification ** //

/**
 * Refines a just-read-from-Newick s_node group with tree structure.
 *
 * This function puts the correct s_node::n_nodes and s_node::n_gen in each
 * s_node which has just been read from Newick.
 *
 * \param node
 *  s_node to refine (root in the first call).
 * \param ind_persp
 *  Number of individuals for nodes without a pre-specified s_node::n_replicas.
 * \param gen_time
 *  Generation time.
 *******************************************************************************/
static void RefineSNodes(s_node * node, int ind_persp, double gen_time);

/**
 * Refines a just-read-from-Newick l_node group with tree structure.
 *
 * This function puts the correct l_node::n_nodes, l_node::n_gen and allocates
 * l_node::g_nodes in each l_node which has just been read from Newick.
 *
 * \param node
 *  l_node to refine (root in the first call).
 * \param n_gleaves
 *  Number of gene tree leaves to allocate the l_node::g_nodes properly.
 * \param ind_persp
 *  Number of individuals for nodes without a pre-specified l_node::n_nodes.
 * \param gen_time
 *  Generation time.
 * \param n_dup
 * Pointer to get the number of duplications.
 * \param n_loss
 * Pointer to get the number of losses.
 * \param n_trans
 * Pointer to get the number of transferences.
 * \param n_gc
 * Pointer to get the number of gene conversions.
 *******************************************************************************/
static void RefineLNodes(l_node * node, int n_gleaves, int ind_persp, double gen_time, int *n_dup, int *n_loss, int *n_trans, int *n_gc);

//\cond DOXYGEN_EXCLUDE
///**
// * Deletes superfluous nodes due to losses.
// *
// * This function deletes superfluous nodes due to losses of a \ref l_tree after
// * their generation by a birth-death process.
// *
// * \param node
// *  Node to delete if it is a superfluous one.
// * \param root
// *  Point to the root of the tree, to change if it is deleted.
// * \param n_deletions
// *  Number of deleted nodes after the work of the function.
// * \param n_leaves
// *  Number of leaves after the work of the function.
// *******************************************************************************/
//\endcond
//static void CleanlossesLNodes(l_node * node, l_node ** root,int * n_deletions, int * n_leaves);

/**
 * Pre-order recursive addition of a constant number of generations and recalculation of times.
 *
 * \param node
 *  Species tree node to modify (root in the first recursion).
 * \param constant
 *  Number of generations to add.
 * \param gen_time
 *  Generation time
 * *******************************************************************************/
static void AddConstantNgenSNodes(s_node *node, double constant, double gen_time);

/**
 * Recursively (pre-order) updates the l_node::time of a bunch/tree of l_nodes.
 *
 * \param node
 *  Imput node.
 * \param gen_time
 *  Global generation time.
 *******************************************************************************/
static void TemporalizeLNodes(l_node *node, double gen_time);

// ** Multilocus Coalescent related functions ** //

/**
 * Calculates input and output probabilities for the different possible number of lineages performing a pre-order traversal.
 *
 * \param node
 *  Current locus tree node (root in the first call).
 * \param Ne
 *  Global population size.
 * \param in_suntree
 *  Logical flag. If =0 probs are not calculated for the current node.
 * \param verbosity
 *  Verbosity level.
 * \param includelosses
 *  Number of g_nodes associated to each lost l_node (0= normal behaviour, 1= simulate lost lineages).
 * \return maximum number of possible gene tree lineages going out from this node.
 *******************************************************************************/
static int CalcProbsNLineagesLTree(l_node *node, int Ne, int in_subtree, int verbosity);

/**
 * Sample a candidate number of output lineages for a given node using the precalculated probabilities,
 *
 * \param node
 *  Node of interest.
 * \param n_gleaves
 *  Number of gene tree leaves of the considered subtree (the full tree can be used as an upper-bound).
 * \param seed
 *  Random number generator's seed.
 *******************************************************************************/
static int SampleNLineagesLNode(l_node *node,int n_gleaves, gsl_rng *seed);

/**
 * Sample locus tree lineage counts based on the precalculated probabilities using a rejection sampling algorithm.
 *
 * \param node
 *  Locus tree node (root in the first call).
 * \param n_gleaves
 *  Number of gene tree leaves.
 * \param Ne
 *  Common population size.
 * \param verbosity
 *  Verbosity level.
 * \param seed
 *  Random number generator's seed.
 *******************************************************************************/
static void SampleNLineagesLTree(l_node *node, int n_gleaves, int Ne, int verbosity, gsl_rng *seed);

/**
 * Sample locus tree lineage counts based on the precalculated probabilities sampling the CDF of a given number of input lineages given the output count and the trees.
 *
 * \param node
 *  Locus tree node (root in the first call).
 * \param n_gleaves
 *  Number of gene tree leaves.
 * \param Ne
 *  Common population size.
 * \param verbosity
 *  Verbosity level.
 * \param seed
 *  Random number generator's seed.
 *******************************************************************************/
static void SampleNLineagesCDFLTree(l_node *node, int n_gleaves, int Ne, int verbosity, gsl_rng *seed);

/**
 * Calculates input and output probabilities for the different possible number of lineages.
 *
 * \param node
 *  l_tree (locus tree).
 * \param n_gleaves
 *  Number of gene tree leaves.
 * \param verbosity
 *  Verbosity level.
 * \param seed
 *  Random number generator's seed.
 *******************************************************************************/
static void SampleNLineagesLTree(l_node *node, int n_gleaves, int Ne, int verbosity, gsl_rng *seed);

// ** Node I/O ** //

/**
 * Writes a given group of s_nodes with tree structure in Newick format in
 * stdout, with branchs given in number of generations.
 *
 * \param root
 *  s_node to print (root in the first call).
 * \param names
 *  Names (name_c *).
 *******************************************************************************/
static void WriteSNodesGen (s_node * root, name_c * names);

/**
 * Writes a given group of s_nodes with tree structure in Newick format in
 * stdout, with branchs given in number of generations and with labels in internal nodes.
 *
 * \param root
 *  s_node to print (root in the first call).
 * \param names
 *  Names (name_c *).
 *******************************************************************************/
static void WriteSNodesGenIntlabel (s_node * root, name_c * names);

/**
 * Writes a given group of s_nodes with tree structure in Newick format in
 * stdout, with branchs given in time units.
 *
 * \param root
 *  s_node to print (root in the first call).
 * \param names
 *  Names (name_c *).
 * \param gen_time
 *  Generation time
 *******************************************************************************/
static void WriteSNodesTime (s_node * root, name_c * names, double gen_time);

/**
 * Writes a given group of s_nodes with tree structure in Newick format in
 * stdout, with branchs given in time units and with labels in internal nodes.
 *
 * \param root
 *  s_node to print (root in the first call).
 * \param names
 *  Names (name_c *).
 * \param gen_time
 *  Generation time
 *******************************************************************************/
static void WriteSNodesTimeIntlabel (s_node * root, name_c * names, double gen_time);

/**
 * Writes a given group of s_nodes with tree structure in Newick format in a
 * file, with branch lengths given in number of generations.
 *
 * \param file
 *  Output opened file.
 * \param root
 *  s_node to print (root in the first call).
 * \param names
 *  Names (name_c *).
 * \note The FILE * should be previously opened and checked to avoid errors.
 *******************************************************************************/
static void WriteSNodesFileGen (FILE * file, s_node * root, name_c * names);

/**
 * Writes a given group of s_nodes with tree structure in Newick format in a
 * file, with branch lengths given in number of generations and with labels in internal nodes.
 *
 * \param file
 *  Output opened file.
 * \param root
 *  s_node to print (root in the first call).
 * \param names
 *  Names (name_c *).
 * \note The FILE * should be previously opened and checked to avoid errors.
 *******************************************************************************/
static void WriteSNodesFileGenIntlabel (FILE * file, s_node * root, name_c * names);

/**
 * Writes a given group of s_nodes with tree structure in Newick format in a
 * file.
 *
 * \param file
 *  Output opened file.
 * \param root
 *  s_node to print (root in the first call).
 * \param names
 *  Names (name_c *).
 * \param gen_time
 *  Generation time
 * \note The FILE * should be previously opened and checked to avoid errors.
 *******************************************************************************/
static void WriteSNodesFileTime (FILE * file, s_node * root, name_c * names, double gen_time);

/**
 * Writes a given group of s_nodes with tree structure in Newick format in a
 * file and with labels in internal nodes.
 *
 * \param file
 *  Output opened file.
 * \param root
 *  s_node to print (root in the first call).
 * \param names
 *  Names (name_c *).
 * \param gen_time
 *  Generation time
 * \note The FILE * should be previously opened and checked to avoid errors.
 *******************************************************************************/
static void WriteSNodesFileTimeIntlabel (FILE * file, s_node * root, name_c * names, double gen_time);

/**
 * Writes a given group of l_nodes with tree structure in Newick format in
 * stdout, with branch length given in number of generations.
 *
 * \param root
 *  l_node to print (root in the first call).
 * \param names
 *  Names (name_c *).
 *******************************************************************************/
static void WriteLNodesGen (l_node * root, name_c * names);

/**
 * Writes a given group of l_nodes with tree structure in Newick format in
 * stdout, with branch length given in number of generations and labels in internal nodes.
 *
 * \param root
 *  l_node to print (root in the first call).
 * \param names
 *  Names (name_c *).
 *******************************************************************************/
static void WriteLNodesGenIntlabel (l_node * root, name_c * names);

/**
 * Writes a given group of l_nodes with tree structure in Newick format in
 * stdout, with branch lengths given in time.
 *
 * \param root
 *  l_node to print (root in the first call).
 * \param names
 *  Names (name_c *).
 * \param gen_time
 *  Generation time
 *******************************************************************************/
static void WriteLNodesTime (l_node * root, name_c * names, double gen_time);

/**
 * Writes a given group of l_nodes with tree structure in Newick format in
 * stdout, with branch lengths given in time.
 *
 * \param root
 *  l_node to print (root in the first call).
 * \param names
 *  Names (name_c *).
 * \param gen_time
 *  Generation time
 *******************************************************************************/
static void WriteLNodesTimeIntlabel (l_node * root, name_c * names, double gen_time);

/**
 * Writes a given group of l_nodes with tree structure in Newick format in a
 * file, with branch lengths given in number of generations.
 *
 * \param file
 *  Output opened file.
 * \param root
 *  l_node to print (root in the first call).
 * \param names
 *  Names (name_c *).
 * \note The FILE * should be previously opened and checked to avoid errors.
 *******************************************************************************/
static void WriteLNodesFileGen (FILE * file, l_node * root, name_c * names);

/**
 * Writes a given group of l_nodes with tree structure in Newick format in a
 * file, with branch lengths given in number of generations and labels in internal nodes.
 *
 * \param file
 *  Output opened file.
 * \param root
 *  l_node to print (root in the first call).
 * \param names
 *  Names (name_c *).
 * \note The FILE * should be previously opened and checked to avoid errors.
 *******************************************************************************/
static void WriteLNodesFileGenIntlabel (FILE * file, l_node * root, name_c * names);

/**
 * Writes a given group of l_nodes with tree structure in Newick format in a
 * file, with branch lengths given in time units.
 *
 * \param file
 *  Output opened file.
 * \param root
 *  l_node to print (root in the first call).
 * \param names
 *  Names (name_c *).
 * \param gen_time
 *  Generation time.
 * \note The FILE * should be previously opened and checked to avoid errors.
 *******************************************************************************/
static void WriteLNodesFileTime (FILE * file, l_node * root, name_c * names, double gen_time);

/**
 * Writes a given group of l_nodes with tree structure in Newick format in a
 * file, with branch lengths given in time units and internal nodes with labels
 *
 * \param file
 *  Output opened file.
 * \param root
 *  l_node to print (root in the first call).
 * \param names
 *  Names (name_c *).
 * \param gen_time
 *  Generation time.
 * \note The FILE * should be previously opened and checked to avoid errors.
 *******************************************************************************/
static void WriteLNodesFileTimeIntlabel (FILE * file, l_node * root, name_c * names, double gen_time);

/**
 * Writes a comma-separated list of daughter lineages in a
 * file.
 *
 * \param file
 *  Output opened file.
 * \param p
 *  l_node to treat (root in the first call).
 * \param names
 *  Names (name_c *).
 * \note The FILE * should be previously opened and checked to avoid errors.
 *******************************************************************************/
static void WriteDaughtersNodesFile(FILE * file,l_node * p, name_c *names);

/**
 * Writes a given group of g_nodes with tree structure in Newick format in stdout.
 *
 * \param root
 *  g_node to print (root in the first call).
 * \param names
 *  Names (name_c *).
 *******************************************************************************/
static void WriteGNodes (g_node * root, name_c * names);

/**
 * Writes a given group of g_nodes with tree structure in Newick format in stdout and internal nodes with labels
 *
 * \param root
 *  g_node to print (root in the first call).
 * \param names
 *  Names (name_c *).
 *******************************************************************************/
static void WriteGNodesIntlabel (g_node * root, name_c * names);

/**
 * Writes a given group of g_nodes with tree structure in Newick format in a
 * string.
 *
 * \param str
 *  String where write the Newick tree.
 * \param root
 *  g_node to print (root in the first call).
 * \param names
 *  Names (name_c *).
 *******************************************************************************/
static void WriteGNodesStr (char * str, g_node * root, name_c * names);

/**
 * Writes a given group of g_nodes with tree structure in Newick format in a
 * string and internal nodes with labels.
 *
 * \param str
 *  String where write the Newick tree.
 * \param root
 *  g_node to print (root in the first call).
 * \param names
 *  Names (name_c *).
 *******************************************************************************/
static void WriteGNodesStrIntlabel (char * str, g_node * root, name_c * names);

/**
 * Writes a given group of g_nodes with tree structure in Newick format in a
 * file.
 *
 * \param file
 *  Output opened file.
 * \param root
 *  g_node to print (root in the first call).
 * \param names
 *  Names (name_c *).
 * \note The FILE * should be previously opened and checked to avoid errors.
 *******************************************************************************/
static void WriteGNodesFile (FILE * file, g_node * root, name_c * names);

/**
 * Writes a given group of g_nodes with tree structure in Newick format in a
 * file and internal nodes with labels.
 *
 * \param file
 *  Output opened file.
 * \param root
 *  g_node to print (root in the first call).
 * \param names
 *  Names (name_c *).
 * \note The FILE * should be previously opened and checked to avoid errors.
 *******************************************************************************/
static void WriteGNodesFileIntlabel (FILE * file, g_node * root, name_c * names);
///@}

/** \name Non-recursive functions **/
///@{
/**
 * Relink the l_node::g_nodes pointers with a newly allocated g_tree memory.
 *
 * \param wl_tree
 *  l_tree to relink
 * \param m_node
 *  New g_tree memory to relink.
 * \note The original gene tree (pointed by l_node::g_nodes) have to be alive, and
 *  both g_trees with the same indices.
 * \return \ref NO_ERROR on OK or an \ref ERRORS "error code" if any error
 *  ocurrs.
 *******************************************************************************/
static long int RelinkLGTrees(l_tree * wl_tree, g_node * m_node);

/**
 * Samples one node from an array, with a probability inversely proportional to the MRCA time between them and one "reference" node.
 *
 * \param l_pointers
 *  Pointer to the array of l_node pointers of candidate nodes.
 * \param n_nodes
 *  Number of candidate nodes (array length)
 * \param t_node
 *  Reference node
 * \param u_num
 *  Previously sampled probability from an U(0-1).
 * \param verbosity
 *  Verbosity level.
 * \return Selected node pointer.
 *******************************************************************************/
static l_node * ChooseLNodePeriod(l_node **l_pointers, int n_nodes, l_node * t_node, double u_num, int verbosity);

/**
 * Returns the position of the next no-blank character in a given string.
 *
 * \param string
 *  char pointer to the string to analyze.
 * \return Position of the first no-blank character.
 *******************************************************************************/
static int firstnoblank(char *string);

static long int GetSnodeParamsFromNexusComments(char *string,s_node *node,int *n_char);

static long int GetLnodeParamsFromNexusComments(char *string,l_node *node,int *n_char);

///@}
///@}
// **** Declarations of public functions **** //

/**
 * \name Public functions
 * Better organized in \ref trees.h
 *******************************************************************************/
///@{

// *** Node memory manage *** //

// ** S_Node creation ** //
s_node * NewSNodes(int n_nodes, int max_children)
{
    s_node * nodes=NULL, *w_node=NULL;
    int i=0,j=0;
    
    // ***
    /// <dl><dt> Function structure </dt><dd>
    
    // **
    /// \ref s_node "S_node" memory allocation
    nodes=calloc(n_nodes,sizeof(s_node));
    ErrorReporter((int long)nodes, NULL);
    
    // **
    /// <dl><dt>\ref s_node "S_node" loop</dt><dd>
    for (i=0; i<n_nodes; ++i)
    {
        // *
        /// Variable initialization
        w_node=nodes+i;
        w_node->index=i;
        /*w_node->sp_index=0;
         w_node->n_lodes=0;
         w_node->n_replicas=0;
         w_node->n_child=0;
         w_node->Ne=0;*/ // Implicit in calloc
        w_node->anc_node=NULL;
        w_node->l_nodes=NULL;
        w_node->gen_length=0.0;
        w_node->n_gen=0.0;
        w_node->time=0;
        w_node->mu_mult=1;
        w_node->gtime_mult=1;
        
        // *
        /// s_node::children memory allocation and initialization</dd></dl></dd></dl>
        w_node->children=calloc(max_children,sizeof(s_node *));
        ErrorReporter((long int) w_node->children, NULL);
        
        for (j=0; j<max_children; ++j)
        {
            *(w_node->children+j)=NULL;
        }
        
    }
    
    return (nodes);
}

// ** L_Node creation ** //
l_node * NewLNodes(int n_nodes, int n_gleaves, int max_children)
{
    l_node * nodes=NULL, *w_node=NULL;
    int i=0,j=0;
    
    // ***
    /// <dl><dt> Function structure </dt><dd>
    
    // **
    /// \ref l_node "L_node" memory allocation
    nodes=calloc(n_nodes,sizeof(l_node));
    ErrorReporter((int long)nodes, NULL);
    
    // **
    /// <dl><dt>\ref l_node "L_node" loop</dt><dd>
    for (i=0; i<n_nodes; ++i)
    {
        // *
        /// Variable initialization
        w_node=nodes+i;
        w_node->index=i;
        /*w_node->sp_index=0;
         w_node->n_nodes=0;
         w_node->n_child=0;
         w_node->Ne=0;
         w_node->kind_node=SP;
         w_node->paralog=0;
        w_node->n_ilin=0;
        w_node->n_olin=0;
        w_node->fmax_nlin=0;*/ // Implicit in calloc
        w_node->anc_node=NULL;
        w_node->lat_node=NULL;
        w_node->conts=NULL;
        w_node->i_probs=NULL;
        w_node->i_combprobs=NULL;
        w_node->o_probs=NULL;
        w_node->gen_length=0.0;
        w_node->n_gen=0.0;
        w_node->time=0.0;
        w_node->mu_mult=1;
        w_node->gtime_mult=1;
        
        // *
        /// l_node::children memory allocation and initialization
        w_node->children=calloc(max_children,sizeof(l_node *));
        ErrorReporter((long int) w_node->children, NULL);
        
        for (j=0; j<max_children; ++j)
        {
            *(w_node->children+j)=NULL;
        }
        
        // *
        /// l_node::g_nodes allocation and initialization (if required)</dd></dl></dd></dl>
        
        if (n_gleaves!=0)
        {
            w_node->g_nodes=calloc(n_gleaves,sizeof(g_node *));
            ErrorReporter((long int)w_node->g_nodes,NULL);
            
            for (j=0; j<n_gleaves; ++j)
            {
                *(w_node->g_nodes+j)=NULL;
            }
        }
        else
            w_node->g_nodes=NULL;
        
        
    }
    
    return (nodes);
}

// ** G_Node creation ** //

g_node * NewGNodes(int n_nodes, int max_children)
{
    g_node * nodes=NULL, *w_node=NULL;
    int i=0,j=0;
    
    // ***
    /// <dl><dt> Function structure </dt><dd>
    
    // **
    /// \ref g_node "G_node" memory allocation
    if (n_nodes==0)
        return (NULL);
    
    nodes=calloc(n_nodes,sizeof(g_node));
    ErrorReporter((int long)nodes,NULL);
    
    // **
    /// <dl><dt>\ref g_node "g_node" loop</dt><dd>
    for (i=0; i<n_nodes; ++i)
    {
        // *
        /// Variable initialization </dd></dl>
        
        w_node=nodes+i;
        w_node->index=i;
        /*w_node->sp_index=0;
         w_node->replica=0;
         w_node->paralog=0;
         w_node->n_child=0 */ // Implicit in calloc
        w_node->children=NULL;
        w_node->anc_node=NULL;
        w_node->contl=NULL;
        w_node->conts=NULL;
        w_node->n_gen=0.0;
        w_node->bl=0.0;
        w_node->gen_length=0.0;
        
        // *
        /// g_node::children memory allocation and initialization</dd></dl>
        w_node->children=calloc(max_children,sizeof(g_node *));
        ErrorReporter((long int) w_node->children,NULL);
        
        for (j=0; j<max_children; ++j)
        {
            *(w_node->children+j)=NULL;
        }
    }
    
    return (nodes);
}

period * NewPeriods(int n_periods, int max_nodes)
{
    period *periods=NULL, *w_period;
    int i=0;
    
    periods=calloc(n_periods, sizeof(struct period));
    ErrorReporter((long int)periods,NULL);
    
    for (i=0; i<n_periods; ++i)
    {
        w_period=periods+i;
        w_period->r_bound=0;
        w_period->l_nodes=calloc(max_nodes,sizeof(l_node *));
        ErrorReporter((long int)w_period->l_nodes,NULL);
        w_period->n_lnodes=0;
    }
    return periods;
}

long double NNexusTrees(FILE *ifile, int *n_trees)
{
    int i=0,in_trees=0,j=0;
    size_t LENGTH=IO_BUFFER;
    char *buffer=NULL, *tok=NULL; //BEL
    
    // ******
    /// <dl><dt> Function structure </dt><dd>
    
    *n_trees=0;
    reallocBuffer(&buffer, &LENGTH, LENGTH);
    ResetBuffer(buffer, LENGTH);
    rewind(ifile);
    
    if (fgets(buffer,(int)LENGTH,ifile)==NULL || strcmp(buffer,"#NEXUS\n")!=0)
    {
        fprintf(stderr,"Improper Nexus trees input file\n");
        return SETTINGS_ERROR;
        
    }
    
    // *****
    /// <dl><dt> Line's loop </dt><dd></dd></dl></dd></dl>
    while (fgets(buffer,(int)LENGTH,ifile)!=NULL && i<=MAX_IT)
    {
        if (in_trees==1)
        {
            tok=strtok(buffer,"\n");
            if (tok!=NULL && (strcmp(tok,"End;")==0 || strcmp(tok,"END;")==0 || strcmp(tok,"end;")==0))
                break;
            else
            {
                if(strcmp(strtok((buffer+firstnoblank(buffer))," "),"tree")==0)
                    ++*n_trees;
                else
                {
                    fprintf(stderr,"Improper Nexus trees input file: Error in line %s\n",buffer);
                    return (IO_ERROR);
                }
                
            }
        }
        else
        {
            for (j=0; j<16;++j) // BEGIN TREES length
            {
                if (*(buffer+j)=='\0')
                    break;
                *(buffer+j)=toupper(*(buffer+j));
            }
            
            if (strcmp(buffer,"BEGIN TREES;\n")==0)
                in_trees=1;
        }
        
        if (*(buffer+LENGTH-2)!=TEST_CHAR) // It reached the buffer's limit
        {
            
#ifdef DBG
            fprintf(stderr,"\n\tWARNING, the buffer used reading the parameters file is not enough. Using a bigger size\n");
            fflush(stderr);
#endif
            reallocBuffer(&buffer,&LENGTH, LENGTH*10);
            ResetBuffer(buffer, LENGTH);
            rewind(ifile);
            *n_trees=0;
            in_trees=0;
            continue;
        }
        if (ferror(ifile)!=0 || feof(ifile)!=0)
        {
            fprintf(stderr,"Improper Nexus trees input file: Error in line %s\n",buffer);
            return (SETTINGS_ERROR);
        }

        ++i;
    }
    if (i>MAX_IT) return LOOP_ERROR; // Avoids ininite loops
    
    free(buffer);
    return NO_ERROR;
}

long double InitNexusParser(FILE *ifile)
{
    int i=0,j=0;
    size_t LENGTH=IO_BUFFER;
    char *buffer=NULL; //BEL
    
    // ******
    /// <dl><dt> Function structure </dt><dd>
    
    reallocBuffer(&buffer, &LENGTH, LENGTH);
    ResetBuffer(buffer, LENGTH);
    
    rewind(ifile);
    
    // *****
    /// <dl><dt> Line's loop </dt><dd></dd></dl></dd></dl>
    while (fgets(buffer,(int)LENGTH,ifile)!=NULL && i<=MAX_IT)
    {
        if (ferror(ifile)!=0 || feof(ifile)!=0)
        {
            fprintf(stderr,"Error in line: %s\n",buffer);
            return (IO_ERROR);
        }
        
        for (j=0; j<LENGTH-1;++j)
        {
            if (*(buffer+j)=='\0')
                break;
            *(buffer+j)=toupper(*(buffer+j));
        }
        
        if (strcmp(buffer,"BEGIN TREES;\n")==0)
        {
            free(buffer);
            return NO_ERROR;
        }
        
        if (*(buffer+LENGTH-2)!=TEST_CHAR) // It reached the buffer's limit
        {
            
#ifdef DBG
            fprintf(stderr,"\n\tWARNING, the buffer used reading the parameters file is not enough. Using a bigger size\n");
            fflush(stderr);
#endif
            reallocBuffer(&buffer, &LENGTH,LENGTH*10);
            ResetBuffer(buffer, LENGTH);
            rewind(ifile);
            continue;
        }
        
        ++i;
    }
    
    free(buffer);
    return IO_ERROR;
}

long double NextNexusTree(FILE *ifile,char **tree_str)
{
    int i=1,j;
    size_t LENGTH=IO_BUFFER;
    size_t BLENGTH=IO_BUFFER;
    char *buffer=NULL, *big_buffer=NULL, *p=NULL;
    
    reallocBuffer(&big_buffer, &BLENGTH, BLENGTH);
    reallocBuffer(&buffer, &LENGTH, LENGTH);
    ResetBuffer(buffer, LENGTH);
                
    if(fgets(buffer,(int)LENGTH,ifile)==NULL || ferror(ifile)!=0 || feof(ifile)!=0)
    {
        fprintf(stderr,"Error getting next NEXUS tree, line: %s\n",buffer);
        return (IO_ERROR);
    }
    
    while(*(buffer+LENGTH-2)!=TEST_CHAR && i<=MAX_IT)
    {
        ++i;
        
#ifdef DBG
        fprintf(stderr,"\n\tWARNING, the buffer used reading the parameters file is not enough. Concatenating more reads\n");
        fflush(stderr);
#endif
        reallocBuffer(&big_buffer, &BLENGTH, IO_BUFFER*i);
        strncat(big_buffer,buffer,LENGTH);
        ResetBuffer(buffer, LENGTH);
        if(fgets(buffer,(int)LENGTH,ifile)==NULL || ferror(ifile)!=0 || feof(ifile)!=0)
        {
            fprintf(stderr,"Error getting next NEXUS tree, line: %s\n",buffer);
            return (IO_ERROR);
        }
    }
    if (i>=MAX_IT)
    {
        free(buffer);
        free(big_buffer);
        return LOOP_ERROR;
    }
    strncat(big_buffer,buffer,LENGTH);
    *tree_str=realloc(*tree_str,LENGTH*i*sizeof(char));
    if (*tree_str==NULL)
        return MEM_ERROR;
    
    j=0;
    p=(big_buffer);
    while(*p!='(' && j<=MAX_IT && j<i*LENGTH)
    {
        p=(big_buffer+j);
        ++j;
    }
    
    strncpy(*tree_str,big_buffer+j-1,LENGTH*i-j);
    *tree_str=strtok(*tree_str,"\n");
    free(buffer);
    free(big_buffer);
    if (j>=MAX_IT || j>i*LENGTH)
    {
        fprintf(stderr,"Error parsing Nexus tree, no parenthesis present!!!\n");
        return LOOP_ERROR;
    }
    
    return NO_ERROR;
}

s_tree * ParseNexusSTree (char * nexus,name_c **names_ptr, int verbosity, double gen_time, int Ne, double mu, int ind_persp)
{
    s_node *current_node=NULL, *anc_node=NULL,*root=NULL;
    s_tree *tree=NULL;
    name_c * names=NULL;
    char code=' ';
    int step=0, in_coment=0;
    int n_char=0, iteration=0, ffree_codename=1, n_leaves=0, n_gleaves=0, n_inodes=0, n_nodes=0, max_children=0, max_lname=0,index=0;
    char *buffer=NULL,*name_buffer=NULL;
    size_t tbuffer=NUM_BUFFER;
    
    // ****
    /// <dl><dt> Function structure </dt><dd>
    
    // ***
    /// Test of the nexus tree string by \ref ChecknexusSTree
    
    ErrorReporter(CheckNexusTree(nexus),": Error in the Nexus species tree file\n");
    
    // ***
    /// Error control
    
    if (gen_time==0)
        ErrorReporter(UNEXPECTED_VALUE,": Inproper generation time");
    
    // ***
    /// Buffer allocation and initialization
    reallocBuffer(&buffer, &tbuffer, tbuffer);
    ResetBuffer(buffer,tbuffer);
    name_buffer=calloc(MAX_NAME,sizeof(char));
    // ***
    /// First read of the Nexus tree. Obtains the number of internal nodes (")"), s_tree::n_leaves ("(" or "," not followed by "(")).
    while (*(nexus+step)!=';')
    {
        switch (in_coment)
        {
            case 0:
                if((*(nexus+step)=='(' ||*(nexus+step)==',' )&&*(nexus+step+1)!='(')
                    ++n_leaves;
                else if(*(nexus+step)==')')
                    ++n_inodes;
                else if(*(nexus+step)=='[')
                {
                    in_coment=1;
                    ++step;
                }
                break;
            case 1:
                switch (*(nexus+step))
            {
                case ']':
                    in_coment=0;
                    break;
            }
                break;
                
            default:
                ErrorReporter(UNEXPECTED_VALUE,NULL);
                break;
        }
        ++step;
        if (step==MAX_IT) ErrorReporter(LOOP_ERROR, ": parsing a Nexus species tree"); // Avoids ininite loops
    }
    n_nodes=n_leaves+n_inodes;
    
    // ***
    /// Name container allocation by \ref NewNames
    *names_ptr=NewNames(n_leaves,0);
    names=*names_ptr;
    
    // ***
    /// <dl><dt>Second read of the nexus tree (main loop)</dt><dd>
    
    step=0; //Reset
    code=*(nexus+step); //First char
    
    while (code!=';') //End of the nexus tree
    {
        switch (code) {
            case '(':
                // **
                /// <dl><dt>New s_node (code="(").</dt><dd>
                if (iteration!=0) //New normal node
                {
                    anc_node=current_node;
                    // *
                    /// New s_node by \ref NewSNodes
                    current_node=NewSNodes(1,max_children);
                    // *
                    /// Points the pointers of this new node and its ancestor
                    *(anc_node->children+anc_node->n_child)=current_node;
                    ++anc_node->n_child;
                    current_node->anc_node=anc_node;
                    // *
                    /// Searches for the maximum number of children in the tree.</dd></dl>
                    if(max_children<anc_node->n_child)
                    {
                        max_children=anc_node->n_child;
                    }
                    
                }
                else //New root node
                {
                    current_node=NewSNodes(1,max_children);
                    root=current_node;
                }
                
                break;
            case ',':
            case ')':
                // **
                /// Goes down one tree node (code=", or )")
                current_node=current_node->anc_node;
                anc_node=current_node->anc_node;
                break;
            case ':':
                // **
                /// <dl><dt>New branch length (code=":").</dt><dd>
                
                // *
                /// Reads all the following integers and . in a buffer
                ResetBuffer(buffer,tbuffer);
                n_char=0;
                ++step;
                code=*(nexus+step);
                while (code<58 && (code>47 || code==46) && n_char<tbuffer && n_char<MAX_IT)
                {
                    *(buffer+n_char)=code;
                    ++n_char;
                    ++step;
                    code=*(nexus+step);
                }
                if (n_char==MAX_IT) ErrorReporter(LOOP_ERROR, ": parsing a branch length"); // Avoids ininite loops
                else if (n_char==tbuffer)
                {
                    step-=n_char+2;
                    reallocBuffer(&buffer, &tbuffer, tbuffer*10);
                }
                else
                {
                    --step;
                    *(buffer+n_char)=0; //Sets the end of the string to avoid problems in the sscanf
                    // *
                    /// Translates the buffer in a float number and asigns it as s_node::gen_length converted in number of generations of the current node</dd></dl>
                    sscanf(buffer,"%lf",&current_node->gen_length);
                }

                break;
            case '[':
                step+=2;
                ErrorReporter(GetSnodeParamsFromNexusComments(nexus,current_node,&step), ": Error parsing Nexus tree comments as branch-specific parameters");
                if (current_node->n_replicas>1 && current_node->n_child!=0)
                {
                    PrintXCharError(nexus, step, "\nNEWICK PARSING ERROR\n", "|<- Fixing the number of nodes/replicates in an internal node is not allowed, please, revisit the manual and your input trees\n");
                    ErrorReporter(SETTINGS_ERROR, NULL);
                }
                if (current_node->n_replicas>1)
                    n_gleaves+=current_node->n_replicas-ind_persp;
                break;
                
            default:
                // **
                /// <dl><dt>New leaf(code!= former ones).</dt><dd>
                
                // *
                /// Reads the name of the leaf in a buffer
                ResetBuffer(name_buffer, MAX_NAME);
                n_char=0;
                while (code!='(' && code!=')' && code!= ',' && code!= ';' && code!= ':' && code!= '[' && n_char<MAX_NAME && n_char<MAX_IT)
                {
                    *(name_buffer+n_char)=code;
                    ++n_char;
                    ++step;
                    code=*(nexus+step);
                }
                if (n_char==MAX_IT) ErrorReporter(LOOP_ERROR,": parsing leaf information"); // Avoids ininite loops
                else if (n_char==MAX_NAME)
                {
                    fprintf(stderr,"Species names are being truncated, since they are bigger than the buffer. You could change this behaviour increasing the buffer size by changing the environment variable SIMPHY_MAXNAME\n"); //\todo Change the name_c implementation to avoid this behaviour.
#ifdef DEBUG
                    fflush(stderr);
#endif
                    while (code!='(' && code!=')' && code!= ',' && code!= ';' && code!= ':' && code!= '[')
                    {
                        ++step;
                        code=*(nexus+step);
                    }
                    if (n_char==MAX_IT) ErrorReporter(LOOP_ERROR,": parsing leaf information"); // Avoids ininite loops
                }
                
                *(name_buffer+n_char)=0; //Sets the end of the string to avoid problems in the sscanf
                
                // *
                /// Searches for the maximum number of characters in a name.
                if (max_lname<n_char)
                    max_lname=n_char;
                
                // *
                /// Asigns the new name to the name_c::names
                strncpy(names->names+(ffree_codename*names->max_lname), name_buffer, names->max_lname);
                --step;
                
                // *
                /// New s_node by \ref NewSNodes
                anc_node=current_node;
                current_node=NewSNodes(1,max_children);
                n_gleaves+=ind_persp;
                // *
                /// Points the pointers of this new node and its ancestor
                *(anc_node->children+anc_node->n_child)=current_node;
                ++anc_node->n_child;
                current_node->anc_node=anc_node;
                // *
                /// Asigns the new name code (position) to the current node (s_node::sp_index)
                current_node->sp_index=ffree_codename;
                ++ffree_codename;
                // *
                /// Searches for the maximum number of children in the tree.</dd></dl></dd></dl>
                if(max_children<anc_node->n_child)
                {
                    max_children=anc_node->n_child;
                }
                
                break;
        }
        
        ++step;
        code=*(nexus+step);
        ++iteration;
        if (iteration==MAX_IT) ErrorReporter(LOOP_ERROR, NULL); // Avoids ininite loops
        
    }
    
    // ***
    /// Refines the readed nodes by \ref RefineSNodes
    
    RefineSNodes(root,ind_persp,gen_time);
    max_lname++; //One extra character for \0.
    
    // ***
    /// Reallocates the name_c memory by \ref ReallocNames
    ReallocNames(names,max_lname);
    
    // ***
    /// Allocates memory for the readed tree and completes it.
    
    tree=NewSTree(0, n_leaves, n_gleaves, max_children,gen_time,Ne,mu); //Nodes have already been allocated.
    tree->root=root;
    tree->n_nodes=n_nodes;
    tree->locus_tree=NULL;
    tree->gene_tree=NULL;
    
    // ***
    /// Fills the node indexes following a post_order.
    PostReorderSNodes(tree->root,&index);
    
    // ***
    /// Checks ultrametricity when measuring the species tree in time units.
    
    if(CheckUltrametricitySTree(tree)==-1)
        ErrorReporter(UNEXPECTED_VALUE,"The species tree is not ultrametric in time units. This doesn't allow to use the full model (transfers)\n");
    
    // ***
    /// Frees dynamic memory (buffers)</dd></dl>
    free(buffer);
    free(name_buffer);
    buffer=NULL;
    
    if (verbosity>2)
    {
        printf("\n\t%d-node species tree correctly built",(n_leaves*2)-1);
        if (verbosity>3)
        {
            printf(": ");
            WriteSNodesGen(root,names);
        }
        printf("\n");
        
    }
    
    return (tree);
}

s_tree * ParseNewickSTree (char * newick,name_c **names_ptr, int verbosity, double gen_time, int Ne, double mu, int ind_persp)
{
    s_node *current_node=NULL, *anc_node=NULL,*root=NULL;
    s_tree *tree=NULL;
    name_c * names=NULL;
    char code=' ';
    int step=0;
    int n_char=0, iteration=0, ffree_codename=1, n_leaves=0, n_gleaves=0, n_inodes=0, n_nodes=0, max_children=0, max_lname=0, n_replica=0,index=0, n_priv_ngleaves=0;
    char *buffer=NULL,*name_buffer=NULL;
    size_t tbuffer=NUM_BUFFER;
    
    /// \attention Changes on buffer managing applied like in Nexus parsers but without posterior testting.
    
    // ****
    /// <dl><dt> Function structure </dt><dd>
    
    // ***
    /// Test of the newick tree string by \ref CheckNewickSTree
    
    ErrorReporter(CheckNewickSTree(newick),": Error in the Newick species tree\n");
    
    // ***
    /// Error control
    
    if (gen_time==0)
        ErrorReporter(UNEXPECTED_VALUE, ": Inproper generation time");
    
    // ***
    /// Buffer allocation and initialization
    
    reallocBuffer(&buffer,&tbuffer,tbuffer);
    ResetBuffer(buffer, tbuffer);
    name_buffer=calloc(MAX_NAME,sizeof(char));
    
    // ***
    /// First read of the Newick tree. Obtains the number of internal nodes (")"), s_tree::n_leaves ("(" or "," not followed by "(") and s_tree::n_gleaves (s_tree::n_leaves + (replicas "/" -1 for each node)).
    while (*(newick+step)!=';')
    {
        if((*(newick+step)=='(' ||*(newick+step)==',' )&&*(newick+step+1)!='(') ++n_leaves;
        if(*(newick+step)==')') ++n_inodes;
        if(*(newick+step)=='/')
        {
            // * Reseting variables * //
            ResetBuffer(buffer,tbuffer);
            n_replica=0;
            n_char=0;
            ++step;
            code=*(newick+step);
            // * Reading replicas digits * //
            while (code<58 && code>47 && n_char<tbuffer && n_char<MAX_IT) // Int numbers or maxbuffer
            {
                *(buffer+n_char)=code;
                ++n_char;
                ++step;
                code=*(newick+step);
            }
            if (n_char==MAX_IT) ErrorReporter(LOOP_ERROR,""); // Avoids ininite loops
            else if (n_char==tbuffer)
            {
                step-=n_char+2;
                reallocBuffer(&buffer, &tbuffer, tbuffer*10);
            }
            else
            {
                --step;
                // * Translating * //
                sscanf(buffer,"%ui",&n_replica);
                n_gleaves+=n_replica;
                ++n_priv_ngleaves;
            }

        }
        ++step;
        if (step==MAX_IT) ErrorReporter(LOOP_ERROR, ": parsing a Newick tree"); // Avoids ininite loops
    }
    n_nodes=n_leaves+n_inodes;
    n_gleaves+=(n_leaves-n_priv_ngleaves)*ind_persp;//Addition of the ngleaves of the nodes using the common number of individuals (ind_persp)
    
    // ***
    /// Name container allocation by \ref NewNames
    *names_ptr=NewNames(n_leaves,0);
    names=*names_ptr;
    
    // ***
    /// <dl><dt>Second read of the Newick tree (main loop)</dt><dd>
    
    step=0; //Reset
    code=*(newick+step); //First char
    
    while (code!=';') //End of the Newick tree
    {
        switch (code) {
            case '(':
                // **
                /// <dl><dt>New s_node (code="(").</dt><dd>
                if (iteration!=0) //New normal node
                {
                    anc_node=current_node;
                    // *
                    /// New s_node by \ref NewSNodes
                    current_node=NewSNodes(1,max_children);
                    // *
                    /// Points the pointers of this new node and its ancestor
                    *(anc_node->children+anc_node->n_child)=current_node;
                    ++anc_node->n_child;
                    current_node->anc_node=anc_node;
                    // *
                    /// Searches for the maximum number of children in the tree.</dd></dl>
                    if(max_children<anc_node->n_child)
                    {
                        max_children=anc_node->n_child;
                    }
                    
                }
                else //New root node
                {
                    current_node=NewSNodes(1,max_children);
                    root=current_node;
                }
                
                break;
            case ',':
            case ')':
                // **
                /// Goes down one tree node (code=", or )")
                current_node=current_node->anc_node;
                anc_node=current_node->anc_node;
                break;
            case ':':
                // **
                /// <dl><dt>New branch length (code=":").</dt><dd>
                
                // *
                /// Reads all the following integers and . in a buffer
                ResetBuffer(buffer, tbuffer);
                n_char=0;
                ++step;
                code=*(newick+step);
                while (code<58 && (code>47 || code==46) && n_char<tbuffer && n_char<MAX_IT)
                {
                    *(buffer+n_char)=code;
                    ++n_char;
                    ++step;
                    code=*(newick+step);
                }
                if (n_char==MAX_IT) ErrorReporter(LOOP_ERROR, ": parsing a branch length"); // Avoids ininite loops
                else if (n_char==tbuffer)
                {
                    step-=n_char+2;
                    reallocBuffer(&buffer, &tbuffer, tbuffer*10);
                }
                else
                {
                    --step;
                    *(buffer+n_char)=0; //Sets the end of the string to avoid problems in the sscanf
                    // *
                    /// Translates the buffer in a float number and asigns it as s_node::gen_length converted in number of generations of the current node</dd></dl>
                    sscanf(buffer,"%lf",&current_node->gen_length);
                }
                
                break;
            case '*':
                // **
                /// <dl><dt>Lineage specific substitution rate multi (code="*").</dt><dd>
                
                // *
                /// Reads all the following integers and . in a buffer
                ResetBuffer(buffer, tbuffer);
                n_char=0;
                ++step;
                code=*(newick+step);
                while (code<58 && (code>47 || code==46)&& n_char<tbuffer && n_char<MAX_IT)
                {
                    *(buffer+n_char)=code;
                    ++n_char;
                    ++step;
                    code=*(newick+step);
                }
                if (n_char==MAX_IT) ErrorReporter(LOOP_ERROR,": parsing a branch specific substitution rate"); // Avoids ininite loops
                else if (n_char==tbuffer)
                {
                    step-=n_char+2;
                    reallocBuffer(&buffer, &tbuffer, tbuffer*10);
                }
                else
                {
                    --step;
                    *(buffer+n_char)=0; //Sets the end of the string to avoid problems in the sscanf
                    // *
                    /// Translates the buffer in a float number and asigns it as s_node::mu_mult </dd></dl>
                    sscanf(buffer,"%lf",&current_node->mu_mult);
                }
                
                break;
            case '~':
                // **
                /// <dl><dt>Lineage specific generation time multi (code="~").</dt><dd>
                
                // *
                /// Reads all the following integers and . in a buffer
                ResetBuffer(buffer, tbuffer);
                n_char=0;
                ++step;
                code=*(newick+step);
                while (code<58 && (code>47 || code==46) && n_char<tbuffer && n_char<MAX_IT)
                {
                    *(buffer+n_char)=code;
                    ++n_char;
                    ++step;
                    code=*(newick+step);
                }
                if (n_char==MAX_IT) ErrorReporter(LOOP_ERROR,": parsing a branch specific generation time"); // Avoids ininite loops
                else if (n_char==tbuffer)
                {
                    step-=n_char+2;
                    reallocBuffer(&buffer, &tbuffer, tbuffer*10);
                }
                else
                {
                    --step;
                    *(buffer+n_char)=0; //Sets the end of the string to avoid problems in the sscanf
                    // *
                    /// Translates the buffer in a float number and asigns it as s_node::gtime_mult </dd></dl>
                    sscanf(buffer,"%lf",&current_node->gtime_mult);
                }
                break;
            case '#':
                // **
                /// <dl><dt>New effective population size (Ne) (code="#").</dt><dd>
                
                // *
                /// Reads all the following integers in a buffer
                ResetBuffer(buffer, tbuffer);
                n_char=0;
                ++step;
                code=*(newick+step);
                while (code<58 && code>47 && n_char<tbuffer && n_char<MAX_IT)
                {
                    *(buffer+n_char)=code;
                    ++n_char;
                    ++step;
                    code=*(newick+step);
                }
                if (n_char==MAX_IT) ErrorReporter(LOOP_ERROR,":  parsing a branch specific population size"); // Avoids ininite loops
                else if (n_char==tbuffer)
                {
                    step-=n_char+2;
                    reallocBuffer(&buffer, &tbuffer, tbuffer*10);
                }
                else
                {
                    --step;
                    *(buffer+n_char)=0; //Sets the end of the string to avoid problems in the sscanf
                    // *
                    /// Translates the buffer in a integer number and asigns it as s_node::Ne of the current node</dd></dl>
                    sscanf(buffer,"%ui",&current_node->Ne);
                    if (current_node->Ne==0)
                        ErrorReporter(SETTINGS_ERROR,": Null branch specific population size is not allowed");
                }
                break;
            case '/':
                // **
                /// <dl><dt>New number of replicates (code="/").</dt><dd>
                
                // *
                /// Reads all the following integers in a buffer
                ResetBuffer(buffer, tbuffer);
                n_replica=0;
                n_char=0;
                ++step;
                code=*(newick+step);
                while (code<58 && code>47 && n_char<tbuffer && n_char<MAX_IT)
                {
                    *(buffer+n_char)=code;
                    ++n_char;
                    ++step;
                    code=*(newick+step);
                }
                if (n_char==MAX_IT) ErrorReporter(LOOP_ERROR,": parsing a branch specific number of replicates"); // Avoids ininite loops
                else if (n_char==tbuffer)
                {
                    step-=n_char+2;
                    reallocBuffer(&buffer, &tbuffer, tbuffer*10);
                }
                else
                {
                    --step;
                    *(buffer+n_char)=0; //Sets the end of the string to avoid problems in the sscanf
                    // *
                    /// Translates the buffer in a integer number and asigns it as s_node::n_nodes of the current node</dd></dl>
                    sscanf(buffer,"%ui",&n_replica);
                    if (n_replica==0)
                        ErrorReporter(SETTINGS_ERROR,": Null branch specific number of replicates is not allowed");
                    current_node->n_replicas=n_replica;
                }
                break;
                
            default:
                // **
                /// <dl><dt>New leaf(code!= former ones).</dt><dd>
                
                // *
                /// Reads the name of the leaf in a buffer
                ResetBuffer(name_buffer, MAX_NAME);
                n_char=0;
                while (code!='(' && code!=')' && code!= ',' && code!= ';' && code!= ':' && code != '~' && n_char<MAX_NAME && n_char<MAX_IT)
                {
                    *(name_buffer+n_char)=code;
                    ++n_char;
                    ++step;
                    code=*(newick+step);
                }
                if (n_char==MAX_IT) ErrorReporter(LOOP_ERROR,": parsing leaf information"); // Avoids ininite loops
                else if (n_char==MAX_NAME)
                {
                    fprintf(stderr,"Species names are being truncated, since they are bigger than the buffer. You could change this behaviour increasing the buffer size by changing the environment variable SIMPHY_MAXNAME\n"); //\todo Change the name_c implementation to avoid this behaviour.
#ifdef DEBUG
                    fflush(stderr);
#endif
                    while (code!='(' && code!=')' && code!= ',' && code!= ';' && code!= ':' && code!= '[' && n_char<MAX_IT)
                    {
                        ++step;
                        code=*(newick+step);
                    }
                    if (n_char==MAX_IT) ErrorReporter(LOOP_ERROR,": parsing leaf information"); // Avoids ininite loops
                }
                *(name_buffer+n_char)=0; //Sets the end of the string to avoid problems in the sscanf
                
                // *
                /// Searches for the maximum number of characters in a name.
                if (max_lname<n_char)
                    max_lname=n_char;
                
                // *
                /// Asigns the new name to the name_c::names
                strncpy(names->names+(ffree_codename*names->max_lname), name_buffer, names->max_lname);
                --step;
                
                // *
                /// New s_node by \ref NewSNodes
                anc_node=current_node;
                current_node=NewSNodes(1,max_children);
                // *
                /// Points the pointers of this new node and its ancestor
                *(anc_node->children+anc_node->n_child)=current_node;
                ++anc_node->n_child;
                current_node->anc_node=anc_node;
                // *
                /// Asigns the new name code (position) to the current node (s_node::sp_index)
                current_node->sp_index=ffree_codename;
                ++ffree_codename;
                // *
                /// Searches for the maximum number of children in the tree.</dd></dl></dd></dl>
                if(max_children<anc_node->n_child)
                {
                    max_children=anc_node->n_child;
                }
                
                break;
        }
        
        ++step;
        code=*(newick+step);
        ++iteration;
        if (iteration==MAX_IT) ErrorReporter(LOOP_ERROR,NULL); // Avoids ininite loops
        
    }
    
    // ***
    /// Refines the readed nodes by \ref RefineSNodes
    
    RefineSNodes(root,ind_persp,gen_time);
    max_lname++; //One extra character for \0.
    
    // ***
    /// Reallocates the name_c memory by \ref ReallocNames
    ReallocNames(names,max_lname);
    
    // ***
    /// Allocates memory for the readed tree and completes it.
    
    tree=NewSTree(0, n_leaves, n_gleaves, max_children,gen_time,Ne,mu); //Nodes have already been allocated.
    tree->root=root;
    tree->n_nodes=n_nodes;
    tree->locus_tree=NULL;
    tree->gene_tree=NULL;
    
    // ***
    /// Fills the node indexes following a post_order.
    PostReorderSNodes(tree->root,&index);
    
    // ***
    /// Checks ultrametricity when measuring the species tree in time units.
    
    if(CheckUltrametricitySTree(tree)==-1)
        ErrorReporter(UNEXPECTED_VALUE,"The species tree is not ultrametric in time units. This doesn't allow to use the full model (transfers)\n");
    // ***
    /// Frees dynamic memory (buffers)</dd></dl>
    free(buffer);
    free(name_buffer);
    buffer=NULL;
    
    if (verbosity>2)
    {
        printf("\n\t%d-node species tree correctly built",(n_leaves*2)-1);
        if (verbosity>3)
        {
            printf(": ");
            WriteSNodesGen(root,names);
        }
        printf("\n");
        
    }
    
    return (tree);
}

l_tree * ParseNexusLTree (char * nexus,name_c **names_ptr, int verbosity, double gen_time, int Ne, double mu, int ind_persp, int *n_dup, int *n_loss, int *n_trans, int *n_gc)
{
    l_node *current_node=NULL, *anc_node=NULL,*root=NULL;
    l_tree *tree=NULL;
    name_c * names=NULL;
    char code=' ';
    int step=0, in_coment=0;
    int n_char=0, iteration=0, ffree_codename=1, n_leaves=0, n_tleaves=0,n_gleaves=0, n_priv_ngleaves=0, n_inodes=0, n_nodes=0, max_children=0, max_lname=0,index=0;
    char *buffer=NULL,*name_buffer=NULL;
    size_t tbuffer=NUM_BUFFER;
    
    // ****
    /// <dl><dt> Function structure </dt><dd>
    
    // ***
    /// Test of the nexus tree string by \ref ChecknexusSTree
    
    ErrorReporter(CheckNexusTree(nexus),": Nexus locus tree error");
    
    // ***
    /// Error control
    
    if (gen_time==0)
        ErrorReporter(UNEXPECTED_VALUE, ": Inproper generation time");
    
    // ***
    /// Buffer allocation and initialization
    
    reallocBuffer(&buffer,&tbuffer,tbuffer);
    ResetBuffer(buffer,tbuffer);
    name_buffer=calloc(MAX_NAME,sizeof(char));
    
    // ***
    /// First read of the Nexus tree. Obtains the number of internal nodes (")"), s_tree::n_leaves ("(" or "," not followed by "(")).
    while (*(nexus+step)!=';')
    {
        switch (in_coment)
        {
            case 0:
                if((*(nexus+step)=='(' ||*(nexus+step)==',' )&&*(nexus+step+1)!='(')
                    ++n_leaves;
                else if(*(nexus+step)==')')
                    ++n_inodes;
                else if(*(nexus+step)=='[')
                {
                    in_coment=1;
                    ++step;
                }
                break;
            case 1:
                switch (*(nexus+step))
            {
                case ']':
                    in_coment=0;
                    break;
            }
                break;
                
            default:
                ErrorReporter(UNEXPECTED_VALUE,": parsing Nexus locus tree");
                break;
        }
        ++step;
        if (step==MAX_IT) ErrorReporter(LOOP_ERROR,NULL); // Avoids ininite loops
    }
    n_nodes=n_leaves+n_inodes;
    n_tleaves=n_leaves;
    
    // ***
    /// Name container allocation by \ref NewNames
    *names_ptr=NewNames(n_leaves,0);
    names=*names_ptr;
    
    // ***
    /// <dl><dt>Second read of the nexus tree (main loop)</dt><dd>
    
    step=0; //Reset
    code=*(nexus+step); //First char
    
    while (code!=';') //End of the nexus tree
    {
        switch (code) {
            case '(':
                // **
                /// <dl><dt>New s_node (code="(").</dt><dd>
                if (iteration!=0) //New normal node
                {
                    anc_node=current_node;
                    // *
                    /// New s_node by \ref NewSNodes
                    current_node=NewLNodes(1,0,max_children);
                    // *
                    /// Points the pointers of this new node and its ancestor
                    *(anc_node->children+anc_node->n_child)=current_node;
                    ++anc_node->n_child;
                    current_node->anc_node=anc_node;
                    // *
                    /// Searches for the maximum number of children in the tree.</dd></dl>
                    if(max_children<anc_node->n_child)
                    {
                        max_children=anc_node->n_child;
                    }
                    
                }
                else //New root node
                {
                    current_node=NewLNodes(1,0,max_children);
                    root=current_node;
                }
                
                break;
            case ',':
            case ')':
                // **
                /// Goes down one tree node (code=", or )")
                current_node=current_node->anc_node;
                anc_node=current_node->anc_node;
                break;
            case ':':
                // **
                /// <dl><dt>New branch length (code=":").</dt><dd>
                
                // *
                /// Reads all the following integers and . in a buffer
                ResetBuffer(buffer, tbuffer);
                n_char=0;
                ++step;
                code=*(nexus+step);
                while (code<58 && (code>47 || code==46) && n_char<tbuffer && n_char<MAX_IT)
                {
                    *(buffer+n_char)=code;
                    ++n_char;
                    ++step;
                    code=*(nexus+step);
                }
                if (n_char==MAX_IT) ErrorReporter(LOOP_ERROR,": parsing a branch length"); // Avoids ininite loops
                else if (n_char==tbuffer)
                {
                    step-=n_char+2;
                    reallocBuffer(&buffer, &tbuffer, tbuffer*10);
                }
                else
                {
                    --step;
                    *(buffer+n_char)=0; //Sets the end of the string to avoid problems in the sscanf
                    // *
                    /// Translates the buffer in a float number and asigns it as s_node::gen_length converted in number of generations of the current node</dd></dl>
                    sscanf(buffer,"%lf",&current_node->gen_length);
                }
                
                break;
            case '[':
                step+=2;
                ErrorReporter(GetLnodeParamsFromNexusComments(nexus,current_node,&step), ": parsing Nexus tree comments");
                if (current_node->n_nodes>0 && current_node->n_nodes!=ind_persp)
                {
                    n_gleaves+=current_node->n_nodes;
                    ++n_priv_ngleaves;
                    if (current_node->n_child!=0)
                    {
                        PrintXCharError(nexus, step, "\nNEWICK PARSING ERROR\n", "|<- Fixing the number of nodes/replicates in an internal node is not allowed, please, revisit the manual and your input trees\n");
                        ErrorReporter(SETTINGS_ERROR,NULL);
                    }
                }
                if (current_node->kind_node==LOSS || current_node->kind_node==RGC || current_node->kind_node==RTRFR)
                    --n_tleaves;
                
                break;
                
            default:
                // **
                /// <dl><dt>New leaf(code!= former ones).</dt><dd>
                
                // *
                /// Reads the name of the leaf in a buffer
                
                ResetBuffer(name_buffer, MAX_NAME);
                n_char=0;
                while (code!='(' && code!=')' && code!= ',' && code!= ';' && code!= ':' && code!= '[' && n_char<MAX_NAME && n_char<MAX_IT)
                {
                    *(name_buffer+n_char)=code;
                    ++n_char;
                    ++step;
                    code=*(nexus+step);
                }
                if (n_char==MAX_IT) ErrorReporter(LOOP_ERROR, ": parsing leaf info"); // Avoids ininite loops
                else if (n_char==MAX_NAME)
                {
                    fprintf(stderr,"\nSpecies names are being truncated, since they are bigger than the buffer. You could change this behaviour increasing the buffer size by changing the environment variable SIMPHY_MAXNAME\n"); //\todo Change the name_c implementation to avoid this behaviour.
#ifdef DEBUG
                    fflush(stderr);
#endif
                    while (code!='(' && code!=')' && code!= ',' && code!= ';' && code!= ':' && code!= '[' && n_char<MAX_IT)
                    {
                        ++step;
                        code=*(nexus+step);
                    }
                    if (n_char==MAX_IT) ErrorReporter(LOOP_ERROR, ": parsing leaf info"); // Avoids ininite loops
                }
                *(name_buffer+n_char)=0; //Sets the end of the string to avoid problems in the sscanf
                
                // *
                /// Searches for the maximum number of characters in a name.
                if (max_lname<n_char)
                    max_lname=n_char;
                
                // *
                /// Asigns the new name to the name_c::names
                strncpy(names->names+(ffree_codename*names->max_lname), name_buffer, names->max_lname);
                --step;
                
                // *
                /// New s_node by \ref NewSNodes
                anc_node=current_node;
                current_node=NewLNodes(1,0,max_children);
                // *
                /// Points the pointers of this new node and its ancestor
                *(anc_node->children+anc_node->n_child)=current_node;
                ++anc_node->n_child;
                current_node->anc_node=anc_node;
                // *
                /// Asigns the new name code (position) to the current node (s_node::sp_index)
                current_node->sp_index=ffree_codename;
                ++ffree_codename;
                // *
                /// Searches for the maximum number of children in the tree.</dd></dl></dd></dl>
                if(max_children<anc_node->n_child)
                {
                    max_children=anc_node->n_child;
                }
                
                break;
        }
        
        ++step;
        code=*(nexus+step);
        ++iteration;
        if (iteration==MAX_IT) ErrorReporter(LOOP_ERROR, NULL); // Avoids ininite loops
        
    }
    n_gleaves+=(n_tleaves-n_priv_ngleaves)*ind_persp;
    // ***
    /// Refines the readed nodes by \ref RefineLNodes
    
    RefineLNodes(root,n_gleaves,ind_persp,gen_time,n_dup,n_loss,n_trans,n_gc);
    max_lname++; //One extra character for \0.
    
    // ***
    /// Reallocates the name_c memory by \ref ReallocNames
    ReallocNames(names,max_lname);
    
    // ***
    /// Allocates memory for the readed tree and completes it.
    
    tree=NewLTree(0, n_leaves, n_gleaves, max_children, gen_time, Ne, mu); //Nodes have already been allocated.
    tree->root=root;
    tree->n_nodes=n_nodes;
    tree->species_tree=NULL;
    tree->gene_tree=NULL;
    
    // ***
    /// Fills the node indexes following a post_order.
    PostReorderLNodes(tree->root,&index);
    
    // ***
    /// Checks ultrametricity when measuring the species tree in time units.
    
    if(CheckUltrametricityLTree(tree)==-1)
        ErrorReporter(UNEXPECTED_VALUE,": locus tree is not ultrametric in time units. This doesn't allow to use the full model (transfers)\n");
    
    // ***
    /// Calculate the count probabilities by \ref RefineLNodes
    if (verbosity>4)
    {
        printf("\n\t\tCalculating lineage count probabilities of entering and leaving each locus tree branch... ");
#ifdef DBG
        fflush(stdout);
#endif
    }
    
    CalcProbsNLineagesLTree(root,Ne, 0, verbosity);//I should try to extend it for using polytomies
    
    if (verbosity>4)
    {
        printf("\n\t\tDone\n");
#ifdef DBG
        fflush(stdout);
#endif
    }
    
    // ***
    /// Frees dynamic memory (buffers)</dd></dl>
    free(buffer);
    free(name_buffer);
    buffer=NULL;
    
    if (verbosity>2)
    {
        printf("\n\t%d-node locus tree correctly built",(n_leaves*2)-1);
        if (verbosity>3)
        {
            printf(": ");
            WriteLNodesGen(root,names);
        }
        printf("\n");
        
    }
    
    return (tree);
}

l_tree * ParseNewickLTree (char * newick,name_c **names_ptr, int verbosity, double gen_time, int Ne, double mu, int ind_persp, int *n_dup, int *n_loss, int *n_trans, int *n_gc)
{
    l_node *current_node=NULL, *anc_node=NULL,*root=NULL;
    l_tree *tree=NULL;
    name_c * names=NULL;
    char code=' ';
    int register step=0;
    int n_char=0, iteration=0, ffree_codename=1, n_leaves=0, n_tleaves=0, n_gleaves=0, n_inodes=0, n_nodes=0, max_children=0, max_lname=0, n_replica=0,index=0, is_dup=0, n_paralog=0, n_priv_ngleaves=0;
    char *buffer=NULL,*name_buffer=NULL;
    size_t tbuffer=NUM_BUFFER;
    
    /// \attention Changes on buffer managing applied like in Nexus parsers but without posterior testting.
    
    // ****
    /// <dl><dt> Function structure </dt><dd>
    
    // ***
    /// Test of the newick tree string by \ref CheckNewickLTree
    
    ErrorReporter(CheckNewickLTree(newick), ": Newick locus tree error");
    
    // ***
    /// Error control
    
    if (gen_time==0)
        ErrorReporter(UNEXPECTED_VALUE, ": Inproper generation time");
    
    // ***
    /// Buffer allocation and initialization
    
    reallocBuffer(&buffer,&tbuffer,tbuffer);
    ResetBuffer(buffer,tbuffer);
    name_buffer=calloc(MAX_NAME,sizeof(char));
    
    // ***
    /// First read of the Newick tree. Obtains the number of internal nodes (")"), s_tree::n_leaves ("(" or "," not followed by "(") and s_tree::n_gleaves (s_tree::n_leaves + (replicas "/" -1 for each node)).
    while (*(newick+step)!=';')
    {
        if((*(newick+step)=='(' ||*(newick+step)==',' )&&*(newick+step+1)!='(') ++n_leaves;
        if(*(newick+step)==')') ++n_inodes;
        if(*(newick+step)=='/')
        {
            // * Reseting variables * //
            ResetBuffer(buffer, tbuffer);
            n_replica=0;
            n_char=0;
            ++step;
            code=*(newick+step);
            // * Reading replicas digits * //
            while (code<58 && code>47 && n_char<tbuffer && n_char<MAX_IT) // Int numbers
            {
                *(buffer+n_char)=code;
                ++n_char;
                ++step;
                code=*(newick+step);
            }
            if (n_char==MAX_IT) ErrorReporter(LOOP_ERROR, ": parsing replicates"); // Avoids ininite loops
            else if (n_char==tbuffer)
            {
                step-=n_char+2;
                reallocBuffer(&buffer, &tbuffer, tbuffer*10);
            }
            else
            {
                --step;
                // * Translating * //
                sscanf(buffer,"%ui",&n_replica);
                n_gleaves+=n_replica;
                ++n_priv_ngleaves;
            }
        }
        ++step;
        if (step==MAX_IT) ErrorReporter(LOOP_ERROR, NULL); // Avoids ininite loops
    }
    n_nodes=n_leaves+n_inodes;
    n_tleaves=n_leaves;
    
    // ***
    /// Name container allocation by \ref NewNames
    *names_ptr=NewNames(n_leaves,0);
    names=*names_ptr;
    
    // ***
    /// <dl><dt>Second read of the Newick tree (main loop)</dt><dd>
    
    step=0; //Reset
    code=*(newick+step); //First char
    
    while (code!=';') //End of the Newick tree
    {
        switch (code) {
            case '(':
                // **
                /// <dl><dt>New l_node (code="(").</dt><dd>
                if (iteration!=0) //New normal node
                {
                    anc_node=current_node;
                    // *
                    /// New l_node by \ref NewLNodes
                    current_node=NewLNodes(1, 0, max_children);
                    // *
                    /// Points the pointers of this new node and its ancestor
                    *(anc_node->children+anc_node->n_child)=current_node;
                    ++anc_node->n_child;
                    current_node->anc_node=anc_node;
                    // *
                    /// Searches for the maximum number of children in the tree.</dd></dl>
                    if(max_children<anc_node->n_child)
                    {
                        max_children=anc_node->n_child;
                    }
                    
                }
                else //New root node
                {
                    current_node=NewLNodes(1,0,max_children);
                    root=current_node;
                }
                
                break;
            case ',':
            case ')':
                // **
                /// Goes down one tree node (code=", or )")
                current_node=current_node->anc_node;
                anc_node=current_node->anc_node;
                break;
            case ':':
                // **
                /// <dl><dt>New branch length (code=":").</dt><dd>
                
                // *
                /// Reads all the following integers and . in a buffer
                ResetBuffer(buffer, tbuffer);
                n_char=0;
                ++step;
                code=*(newick+step);
                while (code<58 && (code>47 || code==46) && n_char<tbuffer && n_char<MAX_IT)
                {
                    *(buffer+n_char)=code;
                    ++n_char;
                    ++step;
                    code=*(newick+step);
                }
                if (n_char==MAX_IT) ErrorReporter(LOOP_ERROR, ": parsing branch length"); // Avoids ininite loops
                else if (n_char==tbuffer)
                {
                    step-=n_char+2;
                    reallocBuffer(&buffer, &tbuffer, tbuffer*10);
                }
                else
                {
                    --step;
                    *(buffer+n_char)=0; //Sets the end of the string to avoid problems in the sscanf
                    // *
                    /// Translates the buffer in a float number and asigns it as s_node::gen_length of the current node</dd></dl>
                    sscanf(buffer,"%lf",&current_node->gen_length);
                }
                
                break;
            case '*':
                // **
                /// <dl><dt>Lineage specific substitution rate multi(code="*").</dt><dd>
                
                // *
                /// Reads all the following integers and . in a buffer
                ResetBuffer(buffer, tbuffer);
                n_char=0;
                ++step;
                code=*(newick+step);
                while (code<58 && (code>47 || code==46) && n_char<tbuffer && n_char<MAX_IT)
                {
                    *(buffer+n_char)=code;
                    ++n_char;
                    ++step;
                    code=*(newick+step);
                }
                if (n_char==MAX_IT) ErrorReporter(LOOP_ERROR, ": parsing branch specific substitution rate"); // Avoids ininite loops
                else if (n_char==tbuffer)
                {
                    step-=n_char+2;
                    reallocBuffer(&buffer, &tbuffer, tbuffer*10);
                }
                else
                {
                    --step;
                    *(buffer+n_char)=0; //Sets the end of the string to avoid problems in the sscanf
                    // *
                    /// Translates the buffer in a float number and asigns it as l_node::mu_mult </dd></dl>
                    sscanf(buffer,"%lf",&current_node->mu_mult);
                }
                break;
                
            case '~':
                // **
                /// <dl><dt>Lineage specific generation time multi (code="~").</dt><dd>
                
                // *
                /// Reads all the following integers and . in a buffer
                ResetBuffer(buffer, tbuffer);
                n_char=0;
                ++step;
                code=*(newick+step);
                while (code<58 && (code>47 || code==46) && n_char<tbuffer && n_char<MAX_IT)
                {
                    *(buffer+n_char)=code;
                    ++n_char;
                    ++step;
                    code=*(newick+step);
                }
                if (n_char==MAX_IT) ErrorReporter(LOOP_ERROR, ": parsing branch specific generation time"); // Avoids ininite loops
                else if (n_char==tbuffer)
                {
                    step-=n_char+2;
                    reallocBuffer(&buffer, &tbuffer, tbuffer*10);
                }
                else
                {
                    --step;
                    *(buffer+n_char)=0; //Sets the end of the string to avoid problems in the sscanf
                    // *
                    /// Translates the buffer in a float number and asigns it as s_node::gtime_mult </dd></dl>
                    sscanf(buffer,"%lf",&current_node->gtime_mult);
                }
                break;
                
            case '#':
                // **
                /// <dl><dt>New effective population size (Ne) (code="#").</dt><dd>
                
                // *
                /// Reads all the following integers in a buffer
                ResetBuffer(buffer, tbuffer);
                n_char=0;
                ++step;
                code=*(newick+step);
                while (code<58 && code>47 && n_char<tbuffer && n_char<MAX_IT)
                {
                    *(buffer+n_char)=code;
                    ++n_char;
                    ++step;
                    code=*(newick+step);
                }
                if (n_char==MAX_IT) ErrorReporter(LOOP_ERROR, ": parsing branch specific population size"); // Avoids ininite loops
                else if (n_char==tbuffer)
                {
                    step-=n_char+2;
                    reallocBuffer(&buffer, &tbuffer, tbuffer*10);
                }
                else
                {
                    --step;
                    *(buffer+n_char)=0; //Sets the end of the string to avoid problems in the sscanf
                    // *
                    /// Translates the buffer in a integer number and asigns it as l_node::Ne of the current node</dd></dl>
                    sscanf(buffer,"%ui",&current_node->Ne);
                    if (current_node->Ne==0)
                        ErrorReporter(SETTINGS_ERROR,": Null branch specific population size is not allowed");
                }
                
                break;
            case '/':
                // **
                /// <dl><dt>New number of replicas (code="/").</dt><dd>
                
                // *
                /// Reads all the following integers in a buffer
                ResetBuffer(buffer, tbuffer);
                n_replica=0;
                n_char=0;
                ++step;
                code=*(newick+step);
                while (code<58 && code>47 && n_char<tbuffer && n_char<MAX_IT)
                {
                    *(buffer+n_char)=code;
                    ++n_char;
                    ++step;
                    code=*(newick+step);
                }
                if (n_char==MAX_IT) ErrorReporter(LOOP_ERROR, ": parsing number of replicates"); // Avoids ininite loops
                else if (n_char==tbuffer)
                {
                    step-=n_char+2;
                    reallocBuffer(&buffer, &tbuffer, tbuffer*10);
                }
                else
                {
                    --step;
                    *(buffer+n_char)=0; //Sets the end of the string to avoid problems in the sscanf
                    // *
                    /// Translates the buffer in a integer number and asigns it as l_node::n_nodes of the current node</dd></dl>
                    sscanf(buffer,"%ui",&n_replica);
                    if (n_replica==0)
                        ErrorReporter(SETTINGS_ERROR, ": Null number of replicates is not allowed");
                    current_node->n_nodes=n_replica;
                }
                break;
            case '%':
                // **
                /// <dl><dt>New duplication info (code="%").</dt><dd>
                
                // *
                /// Reads all the following integers in a buffer
                ResetBuffer(buffer, tbuffer);
                n_replica=0;
                n_char=0;
                ++step;
                code=*(newick+step);
                while (code<58 && code>47 && n_char<tbuffer && n_char<MAX_IT)
                {
                    *(buffer+n_char)=code;
                    ++n_char;
                    ++step;
                    code=*(newick+step);
                }
                if (n_char==MAX_IT) ErrorReporter(LOOP_ERROR, ": parsing node kind"); // Avoids ininite loops
                else if (n_char==tbuffer)
                {
                    step-=n_char+2;
                    reallocBuffer(&buffer, &tbuffer, tbuffer*10);
                }
                else
                {
                    --step;
                    *(buffer+n_char)=0; //Sets the end of the string to avoid problems in the sscanf
                    // *
                    /// Translates the buffer in a integer number and asigns it as l_node::kind_node of the current node</dd></dl>
                    sscanf(buffer,"%ui",&is_dup);
                    current_node->kind_node=is_dup;
                    if (is_dup==LOSS || is_dup==RGC || is_dup==RTRFR)
                    {
                        --n_tleaves;
                    }
                }
                break;
            case '_':
                
                // **
                /// <dl><dt>New paralog info (code="_").</dt><dd>
                
                // *
                /// Reads all the following integers in a buffer
                ResetBuffer(buffer, tbuffer);
                n_paralog=0;
                n_char=0;
                ++step;
                code=*(newick+step);
                while (code<58 && code>47 && n_char<tbuffer && n_char<MAX_IT)
                {
                    *(buffer+n_char)=code;
                    ++n_char;
                    ++step;
                    code=*(newick+step);
                }
                if (n_char==MAX_IT) ErrorReporter(LOOP_ERROR, ": parsing locus id"); // Avoids ininite loops
                else if (n_char==tbuffer)
                {
                    step-=n_char+2;
                    reallocBuffer(&buffer, &tbuffer, tbuffer*10);
                }
                else
                {
                    --step;
                    *(buffer+n_char)=0; //Sets the end of the string to avoid problems in the sscanf
                    // *
                    /// Translates the buffer in a integer number and asigns it as l_node::paralog of the current node</dd></dl>
                    sscanf(buffer,"%ui",&n_paralog);
                    current_node->paralog=n_paralog;
                }
                break;
                
            default:
                // **
                /// <dl><dt>New leaf(code!= former ones).</dt><dd>
                
                // *
                /// Reads the name of the leaf in a buffer
                
                ResetBuffer(name_buffer, MAX_NAME);
                n_char=0;
                while (code!='(' && code!=')' && code!= ',' && code!= ';' && code!= ':' && code!= '_' && code != '~' && n_char<MAX_NAME && n_char<MAX_IT)
                {
                    *(name_buffer+n_char)=code;
                    ++n_char;
                    ++step;
                    code=*(newick+step);
                }
                if (n_char==MAX_IT) ErrorReporter(LOOP_ERROR, ": parsing leaf info"); // Avoids ininite loops
                else if (n_char==MAX_NAME)
                {
                    fprintf(stderr,"Species names are being truncated, since they are bigger than the buffer. You could change this behaviour increasing the buffer size by changing the environment variable SIMPHY_MAXNAME\n"); //\todo Change the name_c implementation to avoid this behaviour.
#ifdef DEBUG
                    fflush(stderr);
#endif
                    while (code!='(' && code!=')' && code!= ',' && code!= ';' && code!= ':' && code!= '[' && n_char<MAX_IT)
                    {
                        ++step;
                        code=*(newick+step);
                    }
                    if (n_char==MAX_IT) ErrorReporter(LOOP_ERROR, ": parsing leaf info"); // Avoids ininite loops

                }
                *(name_buffer+n_char)=0; //Sets the end of the string to avoid problems in the sscanf
                
                // *
                /// Searches for the maximum number of characters in a name.
                if (max_lname<n_char)
                    max_lname=n_char;
                
                // *
                /// Asigns the new name to the name_c::names
                strncpy(names->names+(ffree_codename*names->max_lname), name_buffer, names->max_lname);
                --step;
                
                // *
                /// New l_node by \ref NewLNodes
                anc_node=current_node;
                current_node=NewLNodes(1,0,max_children);
                // *
                /// Points the pointers of this new node and its ancestor
                *(anc_node->children+anc_node->n_child)=current_node;
                ++anc_node->n_child;
                current_node->anc_node=anc_node;
                // *
                /// Asigns the new name code (position) to the current node (l_node::lp_index)
                current_node->sp_index=ffree_codename;
                ++ffree_codename;
                // *
                /// Searches for the maximum number of children in the tree.</dd></dl></dd></dl>
                if(max_children<anc_node->n_child)
                {
                    max_children=anc_node->n_child;
                }
                
                break;
        }
        
        ++step;
        code=*(newick+step);
        ++iteration;
        if (iteration==MAX_IT) ErrorReporter(LOOP_ERROR, NULL); // Avoids ininite loops
        
    }
    
    n_gleaves+=(n_tleaves-n_priv_ngleaves)*ind_persp; //Addition of the ngleaves of the nodes using the common number of individuals (ind_persp)
    
    // ***
    /// Refines the readed nodes by \ref RefineLNodes
    
    RefineLNodes(root,n_gleaves,ind_persp,gen_time,n_dup,n_loss,n_trans,n_gc);
    max_lname++; //One extra character for \0.
    
    // ***
    /// Reallocates the name_c memory by \ref ReallocNames
    ReallocNames(names,max_lname);
    
    // ***
    /// Allocates memory for the readed tree and completes it.
    
    tree=NewLTree(0, n_leaves, n_gleaves, max_children, gen_time, Ne, mu);//Nodes have already been allocated.
    tree->root=root;
    tree->n_nodes=n_nodes;
    tree->species_tree=NULL;
    tree->gene_tree=NULL;
    
    // ***
    /// Fills the node indexes following a post_order.
    PostReorderLNodes(tree->root,&index);
    
    // ***
    /// Checks ultrametricity when measuring the species tree in time units.
    
    if(CheckUltrametricityLTree(tree)==-1)
        ErrorReporter(UNEXPECTED_VALUE,"The species tree is not ultrametric in time units. This doesn't allow to use the full model (transfers)\n");
    
    // ***
    /// Calculate the count probabilities by \ref RefineLNodes
    if (verbosity>4)
    {
        printf("\n\t\tCalculating lineage count probabilities of entering and leaving each locus tree branch... ");
#ifdef DBG
        fflush(stdout);
#endif
    }
    
    CalcProbsNLineagesLTree(root,Ne, 0, verbosity);//I should try to extend it for using polytomies
    
    if (verbosity>4)
    {
        printf("\n\t\tDone\n");
#ifdef DBG
        fflush(stdout);
#endif
    }
    
    // ***
    /// Frees dynamic memory (buffers)</dd></dl>
    free(buffer);
    free(name_buffer);
    buffer=NULL;
    
    if (verbosity>2)
    {
        printf("\n\t\t %d-node locus tree correctly built",(n_leaves*2)-1);
        if (verbosity>3)
        {
            printf(": ");
            WriteLNodesGen(root,names);
        }
        printf("\n");
        
    }
    
    return (tree);
}


// *** Tree memory manage *** //

// ** Tree creation ** //

s_tree * NewSTree (int n_nodes, int n_leaves, int n_gleaves, int max_children, double gen_time, int Ne, double mu)
{
    s_tree * tree=NULL;
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Tree memory allocation
    tree=calloc(1,sizeof(s_tree));
    ErrorReporter((long int)tree, NULL);
    
    // *
    /// Tree initialization
    tree->n_nodes=n_nodes;
    tree->n_leaves=n_leaves;
    tree->n_gleaves=n_gleaves;
    tree->max_children=max_children;
    tree->gen_time=gen_time;
    tree->Ne=Ne;
    tree->mu=mu;
    tree->locus_tree=NULL;
    tree->gene_tree=NULL;
    
    // *
    /// Tree nodes allocation and initialization by \ref NewSNodes </dd></dl>
    
    if (n_nodes>1)
    {
        tree->m_node=NewSNodes(n_nodes,max_children);
        tree->root=NULL;
    }
    else if (n_nodes==1)
    {
        tree->root=NewSNodes(n_nodes, max_children);
        tree->m_node=NULL;
    }
    else
    {
        tree->root=NULL;
        tree->m_node=NULL;
    }
    
    return (tree);
}

long int NewBDSTree (s_tree ** out_tree, int leaves, double time, double b_rate, double d_rate, double gen_time, int Ne, double mu, int ind_per_sp, double outgroup, int complete, int mrca_time, gsl_rng *seed,int out_time ,int labels,int verbosity)
{
    // ***** Declaration and initialization of variables ***** //
    
    int eq_rates=0, j=0,k=0, node_index=0, n_leaves=0,pn_leaves=0,pn_events=0,n_nodes=0, avail_leaves=0, extra_nodes=0, iter=0;
    //unsigned long int stats_leaves=0;
    //double stats_time=0;
    double w_prob=d_rate+b_rate, b_prob=(b_rate/(d_rate+b_rate)), random=0, inv_leaves=0, res1=0, res2=0, *i_nodes=NULL, s_time=0, c_time=0,r_time=0;
    s_node *w_node1=NULL,*w_node2=NULL,*anc_node=NULL, **node_ptrs=NULL;
    
    
    // *******
    /// <dl><dt> Function structure </dt><dd>
    
    // ******
    /// Error control
    
    if (w_prob<=0 || (b_rate<d_rate && leaves!=0))
        return (SETTINGS_ERROR); ///\todo Mejorar el control de errores
    if (gen_time==0)
        return (UNEXPECTED_VALUE);
    if (d_rate==b_rate)
        eq_rates=1;
    
    // ******
    /// Method selection, BDSA or SSA algorithms.
    
    if (leaves>0)
    {
        // *****
        /// <dl><dt>BDSA (fixed number of leaves). Faster than GSA, but without information about extinct lineages.</dt><dd>
        if (verbosity>3)
        {
            printf("\n\tUsing BDSA algorithm (fixed number of leaves) to simulate the trees...");
#ifdef DBG
            fflush(stdout);
#endif
        }
        
        // ****
        /// Dynamic memory allocation and initialization (tree and other variables)
        if (outgroup>0)
        {
            (*out_tree)=NewSTree((((leaves+1)*2)-1), leaves+1, leaves+1, 2,gen_time,Ne, mu);
            
        }
        else
        {
            (*out_tree)=NewSTree(((leaves*2)-1), leaves, leaves, 2,gen_time,Ne, mu);
        }
        
        i_nodes=calloc(leaves-1, sizeof(double));
        ErrorReporter((long int)i_nodes, NULL);
        node_ptrs=calloc((*out_tree)->n_nodes,sizeof(s_node *));
        ErrorReporter((long int)node_ptrs, NULL);
        
        // ****
        /// If the TMRCA is given, it is necessary to sample two BDSA processes with n-x and x leaves respectively x:U[1,n-1]. Else, it is sampled using the inverse of the cdf of the time of the origin of the tree (flat prior), conditional on having n species at the present, Hartmann et al., 2010; Gernhard, 2008.
        if (time==0)
        {
            // ****
            /// Samples the tree origin (From user options or sampled )
            random=gsl_rng_uniform_pos(seed); //If it is 0 we get negative times, although in the paper they say [0,1]
            inv_leaves=1/(double)leaves;
            
            if (eq_rates>0)
                time=1/(b_rate*(pow(random,-inv_leaves)-1));
            else
            {
                res1=pow(random,inv_leaves);
                time=log((1-(res1*d_rate/b_rate))/(1-res1))/(b_rate-d_rate);
            }
            
            // ****
            /// Sample internal node branches using the quantile function of the speciation times, Hartmann et al., 2010; Gernhard, 2008.
            
            for (j=0;j<leaves-1;++j)
            {
                random=gsl_rng_uniform_pos(seed); //If it is 0 we get branch lengths=0, although in the paper they say [0,1]
                if (eq_rates>0)
                    *(i_nodes+j)=(random*time)/(1+(b_rate*time*(1-random)));
                else
                {
                    res1=exp((b_rate-d_rate)*-time);
                    res2=b_rate-(d_rate*res1);
                    *(i_nodes+j)=(log((res2-(d_rate*random*(1-res1)))/(res2-(b_rate*random*(1-res1))))/(b_rate-d_rate));
                }
            }
            
            // ****
            /// Orders the times in ascending order using a quicksort algorithm
            qsort(i_nodes, leaves-1, sizeof(*i_nodes), Compare_DBL);
            
            r_time=*(i_nodes+leaves-2);
            
            if (verbosity>4)
            {
                if (outgroup==0)
                    printf("\n\t\t\tTree TMRCA (sampled) = %.8lf",r_time);
                else
                    printf("\n\t\t\tIngroup TMRCA (sampled) = %.8lf, root time = %.8lf",r_time,r_time+r_time*outgroup);
#ifdef DBG
                fflush(stdout);
#endif
            }
            
            if (verbosity>4)
            {
                printf("\n\t\t\tReconstructing backwards the tree from sampled duplication times...");
                if (verbosity == 6)
                {
                    printf("\n\t\t\t\tDuplication times:");
                    for (j=0; j<leaves-1; ++j)
                    {
                        printf(" %.8lf,",r_time+r_time*outgroup-*(i_nodes+j));
                    }
                }
#ifdef DBG
                fflush(stdout);
#endif
            }
            
            for (j=0;j<(*out_tree)->n_nodes;++j)
            {
                *(node_ptrs+j)=(*out_tree)->m_node+j;
                (*(node_ptrs+j))->n_gen=(r_time+r_time*outgroup)/gen_time;
                (*(node_ptrs+j))->time=r_time+r_time*outgroup;
            }
            
            n_leaves=leaves;
            
            // ****
            /// Reconstructs the tree from the birth-death times
            for (j=0; j<leaves-1;++j)
            {
                // * First offspring node * //
                node_index= (int)gsl_rng_uniform_int(seed,n_leaves);//Samples nodes
                w_node1=*(node_ptrs+node_index);
                --n_leaves;
                *(node_ptrs + node_index) = *(node_ptrs + n_leaves);/* The pointer to this node now points to the last node_ptr, inaccesible due to --n_leaves*/
                
                // * Second offspring node * //
#ifndef NO_VAR
                if (n_leaves>1)
                    node_index= (int)gsl_rng_uniform_int(seed,n_leaves);//Samples nodes
                else
                    node_index=0;
#endif
#ifdef NO_VAR
                node_index=(int)gsl_rng_uniform_int(seed,n_leaves);
#endif
                w_node2= *(node_ptrs + node_index);
                *(node_ptrs + node_index)= (*out_tree)->m_node + j + leaves; //Now the node_ptr that pointed to w_node2 points to the new ancestral node, for using it in the next loop iteration
                
                // ** Ancestral node ** //
                anc_node= (*out_tree)->m_node + j + leaves; //Uses as anc node the next avaliable memory block for internal nodes
                
                // ** Filling working nodes ** //
                // * Pointers * //
                *(anc_node->children)=w_node1;
                *(anc_node->children+1)=w_node2;
                w_node1->anc_node=anc_node;
                w_node2->anc_node=anc_node;
                anc_node->time=r_time+r_time*outgroup-*(i_nodes+j);
                // * Branch lenghts and times * //
                anc_node->n_gen=anc_node->time/gen_time;
                anc_node->n_child=2;
                w_node1->gen_length=w_node1->n_gen-anc_node->n_gen;
                w_node2->gen_length=w_node2->n_gen-anc_node->n_gen;
                
            }
            (*out_tree)->root=anc_node;
            
        }
        else
        {
            // ****
            /// TMRCA given
            if (verbosity>4)
            {
                if (outgroup==0)
                    printf("\n\t\t\tTMRCA (user defined) = %.8lf",time);
                else
                    printf("\n\t\t\tIngroup TMRCA (user defined) = %.8lf, root time = %.8lf",time, time+time*outgroup);
#ifdef DBG
                fflush(stdout);
#endif
            }
            
            // ***
            /// Randomly splitting the number of leaves in two processes
            r_time=time;
            
            n_leaves=(int)gsl_rng_uniform_int(seed,leaves-1)+1; //Number of leaves for the first BDSA process
            
            // ***
            /// Initializating and pointing tree nodes (only retained for leaves
            for (j=0;j<(*out_tree)->n_nodes;++j)
            {
                *(node_ptrs+j)=(*out_tree)->m_node+j;
                (*(node_ptrs+j))->n_gen=(time+time*outgroup)/gen_time;
                (*(node_ptrs+j))->time=time+time*outgroup;
            }
            if (leaves==2)
                k=3;
            else if (n_leaves>1 && n_leaves<leaves-1)
            {
                k=1;
                pn_leaves=0;
            }
            else if (n_leaves==1)
            {
                pn_leaves=1;
                n_leaves=leaves-pn_leaves;
                k=2;
            }
            else
            {
                k=1;
                pn_leaves=0;
            }

            pn_events=0;
            // ****
            /// Loop for the simulated subtrees
            for (;k<=2;++k)
            {
                extra_nodes=pn_leaves==0?0:pn_leaves*2-1; //Number of previously used nodes, number of previous leaves + internal nodes.
                pn_events=pn_leaves==0?0:pn_leaves-1; //Number of previously used internal nodes.
                
                // ****
                /// Sample internal node branches using the quantile function of the speciation times, Hartmann et al., 2010; Gernhard, 2008.
                for (j=0;j<n_leaves-1;++j)
                {
                    random=gsl_rng_uniform_pos(seed); //If it is 0 we get branch lengths=0, although in the paper they say [0,1]
                    if (eq_rates>0)
                        *(i_nodes+pn_events+j)=(random*time)/(1+(b_rate*time*(1-random)));
                    else
                    {
                        res1=exp((b_rate-d_rate)*-time);
                        res2=b_rate-(d_rate*res1);
                        *(i_nodes+pn_events+j)=(log((res2-(d_rate*random*(1-res1)))/(res2-(b_rate*random*(1-res1))))/(b_rate-d_rate));
                    }
                }
                
                // ****
                /// Orders the times in ascending order using a quicksort algorithm
                qsort(i_nodes+pn_events, n_leaves-1, sizeof(*i_nodes), Compare_DBL);
                
                if (verbosity>4 && n_leaves>1)
                {
                    printf("\n\t\t\tReconstructing backwards the tree %d from sampled duplication times...",k);
                    if (verbosity == 6)
                    {
                        printf("\n\t\t\t\tDuplication times:");
                        for (j=0; j<n_leaves-1; ++j)
                        {
                            printf(" %.8lf,",r_time+r_time*outgroup-*(i_nodes+j+pn_events));
                        }
                    }
#ifdef DBG
                    fflush(stdout);
#endif
                }
                
                // ****
                /// Reconstructs the tree from the birth-death times
                avail_leaves=n_leaves;
                for (j=0; j<n_leaves-1;++j)
                {
                    // * First offspring node * //
                    node_index= (int)gsl_rng_uniform_int(seed,avail_leaves)+pn_leaves;//Samples nodes
                    w_node1=*(node_ptrs+node_index);
                    --avail_leaves;
                    *(node_ptrs + node_index) = *(node_ptrs + avail_leaves+pn_leaves);/* The pointer to this node now points to the last node_ptr, inaccesible due to --n_leaves*/
                    
                    // * Second offspring node * //
#ifndef NO_VAR
                    if (avail_leaves>1)
                        node_index= (int)gsl_rng_uniform_int(seed,avail_leaves)+pn_leaves;//Samples nodes
                    else
                        node_index=0+pn_leaves;
#endif
#ifdef NO_VAR
                    node_index=(int)gsl_rng_uniform_int(seed,avail_leaves)+pn_leaves;
#endif
                    w_node2= *(node_ptrs + node_index);
                    *(node_ptrs + node_index)= (*out_tree)->m_node+j+leaves+pn_events;//Now the node_ptr that pointed to w_node2 points to the new ancestral node, for using it in the next loop iteration
                    
                    // ** Ancestral node ** //
                    anc_node= (*out_tree)->m_node+j+leaves+pn_events; //Uses as anc node the next avaliable memory block for internal nodes: Previous used internal nodes (extra_nodes) + currently used internal nodes (j) + displacement (leaves)
                    
                    // ** Filling working nodes ** //
                    // * Pointers * //
                    *(anc_node->children)=w_node1;
                    *(anc_node->children+1)=w_node2;
                    w_node1->anc_node=anc_node;
                    w_node2->anc_node=anc_node;
                    anc_node->time=r_time+r_time*outgroup-*(i_nodes+j+pn_events);
                    // * Branch lenghts and times * //
                    anc_node->n_gen=anc_node->time/gen_time;
                    anc_node->n_child=2;
                    w_node1->gen_length=w_node1->n_gen-anc_node->n_gen;
                    w_node2->gen_length=w_node2->n_gen-anc_node->n_gen;
                    
                }
                
                pn_leaves=n_leaves;
                n_leaves=leaves-pn_leaves;
            }
            
            (*out_tree)->root=(*out_tree)->m_node+leaves*2-2;
            ((*out_tree)->root)->time=r_time*outgroup;
            ((*out_tree)->root)->gen_length=0;
            ((*out_tree)->root)->n_gen=((*out_tree)->root)->time/gen_time;
            ((*out_tree)->root)->n_child=2;
            if (pn_leaves==1)
            {
                *(((*out_tree)->root)->children)=(*out_tree)->m_node+leaves*2-3;
                *(((*out_tree)->root)->children+1)=(*out_tree)->m_node+leaves-1;
            }
            else if (n_leaves==1)
            {
                *(((*out_tree)->root)->children)=(*out_tree)->m_node;
                *(((*out_tree)->root)->children+1)=(*out_tree)->m_node+leaves*2-3;
            }
            else
            {
                *(((*out_tree)->root)->children)=(*out_tree)->m_node+leaves+n_leaves-2;
                *(((*out_tree)->root)->children+1)=(*out_tree)->m_node+leaves*2-3;
            }
            (*(((*out_tree)->root)->children+1))->anc_node=(*out_tree)->root;
            (*(((*out_tree)->root)->children+1))->gen_length=(*(((*out_tree)->root)->children+1))->n_gen-((*out_tree)->root)->n_gen;
            (*(((*out_tree)->root)->children))->anc_node=(*out_tree)->root;
            (*(((*out_tree)->root)->children))->gen_length=(*(((*out_tree)->root)->children))->n_gen-((*out_tree)->root)->n_gen;
            
        }

        // ****
        /// Node info update
        for (j=0;j<leaves;++j)
        {
            w_node1=(*out_tree)->m_node+j;
            w_node1->sp_index=j+1;
            w_node1->n_replicas=ind_per_sp;
        }
        if (outgroup>0)
        {
            w_node1=(*out_tree)->m_node+leaves*2; //Outgroup
            w_node2=(*out_tree)->m_node+leaves*2-1; //New root
            
            w_node1->anc_node=w_node2;
            *w_node2->children=(*out_tree)->root;
            *(w_node2->children+1)=w_node1;
            w_node2->n_gen=0;
            w_node2->time=0;
            w_node2->n_child=2;
            (*out_tree)->root->gen_length=outgroup*r_time/gen_time;
            (*out_tree)->root->anc_node=w_node2;
            
            w_node1->time=r_time+outgroup*r_time;
            w_node1->gen_length=w_node1->time/gen_time;
            w_node1->n_gen=w_node1->gen_length;
            
            w_node1->n_replicas=1;
            
            (*out_tree)->root=w_node2;
            (*out_tree)->n_leaves=leaves+1;
            (*out_tree)->n_gleaves=leaves*ind_per_sp+1;
            (*out_tree)->locus_tree=NULL;
            (*out_tree)->gene_tree=NULL;
        }
        else
        {
            // ****
            /// Tree info update
            (*out_tree)->n_leaves=leaves;
            (*out_tree)->n_gleaves=leaves*ind_per_sp;
            (*out_tree)->locus_tree=NULL;
            (*out_tree)->gene_tree=NULL;
        }

        // ****
        /// Frees allocated memory</dd></dl>
        free (i_nodes);
        i_nodes=NULL;
        free (node_ptrs);
        node_ptrs=NULL;
        
        if (verbosity>4)
        {
            printf("\n\t\t\tDone");
#ifdef DBG
            fflush(stdout);
#endif
        }


    }
    else
    {
        
        // *****
        /// <dl><dt>SSA algorithm (GSA is not required with fixed time)</dt><dd>
        if (verbosity>3)
        {
            printf("\n\tUsing SSA algorithm (fixed generations) to simulate the species tree...");
#ifdef DBG
            fflush(stdout);
#endif
            //            if (verbosity>5)
            //            {
            //                stats_leaves=0;
            //            }
        }
        while (iter<=MAX_IT)
        {
            ++iter;
            
            // ****
            /// Variable initialization
            node_ptrs=calloc(MAX_LEAVES,sizeof(s_node *));
            ErrorReporter((long int)node_ptrs,NULL);
            (*out_tree)=NewSTree(1, 1, 0, 2,gen_time,Ne, mu); // Recursion oriented tree (root).
            
            
            // ****
            /// Iteration-dependent variable initialization
            
            // ****
            /// Start of the tree (one or two lineages)
            w_node1=NewSNodes(1, 2);
            *((*out_tree)->root->children)=w_node1;//This root node is only a way to set the start of the tree. It is deleted at the end.
            w_node1->anc_node=(*out_tree)->root;
            anc_node=w_node1;
            extra_nodes=0;
            c_time=0;
            
            if (mrca_time==0)
            {
                *(node_ptrs)=w_node1; // The actual root of the tree
                n_leaves=1;
                avail_leaves=1;
                n_nodes=1;
                (*out_tree)->root->n_gen=0;
                s_time=-log(gsl_rng_uniform_pos(seed))/(w_prob*avail_leaves);
            }
            else
            {
                w_node1=NewSNodes(1,2);
                w_node2=NewSNodes(1,2);
                *(anc_node->children)=w_node1;
                *(anc_node->children+1)=w_node2;
                w_node1->anc_node=anc_node;
                w_node2->anc_node=anc_node;
                anc_node->n_gen=0;
                anc_node->n_child=2;
                n_leaves=2;
                avail_leaves=2;
                n_nodes=3;
                *(node_ptrs)=w_node1;
                *(node_ptrs+1)=w_node2;
                s_time=-log(gsl_rng_uniform_pos(seed))/(w_prob*avail_leaves);
            }
            
            // ****
            /// <dl><dt>Birth-death process loop (while final time is not exceeded)</dt><dd>
            while (c_time+s_time<=time && avail_leaves!=0)
            {
                // ***
                /// Memory error test (maximum number of leaves)
                if (n_leaves>=MAX_LEAVES)
                {
                    fprintf(stderr,"\nMaximum number of leaves reached. After checking your simulation parameters carefully (mainly species tree simulation birth rate) you can sort this problem out setting the environmental variable SIMPHY_MAXLEAVES with a value higher than %d\n",MAX_LEAVES);
                    return(MEM_ERROR);
                }
                
                // ***
                /// Chooses the leave for the next event
#ifndef NO_VAR
                if(avail_leaves>1)
                    node_index=(int)gsl_rng_uniform_int(seed,avail_leaves);
                else
                    node_index=0;
#endif
#ifdef NO_VAR
                node_index=(int)gsl_rng_uniform_int(seed,avail_leaves);
#endif
                // ***
                /// <dl><dt>Chooses the process</dt><dd>
                if (b_prob>gsl_rng_uniform(seed))
                {
                    // **
                    /// <dl><dt>Birth</dt><dd>
                    
                    // *
                    /// New nodes
                    w_node1=NewSNodes(1,2);
                    w_node2=NewSNodes(1,2);
                    
                    // *
                    /// Tree pointers and branch lengths reconfiguration
                    anc_node=*(node_ptrs+node_index);
                    *(anc_node->children)=w_node1;
                    *(anc_node->children+1)=w_node2;
                    w_node1->anc_node=anc_node;
                    w_node2->anc_node=anc_node;
                    
                    anc_node->time=c_time+s_time;
                    anc_node->n_gen=anc_node->time/gen_time;
                    anc_node->gen_length=anc_node->n_gen-anc_node->anc_node->n_gen;
                    anc_node->n_child=2;
                    
                    // *
                    /// New leaves addition</dl>
                    *(node_ptrs+node_index)=w_node1;
                    *(node_ptrs+avail_leaves)=w_node2;
                    ++avail_leaves;
                    ++n_leaves;
                    n_nodes+=2;
                    
                    if (verbosity==5)
                    {
                        printf("\n\t\t\tNew duplication");
#ifdef DBG
                        fflush(stdout);
#endif
                    }
                    else if (verbosity==6)
                    {
                        printf("\n\t\t\tNew duplication, %.8lf time units",anc_node->time);
#ifdef DBG
                        fflush(stdout);
#endif
                        
                    }
                    
                }
                else
                {
                    // **
                    /// <dl><dt>Death</dt><dd>
                    
                    // *
                    /// It finds the selected node. If complete == 0, it deletes the selected node and his ancestor, reconfiguring the brother node as it had been the only node of this lineage. Else, the selected node is set as non avaliable.</dl></dl>
                    
                    w_node1=*(node_ptrs+node_index);
                    
                    if (avail_leaves==1)
                    {
                        avail_leaves=0;
                        --n_leaves;
                        n_nodes=0;
                    }
                    else
                    {
                        // * Selecting nodes * //
                        anc_node=w_node1->anc_node;
                        if (*(anc_node->children) == w_node1)
                            w_node2=*(anc_node->children+1);
                        else
                            w_node2=*(anc_node->children);
                        
                        if (complete==0)
                        {
                            if (anc_node->anc_node==(*out_tree)->root)
                                avail_leaves=0; // Tree must be reset, since one of the original branches from the MRCA has been lost, and therefore the tMRCA would not fit the asked height.
                            else
                            {
                                // * Reconfiguring pointers * //
                                w_node2->anc_node=anc_node->anc_node;
                                if (*(w_node2->anc_node->children)==anc_node)
                                    *(w_node2->anc_node->children)=w_node2;
                                else
                                    *(w_node2->anc_node->children+1)=w_node2;
                                
                                // * Reconfiguring branches * //
                                if (w_node2->anc_node!=NULL)
                                    w_node2->gen_length=w_node2->n_gen-w_node2->anc_node->n_gen;
                                
                                // * Deleting one available leaf and one internal node * //
                                *(node_ptrs+node_index)=*(node_ptrs+avail_leaves-1);
                                free(w_node1);
                                w_node1=NULL;
                                free(anc_node);
                                anc_node=NULL;
                                --avail_leaves;
                                --n_leaves;
                                n_nodes-=2;
                            }
                        }
                        else
                        {
                            // * Reconfiguring branches * //
                            w_node1->time=c_time+s_time;
                            w_node1->n_gen=w_node1->time/gen_time;
                            w_node1->gen_length=w_node1->n_gen-anc_node->n_gen;
                            
                            // * Setting the leaf as non avaliable * //
                            *(node_ptrs+node_index)=*(node_ptrs+avail_leaves-1);
                            --avail_leaves;
                            
                            w_node1->sp_index=n_leaves-avail_leaves;
                            n_nodes-=2;
                            extra_nodes+=2;
                            
                        }
                    }
                    if (verbosity==5)
                    {
                        printf("\n\t\t\tNew loss");
#ifdef DBG
                        fflush(stdout);
#endif
                    }
                    else if (verbosity==6)
                    {
                        printf("\n\t\t\tNew loss, time %.8lf",c_time+s_time);
#ifdef DBG
                        fflush(stdout);
#endif
                    }
                    
                }
                
                // ***
                /// New step time </dl>
                c_time+=s_time;
                s_time=-log(gsl_rng_uniform_pos(seed))/(w_prob*avail_leaves);
                
            }
            // ****
            /// If the tree ends without leaves, it is discarted and the algorithm retries its construction.
            if (avail_leaves<2)
            {
                FreeSTree(out_tree);
                if (verbosity>2)
                {
                    printf("\n\t\t\tSpecies tree with less than 2 leaves or mrca smaller than expected, restart of the simulation of this tree. Try %d of %d",iter,MAX_IT);
                    if (verbosity>4)
                        printf("\n");
#ifdef DBG
                    fflush(stdout);
#endif
                }
                continue;
            }
            else
            {
                iter=0;
                break;
            }
            
        }
        
        // ****
        /// Addition of the last period without any event.
        for (j=0;j<avail_leaves;++j)
        {
            w_node1=*(node_ptrs+j);
            w_node1->time=time;
            w_node1->n_gen=time/gen_time;
            w_node1->gen_length=w_node1->n_gen-w_node1->anc_node->n_gen;
            w_node1->n_replicas=ind_per_sp;
            w_node1->sp_index=j+n_leaves+1-avail_leaves;
        }
        
        // ****
        /// Deletion of the dummie root of the tree if there is no outgroup addition, and reescaling and outgroup addition otherwise
        if (outgroup>0)
        {
            outgroup=time*outgroup; //From a deviation to the half of the tree heigth to a real internal branch length
            w_node1=NewSNodes(1,2);
            *((*out_tree)->root->children+1)=w_node1;
            w_node1->anc_node=(*out_tree)->root;
            (*out_tree)->root->n_gen=0;
            (*out_tree)->root->n_child=2;
            (*(*out_tree)->root->children)->gen_length=outgroup/gen_time;
            w_node1->time=time+outgroup;
            w_node1->n_gen=w_node1->time/gen_time;
            w_node1->gen_length=w_node1->n_gen;
            w_node1->n_replicas=1;
            
            AddConstantNgenSNodes((*(*out_tree)->root->children), outgroup/gen_time,gen_time);
  
            (*out_tree)->n_nodes=n_nodes+extra_nodes+2;
            (*out_tree)->n_leaves=avail_leaves+(extra_nodes/2)+1;
            (*out_tree)->n_gleaves=(avail_leaves+(extra_nodes/2))*ind_per_sp+1;
            (*out_tree)->locus_tree=NULL;
            (*out_tree)->gene_tree=NULL;

        }
        else
        {
            w_node1=*((*out_tree)->root->children);
            w_node1->anc_node=NULL;
            free((*out_tree)->root);
            (*out_tree)->root=w_node1;
            
            // ****
            /// Tree info fill
            (*out_tree)->n_nodes=n_nodes+extra_nodes;
            (*out_tree)->n_leaves=avail_leaves+(extra_nodes/2);
            (*out_tree)->n_gleaves=(avail_leaves+(extra_nodes/2)*ind_per_sp);
            (*out_tree)->locus_tree=NULL;
            (*out_tree)->gene_tree=NULL;
        }

        // *****
        /// Frees dynamic memory</dl></dl>
        free(node_ptrs);
        node_ptrs=NULL;
        
        if (iter>MAX_IT)
        {
            if (verbosity>3)
            {
                printf("\n\t\tMax iteration reached, it is imposible to simulate the species tree using the current parameters\n");
#ifdef DBG
                fflush(stdout);
#endif
            }
            return (LOOP_ERROR);
        }
    }
    switch (verbosity)
    {
        case 0:
        case 1:
        case 2:
            break;
        case 3:
            printf("\n\t%d-node species tree correctly simulated\n",(*out_tree)->n_nodes);
            break;
        default:
            printf("\n\t\t%d-node species tree correctly simulated: ",(*out_tree)->n_nodes);
            WriteSTree((*out_tree),(name_c *)NULL,out_time,labels);
            printf("\n");
            break;
    }

#ifdef DBG
        fflush(stdout);
#endif
    
    // ******
    /// Return
    return (NO_ERROR);
}

long int SimBDLTree(s_tree *wsp_tree,l_tree **wlocus_tree, l_node **node_ptrs, double b_rate,double d_rate,gsl_rng *seed, int min_lleaves, int min_lsleaves, int verbosity, int *st_losses, int *st_dups, int *st_leaves, int *st_gleaves)
{
    
    // *******
    /// <dl><dt>Declarations</dt><dd>
    
    // ******
    ///Structures
    l_node dummy_node,*aux_lnode=NULL,*w_lnode3=NULL,*anc_lnode=NULL,*w_lnode=NULL,*w_lnode2=NULL;
    s_node *w_snode=NULL;
    
    // ******
    /// B-D process related variables
    int n_leaves=0,lt_true_leaves=0,lt_diffs_true_leaves=0,diffs_true_leaves=0,tn_nodes=0,extra_nodes=0,n_losses=0;
    int next_paralog=0,node_index=0,avail_leaves=0;
    int n_nodes=0;
    double w_prob=d_rate+b_rate,b_prob=(b_rate/(d_rate+b_rate)),current_ngen=0,max_ngen=0,sampled_ngen=0, gen_time=wsp_tree->gen_time;
    
    // ******
    /// Loop related variables</dd></dl>
    int i=0,j=0,k=0,ltree_iter=0,bd_iter=0,l_tree_retries=0, maxnleaves_reached=0,done=0;
    
    while (done!=1 && ltree_iter<=MAX_IT)
    {
        // ******
        /// Locus tree re/initialization
        if (*wlocus_tree!=NULL)
            FreeLTree(wlocus_tree);
        *wlocus_tree=NewLTree(1, 0, 0, wsp_tree->max_children,wsp_tree->gen_time,wsp_tree->Ne,wsp_tree->mu); //Tree with only a root node and without g_node pointers
        
        // ******
        /// Species tree reset by \ref ResetSTreeSimL
        
        if(ResetSTreeSimL(wsp_tree)!= NO_ERROR)
            return MEM_ERROR;
        
        // ******
        /// Reseting of B-D process variables and linking the root of both trees
        wsp_tree->root->l_nodes=(*wlocus_tree)->root;
        wsp_tree->root->n_lnodes=1;
        (*wlocus_tree)->root->n_gen=0;
        (*wlocus_tree)->root->time=0;
        (*wlocus_tree)->root->conts=wsp_tree->root;
        next_paralog=1;
        n_leaves=0;
        lt_true_leaves=0;
        lt_diffs_true_leaves=0;
        tn_nodes=1;
        extra_nodes=0;
        n_losses=0;
        *st_leaves=0;
        *st_gleaves=0;
        *st_dups=0;
        *st_losses=0;
        ++ltree_iter;
        maxnleaves_reached=0;
        
        // ******
        /// <dl><dt>Pre-order s_node loop, to do a SSA (Simple sampling approach) algorithm along each branch of the species tree (GSA (General Sampling approach) is not required with fixed time)</dt><dd>
        for (i=1;i<wsp_tree->n_nodes;++i) //The root is not iterated.
        {
            
            if (verbosity>4)
            {
                printf("\n\t\tBranch %u of %u...",i,wsp_tree->n_nodes-1);
#ifdef DBG
                fflush(stdout);
#endif
            }
            
            w_snode=wsp_tree->m_node+i;
            
            if (w_snode->anc_node->n_lnodes==0) //There is no l_nodes in the ancestor of this s_node (lost lineage)
                continue;
            
            dummy_node.lat_node=w_snode->anc_node->l_nodes; //dummy node to allow the loop to start in the first iteration, due to the use of different structures (first iteration: s_node->l_nodes, rest: l_node->lat_node)
            w_lnode3=&dummy_node;
            diffs_true_leaves=0;
            bd_iter=0;
            
            // *****
            /// <dl><dt>Loop of associated \ref s_node::l_nodes "l_nodes" to the working s_node, i.e paralogs</dt><dd>
            
            
            for (j=0;j<w_snode->anc_node->n_lnodes;++j)
            {
                if (verbosity>4)
                {
                    printf("\n\t\t\tLocus tree node %u of %u...",j+1,w_snode->anc_node->n_lnodes);
#ifdef DBG
                    fflush(stdout);
#endif
                }
                
                // ****
                /// Waiting-time sampling
                w_lnode3=w_lnode3->lat_node;
                anc_lnode=w_lnode3;
                
                current_ngen=w_snode->anc_node->n_gen;
                max_ngen=w_snode->n_gen;
                sampled_ngen= -log(gsl_rng_uniform_pos(seed))/(w_prob);
                
                if (sampled_ngen>=max_ngen-current_ngen) //No events
                {
                    if (verbosity>4)
                    {
                        printf("\n\t\t\tThere is neither birth nor death event");
#ifdef DBG
                        fflush(stdout);
#endif
                    }
                    
                    // ****
                    /// New l-node allocation and configuration if there is neither birth nor death events along this branch
                    
                    w_lnode=NewLNodes(1, 0, wsp_tree->max_children); //New node
                    //Locus tree pointers
                    w_lnode->anc_node=anc_lnode;
                    *(anc_lnode->children+anc_lnode->n_child)=w_lnode;
                    //anc_lnode info
                    anc_lnode->n_child++;
                    
                    //S_tree/l_tree pointers
                    if (w_snode->n_lnodes==0)
                    {
                        w_snode->l_nodes=w_lnode;
                    }
                    else
                    {
                        w_lnode2=w_snode->l_nodes;
                        for (k=1; k<w_snode->n_lnodes; ++k)
                        {
                            w_lnode2=w_lnode2->lat_node;
                        }
                        w_lnode2->lat_node=w_lnode;
                    }
                    w_lnode->conts=w_snode;
                    
                    // Branches
                    w_lnode->n_gen=w_snode->n_gen;
                    w_lnode->gen_length=w_snode->gen_length;
                    w_lnode->time=w_snode->time;
                    //w_lnode info
                    w_lnode->sp_index=w_snode->sp_index;
                    w_lnode->Ne=w_snode->Ne;
                    w_lnode->mu_mult=w_snode->mu_mult;
                    w_lnode->gtime_mult=w_snode->gtime_mult;
                    w_lnode->paralog=w_lnode->anc_node->paralog;
                    w_lnode->n_nodes=w_snode->n_replicas;
                    //S_Node info
                    w_snode->n_lnodes++;
                    
                    ++tn_nodes;
                    if (w_snode->n_child==0)
                    {
                        ++(*st_leaves);
                        ++lt_true_leaves;
                        diffs_true_leaves|=1;
                        (*st_gleaves)+=w_lnode->n_nodes;
                    }
                    
                    
                }
                else
                {
                    if (verbosity>4)
                    {
                        printf("\n\t\t\t\tBirth-death process");
#ifdef DBG
                        fflush(stdout);
#endif
                    }
                    
                    // ****
                    /// <dl><dt>Else, Birth-death process along this branch</dt><dd>
                    
                    // ****
                    /// Iteration-dependent variable initialization
                    
                    w_lnode=NewLNodes(1,0,wsp_tree->max_children); //New node
                    
                    //Pointers
                    w_lnode->anc_node=anc_lnode;
                    *(anc_lnode->children+anc_lnode->n_child)=w_lnode;
                    anc_lnode->n_child++;
                    w_lnode->conts=w_snode;
                    //Info
                    w_lnode->paralog=w_lnode->anc_node->paralog;
                    
                    extra_nodes=0;
                    *(node_ptrs)=w_lnode;
                    n_leaves=1;
                    avail_leaves=1;
                    n_nodes=1;
                    anc_lnode=w_lnode;
                    n_losses=0;
                    bd_iter=0;
                    
                    // ****
                    /// <dl><dt>Birth-death process loop (while final time is not exceeded)</dt><dd>
                    while(current_ngen+sampled_ngen<max_ngen && avail_leaves!=0 && bd_iter<=MAX_IT && bd_iter<MAX_LEAVES-1)
                    {
                        ++bd_iter;
                        printf("iter %d\n",bd_iter);
                        fflush(stdout);
                        // ***
                        /// Chooses the leave for the next event
                        if(avail_leaves>1)
                            node_index=(int)gsl_rng_uniform_int(seed,avail_leaves);
                        else
                            node_index=0;
                        
                        // ***
                        /// <dl><dt>Chooses the process</dt><dd>
                        if (b_prob>gsl_rng_uniform(seed))
                        {
                            // **
                            /// <dl><dt>Birth</dt><dd>
                            
                            // *
                            /// New nodes allocation and reconfiguration of tree pointers/info</dd></dl>
                            w_lnode=NewLNodes(1,0,wsp_tree->max_children);
                            w_lnode2=NewLNodes(1,0,wsp_tree->max_children);
                            
                            //Tree pointers and branch lengths reconfiguration
                            anc_lnode=*(node_ptrs+node_index);
                            *(anc_lnode->children)=w_lnode;
                            *(anc_lnode->children+1)=w_lnode2;
                            w_lnode->anc_node=anc_lnode;
                            w_lnode2->anc_node=anc_lnode;
                            
                            anc_lnode->n_gen=current_ngen+sampled_ngen;
                            anc_lnode->gen_length=anc_lnode->n_gen-anc_lnode->anc_node->n_gen;
                            anc_lnode->time=anc_lnode->anc_node->time+anc_lnode->gen_length*w_snode->gtime_mult*gen_time;
                            
                            //Ancestor node info
                            anc_lnode->kind_node=DUP; //Duplication
                            anc_lnode->n_child=2;
                            anc_lnode->sp_index=w_snode->sp_index;
                            anc_lnode->Ne=w_snode->Ne;
                            anc_lnode->mu_mult=w_snode->mu_mult;
                            anc_lnode->gtime_mult=w_snode->gtime_mult;
                            anc_lnode->conts=w_snode;
                            
                            //New nodes info
                            w_lnode->paralog=anc_lnode->paralog;
                            w_lnode2->paralog=next_paralog;
                            ++next_paralog;
                            
                            //New leaves addition
                            *(node_ptrs+node_index)=w_lnode;
                            *(node_ptrs+avail_leaves)=w_lnode2;
                            ++avail_leaves;
                            ++n_leaves;
                            n_nodes+=2;
                            ++(*st_dups);
                            
                            if (verbosity==5)
                            {
                                printf("\n\t\t\tNew duplication");
#ifdef DBG
                                fflush(stdout);
#endif
                            }
                            else if (verbosity==6)
                            {
                                printf("\n\t\t\tNew duplication, generation %lf",anc_lnode->n_gen);
#ifdef DBG
                                fflush(stdout);
#endif
                                
                            }
                            
                        }
                        else
                        {
                            // **
                            /// <dl><dt>Death</dt><dd>
                            
                            // *
                            /// Set the selected node as non avaliable, filling its info.</dd></dl></dd></dl>
                            
                            w_lnode=*(node_ptrs+node_index);
                            anc_lnode=w_lnode->anc_node;
                            
                            //Reconfiguring branches
                            w_lnode->n_gen=current_ngen+sampled_ngen;
                            w_lnode->gen_length=w_lnode->n_gen-anc_lnode->n_gen;
                            w_lnode->time=w_lnode->anc_node->time+w_lnode->gen_length*w_snode->gtime_mult*gen_time;
                            
                            //Setting the leaf as non avaliable
                            *(node_ptrs+node_index)=*(node_ptrs+avail_leaves-1);
                            --avail_leaves;
                            --n_leaves;
                            
                            w_lnode->sp_index=w_snode->sp_index;
                            w_lnode->kind_node=LOSS;
                            w_lnode->Ne=w_snode->Ne;
                            w_lnode->mu_mult=w_snode->mu_mult;
                            w_lnode->gtime_mult=w_snode->gtime_mult;
                            w_lnode->n_nodes=0;
                            w_lnode->conts=w_snode;
                            n_nodes-=2; //2 true nodes lost (w_lnode an its ancestor) becoming extra_nodes (they have allocated memory, and are part of the l_tree, but they will not be represented in the g_tree)
                            extra_nodes+=2;
                            ++n_losses;
                            ++(*st_losses);
                            
                            if (verbosity==5)
                            {
                                printf("\n\t\t\t\t\tNew loss");
#ifdef DBG
                                fflush(stdout);
#endif
                            }
                            else if (verbosity==6)
                            {
                                printf("\n\t\t\t\t\tNew loss, generation %lf",w_lnode->n_gen);
#ifdef DBG
                                fflush(stdout);
#endif
                            }
                            
                            
                        }
                        
                        // ***
                        /// New step time </dd></dl>
                        current_ngen+=sampled_ngen;
                        sampled_ngen= -log(gsl_rng_uniform_pos(seed))/(w_prob*avail_leaves);
                    }
                    
                    if (bd_iter+1>=MAX_LEAVES)
                    {
                        maxnleaves_reached=1;
                        break;
                    }
                    
                    // ****
                    /// Addition of the last period (time between last event and the branch length) without any event and asociation of all the resulting l_nodes (tips) to this s_node, to allow the next step of the by branch SSA algorithm.</dd></dl></dd></dl></dd></dl>
                    
                    //l variable will be the number of l_nodes which has to be transversed through (l_node::lat_nodes) l=0 is different, as it is the first l_node pointed by s_node->l_nodes.
                    
                    //If there is no l_nodes associated to this s_node, the first one will be associated, and the next loop will start with l=1 (one association done). If there is any associated l_node, aux_lnode is goin to point to the last associated l_node, and the association loop will start with l=0. This strange set of code is due to the l_node::lat_node trick to maintain the memory efficiency.
                    if (avail_leaves!=0)
                    {
                        if (w_snode->n_lnodes==0)//Points the first l_node from s_node::l_node if there is no l_node associated to this s_node, and also points this node from aux_lnode to work with it in the association loop
                        {
                            w_lnode=*(node_ptrs);
                            w_lnode->n_gen=w_snode->n_gen;
                            w_lnode->gen_length=w_lnode->n_gen-w_lnode->anc_node->n_gen;
                            w_lnode->time=w_snode->time;
                            w_lnode->sp_index=w_snode->sp_index;
                            w_lnode->Ne=w_snode->Ne;
                            w_lnode->mu_mult=w_snode->mu_mult;
                            w_lnode->gtime_mult=w_snode->gtime_mult;
                            w_lnode->n_nodes=w_snode->n_replicas;
                            w_lnode->conts=w_snode;
                            w_snode->l_nodes=w_lnode;
                            aux_lnode=w_lnode;
                            ++w_snode->n_lnodes;
                            k=1;
                        }
                        else //Sets aux_lnode as the last associated l_node
                        {
                            aux_lnode=w_snode->l_nodes;
                            for (k=1;k<w_snode->n_lnodes;++k)
                            {
                                aux_lnode=aux_lnode->lat_node;
                            }
                            k=0;
                        }
                    }
                    
                    for (;k<avail_leaves;++k) //Association loop
                    {
                        w_lnode=*(node_ptrs+k);
                        w_lnode->n_gen=w_snode->n_gen;
                        w_lnode->gen_length=w_lnode->n_gen-w_lnode->anc_node->n_gen;
                        w_lnode->time=w_snode->time;
                        w_lnode->sp_index=w_snode->sp_index;
                        w_lnode->Ne=w_snode->Ne;
                        w_lnode->mu_mult=w_snode->mu_mult;
                        w_lnode->gtime_mult=w_snode->gtime_mult;
                        w_lnode->n_nodes=w_snode->n_replicas;
                        w_lnode->conts=w_snode;
                        aux_lnode->lat_node=w_lnode;
                        aux_lnode=w_lnode;
                        ++w_snode->n_lnodes;
                    }
                    
                    if (w_snode->n_child==0) //This s_node is a tip
                    {
                        (*st_leaves)+=n_leaves;
                        lt_true_leaves+=n_leaves;
                        diffs_true_leaves|=1;
                        (*st_gleaves)+=n_leaves*w_snode->n_replicas;
                    }
                    
                    (*st_leaves)+=n_losses;
                    tn_nodes+=n_nodes+extra_nodes; //Takes into account the number of nodes
                    
                }
                if (maxnleaves_reached==1)
                {
                    break;
                }
                
            }
            if (maxnleaves_reached==1)
            {
                break;
            }
            else
            {
                lt_diffs_true_leaves+=diffs_true_leaves;
                if (verbosity>4)
                {
                    printf("\n\t\t\tDone");
#ifdef DBG
                    fflush(stdout);
#endif
                }
            }
            
        }
        
        if (lt_true_leaves<min_lleaves || maxnleaves_reached==1 || lt_diffs_true_leaves<min_lsleaves)
        {
            if (verbosity>2)
            {
                if (maxnleaves_reached==1)
                {
                    printf("\n\t\tThe locus tree simulation reached the maximum number of locus tree lineages (%d) inside a species tree branch. After checking your simulation parameters carefully (mainly the duplication rate) you can sort this problem out setting the environmental variable SIMPHY_MAXLEAVES with a value higher than Try %d of %d",MAX_LEAVES,l_tree_retries,MAX_IT);
                }
                else if (lt_true_leaves<min_lleaves)
                {
                    printf("\n\t\tLocus tree with %u leaves, less than the minimum %u , restart of the simulation of this tree. Try %d of %d",lt_true_leaves,min_lleaves,l_tree_retries,MAX_IT);
                }
                else
                {
                    printf("\n\t\tLocus tree with %u leaves from different species, less than the minimum %u , restart of the simulation of this tree. Try %d of %d",lt_diffs_true_leaves,min_lsleaves,l_tree_retries,MAX_IT);
                }
                if (verbosity>4)
                {
                    printf("\n");
                }
#ifdef DBG
                fflush(stdout);
#endif
            }
            
            
            ++l_tree_retries;
        }
        else
        {
            done=1;
        }
    }
    
    
    if (ltree_iter>MAX_IT)
        return (LOOP_ERROR);
    else
    {
        // ******
        /// Filling information of the just simulated l_tree
        (*wlocus_tree)->n_leaves=*st_leaves;
        (*wlocus_tree)->n_nodes=tn_nodes;
        (*wlocus_tree)->n_gleaves=*st_gleaves;
        (*wlocus_tree)->species_tree=wsp_tree;
        wsp_tree->locus_tree=*wlocus_tree;
        ErrorReporter(CollapseLTree(*wlocus_tree,1,0,0), NULL);
        return (NO_ERROR);
    }
    
}

long int SimBDLHTree(s_tree *wsp_tree,l_tree **wlocus_tree, l_node **node_ptrs, double b_rate,double d_rate, double h_rate, double gc_rate, int t_kind, gsl_rng *seed, int min_lleaves, int min_lsleaves, int verbosity, int *st_losses, int *st_dups, int *st_transfr, int *st_gc, int *st_leaves, int *st_gleaves, int *st_ntrials)
{
    // *******
    /// <dl><dt>Declarations</dt><dd>
    
    // ******
    ///Structures
    l_node dummy_node,*aux_lnode=NULL,*anc_lnode=NULL,*w_lnode=NULL,*w_lnode2=NULL,*w_lnode3=NULL;
    l_node **avail_receptors=NULL;
    s_node *w_snode=NULL;
    period *periods=NULL, *w_period=NULL, *w_period2=NULL;
    
    // ******
    /// B-D + Poisson process related variables
    int n_leaves=0,slice_leaves=0,lt_true_leaves=0,lt_diffs_true_leaves=0,diffs_true_leaves=0,tn_nodes=0,extra_nodes=0,n_losses=0, n_periods=0,n_ltransf=0, n_lgc=0, n_avail_receptors=0, t_event=0, expected_nleaves=0;
    int next_paralog=0,node_index=0,avail_leaves=0,n_transfer=0,n_gc=0;
    int n_nodes=0;
    double w_prob=d_rate+b_rate+h_rate+gc_rate, dtgc_prob=((d_rate+h_rate+gc_rate)/(b_rate+d_rate+h_rate+gc_rate)), tgc_prob=((h_rate+gc_rate)/(b_rate+d_rate+h_rate+gc_rate)), gc_prob=(gc_rate/(b_rate+d_rate+h_rate+gc_rate)),current_ngen=0,max_ngen=0,sampled_ngen=0, rnumber=0, max_time=0, gen_time=wsp_tree->gen_time;
    
    // ******
    /// Loop related variables</dd></dl>
    int bd_iter=0,l_tree_retries=0, maxnleaves_reached=0,done=0, discarded_t=0;
    int i=0,j=0,k=0;
    
    while (done!=1 && *st_ntrials<=MAX_IT)
    {
        // ******
        /// Locus tree re/initialization
        if (*wlocus_tree!=NULL)
            FreeLTree(wlocus_tree);
        *wlocus_tree=NewLTree(1, 0, 0, wsp_tree->max_children,wsp_tree->gen_time,wsp_tree->Ne,wsp_tree->mu); //Tree with only a root node and without g_node pointers
        
        // ******
        /// Species tree reset by \ref ResetSTreeSimL
        
        if(ResetSTreeSimL(wsp_tree)!= NO_ERROR)
            return MEM_ERROR;
        expected_nleaves=ExpectedPrunedLtreeNleavesSNodes(wsp_tree,b_rate,d_rate);
        if (expected_nleaves>=MAX_LEAVES)
        {
            fprintf(stderr,"\n\nLocus tree birth rate is too hight for the rest of the settings, with an expected number of locus tree leaves bigger than MAX_LEAVES. You could increase this limit changing the environmental variable SIMPHY_MAXLEAVES\n");
            
#ifdef DBG
            fflush(stderr);
#endif
            return SETTINGS_ERROR;
        }
        
        // ******
        /// Reseting of B-D process variables and linking the root of both trees
        wsp_tree->root->l_nodes=(*wlocus_tree)->root;
        wsp_tree->root->n_lnodes=1;
        (*wlocus_tree)->root->n_gen=0;
        (*wlocus_tree)->root->time=0;
        (*wlocus_tree)->root->conts=wsp_tree->root;
        next_paralog=1;
        n_leaves=0;
        lt_true_leaves=0;
        lt_diffs_true_leaves=0;
        tn_nodes=1;
        extra_nodes=0;
        n_losses=0;
        n_ltransf=0;
        n_lgc=0;
        *st_leaves=0;
        *st_gleaves=0;
        *st_dups=0;
        *st_losses=0;
        *st_transfr=0;
        *st_gc=0;
        ++(*st_ntrials);
        maxnleaves_reached=0;
        max_time=0;
        
        // ******
        /// <dl><dt>Pre-order s_node loop, to do a SSA (Simple sampling approach) algorithm along each branch of the species tree (GSA (General Sampling approach) is not required with fixed time)</dt><dd>
        for (i=1;i<wsp_tree->n_nodes;++i) //The root is not iterated.
        {
            
            if (verbosity>4)
            {
                printf("\n\t\tBranch %u of %u...",i,wsp_tree->n_nodes-1);
#ifdef DBG
                fflush(stdout);
#endif
            }
            
            w_snode=wsp_tree->m_node+i;
            
            if (w_snode->anc_node->n_lnodes==0) //There is no l_nodes in the ancestor of this s_node (lost lineage)
                continue;
            
            dummy_node.lat_node=w_snode->anc_node->l_nodes; //dummy node to allow the loop to start in the first iteration, due to the use of different structures (first iteration: s_node->l_nodes, rest: l_node->lat_node)
            w_lnode3=&dummy_node;
            diffs_true_leaves=0;
            slice_leaves=0;
            
            // *****
            /// <dl><dt>Loop of associated \ref s_node::l_nodes "l_nodes" to the working s_node, i.e paralogs</dt><dd>
            
            
            for (j=0;j<w_snode->anc_node->n_lnodes;++j)
            {
                if (verbosity>4)
                {
                    printf("\n\t\t\tLocus tree node %u of %u...",j+1,w_snode->anc_node->n_lnodes);
#ifdef DBG
                    fflush(stdout);
#endif
                }

                // ****
                /// Waiting-time sampling
                w_lnode3=w_lnode3->lat_node;
                anc_lnode=w_lnode3;
                
                current_ngen=w_snode->anc_node->n_gen;
                max_ngen=w_snode->n_gen;
                sampled_ngen= -log(gsl_rng_uniform_pos(seed))/(w_prob);
                
                if (sampled_ngen>=max_ngen-current_ngen) //No events
                {
                    if (verbosity>4)
                    {
                        printf("\n\t\t\t\tThere are no events");
#ifdef DBG
                        fflush(stdout);
#endif
                    }
                    
                    ++slice_leaves;
                    // ****
                    /// New l-node allocation and configuration if there are no events along this branch
                    
                    w_lnode=NewLNodes(1, 0, wsp_tree->max_children); //New node
                    //Locus tree pointers
                    w_lnode->anc_node=anc_lnode;
                    *(anc_lnode->children+anc_lnode->n_child)=w_lnode;
                    //anc_lnode info
                    anc_lnode->n_child++;
                    
                    //S_tree/l_tree pointers
                    if (w_snode->n_lnodes==0)
                    {
                        w_snode->l_nodes=w_lnode;
                    }
                    else
                    {
                        w_lnode2=w_snode->l_nodes;
                        for (k=1; k<w_snode->n_lnodes; ++k)
                        {
                            w_lnode2=w_lnode2->lat_node;
                        }
                        w_lnode2->lat_node=w_lnode;
                    }
                    w_lnode->conts=w_snode;
                    
                    // Branches
                    w_lnode->n_gen=w_snode->n_gen;
                    w_lnode->gen_length=w_snode->gen_length;
                    
                    //w_lnode info
                    w_lnode->sp_index=w_snode->sp_index;
                    w_lnode->Ne=w_snode->Ne;
                    w_lnode->mu_mult=w_snode->mu_mult;
                    w_lnode->gtime_mult=w_snode->gtime_mult;
                    w_lnode->paralog=w_lnode->anc_node->paralog;
                    w_lnode->n_nodes=w_snode->n_replicas;
                    w_lnode->n_ilin=w_snode->n_replicas;
                    w_lnode->time=w_snode->time;
                    if (max_time<w_lnode->time)
                        max_time=w_lnode->time;
                    //S_Node info
                    w_snode->n_lnodes++;
                    
                    ++tn_nodes;
                    if (w_snode->n_child==0)
                    {
                        ++(*st_leaves);
                        ++lt_true_leaves;
                        diffs_true_leaves|=1;
                        (*st_gleaves)+=w_lnode->n_nodes;
                    }
                    
                    
                }
                else
                {
                    if (verbosity>4)
                    {
                        printf("\n\t\t\t\tBirth-death process");
#ifdef DBG
                        fflush(stdout);
#endif
                    }
                    
                    // ****
                    /// <dl><dt>Else, Birth-death process along this branch</dt><dd>
                    
                    // ****
                    /// Iteration-dependent variable initialization
                    
                    w_lnode=NewLNodes(1,0,wsp_tree->max_children); //New node
                    
                    //Pointers
                    w_lnode->anc_node=anc_lnode;
                    *(anc_lnode->children+anc_lnode->n_child)=w_lnode;
                    anc_lnode->n_child++;
                    w_lnode->conts=w_snode;
                    //Info
                    w_lnode->paralog=w_lnode->anc_node->paralog;
                    
                    extra_nodes=0;
                    *(node_ptrs)=w_lnode;
                    n_leaves=1;
                    avail_leaves=1;
                    n_nodes=1;
                    anc_lnode=w_lnode;
                    n_losses=0;
                    bd_iter=0;
                    
                    // ****
                    /// <dl><dt>Birth-death process loop (while final time is not exceeded)</dt><dd>
                    while(current_ngen+sampled_ngen<max_ngen && avail_leaves!=0 && bd_iter<=MAX_IT && slice_leaves<MAX_LEAVES-1)
                    {
                        ++bd_iter;
                        // ***
                        /// Chooses the leave for the next event
                        if(avail_leaves>1)
                            node_index=(int)gsl_rng_uniform_int(seed,avail_leaves);
                        else
                            node_index=0;
                        
                        // ***
                        /// <dl><dt>Chooses the process</dt><dd>
                        rnumber=gsl_rng_uniform(seed);
                        if (rnumber>dtgc_prob)
                        {
                            // **
                            /// <dl><dt>Birth</dt><dd>
                            
                            ++slice_leaves;
                            // *
                            /// New nodes allocation and reconfiguration of tree pointers/info</dd></dl>
                            w_lnode=NewLNodes(1,0,wsp_tree->max_children);
                            w_lnode2=NewLNodes(1,0,wsp_tree->max_children);
                            
                            //Tree pointers and branch lengths reconfiguration
                            anc_lnode=*(node_ptrs+node_index);
                            *(anc_lnode->children)=w_lnode;
                            *(anc_lnode->children+1)=w_lnode2;
                            w_lnode->anc_node=anc_lnode;
                            w_lnode2->anc_node=anc_lnode;
                            
                            anc_lnode->n_gen=current_ngen+sampled_ngen;
                            anc_lnode->gen_length=anc_lnode->n_gen-anc_lnode->anc_node->n_gen;
                            anc_lnode->time=anc_lnode->anc_node->time+anc_lnode->gen_length*w_snode->gtime_mult*gen_time;
                            //Ancestor node info
                            anc_lnode->kind_node=DUP; //Duplication
                            anc_lnode->n_child=2;
                            anc_lnode->sp_index=w_snode->sp_index;
                            anc_lnode->Ne=w_snode->Ne;
                            anc_lnode->mu_mult=w_snode->mu_mult;
                            anc_lnode->gtime_mult=w_snode->gtime_mult;
                            anc_lnode->conts=w_snode;
                            
                            //New nodes info
                            w_lnode->paralog=anc_lnode->paralog;
                            w_lnode2->paralog=next_paralog;
                            w_lnode2->fmax_nlin=1;
                            w_lnode2->n_olin=1;
                            
                            ++next_paralog;
                            
                            //New leaves addition
                            *(node_ptrs+node_index)=w_lnode;
                            *(node_ptrs+avail_leaves)=w_lnode2;
                            ++avail_leaves;
                            ++n_leaves;
                            n_nodes+=2;
                            ++(*st_dups);
                            
                            if (verbosity==5)
                            {
                                printf("\n\t\t\t\t\tNew duplication");
#ifdef DBG
                                fflush(stdout);
#endif
                            }
                            else if (verbosity==6)
                            {
                                printf("\n\t\t\t\t\tNew duplication, generation %lf",anc_lnode->n_gen);
#ifdef DBG
                                fflush(stdout);
#endif
                                
                            }
                        }
                        else if (rnumber>tgc_prob)
                        {
                            // **
                            /// <dl><dt>Death</dt><dd>
                            
                            // *
                            /// Set the selected node as non avaliable, filling its info.</dd></dl></dd></dl>
                            --slice_leaves;
                            w_lnode=*(node_ptrs+node_index);
                            anc_lnode=w_lnode->anc_node;
                            
                            //Reconfiguring branches
                            w_lnode->n_gen=current_ngen+sampled_ngen;
                            w_lnode->gen_length=w_lnode->n_gen-anc_lnode->n_gen;
                            w_lnode->time=w_lnode->anc_node->time+w_lnode->gen_length*w_snode->gtime_mult*gen_time;
                            
                            //Setting the leaf as non avaliable
                            *(node_ptrs+node_index)=*(node_ptrs+avail_leaves-1);
                            --avail_leaves;
                            --n_leaves;
                            
                            w_lnode->sp_index=w_snode->sp_index;
                            w_lnode->kind_node=LOSS;
                            w_lnode->Ne=w_snode->Ne;
                            w_lnode->mu_mult=w_snode->mu_mult;
                            w_lnode->gtime_mult=w_snode->gtime_mult;
                            w_lnode->n_nodes=0;
                            w_lnode->n_ilin=0;
                            w_lnode->n_olin=0;
                            w_lnode->conts=w_snode;
                            n_nodes-=2; //2 true nodes lost (w_lnode an its ancestor) becoming extra_nodes (they have allocated memory, and are part of the l_tree, but they will not be represented in the g_tree)
                            extra_nodes+=2;
                            ++n_losses;
                            ++(*st_losses);
                            
                            if (verbosity==5)
                            {
                                printf("\n\t\t\t\t\tNew loss");
#ifdef DBG
                                fflush(stdout);
#endif
                            }
                            else if (verbosity==6)
                            {
                                printf("\n\t\t\t\t\tNew loss, generation %lf",w_lnode->n_gen);
#ifdef DBG
                                fflush(stdout);
#endif
                            }
                        }
                        else if (rnumber>gc_prob)
                        {
                            // **
                            /// <dl><dt>Transfer (donnor)</dt><dd>
                            
                            // *
                            /// New nodes allocation and reconfiguration of tree pointers/info</dd></dl>
                            w_lnode=NewLNodes(1,0,wsp_tree->max_children);
                            
                            //Tree pointers and branch lengths reconfiguration
                            anc_lnode=*(node_ptrs+node_index);
                            *(anc_lnode->children)=w_lnode;
                            w_lnode->anc_node=anc_lnode;
                            
                            anc_lnode->n_gen=current_ngen+sampled_ngen;
                            anc_lnode->gen_length=anc_lnode->n_gen-anc_lnode->anc_node->n_gen;
                            anc_lnode->time=anc_lnode->anc_node->time+anc_lnode->gen_length*w_snode->gtime_mult*gen_time;
                            
                            //Ancestor node info
                            anc_lnode->kind_node=TRFR; //Transfer (donnor)
                            anc_lnode->n_child=1;
                            anc_lnode->sp_index=w_snode->sp_index;
                            anc_lnode->Ne=w_snode->Ne;
                            anc_lnode->mu_mult=w_snode->mu_mult;
                            anc_lnode->gtime_mult=w_snode->gtime_mult;
                            anc_lnode->conts=w_snode;
                            
                            //New nodes info
                            w_lnode->paralog=anc_lnode->paralog;
                            
                            //New leaves addition
                            *(node_ptrs+node_index)=w_lnode;
                            ++n_nodes;
                            ++(*st_transfr);
                            
                            if (verbosity==5)
                            {
                                printf("\n\t\t\t\t\tNew transfer donnor");
#ifdef DBG
                                fflush(stdout);
#endif
                            }
                            else if (verbosity==6)
                            {
                                printf("\n\t\t\t\t\tNew transfer donnor, generation %lf",anc_lnode->n_gen);
#ifdef DBG
                                fflush(stdout);
#endif
                                
                            }
                        }
                        else
                        {
                            // **
                            /// <dl><dt>Gene conversion (donnor)</dt><dd>
                            
                            // *
                            /// New nodes allocation and reconfiguration of tree pointers/info</dd></dl>
                            w_lnode=NewLNodes(1,0,wsp_tree->max_children);
                            
                            //Tree pointers and branch lengths reconfiguration
                            anc_lnode=*(node_ptrs+node_index);
                            *(anc_lnode->children)=w_lnode;
                            w_lnode->anc_node=anc_lnode;
                            
                            anc_lnode->n_gen=current_ngen+sampled_ngen;
                            anc_lnode->gen_length=anc_lnode->n_gen-anc_lnode->anc_node->n_gen;
                            anc_lnode->time=anc_lnode->anc_node->time+anc_lnode->gen_length*w_snode->gtime_mult*gen_time;
                            
                            //Ancestor node info
                            anc_lnode->kind_node=GC; //Gene conversion (donnor)
                            anc_lnode->n_child=1;
                            anc_lnode->sp_index=w_snode->sp_index;
                            anc_lnode->Ne=w_snode->Ne;
                            anc_lnode->mu_mult=w_snode->mu_mult;
                            anc_lnode->gtime_mult=w_snode->gtime_mult;
                            anc_lnode->conts=w_snode;
                            
                            //New nodes info
                            w_lnode->paralog=anc_lnode->paralog;
                            
                            //New leaves addition
                            *(node_ptrs+node_index)=w_lnode;
                            ++n_nodes;
                            ++(*st_gc);
                            
                            if (verbosity==5)
                            {
                                printf("\n\t\t\t\t\tNew gene conversion donnor");
#ifdef DBG
                                fflush(stdout);
#endif
                            }
                            else if (verbosity==6)
                            {
                                printf("\n\t\t\t\t\tNew gene conversion donnor, generation %lf",anc_lnode->n_gen);
#ifdef DBG
                                fflush(stdout);
#endif
                                
                            }
                        }
                        
                        // ***
                        /// New step time </dd></dl>
                        current_ngen+=sampled_ngen;
                        sampled_ngen= -log(gsl_rng_uniform_pos(seed))/(w_prob*avail_leaves);
                    }
                    
                    if (slice_leaves+1>=MAX_LEAVES)
                    {
                        maxnleaves_reached=1;
                        break;
                    }
                    
                    // ****
                    /// Addition of the last period (time between last event and the branch length) without any event and asociation of all the resulting l_nodes (tips) to this s_node, to allow the next step of the by branch SSA algorithm.</dd></dl></dd></dl></dd></dl>
                    
                    //l variable will be the number of l_nodes which has to be transversed through (l_node::lat_nodes) l=0 is different, as it is the first l_node pointed by s_node->l_nodes.
                    
                    //If there is no l_nodes associated to this s_node, the first one will be associated, and the next loop will start with l=1 (one association done). If there is any associated l_node, aux_lnode is goin to point to the last associated l_node, and the association loop will start with l=0. This strange set of code is due to the l_node::lat_node trick to maintain the memory efficiency.
                    if (avail_leaves!=0)
                    {
                        if (w_snode->n_lnodes==0)//Points the first l_node from s_node::l_node if there is no l_node associated to this s_node, and also points this node from aux_lnode to work with it in the association loop
                        {
                            w_lnode=*(node_ptrs);
                            w_lnode->n_gen=w_snode->n_gen;
                            w_lnode->gen_length=w_lnode->n_gen-w_lnode->anc_node->n_gen;
                            w_lnode->time=w_snode->time;
                            w_lnode->sp_index=w_snode->sp_index;
                            w_lnode->Ne=w_snode->Ne;
                            w_lnode->mu_mult=w_snode->mu_mult;
                            w_lnode->gtime_mult=w_snode->gtime_mult;
                            w_lnode->n_nodes=w_snode->n_replicas;
                            w_lnode->n_ilin=w_lnode->n_nodes;
                            w_lnode->conts=w_snode;
                            if (max_time<w_lnode->time)
                                max_time=w_lnode->time;
                            w_snode->l_nodes=w_lnode;
                            aux_lnode=w_lnode;
                            ++w_snode->n_lnodes;
                            k=1;
                        }
                        else //Sets aux_lnode as the last associated l_node
                        {
                            aux_lnode=w_snode->l_nodes;
                            for (k=1;k<w_snode->n_lnodes;++k)
                            {
                                aux_lnode=aux_lnode->lat_node;
                            }
                            k=0;
                        }
                    }
                    
                    for (;k<avail_leaves;++k) //Association loop
                    {
                        w_lnode=*(node_ptrs+k);
                        w_lnode->n_gen=w_snode->n_gen;
                        w_lnode->gen_length=w_lnode->n_gen-w_lnode->anc_node->n_gen;
                        w_lnode->time=w_snode->time;
                        w_lnode->sp_index=w_snode->sp_index;
                        w_lnode->Ne=w_snode->Ne;
                        w_lnode->mu_mult=w_snode->mu_mult;
                        w_lnode->gtime_mult=w_snode->gtime_mult;
                        w_lnode->n_nodes=w_snode->n_replicas;
                        w_lnode->n_ilin=w_lnode->n_nodes;
                        w_lnode->conts=w_snode;
                        if (max_time<w_lnode->time)
                            max_time=w_lnode->time;
                        aux_lnode->lat_node=w_lnode;
                        aux_lnode=w_lnode;
                        ++w_snode->n_lnodes;
                    }
                    
                    if (w_snode->n_child==0) //This s_node is a tip
                    {
                        (*st_leaves)+=n_leaves;
                        lt_true_leaves+=n_leaves;
                        diffs_true_leaves|=1;
                        (*st_gleaves)+=n_leaves*w_snode->n_replicas;
                    }
                    
                    (*st_leaves)+=n_losses;
                    tn_nodes+=n_nodes+extra_nodes; //Takes into account the number of nodes
                    
                }
                if (maxnleaves_reached==1)
                {
                    break;
                }
                
            }
            if (maxnleaves_reached==1)
            {
                break;
            }
            else
            {
                lt_diffs_true_leaves+=diffs_true_leaves;
                if (verbosity>4)
                {
                    printf("\n\t\tDone");
#ifdef DBG
                    fflush(stdout);
#endif
                }
            }
            
        }
        
        if (lt_true_leaves<min_lleaves || maxnleaves_reached==1 || lt_diffs_true_leaves<min_lsleaves || *st_gleaves<=1)
        {
            if (verbosity>2)
            {
                if (maxnleaves_reached==1)
                {
                    printf("\n\t\tThe locus tree simulation reached the maximum number of locus tree lineages (%d) inside a species tree branch. After checking your simulation parameters carefully (mainly the duplication rate) you can sort this problem out setting the environmental variable SIMPHY_MAXLEAVES with a value higher than %d. Try %d of %d",MAX_LEAVES,MAX_LEAVES,l_tree_retries,MAX_IT);
                }
                else if (lt_true_leaves<min_lleaves)
                {
                    printf("\n\t\tLocus tree with %u leaves, less than the minimum %u , restart of the simulation of this tree. Try %d of %d",lt_true_leaves,min_lleaves,l_tree_retries,MAX_IT);
                }
                else if (lt_diffs_true_leaves<min_lsleaves)
                {
                    printf("\n\t\tLocus tree with %u leaves from different species, less than the minimum %u , restart of the simulation of this tree. Try %d of %d",lt_diffs_true_leaves,min_lsleaves,l_tree_retries,MAX_IT);
                }
                else
                    printf("\n\t\tLocus tree with %u associated gene tree leaves, less than 2, restart of the simulation of this tree. Try %d of %d",*st_gleaves,l_tree_retries,MAX_IT);
                    
                if (verbosity>4)
                {
                    printf("\n");
                }
#ifdef DBG
                fflush(stdout);
#endif
            }
            
            
            ++l_tree_retries;
        }
        else
        {
            
            done=1;
        }
    }
    
    
    if (*st_ntrials>MAX_IT)
        return (LOOP_ERROR);
    else
    {
        
        (*wlocus_tree)->n_leaves=*st_leaves;
        (*wlocus_tree)->n_nodes=tn_nodes;
        (*wlocus_tree)->n_gleaves=*st_gleaves;
        (*wlocus_tree)->species_tree=wsp_tree;
        
        if (*st_transfr+*st_gc>0)
        {
            
            if (verbosity>4)
            {
                printf("\n\n\t\tCompleting transfers and/or gene conversions... ");
#ifdef DBG
                fflush(stdout);
#endif
            }
 
            // Tree realocation (adding necessary memory)
            ErrorReporter(CollapseResizeLTree(*wlocus_tree,tn_nodes+*st_transfr+*st_gc,*st_leaves,1,0,0),NULL); //The root of this tree is going to be badly set due to the extra still not used nodes.
            
            //Variable initialization
            n_periods=tn_nodes-*st_leaves+*st_losses+1; //Maximum number of periods of an ultrametric tree (without the root) =Internal nodes + losses (tip_dates). I add a dummy one, with r_bound==0 to avoid some pointer problems
            
            //Here I'm not using NewPeriods to perform the initialization in a more efficient way, saving one extra loop
            periods=calloc(n_periods, sizeof(struct period));
            if (periods==NULL)
                return MEM_ERROR;
            periods->r_bound=max_time;
            periods->l_nodes=NULL;
            periods->n_lnodes=0;
            j=1;
            
            avail_receptors=calloc(*st_leaves,sizeof(l_node *));
            if (avail_receptors==NULL)
                return MEM_ERROR;
            
            if (verbosity>4)
            {
                printf("\n\t\t\tInitializing periods... ");
#ifdef DBG
                fflush(stdout);
#endif
            }

            for (i=0; i<tn_nodes; ++i)
            {
                w_lnode=((*wlocus_tree)->m_node+i);
                if (w_lnode->kind_node==SP && w_lnode->n_child==0)
                    continue;
                w_period=periods+j;
                w_period->l_nodes=calloc(*st_leaves, sizeof(l_node *)); //The maximum number of lineages in a period is the number of leaves (star tree)
                if (w_period->l_nodes==NULL)
                    return MEM_ERROR;
                w_period->n_lnodes=w_lnode->n_child;
                for (k=0; k<w_lnode->n_child; ++k)
                {
                    *(w_period->l_nodes+k)=*(w_lnode->children+k);
                }
                w_period->r_bound=w_lnode->time;
                ++j;
            }
            if (verbosity>4)
            {
                printf("Done");
#ifdef DBG
                fflush(stdout);
#endif
            }
            
            if (verbosity>4)
            {
                printf("\n\t\t\tOrdering periods... ");
#ifdef DBG
                fflush(stdout);
#endif
            }
            
            qsort(periods, n_periods, sizeof(struct period), Compare_periods);
            
            if (verbosity>4)
            {
                printf("Done");
#ifdef DBG
                fflush(stdout);
#endif
            }
            
            if (verbosity>4)
            {
                printf("\n\t\t\tTraversing periods and performing transferences... ");
#ifdef DBG
                fflush(stdout);
#endif
            }
            for (i=0; i<n_periods; ++i)
            {
                w_period=periods+i;
                w_period2=periods+i+1;
                n_avail_receptors=0;
                discarded_t=0;
                w_lnode2=NULL;
                
                if (verbosity>5)
                {
                    printf("\n\t\t\t\tPeriod %d, lower bound %lf, n_lineages %d",i, w_period->r_bound, w_period->n_lnodes);
#ifdef DBG
                    fflush(stdout);
#endif
                }
                
                for (j=0; j<w_period->n_lnodes; ++j)
                {
                    w_lnode=*(w_period->l_nodes+j);

                    if (w_lnode->time>w_period2->r_bound)
                    {
                        if (verbosity>5)
                        {
                            printf("\n\t\t\t\t\tCopying node %d to the previous period",w_lnode->index);
#ifdef DBG
                            fflush(stdout);
#endif
                        }
                        *(w_period2->l_nodes+w_period2->n_lnodes)=w_lnode;
                        w_period2->n_lnodes+=1;
                    }
                    else if (w_period->n_lnodes>1 && ((w_lnode->kind_node==TRFR || w_lnode->kind_node==GC) && w_lnode->time==w_period2->r_bound))
                    {
                        w_snode=w_lnode->conts;
                        w_lnode2=w_lnode;
                        t_event=w_lnode->kind_node;
                    }
                    else if (w_period->n_lnodes==1 && ((w_lnode->kind_node==TRFR || w_lnode->kind_node==GC) && w_lnode->time==w_period2->r_bound))
                    {
                        w_lnode2=w_lnode;
                        discarded_t=1;
                    }
                }
                
                if (w_lnode2!=NULL)
                {
                    switch (verbosity)
                    {
                        case 0:
                        case 1:
                        case 2:
                        case 3:
                        case 4:
                        case 5:
                            break;
                        default:
                            printf("\n\t\t\t\t\t%s candidate node %d:\n\t\t\t\t\t\tCandidate receptors:",w_lnode2->kind_node==GC?"Gene conversion":"Transfer",w_lnode2->index);
#ifdef DBG
                            fflush(stdout);
#endif
                        break;
                    }
                    
                    switch (discarded_t)
                    {
                        case 0:
                            for (j=0; j<w_period->n_lnodes; ++j)
                            {
                                w_lnode=*(w_period->l_nodes+j);
                                switch (t_event)
                                {
                                    case TRFR:
                                        if (w_lnode->conts!=w_snode)
                                        {
                                            *(avail_receptors+n_avail_receptors)=w_lnode;
                                            ++n_avail_receptors;
                                            switch (verbosity)
                                            {
                                                case 0:
                                                case 1:
                                                case 2:
                                                case 3:
                                                case 4:
                                                case 5:
                                                    break;
                                                default:
                                                    printf(" %u,", w_lnode->index);
#ifdef DBG
                                                    fflush(stdout);
#endif
                                                    break;
                                            }
                                        }
                                        break;
                                    case GC:
                                        if (w_lnode->conts==w_snode && w_lnode!=w_lnode2) //w_lnode2=GC node
                                        {
                                            *(avail_receptors+n_avail_receptors)=w_lnode;
                                            ++n_avail_receptors;
                                            switch (verbosity)
                                            {
                                                case 0:
                                                case 1:
                                                case 2:
                                                case 3:
                                                case 4:
                                                case 5:
                                                    break;
                                                default:
                                                    printf(" %u,", w_lnode->index);
#ifdef DBG
                                                    fflush(stdout);
#endif
                                                    break;
                                            }
                                        }
                                        break;
                                    default:
                                        ErrorReporter(UNEXPECTED_VALUE,NULL);
                                        break;
                                }
                            }
                            if (n_avail_receptors>0)
                            {
                                switch (verbosity)
                                {
                                    case 0:
                                    case 1:
                                    case 2:
                                    case 3:
                                    case 4:
                                    case 5:
                                        break;
                                    default:
                                        printf("\n\t\t\t\t\t\t%s selection between %u candidates", t_kind==1?"Inversely proportional to MRCA":"Random" ,n_avail_receptors);
#ifdef DBG
                                        fflush(stdout);
#endif
                                        break;
                                }
                                switch (t_kind)
                                {
                                    case 1:
                                        w_lnode=w_lnode2;
                                        w_lnode2=ChooseLNodePeriod(avail_receptors,n_avail_receptors,w_lnode2,gsl_rng_uniform_pos(seed),verbosity);
                                        break;
                                    case 0:
                                        w_lnode=w_lnode2;
                                        w_lnode2=*(avail_receptors+gsl_rng_uniform_int(seed,n_avail_receptors));
                                        break;
                                    default:
                                        ErrorReporter(SETTINGS_ERROR,NULL);
                                        break;
                                        
                                }

                                //w_lnode -> donnor
                                //w_lnode2 -> receptor after reception
                                //w_lnode3 -> receptor (lost)
                                
                                w_lnode3=(*wlocus_tree)->m_node+tn_nodes+n_transfer+n_gc; //This will be the receptor
                                anc_lnode=w_lnode2->anc_node;
                                w_lnode3->anc_node=anc_lnode;
                                w_lnode3->paralog=w_lnode2->paralog;
                                w_lnode3->conts=w_lnode2->conts;
                                w_lnode3->mu_mult=w_lnode2->mu_mult;
                                w_lnode3->gtime_mult=w_lnode2->gtime_mult;
                                w_lnode3->Ne=w_lnode2->Ne;
                                w_lnode3->time=w_lnode->time;
                                w_lnode3->n_gen=anc_lnode->n_gen+ (w_lnode3->time-anc_lnode->time)/w_lnode3->gtime_mult/gen_time; //gen_lenth=Time_length * (1/gen_time*gtime_mult)
                                w_lnode3->gen_length=w_lnode3->n_gen-anc_lnode->n_gen;
                                w_lnode3->sp_index=w_lnode3->conts->sp_index;
                                //w_lnode3->n_child=0; Implicit.
                                //w_lnode3->n_nodes=0;
                                //w_lnode3->fmax_nlin=0;
                                //w_lnode3->n_ilin=0;
                                //w_lnode3->n_olin=0;
                                
                                for (k=0; k<anc_lnode->n_child; ++k)
                                {
                                    if (*(anc_lnode->children+k)==w_lnode2)
                                        *(anc_lnode->children+k)=w_lnode3;
                                }
                                
                                w_lnode2->gen_length=w_lnode2->n_gen-w_lnode3->n_gen;
                                w_lnode2->anc_node=w_lnode;
                                w_lnode2->fmax_nlin=1;
                                w_lnode2->n_olin=1;
                                
                                
                                *(w_lnode->children+1)=w_lnode2; //The second node in a TRFR node is allways the transfered lineage (bound!!).
                                w_lnode->n_child=2;
                                
                                switch (t_event)
                                {
                                    case TRFR:
                                        w_lnode3->kind_node=RTRFR;
                                        ++n_transfer;
                                        switch (verbosity)
                                        {
                                            case 0:
                                            case 1:
                                            case 2:
                                            case 3:
                                            case 4:
                                                break;
                                            case 5:
                                                printf("\n\t\t\t\t\t Transfer");
                                                break;
                                            default:
                                                printf("\n\t\t\t\t\tTransfer from %d to %d, time %lf",w_lnode->index, w_lnode2->index, w_period2->r_bound);
                                                break;
                                        }
                                        
#ifdef DBG
                                        fflush(stdout);
#endif
                                        break;
                                    case GC:
                                        w_lnode3->kind_node=RGC;
                                        ++n_gc;
                                        switch (verbosity)
                                        {
                                            case 0:
                                            case 1:
                                            case 2:
                                            case 3:
                                            case 4:
                                                break;
                                            case 5:
                                                printf("\n\t\t\t\t\t Gene conversion");
                                                break;
                                            default:
                                                printf("\n\t\t\t\t\tGene conversion from %d to %d, time %lf",w_lnode->index, w_lnode2->index, w_period2->r_bound);
                                                break;
                                        }
                                        
#ifdef DBG
                                        fflush(stdout);
#endif
                                        break;
                                        
                                    default:
                                        ErrorReporter(UNEXPECTED_VALUE,NULL);
                                        break;
                                }
                                
                                break;
                            }
                            //else go to default
                            
                        default:
                            switch (verbosity)
                            {
                                case 0:
                                case 1:
                                case 2:
                                case 3:
                                case 4:
                                    break;
                                case 5:
                                    printf("\n\t\t\t\t\t\t Discarded %s",w_lnode2->kind_node==TRFR?"transfer":"gene conversion");
                                    break;
                                default:
                                    printf("\n\t\t\t\t\t\tImpossible to perform the sampled %s from the node %u, as there is no candidate receptors. Discarded transfer\n",w_lnode2->kind_node==TRFR?"transfer":"gene conversion",w_lnode2->index);
                                    break;
                            }
                            
    #ifdef DBG
                            fflush(stdout);
    #endif
                            //w_lnode2 -> candidate transfer node
                            
                            (*(w_lnode2->children))->anc_node=w_lnode2->anc_node; //connecting daughter with grandmather
                            for(k=0;k<w_lnode2->anc_node->n_child;++k)
                            {
                                if (*(w_lnode2->anc_node->children+k)==w_lnode2)
                                {
                                    *(w_lnode2->anc_node->children+k)=*(w_lnode2->children); //connecting grandmather with daughter
                                    break;
                                }
                            }
                            (*(w_lnode2->children))->gen_length=(*(w_lnode2->children))->n_gen-w_lnode2->anc_node->n_gen;
                            if (w_lnode2->kind_node==TRFR)
                                ++n_ltransf;
                            else
                                ++n_lgc;
                            
                            break;
                    }
                }
                
            }
            if (verbosity>4)
            {
                printf("\n\t\t\tDone");
#ifdef DBG
                fflush(stdout);
#endif
            }

            FreePeriods(periods, n_periods);
            free(avail_receptors);
            
            (*wlocus_tree)->root=((*wlocus_tree)->m_node+tn_nodes-1);
        }
        // ******
        /// Filling information of the just simulated l_tree
        *st_transfr-=n_ltransf;
        *st_gc-=n_lgc;
        wsp_tree->locus_tree=*wlocus_tree;
        
        if (verbosity>4)
        {
            printf("\n\t\tDone");
#ifdef DBG
            fflush(stdout);
#endif
        }
        
        if (*st_dups>0 || *st_transfr>0 || *st_gc>0)
        {
            ErrorReporter(CollapseResizeLTree(*wlocus_tree,tn_nodes+*st_transfr-n_ltransf+*st_gc-n_lgc,*st_leaves+*st_transfr+*st_gc,1,0,1), NULL);
            if (verbosity>4)
            {
                printf("\n\t\tCalculating lineage count probabilities of entering and leaving each locus tree branch... ");
#ifdef DBG
                fflush(stdout);
#endif
            }
            
            // ****
            /// Calculation of the probability of each possible number of lineages entering and exiting every locus tree branch
            
            CalcProbsNLineagesLTree((*wlocus_tree)->root, (*wlocus_tree)->Ne, 0, verbosity);//I should try to extend it for using polytomies
            
            if (verbosity>4)
            {
                printf("\n\t\tDone\n");
#ifdef DBG
                fflush(stdout);
#endif
            }
        }
        else
            ErrorReporter(CollapseResizeLTree(*wlocus_tree,tn_nodes+*st_transfr-n_ltransf+*st_gc-n_lgc,*st_leaves+*st_transfr+*st_gc,1,0,0), NULL);
        
        return (NO_ERROR);
    }
    
}

long int SimMSCGTree(l_tree *wlocus_tree, g_tree **gene_tree, name_c * names, float epsilon_brent, gsl_rng *seed, int *n_lcoals,int simlosses, int verbosity, double gen_time, int collapse)
{
    // *******
    /// <dl><dt>Variable declaration</dt><dd>
    
    // ******
    /// Structures
    l_node *w_lnode=NULL, *anc_lnode=NULL;
    g_node *w_gnodes=NULL,*off1=NULL,*off2=NULL,*anc_gnode=NULL;
    g_node **w_gnodes_ptr=NULL;
    
    // ******
    /// IO related
    char  *iobuffer=NULL;
    
    // ******
    /// Values
    int next_avail_inode=0,max_coals=0,n_anc_nodes=0,p_Ne=0,k=0;
    double p_mu=0;
    double sampled_ngen=0, min_ngen=0, current_ngen=0;
    int avail_leaves=0,node_index=0;
    
    // ******
    /// Loop related variables</dd></dl>
    int i=0,j=0;
    
    // ******
    /// Initialization of other general variables</dd></dl>
    if (verbosity>4)
    {
        iobuffer=calloc(count_intdigits((long)UINT_MAX, 0), sizeof(char));
        if (iobuffer==NULL)
            return MEM_ERROR;
    }
    
    // ****
    /// Association between locus tree and gene tree and gene tree reset
    if (verbosity>4)
        printf("\n\t\t\tAssociation and initialization of gene tree and locus tree structures...");
#ifdef DBG
    fflush(stdout);
#endif
    MatchTreesMSC(wlocus_tree,*gene_tree,1,simlosses);
    if (verbosity>4)
        printf("\n\t\t\tDone\n\t\t\tSampling the coalescent history along each branch of the locus tree...");
#ifdef DBG
    fflush(stdout);
#endif
    
    w_gnodes=(*gene_tree)->m_node;
    next_avail_inode=wlocus_tree->n_gleaves;
    *n_lcoals=0;
    
    // ****
    /// <dl><dt>Post-order iteration loop over l_nodes to perform a coalescent simulation along each l_tree branch</dt><dd>
    
    for(i=0; i<wlocus_tree->n_nodes; ++i) // Nodes has been pre-ordered in a post-order.
    {
        // ***
        /// Reset of iteration-dependent variables related with the current l_node
        sampled_ngen=0.0;
        node_index=0;
        n_anc_nodes=0;
        anc_lnode=NULL;
        off1=NULL;
        off2=NULL;
        
        // ** Taking into account a new locus tree node ** //
        w_lnode=(wlocus_tree->m_node+i);
        
        if (w_lnode->n_nodes==0)
            continue;
        
        anc_lnode=w_lnode->anc_node;
        k=w_lnode->n_nodes;
        avail_leaves=k;
        w_gnodes_ptr=w_lnode->g_nodes;
        
        if(w_lnode->Ne!=0) //Node with a private Ne.
            p_Ne=w_lnode->Ne;
        else
            p_Ne=wlocus_tree->Ne; //Global Ne.
        
        p_mu=wlocus_tree->mu*w_lnode->mu_mult;
        
        // ** Constraining coalescent time and number of events ** //
        current_ngen=w_lnode->n_gen;
        min_ngen=current_ngen-w_lnode->gen_length;
        max_coals=w_lnode->n_nodes-1;//Max num of coalescent events = number of nodes to coalesce -1
        
        if (verbosity==5)
        {
            if (w_lnode->conts==NULL) //There is no species tree
            {
                if(names!=NULL)
                    printf("\n\t\tLocus tree node %s_%d... ",(w_lnode->kind_node==SP && w_lnode->n_child==0)?(names->names+(w_lnode->sp_index*names->max_lname)):"Internal node",w_lnode->paralog);
                else
                {
                    sprintf(iobuffer, "%d",w_lnode->sp_index);
                    printf("\n\t\tLocus tree node %s_%d... ",(w_lnode->kind_node==SP && w_lnode->n_child==0)?iobuffer:"Internal node",w_lnode->paralog);
                }
            }
            else
            {
                if(names!=NULL)
                    printf("\n\t\tSpecies tree node %s, locus tree node %s_%d... ",(names->names+(w_lnode->conts->sp_index*names->max_lname)),(w_lnode->kind_node==SP && w_lnode->n_child==0)?(names->names+(w_lnode->sp_index*names->max_lname)):"Internal node",w_lnode->paralog);
                else
                {
                    sprintf(iobuffer, "%d",w_lnode->conts->sp_index);
                    printf("\n\t\tSpecies tree node %s, locus tree node %s_%d... ",w_lnode->conts->sp_index==0?"Internal node":iobuffer,(w_lnode->kind_node==SP && w_lnode->n_child==0)?iobuffer:"Internal node",w_lnode->paralog);
                }
            }
            
            
#ifdef DBG
            fflush(stdout);
#endif
        }
        else if (verbosity==6)
        {
            if (w_lnode->conts==NULL)
            {
                if(names!=NULL)
                    printf("\n\t\tLocus tree node %s_%d (index %d, sp_index %d, n_nodes %d) ",(w_lnode->kind_node==SP && w_lnode->n_child==0)?(names->names+(w_lnode->sp_index*names->max_lname)):"Internal node",w_lnode->paralog,w_lnode->index, w_lnode->sp_index, k);
                else
                {
                    sprintf(iobuffer, "%d",w_lnode->sp_index);
                    printf("\n\t\tLocus tree node %s_%d (index %d, sp_index %d, n_nodes %d) ",(w_lnode->kind_node==SP && w_lnode->n_child==0)?iobuffer:"Internal node",w_lnode->paralog,w_lnode->index, w_lnode->sp_index, k);
                }
            }
            else
            {
                if(names!=NULL)
                    printf("\n\t\tSpecies tree node %s, locus tree node %s_%d (index %d, sp_index %d, n_nodes %d) ",(names->names+(w_lnode->conts->sp_index*names->max_lname)),(w_lnode->kind_node==SP && w_lnode->n_child==0)?(names->names+(w_lnode->sp_index*names->max_lname)):"Internal node",w_lnode->paralog,w_lnode->index, w_lnode->sp_index, k);
                else
                {
                    sprintf(iobuffer, "%d",w_lnode->conts->sp_index);
                    printf("\n\t\tSpecies tree node %s, locus tree node %s_%d (index %d, sp_index %d, n_nodes %d) ",w_lnode->conts->sp_index==0?"Internal node":iobuffer,(w_lnode->kind_node==SP && w_lnode->n_child==0)?iobuffer:"Internal node",w_lnode->paralog,w_lnode->index, w_lnode->sp_index, k);
                }
            }
            
#ifdef DBG
            fflush(stdout);
#endif
        }
        
        // ***
        /// <dl><dt> Coalescent simulation loop. Each iteration means one coalescence. Ends if the end of the branch has reached or there is no more posible coalescences (only 1 active g_node left)</dt><dd>
        
        for (j=0; j<max_coals;++j)
        {
            // **
            /// Waiting time (coalescent time) sample. Regular coalescent process:  \f$ Time\sim exp({k\choose 2}) \f$.

            sampled_ngen=((-log(gsl_rng_uniform_pos(seed)))*(2*p_Ne))/(k*(k-1));// - (1/lambda) * ln uniform random variable. Lambda = parameter (k over 2). K, number of nodes avaliable to coalesce.
            
            if(sampled_ngen<=0)
                return UNEXPECTED_VALUE;
            
            current_ngen-=sampled_ngen;
            --k;
            
            //End of the branch reached. There is no restriction if anc_resnode==NULL (root)
            if (current_ngen<min_ngen && anc_lnode!=NULL)
            {
                
                (*n_lcoals)+=(avail_leaves-1);
                
                if (verbosity==5)
                {
                    printf("\n\t\t\t %d extra lineages going deeper",avail_leaves-1);
#ifdef DBG
                    fflush(stdout);
#endif
                }
                if (verbosity==6)
                {
                    printf("\n\t\t\t %d extra lineages going deeper, min coalescent age %f, proposed coalescent time %f",avail_leaves-1,min_ngen,current_ngen);
#ifdef DBG
                    fflush(stdout);
#endif
                }
                break;
                
            }

            // **
            /// Randomly choice of offspring nodes, ancestral node creation and configuration </dd></dl>
            
            //Randomly choice of the first offspring node
            node_index= gsl_rng_uniform(seed)*avail_leaves; //Samples nodes
            off1= *(w_gnodes_ptr + node_index);
            --avail_leaves; //Uses one leaf
            *(w_gnodes_ptr + node_index)= *(w_gnodes_ptr + avail_leaves);// The pointer to this node now points to the last w_gnodes_ptr, inaccesible due to --avail_leaves
            
            //Randomly choice of the second offspring node
            if(avail_leaves>1)
                node_index= gsl_rng_uniform(seed)*avail_leaves;//Samples nodes
            else
                node_index=0;
            
            off2= *(w_gnodes_ptr + node_index);
            *(w_gnodes_ptr + node_index)= w_gnodes + next_avail_inode; //Now the restriction that pointed to off2 points to the new ancestral node, for using it in the next loop iteration (and/or next coalescent process)
            
            
            // Creation of the new ancestral node and node configuration
            anc_gnode= w_gnodes + next_avail_inode; //Uses as anc node the next avaliable memory block into "w_gnodes" (gene nodes)
            
            if (off1->sp_index == off2->sp_index)
            {
                anc_gnode->sp_index=off1->sp_index; //Intraspecific node
            }
            
            ++next_avail_inode; //One node used
            --w_lnode->n_nodes; //One node used
            
            //Pointers
            *(anc_gnode->children+anc_gnode->n_child)=off1;
            *(anc_gnode->children+anc_gnode->n_child+1)=off2;
            off1->anc_node=anc_gnode;
            off2->anc_node=anc_gnode;
            anc_gnode->contl=w_lnode;
            anc_gnode->conts=w_lnode->conts;
            anc_gnode->n_child+=2;
            
            //Branch lengths and times
            anc_gnode->n_gen=current_ngen;
            off1->gen_length= off1->n_gen - anc_gnode->n_gen;
            off1->bl+=((off1->n_gen>w_lnode->n_gen?w_lnode->n_gen:off1->n_gen)-anc_gnode->n_gen)*p_mu; //The difference between the newly created anc_gnode (coalescence) and either the beginning of the locus branch (off1 n_gen comes from shallower branches, so you only have to add this last step) or the node n_gen (node present in this l_branch).
            off2->gen_length= off2->n_gen - anc_gnode->n_gen;
            off2->bl+=((off2->n_gen>w_lnode->n_gen?w_lnode->n_gen:off2->n_gen)-anc_gnode->n_gen)*p_mu;
            
            //Info
            anc_gnode->paralog=w_lnode->paralog;
            
            if (verbosity==5)
            {

                printf("\n\t\t\tNew coalescence\t\t    ");
#ifdef DBG
                fflush(stdout);
#endif
                
                
            }
            else if (verbosity==6)
            {
                if (w_lnode->anc_node==NULL)
                {
                    printf("\n\t\t\tNew coalescence, nodes %u and %u, coalescent time %f",off1->index,off2->index,current_ngen);
#ifdef DBG
                    fflush(stdout);
#endif
                }
                else
                {
                    printf("\n\t\t\tNew coalescence, nodes %u and %u, min coalescent age %f, coalescent time %f",off1->index, off2->index,min_ngen,current_ngen);
#ifdef DBG
                    fflush(stdout);
#endif
                }
            }
            
        }
        
        // ***
        /// Configuration of the current l_node ancestor with the new g_nodes due to the current coalescent process</dd></dl>
        
        if (anc_lnode != NULL)
        {
            n_anc_nodes=anc_lnode->n_nodes;
            for (j=0;j<w_lnode->n_nodes;++j) //Uncoalesced nodes in the current l_node
            {
                *(anc_lnode->g_nodes+n_anc_nodes+j)= *(w_gnodes_ptr+j); //Adding into the ancestor g_nodes
                off1=*(w_gnodes_ptr+j);
                off1->bl+=(off1->n_gen>w_lnode->n_gen?w_lnode->gen_length:off1->n_gen-min_ngen)*p_mu; //Addition of the branch length of this period, which can't be added at the end of the g_node branch due to the lineage especific substitution rates.
                ++anc_lnode->n_nodes;
            }
            
            if (verbosity==6)
            {
                printf("\n\t\t\tConfiguring ancestor restrictions: n_nodes=%d",anc_lnode->n_nodes);
#ifdef DBG
                fflush(stdout);
#endif
            }
        }
        else //Last
        {
            (*gene_tree)->root=anc_gnode; //Saves the root of the gene tree
        }
        
        if (verbosity>4)
        {
            printf("\n\t\tDone");
#ifdef DBG
            fflush(stdout);
#endif
        }
        
    }
    
    if (verbosity>4)
        free(iobuffer);

    switch (collapse)
    {
        case 1:
            CollapseGTree(*gene_tree, 1, 1);
            break;
    }

    return(NO_ERROR);
}

long int SimMLCGTree(l_tree *wlocus_tree, g_tree **gene_tree, name_c * names, float epsilon_brent,gsl_rng *seed, int *tn_lcoals,int verbosity, double gen_time, int collapse)
{
    // *******
    /// <dl><dt>Variable declaration</dt><dd>
    
    // ******
    /// Structures
    l_node *w_lnode=NULL, *anc_lnode=NULL;
    g_node *w_gnodes=NULL,*off1=NULL,*off2=NULL,*anc_gnode=NULL;
    g_node **w_gnodes_ptr=NULL;
    
    // ******
    /// IO related
    char  *iobuffer=NULL;
    char b_subtree[19]=" (bounded subtree)";
    
    // ******
    /// Values
    int next_avail_inode=0,n_coals=0,n_anc_nodes=0,p_Ne=0,no_counts=0;
    double p_mu=0;
    double sampled_ngen=0, min_ngen=0, current_ngen=0, *sampled_ngens=NULL;
    int avail_leaves=0,node_index=0;
    
    // ******
    /// Loop related variables</dd></dl>
    int i=0,j=0,pn_nodes=0;
    
    // ******
    /// Initialization of other general variables</dd></dl>
    if (verbosity>4)
    {
        iobuffer=calloc(count_intdigits((long)UINT_MAX, 0), sizeof(char));
        if (iobuffer==NULL)
            return MEM_ERROR;
    }
    sampled_ngens=calloc(wlocus_tree->n_gleaves, sizeof(double));
    if (sampled_ngens==NULL)
        return MEM_ERROR;
    
    // ****
    /// Association between locus tree and gene tree and gene tree reset
    if (verbosity>4)
        printf("\n\t\t\tAssociation and initialization of gene tree and locus tree structures...");
#ifdef DBG
    fflush(stdout);
#endif
    MatchTreesMLC(wlocus_tree,*gene_tree,1);
    if (verbosity>4)
        printf("Done\n\t\t\tSampling the number of lineages using the probabilities of lineage counts...");
#ifdef DBG
    fflush(stdout);
#endif
    // ****
    /// Sampling of the number of lineages entering and exiting every locus tree branch
    SampleNLineagesCDFLTree(wlocus_tree->root, wlocus_tree->n_gleaves, wlocus_tree->Ne, verbosity, seed); //I should try to extend it for using polytomies
    if (verbosity>4)
        printf("\n\t\t\tDone\n\t\t\tSampling the coalescent history along each branch of the locus tree, using either the sampled counts or the common multispecies coalescent simulation strategy depending on the presence of bounds...");
#ifdef DBG
    fflush(stdout);
#endif
    w_gnodes=(*gene_tree)->m_node;
    next_avail_inode=wlocus_tree->n_gleaves;
    *tn_lcoals=0;
    
    for (i=0; i<wlocus_tree->n_nodes; ++i)
    {
        // ***
        /// Reset of iteration-dependent variables related with the current l_node
        sampled_ngen=0.0;
        node_index=0;
        n_anc_nodes=0;
        anc_lnode=NULL;
        off1=NULL;
        off2=NULL;
        
        // ** Taking into account a new locus tree node ** //
        w_lnode=(wlocus_tree->m_node+i);
        
        switch (w_lnode->n_nodes)
        {
            case 0:
                pn_nodes=0;
                continue;
                break;
        }
        switch (w_lnode->n_olin)
        {
            case 0:
                n_coals=w_lnode->n_nodes-1;
                w_lnode->n_ilin=w_lnode->n_nodes;
                no_counts=1;
                break;
            default:
                no_counts=0;
                n_coals=w_lnode->n_nodes-w_lnode->n_olin;
                *tn_lcoals+=w_lnode->n_olin-1;
                break;
        }
        
        anc_lnode=w_lnode->anc_node;
        avail_leaves=w_lnode->n_nodes;
        
        if (pn_nodes<w_lnode->n_nodes)
        {
            w_gnodes_ptr=realloc(w_gnodes_ptr, sizeof(g_node *)*w_lnode->n_nodes);
            if (w_gnodes_ptr==NULL)
                return MEM_ERROR;
        }
        
        pn_nodes=w_lnode->n_nodes;
        
        memmove(w_gnodes_ptr,w_lnode->g_nodes,sizeof(g_node *)*w_lnode->n_nodes);
        
        if(w_lnode->Ne!=0) //Node with a private Ne.
            p_Ne=w_lnode->Ne;
        else
            p_Ne=wlocus_tree->Ne; //Global Ne.
        
        p_mu=wlocus_tree->mu*w_lnode->mu_mult;
        
        // ** Constraining coalescent time and number of events ** //
        current_ngen=w_lnode->n_gen;

        min_ngen=current_ngen-w_lnode->gen_length;
        
        if (verbosity==5)
        {
            if (w_lnode->conts==NULL) //There is no species tree
            {
                if(names!=NULL)
                    printf("\n\t\t\t\tLocus tree node %s_%d%s... ",(w_lnode->kind_node==SP && w_lnode->n_child==0)?(names->names+(w_lnode->sp_index*names->max_lname)):"Internal node",w_lnode->paralog,no_counts==0?b_subtree:"");
                else
                {
                    sprintf(iobuffer, "%d",w_lnode->sp_index);
                    printf("\n\t\t\t\tLocus tree node %s_%d%s... ",(w_lnode->kind_node==SP && w_lnode->n_child==0)?iobuffer:"Internal node",w_lnode->paralog,no_counts==0?b_subtree:"");
                }
            }
            else
            {
                if(names!=NULL)
                    printf("\n\t\t\t\tSpecies tree node %s, locus tree node %s_%d%s... ",(names->names+(w_lnode->conts->sp_index*names->max_lname)),(w_lnode->kind_node==SP && w_lnode->n_child==0)?(names->names+(w_lnode->sp_index*names->max_lname)):"Internal node",w_lnode->paralog,no_counts==0?b_subtree:"");
                else
                {
                    sprintf(iobuffer, "%d",w_lnode->conts->sp_index);
                    printf("\n\t\t\t\tSpecies tree node %s, locus tree node %s_%d%s... ",w_lnode->conts->sp_index==0?"Internal node":iobuffer,(w_lnode->kind_node==SP && w_lnode->n_child==0)?iobuffer:"Internal node",w_lnode->paralog,no_counts==0?b_subtree:"");
                }
            }
            
            
#ifdef DBG
            fflush(stdout);
#endif
        }
        else if (verbosity==6)
        {
            if (w_lnode->conts==NULL)
            {
                if(names!=NULL)
                    printf("\n\t\t\t\tLocus tree node %s_%d%s (index %d, sp_index %d, n_nodes %d) ",(w_lnode->kind_node==SP && w_lnode->n_child==0)?(names->names+(w_lnode->sp_index*names->max_lname)):"Internal node",w_lnode->paralog,no_counts==0?b_subtree:"",w_lnode->index, w_lnode->sp_index, w_lnode->n_nodes);
                else
                {
                    sprintf(iobuffer, "%d",w_lnode->sp_index);
                    printf("\n\t\t\t\tLocus tree node %s_%d%s (index %d, sp_index %d, n_nodes %d) ",(w_lnode->kind_node==SP && w_lnode->n_child==0)?iobuffer:"Internal node",w_lnode->paralog,no_counts==0?b_subtree:"",w_lnode->index, w_lnode->sp_index, w_lnode->n_nodes);
                }
            }
            else
            {
                if(names!=NULL)
                    printf("\n\t\t\t\tSpecies tree node %s, locus tree node %s_%d%s (index %d, sp_index %d, n_nodes %d) ",(names->names+(w_lnode->conts->sp_index*names->max_lname)),(w_lnode->kind_node==SP && w_lnode->n_child==0)?(names->names+(w_lnode->sp_index*names->max_lname)):"Internal node",w_lnode->paralog,no_counts==0?b_subtree:"",w_lnode->index, w_lnode->sp_index, w_lnode->n_nodes);
                else
                {
                    sprintf(iobuffer, "%d",w_lnode->conts->sp_index);
                    printf("\n\t\t\t\tSpecies tree node %s, locus tree node %s_%d%s (index %d, sp_index %d, n_nodes %d) ",w_lnode->conts->sp_index==0?"Internal node":iobuffer,(w_lnode->kind_node==SP && w_lnode->n_child==0)?iobuffer:"Internal node",w_lnode->paralog,no_counts==0?b_subtree:"",w_lnode->index, w_lnode->sp_index, w_lnode->n_nodes);
                }
            }
            
#ifdef DBG
            fflush(stdout);
#endif
        }
        
        
        // ***
        /// <dl><dt> Coalescent simulation loop. Each iteration means one coalescence. It ends either if the end of the branch has been reached, or there is no more posible coalescences (only 1 active g_node left) </dt><dd>
        // **
        /// Waiting time (coalescent time) sample, depending on the sampled counts (regular coalescence at the root)
        if (n_coals>0)
        {
            if (w_lnode->gen_length!=0)
            {
                switch (no_counts)
                {
                    case 1:
                        j=0;
                        *(sampled_ngens+0)=((-log(gsl_rng_uniform_pos(seed)))*(2*p_Ne))/(avail_leaves*(avail_leaves-1));// - (1/lambda) * ln uniform random variable. Lambda = parameter (k over 2). K, number of nodes avaliable to coalesce.
                        while (j<n_coals && current_ngen-*(sampled_ngens+j)>min_ngen)
                        {
                            current_ngen-=*(sampled_ngens+j);
                            --avail_leaves;
                            ++j;
                            *(sampled_ngens+j)=((-log(gsl_rng_uniform_pos(seed)))*(2*p_Ne))/(avail_leaves*(avail_leaves-1));// - (1/lambda) * ln uniform random variable. Lambda = parameter (k over 2). K, number of nodes avaliable to coalesce.
                        }
                        n_coals=j;
                        break;
                        
                    default:
                        for (j=0;j<n_coals;++j)
                        {
                            *(sampled_ngens+j)=SampleCoalTimeMLCFromXtoYLineages(avail_leaves, w_lnode->n_olin, current_ngen-min_ngen,p_Ne, epsilon_brent, gsl_rng_uniform_pos(seed), verbosity);
                            current_ngen-=*(sampled_ngens+j);
                            if (current_ngen-epsilon_brent<=min_ngen)
                                break;
                            --avail_leaves;
                        }
                        if(*(sampled_ngens+n_coals-1)<=0 || current_ngen-epsilon_brent<=min_ngen)
                        {
#ifdef __MPFR_H
                            current_ngen=w_lnode->n_gen;
                            avail_leaves=w_lnode->n_ilin;
                            for (j=0;j<n_coals;++j)
                            {
                                *(sampled_ngens+j)=MPFRSampleCoalTimeMLCFromXtoYLineages(avail_leaves, w_lnode->n_olin, current_ngen-min_ngen,p_Ne, epsilon_brent, gsl_rng_uniform_pos(seed), verbosity,512);
                                current_ngen-=*(sampled_ngens+j);
                                --avail_leaves;
                            }
#endif
                        }
                        break;
                }
                if(n_coals>1 && (*(sampled_ngens+n_coals-1)<=0 || current_ngen-epsilon_brent<=min_ngen))
                    return UNEXPECTED_VALUE;

            }
            else
            {
                for (j=0;j<n_coals;++j)
                {
                    *(sampled_ngens+j)=((-log(gsl_rng_uniform_pos(seed)))*(2*p_Ne))/(avail_leaves*(avail_leaves-1));// - (1/lambda) * ln uniform random variable. Lambda = parameter (k over 2). K, number of nodes avaliable to coalesce.
                    current_ngen-=*(sampled_ngens+j);
                    --avail_leaves;
                }
            }
            
            current_ngen=w_lnode->n_gen;
            avail_leaves=w_lnode->n_nodes;
            min_ngen=current_ngen-w_lnode->gen_length;
        }

        
        for (j=0; j<n_coals;++j)
        {
            // **
            /// Next sampled value
            
            sampled_ngen=*(sampled_ngens+j);
            current_ngen-=sampled_ngen;
            
            // **
            /// Randomly choice of offspring nodes, ancestral node creation and configuration </dd></dl>
            
            //Randomly choice of the first offspring node
            node_index= gsl_rng_uniform(seed)*avail_leaves; //Samples nodes
            off1= *(w_gnodes_ptr + node_index);
            --avail_leaves; //Uses one leaf
            *(w_gnodes_ptr + node_index)= *(w_gnodes_ptr + avail_leaves);// The pointer to this node now points to the last w_gnodes_ptr, inaccesible due to --avail_leaves
            
            //Randomly choice of the second offspring node
            if(avail_leaves>1)
                node_index= gsl_rng_uniform(seed)*avail_leaves;//Samples nodes
            else
                node_index=0;
            
            off2= *(w_gnodes_ptr + node_index);
            *(w_gnodes_ptr + node_index)= w_gnodes + next_avail_inode; //Now the restriction that pointed to off2 points to the new ancestral node, for using it in the next loop iteration (and/or next coalescent process)
            
            
            // Creation of the new ancestral node and node configuration
            anc_gnode= w_gnodes + next_avail_inode; //Uses as anc node the next avaliable memory block into "w_gnodes" (gene nodes)
            
            if (off1->sp_index == off2->sp_index)
            {
                anc_gnode->sp_index=off1->sp_index; //Intraspecific node
            }
            
            ++next_avail_inode; //One node used
            --w_lnode->n_nodes; //One node used
            
            //Pointers
            *(anc_gnode->children+anc_gnode->n_child)=off1;
            *(anc_gnode->children+anc_gnode->n_child+1)=off2;
            off1->anc_node=anc_gnode;
            off2->anc_node=anc_gnode;
            anc_gnode->contl=w_lnode;
            anc_gnode->conts=w_lnode->conts;
            anc_gnode->n_child+=2;
            
            //Branch lengths and times
            anc_gnode->n_gen=current_ngen;
            off1->gen_length= off1->n_gen - anc_gnode->n_gen;
            off1->bl+=((off1->n_gen>w_lnode->n_gen?w_lnode->n_gen:off1->n_gen)-anc_gnode->n_gen)*p_mu;
            off2->gen_length= off2->n_gen - anc_gnode->n_gen;
            off2->bl+=((off2->n_gen>w_lnode->n_gen?w_lnode->n_gen:off2->n_gen)-anc_gnode->n_gen)*p_mu;
            
            //Info
            anc_gnode->paralog=w_lnode->paralog;
            
            if (verbosity==5)
            {
                
                printf("\n\t\t\t\t\tNew coalescence\t\t    ");
#ifdef DBG
                fflush(stdout);
#endif
                
                
            }
            else if (verbosity==6)
            {
                if (w_lnode->anc_node==NULL)
                {
                    printf("\n\t\t\t\t\tNew coalescence, nodes %u and %u, coalescent time %f",off1->index,off2->index,current_ngen);
#ifdef DBG
                    fflush(stdout);
#endif
                }
                else
                {
                    printf("\n\t\t\t\t\tNew coalescence, nodes %u and %u, min coalescent age %f, coalescent time %f",off1->index, off2->index,min_ngen,current_ngen);
#ifdef DBG
                    fflush(stdout);
#endif
                }
            }
            
        }
        
        if(w_lnode->n_nodes>1)
        {
            if (verbosity==6 && no_counts==1)
            {
                printf("\n\t\t\t\t\t %d extra lineages going deeper, min coalescent age %f, proposed coalescent time %f",avail_leaves-1,min_ngen,current_ngen-*(sampled_ngens+j));
#ifdef DBG
                fflush(stdout);
#endif
            }
            else if (verbosity>4)
            {
                printf("\n\t\t\t\t\t %d extra lineages going deeper",avail_leaves-1);
#ifdef DBG
                fflush(stdout);
#endif
            }

        }
        if (no_counts==1)
            (*tn_lcoals)+=(avail_leaves-1);
        
        // ***
        /// Configuration of the current l_node ancestor with the new g_nodes due to the current coalescent process</dd></dl>
        
        if (anc_lnode != NULL)
        {
            n_anc_nodes=anc_lnode->n_nodes;
            for (j=0;j<w_lnode->n_nodes;++j) //Uncoalesced nodes in the current l_node
            {
                *(anc_lnode->g_nodes+n_anc_nodes+j)= *(w_gnodes_ptr+j); //Adding into the ancestor g_nodes
                off1=*(w_gnodes_ptr+j);
                off1->bl+=(off1->n_gen>w_lnode->n_gen?w_lnode->gen_length:off1->n_gen-min_ngen)*p_mu; //Addition of the branch length of this period, which can't be added at the end of the g_node branch due to the lineage especific substitution rates.
                ++anc_lnode->n_nodes;
            }
        }
        else //Last
        {
            (*gene_tree)->root=anc_gnode; //Saves the root of the gene tree
        }
        
        if (verbosity>4)
        {
            printf("\n\t\t\t\tDone");
#ifdef DBG
            fflush(stdout);
#endif
        }
        
    }
    
    if (verbosity>4)
    {
        free(iobuffer);
        printf("\n\t\t\tDone");
#ifdef DBG
        fflush(stdout);
#endif
    }
    
    free(w_gnodes_ptr);
    free(sampled_ngens);
    
    switch (collapse)
    {
        case 1:
            CollapseGTree(*gene_tree, 1, 1);
            break;
    }
    
    return(NO_ERROR);
}

int CalcProbsNLineagesLTree(l_node *node, int Ne, int in_subtree ,int verbosity)
{
    // **
    /// Variable declaration
    int *max_nlin=NULL, i=0,j=0,k=0, fmax_nlin=0,combfmax_nlin=1,return_nlin=0, nbranch_wlin=0;
    //Kahan summation
    double comp=0,comp_opsum=0,o_prob_sum=0,bkp_prob=0;
    
    // *
    /// Obtaining number of lineages and input probabilities using combinatorics and data of output lineages/probabilities of children
    switch (node->n_child)
    {
        default:
            max_nlin=calloc(node->n_child,sizeof(int));
            ErrorReporter((long int)max_nlin,NULL);
            
            for (i=0; i<node->n_child; ++i)
            {
                if (i==1 && (node->kind_node==DUP || node->kind_node==TRFR || node->kind_node==GC))
                {
                    *(max_nlin+i)=CalcProbsNLineagesLTree(*(node->children+i),Ne,1,verbosity);
                }
                else
                    *(max_nlin+i)=CalcProbsNLineagesLTree(*(node->children+i),Ne,in_subtree,verbosity);
                fmax_nlin+=*(max_nlin+i);
                combfmax_nlin*=*(max_nlin+i);
                switch (*(max_nlin+i))
                {
                    case 0:
                        break;
                    default:
                        ++nbranch_wlin;
                        break;
                }
            }
            switch (in_subtree)
            {
                case 1:
                    node->i_probs=calloc(fmax_nlin+1, sizeof(double));
                    node->i_combprobs=calloc(combfmax_nlin, sizeof(double));
                    node->o_probs=calloc(fmax_nlin+1, sizeof(double));
                    ErrorReporter((long int)node->i_probs, "");
                    ErrorReporter((long int)node->i_combprobs, "");
                    ErrorReporter((long int)node->o_probs, "");
                    switch (*(max_nlin+1))
                    {
                        case 1:
                            bkp_prob=*((*(node->children+1))->o_probs+1);
                            *((*(node->children+1))->o_probs+1)=1;
                            break;
                    }
                    switch (nbranch_wlin)
                    {
                        case 2:
                            k=0;
                            for (i=2; i<=fmax_nlin; ++i)
                                for (j=1; j<i; ++j)
                                {
                                    if (j<=*(max_nlin+0) && (i-j)<=*(max_nlin+1))
                                    {
                                        *(node->i_combprobs+k)= *((*(node->children+0))->o_probs+j) * *((*(node->children+1))->o_probs+(i-j));
                                        *(node->i_probs+i)+=*(node->i_combprobs+k);
                                        ++k;
                                    }
                                }
                            break;
                        case 1:
                            for (j=0; j<node->n_child; ++j)
                            {
                                if (*(max_nlin+j)!=0)
                                {
                                    for (i=1; i<=fmax_nlin; ++i)
                                    {
                                        *(node->i_probs+i)=*((*(node->children+j))->o_probs+i);
                                    }
                                    break;
                                }
                            }
                            break;
                        case 0:
                            fmax_nlin=0;
                            *(node->i_probs)=1;
                            *(node->o_probs)=1;
                            break;
                        default:
                            ErrorReporter(UNEXPECTED_VALUE,"Not implemented yet\n");
                            break;
                    }
                    break;
                    
                default:
                    break;
            }
            
            break;
            
        case 0:
            fmax_nlin=node->n_ilin;//Either lost nodes (RGC/RTRFR/LOSS) or SP leaves. Given n_ilin, previously applied by MatchTreesMLC
            node->i_probs=calloc(fmax_nlin+1, sizeof(double));
            node->o_probs=calloc(fmax_nlin+1, sizeof(double));
            ErrorReporter((long int)node->i_probs, "");
            ErrorReporter((long int)node->o_probs, "");
            *(node->i_probs+fmax_nlin)=1;
            
            break;
    }
    
    return_nlin=fmax_nlin;

    // *
    /// Calculating the output probabilities using \ref LogscaleProbCoalFromXtoYLineages for regular nodes inside bounded subtrees or setting it to fixed values (with the corresponding n_olin) if necessary (1 if bounded node, 0 if either no lineages present or unbounded process)
    switch (in_subtree)
    {
        case 1:
            if (node->anc_node==NULL || fmax_nlin==0) //Without any node
            {
                node->n_olin=0;
                *(node->o_probs)=1;
                return_nlin=0;
            }
            else if (fmax_nlin==1) //Just one input node
            {
                node->n_olin=1;
                *(node->o_probs+1)=1;
            }
            else //Normal
            {
                comp_opsum=0;
                for (i=1; i<=fmax_nlin;++i)
                {
                    comp=0;
                    for (j=i; j<=fmax_nlin;++j)
                        SKahanSum(node->o_probs+i,KahanLogscaleProbCoalFromXtoYLineages(j,i,node->gen_length,node->Ne==0?Ne:node->Ne)* *(node->i_probs+j),&comp);
                    SKahanSum(&o_prob_sum, *(node->o_probs+i), &comp_opsum);
                    if (o_prob_sum>=1)
                        break;
                }
                if(node->anc_node!=NULL && ((node->anc_node->kind_node==DUP || node->anc_node->kind_node==TRFR || node->anc_node->kind_node==GC) && *(node->anc_node->children+1)==node))
                {
                    node->n_olin=1;
                    return_nlin=1;
                }
            }
            
            if (verbosity>4)
            {
                if (verbosity==5)
                {
                    printf("\n\t\t\t\tNode %u, input probabilities sum %lf, output probabilities sum %lf",node->index,VKahanSum(node->i_probs,fmax_nlin+1),VKahanSum(node->o_probs,fmax_nlin+1));
                }
                else
                {
                    printf("\n\t\t\t\tNode %u, kind node %u, n_children %u, input probabilities sum %lf, input probabilities: ", node->index, node->kind_node, node->n_child,VKahanSum(node->i_probs,fmax_nlin+1));
                    for (i=0; i<=fmax_nlin; ++i)
                    {
                        printf("%u:%lf,",i,*(node->i_probs+i));
                    }
                    printf("; Output probabilities sum %lf, output probabilities: ", VKahanSum(node->o_probs,return_nlin+1));
                    for (i=0; i<=return_nlin; ++i)
                    {
                        printf("%u:%lf,",i,*(node->o_probs+i));
                    }
                    
                }
#ifdef DBG
                fflush(stdout);
#endif
            }
            
            if (node->n_child>0 && *(max_nlin+1)==1)
                *((*(node->children+1))->o_probs+1)=bkp_prob;
            break;
            
        default:
            if (verbosity>4)
            {
                if (verbosity==5)
                {
                    printf("\n\t\t\t\tNode %u, unbounded sampling (multispecies coalescent)",node->index);
                }
                else
                {
                    printf("\n\t\t\t\tNode %u, kind node %u, n_children %u, unbounded sampling (multispecies coalescent)", node->index, node->kind_node, node->n_child);
                }
#ifdef DBG
                fflush(stdout);
#endif
            }
            break;
    }

    
    if (max_nlin!=NULL)
        free(max_nlin);
    
    node->fmax_nlin=return_nlin;
    
    return return_nlin;
}

inline int SampleNLineagesLNode(l_node *node, int n_gleaves, gsl_rng *seed)
{
    int result=0;
    switch (node->n_olin)
    {
        case 0:
            result=(int)SampleNDoubles(n_gleaves+1,node->o_probs,seed);
            break;
            
        default:
            result=node->n_olin;
            break;
    }
    return result;
}

void SampleNLineagesLTree(l_node *node, int n_gleaves, int Ne, int verbosity, gsl_rng *seed)
{
    int i=0,n_branch_wlin=0, wn_olin0=0, wn_olin1=0, reject=1, child,n_it=0;
    double p=0,u=0;
    
    if (node->n_child!=0)
    {
        for (i=0; i<node->n_child;++i)
        {
            if (*((*(node->children+i))->o_probs+0)==0) //Or n_olin???
            {
                ++n_branch_wlin;
                child=i;
            }
        }
        
        switch (n_branch_wlin)
        {
            case 2:
                wn_olin0=SampleNLineagesLNode(*(node->children+0), n_gleaves, seed);
                wn_olin1=SampleNLineagesLNode(*(node->children+1), n_gleaves, seed);
                switch (node->n_olin)
                {
                    case 0:
                        switch (verbosity)
                        {
                            case 0:
                            case 1:
                            case 2:
                            case 3:
                            case 4:
                            case 5:
                                break;
                            default:
                                printf("\n\t\t\t\tNode %u: unrestricted, 2 children (node %u and %u): sampled %u sampled input lineages, output children lineages %u %u", node->index, (*(node->children+0))->index, (*(node->children+1))->index, wn_olin0+wn_olin1, wn_olin0,wn_olin1);
                                break;
                        }
#ifdef DBG
                        fflush(stdout);
#endif
                        break;
                        
                    default: //Rejection sampling when it has a concrete number of starting lineages.
                        reject=1;
                        n_it=0;
                        
                        switch (verbosity)
                        {
                            case 0:
                            case 1:
                            case 2:
                            case 3:
                            case 4:
                                break;
                            default:
                                printf("\n\t\t\t\tNode %u: restricted (%u output lineages), 2 children (node %u and %u): Rejection sampling... ", node->index, node->n_olin, (*(node->children+0))->index, (*(node->children+1))->index);
                                break;
                        }
#ifdef DBG
                        fflush(stdout);
#endif
                        
                        while (reject==1 && n_it<MAX_IT)
                        {
                            ++n_it;
                            p=KahanLogscaleProbCoalFromXtoYLineages(wn_olin0+wn_olin1, node->n_olin, node->gen_length, node->Ne==0?Ne:node->Ne);
                            u=gsl_rng_uniform(seed);
                            switch (verbosity)
                            {
                                case 0:
                                case 1:
                                case 2:
                                case 3:
                                case 4:
                                case 5:
                                    break;
                                default:
                                    printf("\n\t\t\t\t\tCandidate probability %lf, calculated probability %lf, candidate number of lineages of child1 %u, child2 %u",u,p,wn_olin0,wn_olin1);
                                    break;
                            }
#ifdef DBG
                            fflush(stdout);
#endif
                            
                            if (u<p)
                                reject=0;
                            else
                            {
                                wn_olin0=SampleNLineagesLNode(*(node->children+0), n_gleaves, seed);
                                wn_olin1=SampleNLineagesLNode(*(node->children+1), n_gleaves, seed);
                            }
                        }
                        
                        if (reject==1)
                            ErrorReporter(LOOP_ERROR,": reached MAX_IT iterations using rejection sampling to sample lineage counts. This problem is usually related to really unlikely scenarios with unrealistic short branches and big population sizes. Please, check your simulation settings.");
                        
                        switch (verbosity)
                        {
                            case 0:
                            case 1:
                            case 2:
                            case 3:
                            case 4:
                                break;
                            case 5:
                                printf("%u input lineages, from output children %u and %u", wn_olin0+wn_olin1,wn_olin0,wn_olin1);
                                break;
                            default:
                                printf("\n\t\t\t\t\tDone: %u input lineages, from outputs child1 %u, child2 %u", wn_olin0+wn_olin1,wn_olin0,wn_olin1);
                                break;
                        }
#ifdef DBG
                        fflush(stdout);
#endif
                        break;
                }
                node->n_ilin=wn_olin0+wn_olin1;
                (*(node->children+0))->n_olin=wn_olin0;
                (*(node->children+1))->n_olin=wn_olin1;
                SampleNLineagesLTree(*(node->children+0),n_gleaves,Ne,verbosity,seed);
                SampleNLineagesLTree(*(node->children+1),n_gleaves,Ne,verbosity,seed);
                break;
            case 1:
                wn_olin0=SampleNLineagesLNode(*(node->children+child), n_gleaves, seed);
                switch (node->n_olin)
                {
                    case 0:
                        switch (verbosity)
                        {
                            case 0:
                            case 1:
                            case 2:
                            case 3:
                            case 4:
                            case 5:
                                break;
                            default:
                                printf("\n\t\t\t\tNode %u: unrestricted, 1 children (node %u): sampled %u sampled input lineages", node->index, (*(node->children+0))->index, wn_olin0);
                                break;
                        }
#ifdef DBG
                        fflush(stdout);
#endif
                        break;
                    
                    default: //Rejection sampling when it has a concrete number of starting lineages.
                        switch (verbosity)
                        {
                            case 0:
                            case 1:
                            case 2:
                            case 3:
                            case 4:
                                break;
                            default:
                                printf("\n\t\t\t\tNode %u: restricted (%u output lineages), 1 children (node %u): Rejection sampling... ", node->index, node->n_olin, (*(node->children+0))->index);
                                break;
                        }
#ifdef DBG
                        fflush(stdout);
#endif
                        reject=1;
                        while (reject==1)
                        {
                            p=KahanLogscaleProbCoalFromXtoYLineages(wn_olin0, node->n_olin, node->gen_length, node->Ne==0?Ne:node->Ne);
                            u=gsl_rng_uniform(seed);
                            switch (verbosity)
                            {
                                case 0:
                                case 1:
                                case 2:
                                case 3:
                                case 4:
                                case 5:
                                    break;
                                default:
                                    printf("\n\t\t\t\t\tCandidate probability %lf, calculated probability %lf, candidate number of lineages of child1 %u",u,p,wn_olin0);
                                    break;
                            }
#ifdef DBG
                            fflush(stdout);
#endif
                            if (u<p)
                                reject=0;
                            else
                            {
                                wn_olin0=SampleNLineagesLNode(*(node->children+child), n_gleaves, seed);
                            }
                        }
                        switch (verbosity)
                        {
                            case 0:
                            case 1:
                            case 2:
                            case 3:
                            case 4:
                                break;
                            case 5:
                                printf("%u input lineages", wn_olin0);
                                break;
                            default:
                                printf("\n\t\t\t\t\tDone: %u input lineages",wn_olin0);
                                break;
                        }
#ifdef DBG
                        fflush(stdout);
#endif
                        break;
                        
                }
                node->n_ilin=wn_olin0;
                (*(node->children+child))->n_olin=wn_olin0;
                SampleNLineagesLTree(*(node->children+child),n_gleaves,Ne,verbosity,seed);
                break;
            case 0:
                if (verbosity>4)
                {
                    printf("\n\t\t\t\tNode %u: %u output lineages and %u input lineages, 0 children with lineages, previously set lineage counts", node->index, node->n_olin, node->n_ilin);
                }
#ifdef DBG
                fflush(stdout);
#endif
                //The n_lineages have to be previosly set here.
                break;
            default:
                ErrorReporter(UNEXPECTED_VALUE,": not implemented yet\n");
                break;
        }
    }
    else
    {
        //The n_lineages have to be previosly set here.
        if (verbosity>4)
        {
            printf("\n\t\t\t\tNode %u: %u output lineages and %u input lineages (leaf), previously set lineage counts", node->index, node->n_olin, node->n_ilin);
        }
#ifdef DBG
        fflush(stdout);
#endif
    }
}

void SampleNLineagesCDFLTree(l_node *node, int n_gleaves, int Ne, int verbosity, gsl_rng *seed)
{
    int i=0,j=0,prob_padding=0,prob_combs=0,*mapcomb_lin0=NULL,n_branch_wlin=0, wn_olin0=0, wn_olin1=0, child=0;
    double p=0,u=0;
    
    int kt=0,ktmax=0;
#ifdef DBG
    int ktdbug=0;
#endif
    double sum=0,cden=0,comp=0;
    
    switch(node->n_child)
    {
        default:
            switch (node->n_olin)
            {
                case 0:
                    for (i=0; i<node->n_child;++i)
                        SampleNLineagesCDFLTree(*(node->children+i),n_gleaves,Ne,verbosity,seed);
                    break;
                    
                default:
                    u=gsl_rng_uniform(seed);
                    ktmax=(*(node->children+0))->fmax_nlin+(*(node->children+1))->fmax_nlin;
                    cden=*(node->o_probs+node->n_olin);
                    comp=0;
#ifndef DBG
                    kt=node->n_olin-1;
                    if (node->n_olin<ktmax)
                        do
                        {
                            ++kt;
                            SKahanSum(&sum,KahanLogscaleProbCoalFromXtoYLineages(kt, node->n_olin, node->gen_length, node->Ne!=0?node->Ne:Ne)**(node->i_probs+kt), &comp);
                            p=sum/cden;
                        } while (p<u && kt<ktmax);
                    else if (node->n_olin==ktmax)
                        kt=ktmax;
                    else
                        ErrorReporter(UNEXPECTED_VALUE, "");
                        
#endif
                    
#ifdef DBG
                    ktdbug=kt=node->n_olin;
                    if (node->n_olin<ktmax)
                    {
                        do
                        {
                            SKahanSum(&sum,KahanLogscaleProbCoalFromXtoYLineages(ktdbug, node->n_olin, node->gen_length, node->Ne!=0?node->Ne:Ne)**(node->i_probs+ktdbug), &comp);
                            p=sum/cden;
                            ++ktdbug;
                            if (p<u)
                                ++kt;
                        } while (ktdbug<=ktmax);
                        if (sum/cden<0.9)
                        {
                            ErrorReporter(UNEXPECTED_VALUE, ": The CDF of the input number of lineages given the output number of lineages is not working properly. The simulator is generating a really unrealistic biological scenario. Please, check your input parameters.");
                        }
                        
                        printf("CDF integrates to %lf",sum/cden);
                    }
                    else if (node->n_olin==ktmax)
                        kt=ktmax;
                    else
                        ErrorReporter(UNEXPECTED_VALUE, "");

#endif
                    for (i=0; i<node->n_child;++i)
                    {
                        if (*((*(node->children+i))->o_probs+0)==0)
                        {
                            ++n_branch_wlin;
                            child=i;
                        }
                    }

                    switch (n_branch_wlin)
                    {
                        case 2:
                            node->n_ilin=kt;
                            prob_padding=0;
                            prob_combs=0;
                            mapcomb_lin0=calloc(kt, sizeof(int));
                            ErrorReporter((long int)mapcomb_lin0, "");
                            for (i=2; i<kt; ++i)
                                for (j=1; j<i; ++j)
                                    if (j<=(*(node->children+0))->fmax_nlin && (i-j)<=(*(node->children+1))->fmax_nlin)
                                        ++prob_padding;
                            for (j=1; j<kt; ++j)
                                if (j<=(*(node->children+0))->fmax_nlin && (kt-j)<=(*(node->children+1))->fmax_nlin)
                                {
                                    *(mapcomb_lin0+prob_combs)=j;
                                    ++prob_combs;
                                }
                            
                            wn_olin0=*(mapcomb_lin0+(int)SampleNDoubles(prob_combs, node->i_combprobs+prob_padding, seed));
                            wn_olin1=kt-wn_olin0;
                            (*(node->children+0))->n_olin=wn_olin0;
                            (*(node->children+1))->n_olin=wn_olin1;
                            if (verbosity>4)
                            {
                                printf("\n\t\t\t\tNode %u: %u output lineages and %u input lineages, 2 children with lineages, node %u %u lineages, node %u %u lineages", node->index, node->n_olin, node->n_ilin, (*(node->children+0))->index,(*(node->children+0))->n_olin,(*(node->children+1))->index,(*(node->children+1))->n_olin);
                            }
#ifdef DBG
                            fflush(stdout);
#endif
                            free(mapcomb_lin0);
                            SampleNLineagesCDFLTree(*(node->children+0),n_gleaves,Ne,verbosity,seed);
                            SampleNLineagesCDFLTree(*(node->children+1),n_gleaves,Ne,verbosity,seed);
                            break;
                        case 1:
                            node->n_ilin=kt;
                            (*(node->children+child))->n_olin=kt;
                            if (verbosity>4)
                            {
                                printf("\n\t\t\t\tNode %u: %u output lineages and %u input lineages, 1 children node %u with lineages", node->index, node->n_olin, node->n_ilin, (*(node->children+child))->index);
                            }
#ifdef DBG
                            fflush(stdout);
#endif
                            SampleNLineagesCDFLTree(*(node->children+child),n_gleaves,Ne,verbosity,seed);
                            break;
                        case 0:
                            if (verbosity>4)
                            {
                                printf("\n\t\t\t\tNode %u: %u output lineages and %u input lineages, 0 children with lineages, previously set lineage counts", node->index, node->n_olin, node->n_ilin);
                            }
    #ifdef DBG
                            fflush(stdout);
    #endif
                            //The n_lineages have to be previosly set here.
                            break;
                        default:
                            ErrorReporter(UNEXPECTED_VALUE,": not implemented yet\n");
                            break;
                    }

                    break;
            }
            break;
        case 0://The n_lineages have to be previosly set here.
            if (verbosity>4 && node->n_olin>1)
            {
                printf("\n\t\t\t\tLeaf node %u: %u preset input lineages (leaf), and %u output lineages", node->index, node->n_ilin, node->n_olin);
            }
    #ifdef DBG
            fflush(stdout);
    #endif
            break;
        
    }
}


l_tree * NewLTree (int n_nodes, int n_leaves, int n_gleaves, int max_children, double gen_time, int Ne, double mu)
{
    l_tree * tree=NULL;
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Tree memory allocation
    tree=calloc(1,sizeof(l_tree));
    ErrorReporter((long int)tree,NULL);
    
    // *
    /// Tree initialization
    tree->n_nodes=n_nodes;
    tree->n_leaves=n_leaves;
    tree->n_gleaves=n_gleaves;
    tree->max_children=max_children;
    tree->gen_time=gen_time;
    tree->Ne=Ne;
    tree->mu=mu;
    tree->species_tree=NULL;
    tree->gene_tree=NULL;
    
    // *
    /// Tree nodes allocation and initialization by \ref NewLNodes </dd></dl>
    
    if (n_nodes>1)
    {
        tree->m_node=NewLNodes(n_nodes,n_gleaves,max_children);
        tree->root=NULL;
    }
    else if (n_nodes==1)
    {
        tree->root=NewLNodes(n_nodes, n_gleaves, max_children);
        tree->m_node=NULL;
    }
    else
    {
        tree->root=NULL;
        tree->m_node=NULL;
    }
    
    return (tree);
}

g_tree * NewGTree (int n_nodes, int max_children, double gen_time)
{
    g_tree * tree=NULL;
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Tree memory allocation
    tree=calloc(1,sizeof(g_tree));
    ErrorReporter((long int) tree,NULL);
    
    // *
    /// Tree initialization
    tree->n_nodes=n_nodes;
    tree->max_children=max_children;
    tree->species_tree=NULL;
    tree->locus_tree=NULL;
    tree->gen_time=gen_time;
    
    // *
    /// Tree nodes allocation and initialization by \ref NewGNodes </dd></dl>
    if (n_nodes<2)
    {
        tree->root=NewGNodes(n_nodes,max_children);
        tree->m_node=NULL;
    }
    else
    {
        tree->m_node=NewGNodes(n_nodes,max_children);
        tree->root=NULL;
    }
    
    
    return (tree);
}

// ** Tree copy ** //

long int CopySTree (s_tree ** out_tree_ptr,s_tree * in_tree,int tree_struct, int l_nodes_ptr)
{
    s_tree * out_tree=NULL;
    s_node * w_output=NULL, * w_input=NULL;
    int i=0, j=0,post_order=0;
    
    // ****
    /// <dl><dt> Function structure </dt><dd>
    
    // ***
    /// <dl><dt>Input tree memory control</dt><dd>
    
    // **
    /// Presence of root node and correct configuration in pre or post-order
    if (in_tree->root == NULL || ((in_tree->root->index != in_tree->n_nodes-1) && (in_tree->root->index!=0)))
        return MEM_ERROR;
    else if (in_tree->root->index == in_tree->n_nodes-1)
    {
        post_order=1;
    }
    
    // **
    /// Input tree is collapsed using a post-order if it has sparse nodes</dd></dl>
    if (in_tree->m_node == NULL)
    {
        CollapseSTree(in_tree,post_order);
    }
    
    // ***
    /// <dl><dt>Output tree memory control</dt><dd>
    
    // **
    /// Allocation if *out_tree_ptr==NULL
    if (*out_tree_ptr==NULL)
    {
        *out_tree_ptr=NewSTree(in_tree->n_nodes, in_tree->n_leaves,in_tree->n_gleaves, in_tree->max_children,in_tree->gen_time, in_tree->Ne, in_tree->mu);
    }
    // **
    /// Reallocation (by \ref FreeSTree and \ref NewSTree) if the output tree pointer has a sparse tree, and s_tree::n_nodes or s_tree::max_children are different. </dd></dl>
    else if ((*out_tree_ptr)->m_node==NULL || (*out_tree_ptr)->n_nodes!=in_tree->n_nodes || (*out_tree_ptr)->max_children!=in_tree->max_children)
    {
        FreeSTree(out_tree_ptr);
        *out_tree_ptr=NewSTree(in_tree->n_nodes, in_tree->n_leaves,in_tree->n_gleaves, in_tree->max_children,in_tree->gen_time, in_tree->Ne, in_tree->mu);
    }
    
    out_tree=*out_tree_ptr;
    
    // ***
    /// <dl><dt>Input tree nodes loop</dt><dd>
    for (i=0;i<in_tree->n_nodes;++i)
    {
        w_output=out_tree->m_node+i;
        w_input=in_tree->m_node+i;
        
        // **
        /// Copy of values
        w_output->index=w_output->index;
        w_output->sp_index=w_input->sp_index;
        w_output->n_child=w_input->n_child;
        w_output->n_lnodes=w_input->n_lnodes;
        w_output->n_replicas=w_input->n_replicas;
        w_output->Ne=w_input->Ne;
        w_output->n_gen=w_input->n_gen;
        w_output->time=w_input->time;
        w_output->gen_length=w_input->gen_length;
        w_output->mu_mult=w_input->mu_mult;
        w_output->gtime_mult=w_input->gtime_mult;
        
        // **
        /// <dl><dt>Copy of tree struct pointers</dt><dd>
        if (tree_struct>0)
        {
            if (w_output->children==NULL)
                return MEM_ERROR;
            // *
            /// Child loop
            for (j=0;j<w_input->n_child;++j)
            {
                *(w_output->children+j)=(out_tree->m_node+(*(w_input->children+j))->index);
            }
            // *
            /// Ancestor </dd></dl>
            if (w_input->anc_node != NULL)
            {
                w_output->anc_node=(out_tree->m_node+w_input->anc_node->index);
            }
            else
            {
                w_output->anc_node=NULL;
            }
        }
        
        // **
        /// Copy of s_node::l_nodes pointers </dd></dl>
        if (l_nodes_ptr>0)
        {
            w_output->l_nodes=w_input->l_nodes;
        }
    }
    
    // ***
    /// Tree values </dd></dl>
    out_tree->n_leaves=in_tree->n_leaves;
    out_tree->n_gleaves=in_tree->n_gleaves;
    out_tree->gen_time=in_tree->gen_time;
    out_tree->Ne=in_tree->Ne;
    out_tree->mu=in_tree->mu;
    out_tree->locus_tree=in_tree->locus_tree;
    out_tree->gene_tree=in_tree->gene_tree;
    if (in_tree->root != NULL)
        out_tree->root=out_tree->m_node+(in_tree->root->index);
    
    return NO_ERROR;
}

long int CopyLTree (l_tree **out_tree_ptr, l_tree *in_tree, int tree_struct, int l_nodes_ptr, int g_nodes_ptr, int probs, int relink, int preserve)
{
    l_tree * out_tree=NULL;
    l_node * w_output=NULL, * w_input=NULL;
    int i=0, j=0, post_order=0, tfmax_nlin=0;
    
    // ****
    /// <dl><dt> Function structure </dt><dd>
    
    // ***
    /// <dl><dt>Input tree memory control</dt><dd>
    
    // **
    /// Presence of root node and correct configuration in pre or post-order
    if (in_tree->root == NULL || ((in_tree->root->index != in_tree->n_nodes-1) && (in_tree->root->index!=0)))
        return MEM_ERROR;
    else if (in_tree->root->index == in_tree->n_nodes-1)
    {
        post_order=1;
    }
    
    // **
    /// Input tree is collapsed using a post-order if it has sparse nodes</dd></dl>
    if (in_tree->m_node == NULL)
    {
        CollapseLTree(in_tree,post_order,1,probs);
    }
    
    // ***
    /// <dl><dt>Output tree memory control</dt><dd>
    
    // **
    /// Allocation if *out_tree_ptr==NULL
    if (*out_tree_ptr==NULL)
    {
        *out_tree_ptr=NewLTree(in_tree->n_nodes, in_tree->n_leaves,in_tree->n_gleaves, in_tree->max_children, in_tree->gen_time, in_tree->Ne, in_tree->mu);
    }
    // **
    /// Reallocation (by \ref FreeSTree and \ref NewLTree) if the output tree pointer has a sparse tree, and l_tree::n_nodes or l_tree::max_children are different. </dd></dl>
    else if ((*out_tree_ptr)->m_node==NULL || (*out_tree_ptr)->n_nodes!=in_tree->n_nodes || (*out_tree_ptr)->max_children!=in_tree->max_children)
    {
        FreeLTree(out_tree_ptr);
        *out_tree_ptr=NewLTree(in_tree->n_nodes, in_tree->n_leaves,in_tree->n_gleaves, in_tree->max_children, in_tree->gen_time, in_tree->Ne, in_tree->mu);
    }
    
    out_tree=*out_tree_ptr;
    
    if(in_tree->n_nodes<2) //HARDCODED very rare special case, when the tree only has one node (root and leave at the same time). This is not a true tree, but it allows to generate a gene tree from replicas of one species.
    {
        w_output=out_tree->root;
        w_input=in_tree->root;
        w_output->index=w_input->index;
        w_output->sp_index=w_input->sp_index;
        w_output->kind_node=w_input->kind_node;
        w_output->paralog=w_input->paralog;
        w_output->n_child=w_input->n_child;
        w_output->n_nodes=w_input->n_nodes;
        w_output->Ne=w_input->Ne;
        w_output->n_gen=w_input->n_gen;
        w_output->time=w_input->time;
        w_output->gen_length=w_input->gen_length;
        w_output->conts=w_input->conts;
        w_output->mu_mult=w_input->mu_mult;
        w_output->gtime_mult=w_input->gtime_mult;
        w_output->n_ilin=w_input->n_ilin;
        w_output->n_olin=w_input->n_olin;
        w_output->fmax_nlin=w_input->fmax_nlin;
        
        switch (probs)
        {
            case 1:
                switch (preserve)
                {
                    case 1:
                        if (w_input->i_probs!=NULL)
                        {
                            w_output->i_probs=calloc(w_input->fmax_nlin+1,sizeof(double));
                            if (w_output->i_probs==NULL)
                                return MEM_ERROR;
                            for (j=0;j<=w_input->fmax_nlin;++j)
                            {
                                *(w_output->i_probs+j)=*(w_input->i_probs+j);
                            }
                        }
                        if (w_input->i_combprobs!=NULL)
                        {
                            tfmax_nlin=1;
                            for (j=0;j<w_input->n_child;++j)
                                tfmax_nlin*=(*(w_input->children+j))->fmax_nlin;
                            w_output->i_combprobs=calloc(tfmax_nlin,sizeof(double));
                            if (w_output->i_combprobs==NULL)
                                return MEM_ERROR;
                            for (j=0;j<tfmax_nlin;++j)
                            {
                                *(w_output->i_combprobs+j)=*(w_input->i_combprobs+j);
                            }
                        }
                        if (w_input->o_probs!=NULL)
                        {
                            w_output->o_probs=calloc(in_tree->n_gleaves+1,sizeof(double));
                            if (w_output->o_probs==NULL)
                                return MEM_ERROR;
                            for (j=0;j<=w_input->fmax_nlin;++j)
                            {
                                *(w_output->o_probs+j)=*(w_input->o_probs+j);
                            }
                        }
                        break;
                        
                    default:
                        w_output->i_probs=w_input->i_probs;
                        w_output->i_combprobs=w_input->i_combprobs;
                        w_output->o_probs=w_input->o_probs;
                        w_input->i_probs=NULL;
                        w_input->i_combprobs=NULL;
                        w_input->o_probs=NULL;
                        break;
                }
                break;
            default:
                w_output->i_probs=NULL;
                w_output->i_combprobs=NULL;
                w_output->o_probs=NULL;
                break;
        }
        
        if (tree_struct>0)
        {
            if (w_output->children==NULL)
                return MEM_ERROR;
            
            for (j=0;j<w_input->n_child;++j)
            {
                *(w_output->children+j)=(out_tree->m_node+(*(w_input->children+j))->index);
            }
            w_output->anc_node=NULL;
        }
        
        if (g_nodes_ptr>0)
        {
            for (j=0; j>w_input->n_nodes; ++j)
            {
                *(w_output->children+j)=out_tree->m_node+((*(w_input->children+j))->index);
            }
        }
        
        if (relink==1 && w_output->conts!=NULL && w_output->conts->l_nodes->index == i)
            w_output->conts->l_nodes=w_output;
        
        i=1;
    }
    else
    {
        i=0;
    }

    // ***
    /// <dl><dt>Input tree nodes loop</dt><dd>
    for (;i<in_tree->n_nodes;++i)
    {
        w_output=out_tree->m_node+i;
        w_input=in_tree->m_node+i;
        
        // **
        /// Copy of values
        w_output->index=w_input->index;
        w_output->sp_index=w_input->sp_index;
        w_output->kind_node=w_input->kind_node;
        w_output->paralog=w_input->paralog;
        w_output->n_child=w_input->n_child;
        w_output->n_nodes=w_input->n_nodes;
        w_output->Ne=w_input->Ne;
        w_output->n_gen=w_input->n_gen;
        w_output->time=w_input->time;
        w_output->gen_length=w_input->gen_length;
        w_output->conts=w_input->conts;
        w_output->mu_mult=w_input->mu_mult;
        w_output->gtime_mult=w_input->gtime_mult;
        w_output->n_ilin=w_input->n_ilin;
        w_output->n_olin=w_input->n_olin;
        w_output->fmax_nlin=w_input->fmax_nlin;
        
        // **
        /// Copy of input and output probabilities of numbers of lineages
        switch (probs)
        {
            case 1:
                switch (preserve)
            {
                case 1:
                    if (w_input->i_probs!=NULL)
                    {
                        w_output->i_probs=calloc(w_input->fmax_nlin+1,sizeof(double));
                        if (w_output->i_probs==NULL)
                            return MEM_ERROR;
                        for (j=0;j<=w_input->fmax_nlin;++j)
                        {
                            *(w_output->i_probs+j)=*(w_input->i_probs+j);
                        }
                    }
                    if (w_input->i_combprobs!=NULL)
                    {
                        tfmax_nlin=1;
                        for (j=0;j<w_input->n_child;++j)
                            tfmax_nlin*=(*(w_input->children+j))->fmax_nlin;
                        w_output->i_combprobs=calloc(tfmax_nlin,sizeof(double));
                        if (w_output->i_combprobs==NULL)
                            return MEM_ERROR;
                        for (j=0;j<tfmax_nlin;++j)
                        {
                            *(w_output->i_combprobs+j)=*(w_input->i_combprobs+j);
                        }
                    }
                    if (w_input->o_probs!=NULL)
                    {
                        w_output->o_probs=calloc(in_tree->n_gleaves+1,sizeof(double));
                        if (w_output->o_probs==NULL)
                            return MEM_ERROR;
                        for (j=0;j<=w_input->fmax_nlin;++j)
                        {
                            *(w_output->o_probs+j)=*(w_input->o_probs+j);
                        }
                    }
                    break;
                    
                default:
                    w_output->i_probs=w_input->i_probs;
                    w_output->i_combprobs=w_input->i_combprobs;
                    w_output->o_probs=w_input->o_probs;
                    w_input->i_probs=NULL;
                    w_input->i_combprobs=NULL;
                    w_input->o_probs=NULL;
                    break;
            }
                break;
            default:
                w_output->i_probs=NULL;
                w_output->i_combprobs=NULL;
                w_output->o_probs=NULL;
                break;
        }
        
        // **
        /// <dl><dt>Copy of tree struct pointers</dt><dd>
        if (tree_struct>0)
        {
            if (w_output->children==NULL)
                return MEM_ERROR;
            
            // *
            /// Child loop
            for (j=0;j<w_input->n_child;++j)
            {
                *(w_output->children+j)=(out_tree->m_node+(*(w_input->children+j))->index);
            }
            // *
            /// Ancestor </dd></dl>
            if (w_input->anc_node != NULL)
            {
                w_output->anc_node=(out_tree->m_node+w_input->anc_node->index);
            }
            else
            {
                w_output->anc_node=NULL;
            }
        }
        
        // **
        /// Copy of l_node::lat_node pointers
        if (l_nodes_ptr>0 && w_input->lat_node!=NULL)
        {
            w_output->lat_node=out_tree->m_node+(w_input->lat_node->index);
        }
        
        // **
        /// Copy of l_node::g_nodes pointers 
        if (g_nodes_ptr>0)
        {
            for (j=0; j>w_input->n_nodes; ++j)
            {
                *(w_output->children+j)=out_tree->m_node+((*(w_input->children+j))->index);
            }
        }
        
        // **
        /// Relinkage of \ref s_node::l_nodes when it is necessary</dd></dl>
        if (relink==1 && w_output->conts!=NULL && w_output->conts->l_nodes->index == i)
            w_output->conts->l_nodes=w_output;
    }
    
    // ***
    /// Tree values </dd></dl>
    out_tree->n_leaves=in_tree->n_leaves;
    out_tree->n_gleaves=in_tree->n_gleaves;
    out_tree->gen_time=in_tree->gen_time;
    out_tree->Ne=in_tree->Ne;
    out_tree->mu=in_tree->mu;
    out_tree->species_tree=in_tree->species_tree;
    out_tree->gene_tree=in_tree->gene_tree;
    
    if (in_tree->root != NULL && in_tree->n_nodes>1)
        out_tree->root=out_tree->m_node+(in_tree->root->index);
    
    return NO_ERROR;
}

long int CopyStoLTree(s_tree *sp_tree, l_tree *locus_tree, int reindex)
{
    int i=0,j=0;
    s_node *w_snode=NULL;
    l_node *w_lnode=NULL;
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Input tree memory control
    
    // *
    /// Tree equivalence 
    if (sp_tree->max_children!=locus_tree->max_children||sp_tree->n_gleaves!=locus_tree->n_gleaves||sp_tree->n_leaves!=locus_tree->n_leaves||sp_tree->n_nodes!=locus_tree->n_nodes)
        return MEM_ERROR;
    if (sp_tree->root==sp_tree->m_node)//Pre-order
    {
        ReindexSTree(sp_tree, 1);
        CopyStoLNodes(locus_tree->m_node,sp_tree->root,&i);///Do it using recursion
        switch (reindex)
        {
            case 1:
                ReindexSTree(sp_tree, 0);
                break;
        }
        
    }
    else
    {
        // *
        /// Tree copy loop
        for (i=0;i<sp_tree->n_nodes;++i)
        {
            w_snode=sp_tree->m_node+i;
            w_lnode=locus_tree->m_node+i;
            w_lnode->sp_index=w_snode->sp_index;
            w_lnode->n_child=w_snode->n_child;
            w_lnode->n_nodes=w_snode->n_replicas;
            w_lnode->Ne=w_snode->Ne;
            w_lnode->n_gen=w_snode->n_gen;
            w_lnode->time=w_snode->time;
            w_lnode->gen_length=w_snode->gen_length;
            w_lnode->conts=w_snode;
            w_lnode->mu_mult=w_snode->mu_mult;
            w_lnode->gtime_mult=w_snode->gtime_mult;
            w_snode->n_lnodes=1;
            w_snode->l_nodes=w_lnode;
            
            for (j=0;j<w_lnode->n_child; ++j)
            {
                *(w_lnode->children+j)=(locus_tree->m_node+((*(w_snode->children+j))->index));
            }
            
            if (w_snode->anc_node != NULL)
                w_lnode->anc_node=(locus_tree->m_node+(w_snode->anc_node->index));
            
        }
    }
    
    // *
    /// Tree info copy </dd></dl>
    locus_tree->root=sp_tree->root->l_nodes;
    locus_tree->gen_time=sp_tree->gen_time;
    locus_tree->Ne=sp_tree->Ne;
    locus_tree->mu=sp_tree->mu;
    locus_tree->species_tree=sp_tree;
    sp_tree->locus_tree=locus_tree;
    
    return NO_ERROR;
}


//\cond DOXYGEN_EXCLUDE
//// ** Tree edition ** //
//
//long int CleanlossesLTree(l_tree *locus_tree)
//{
//    int n_deletions=0;
//    int n_leaves=0;
//    if (locus_tree->m_node!=NULL || locus_tree->root == NULL)
//        return MEM_ERROR;
//    
//    CleanlossesLNodes(locus_tree->root,&locus_tree->root,&n_deletions,&n_leaves);
//    
//    locus_tree->n_nodes-=n_deletions;
//    locus_tree->n_leaves=n_leaves;
//    
//    return NO_ERROR;
//}
//
//\endcond
// ** Tree reset ** //

long int ResetSTreeSimL (s_tree *tree)
{
    int i=0;
    s_node * w_node=NULL;
    
    if (tree->m_node==NULL)
        return MEM_ERROR;
    
    for (i=0;i<tree->n_nodes;++i)
    {
        w_node=tree->m_node+i;
        w_node->index=i;
        w_node->l_nodes=NULL;
        w_node->n_lnodes=0;
    }
    tree->locus_tree=NULL;
    
    return NO_ERROR;
}

long int ResetGTree (g_tree *tree)
{
    int i=0,j=0;
    g_node * w_node=NULL;
    
    if (tree->m_node==NULL)
        return MEM_ERROR;
    
    for (i=0;i<tree->n_nodes;++i)
    {
        w_node=tree->m_node+i;
        
        w_node->sp_index=0;
        w_node->replica=0;
        w_node->paralog=0;
        w_node->n_gen=0;
        w_node->bl=0;
        w_node->gen_length=0.0;
        for (j=0;j<tree->max_children;++j)
        {
            *(w_node->children+j)=NULL;
        }
        w_node->anc_node=NULL;
        
    }
    
    return NO_ERROR;
}

// ** Tree deletion ** //

void FreeSTree(s_tree ** tree)
{
    register int i;
    s_node *w_node=NULL;
    
    // ***
    /// <dl><dt> Function structure </dt><dd>
    
    // **
    /// Frees s_nodes by \ref FreeSNodes (recursive) if the nodes are note allocated as an array.
    if ((*tree)->m_node==NULL)
    {
        FreeSNodes((*tree)->root);
    }
    else
    {
        // **
        /// <dl><dt>Iterative node loop (array of nodes)</dt><dd>
        for (i=0; i<(*tree)->n_nodes;++i)
        {
            w_node=(*tree)->m_node+i;
            
            // *
            /// Frees s_node::children </dd></dl>
            if (w_node->children !=NULL)
            {
                free(w_node->children);
                w_node->children=NULL;
            }
        }
        // **
        /// Frees (*tree) nodes
        free((*tree)->m_node);
        (*tree)->m_node=NULL;
        
    }
    
    // **
    /// Frees (*tree) </dd></dl>
    free((*tree));
    (*tree)=NULL;
}

void FreeLTree(l_tree ** tree)
{
    register int i;
    l_node *w_node=NULL;
    
    // ***
    /// <dl><dt> Function structure </dt><dd>
    
    // **
    /// Frees l_nodes by \ref FreeLNodes (recursive) if the nodes are note allocated as an array.
    if ((*tree)->m_node==NULL)
    {
        FreeLNodes((*tree)->root,1);
    }
    else
    {
        // **
        /// <dl><dt>Iterative node loop (array of nodes)</dt><dd>
        for (i=0; i<(*tree)->n_nodes;++i)
        {
            w_node=(*tree)->m_node+i;
            
            // *
            /// Frees l_node::g_nodes, l_node::children, l_node::n_ilin and l_node::n_olin </dd></dl>
            if (w_node->g_nodes !=NULL)
            {
                free(w_node->g_nodes);
                w_node->g_nodes=NULL;
            }
            if (w_node->children !=NULL)
            {
                free(w_node->children);
                w_node->children=NULL;
            }
            if (w_node->i_probs !=NULL)
            {
                free(w_node->i_probs);
                w_node->i_probs=NULL;
            }
            if (w_node->i_combprobs !=NULL)
            {
                free(w_node->i_combprobs);
                w_node->i_combprobs=NULL;
            }
            if (w_node->o_probs !=NULL)
            {
                free(w_node->o_probs);
                w_node->o_probs=NULL;
            }
        }
        // **
        /// Frees tree nodes
        free((*tree)->m_node);
        (*tree)->m_node=NULL;
        
    }
    
    // **
    /// Frees tree </dd></dl>
    free(*tree);
    *tree=NULL;
}

void FreeGTree(g_tree ** tree, int complete)
{
    int i;
    g_node * w_node=NULL;
    
    // ***
    /// <dl><dt> Function structure </dt><dd>
    
    // **
    /// Frees g_nodes by \ref FreeGNodes (recursive) if the nodes are note allocated as an array.
    if ((*tree)->m_node==NULL)
    {
        FreeGNodes((*tree)->root,complete);
    }
    else
    {
        // **
        /// <dl><dt>Iterative node loop (array of g_nodes)</dt><dd>
        for (i=0; i<(*tree)->n_nodes;++i)
        {
            w_node=(*tree)->m_node+i;
            
            // *
            /// Frees g_node::children</dd></dl>
            if (w_node->children !=NULL)
            {
                free(w_node->children);
                w_node->children=NULL;
            }
        }
        
        // **
        /// Frees (*tree) nodes
        free((*tree)->m_node);
        (*tree)->m_node=NULL;
    }
    
    // **
    /// Frees the (*tree) or resets the root (complete argument)</dd></dl>
    if (complete==1)
    {
        free((*tree));
        (*tree)=NULL;
    }
    else
    {
        (*tree)->n_nodes=1;
        (*tree)->root->index=0;
        (*tree)->root->sp_index=0;
        (*tree)->root->replica=0;
        (*tree)->root->paralog=0;
        (*tree)->root->n_gen=0;
        (*tree)->root->n_child=0;
        (*tree)->root->bl=0;
        (*tree)->root->gen_length=0.0;
        (*tree)->root->contl=NULL;
        (*tree)->root->conts=NULL;
        (*tree)->root->anc_node=NULL;
    }
}

void FreePeriods(period * periods, int n_periods)
{
    int i=0;
    
    for (i=0; i<n_periods; ++i)
    {
        if ((periods+i)->l_nodes!=NULL)
            free((periods+i)->l_nodes);
    }
    
    free(periods);
}

// *** Names memory manage *** //

name_c * NewNames (int n_names, int max_lname)
{
    name_c * names=NULL;
    
    // *
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Control of max_lname variable
    if (max_lname==0)
        max_lname=MAX_NAME;
    
    // *
    /// Allocation of name_c memory
    names=calloc(1,sizeof(name_c));
    ErrorReporter((long int) names,NULL);
    
    // *
    /// Inizialitation of name_c
    names->n_names=n_names;
    names->max_lname=max_lname;
    
    // *
    /// Allocation of the string of names (name_c::names)
    names->names=calloc((n_names+1)*max_lname,sizeof(char));
    ErrorReporter((long int)names->names,NULL);
    
    // *
    /// Inizialitation of the string </dd></dl>
    strncpy(names->names,"Internal node",max_lname);
    
    return (names);
    
}

void ReallocNames (name_c * names,int max_lname)
{
    char * new_names=NULL;
    int register i=0;
    int min_max_lname=0;
    
    // *
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Control of max_lname variable
    min_max_lname=(int)strlen(names->names)+1; // strlen is measuring only the first name, usually "Internal node"". +1 due to \0.
    
    if (max_lname<min_max_lname)
        max_lname=min_max_lname;
    
    // *
    /// Allocation of the new string of names
    new_names=calloc((names->n_names+1)*max_lname,sizeof(char));
    ErrorReporter((long int)new_names,NULL);
    
    // *
    /// Loop for the names copy
    for (i=0;i<=names->n_names;++i)
    {
        strncpy((new_names+i*max_lname),(names->names+i*names->max_lname),max_lname);
    }
    
    // *
    /// Frees the former string (name_c::names)</dd></dl>
    free(names->names);
    
    // *
    /// Copy of new values into name_c (name_c::names and name_c::max_lname)
    names->names=new_names;
    names->max_lname=max_lname;
    
}

void FreeNames(name_c **names)
{
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Frees name_c::names
    free ((*names)->names);
    (*names)->names=NULL;
    //*
    /// Frees name_c </dd></dl>
    free (*names);
    (*names)=NULL;
}

// *** Tree data modification *** //

// ** G and L tree reassociation ** //

long int MatchTreesMSC(l_tree *locus_tree, g_tree *gene_tree, int reset_gtree, int includelosses)
{
    int i=0;
    int j=0,next_leaf=0;
    l_node *w_lnode=NULL;
    g_node *w_gnode=NULL;
    
    // ****
    /// <dl><dt> Function structure </dt><dd>
    
    // ***
    ///  Errors control
    if (locus_tree==NULL || gene_tree==NULL || gene_tree->m_node==NULL)
    {
        return (MEM_ERROR);
    }
    else if (locus_tree->m_node==NULL)
    {
        if (locus_tree->root==NULL)
            return MEM_ERROR;
        
        CollapseLTree(locus_tree,1,0,0); //Post-order
    }
    if (includelosses>1)
        return UNEXPECTED_VALUE;
    
    // ***
    /// <dl><dt>Gene tree reset if it is required</dt><dd>
    
    // **
    /// <dl><dt>Gene tree loop</dt><dd>
    switch (reset_gtree)
    {
        case 1:
            for (i=0;i<gene_tree->n_nodes;++i)
            {
                // *
                /// Gene node info and pointers reset</dd></dl></dd></dl>
                w_gnode=gene_tree->m_node+i;
                w_gnode->index=i;
                w_gnode->sp_index=0;
                w_gnode->paralog=0;
                w_gnode->replica=0;
                w_gnode->conts=NULL;
                w_gnode->contl=NULL;
                w_gnode->n_child=0;
                w_gnode->anc_node=NULL;
                w_gnode->n_gen=0.0;
                w_gnode->bl=0.0;
                w_gnode->gen_length=0.0;
            }
            
            // **
            /// Gene tree root deletion</dd></dl>
            gene_tree->root=NULL;
            break;

    }
    
    // ***
    /// <dl><dt>Locus tree loop</dt><dd>
    for (i=0;i<locus_tree->n_nodes;++i)
    {
        w_lnode=locus_tree->m_node+i;
        
        // **
        /// <dl><dt>Leaf nodes association </dt><dd>
        switch (w_lnode->n_child)
        {
            case 0:
                // *
                /// The l_node::n_nodes is reset using l_node::conts s_node::n_replicas
                switch (w_lnode->kind_node)
                {
                    case LOSS:
                    case RTRFR:
                    case RGC:
                        w_lnode->n_nodes=includelosses;
                        break;
                    default:
                        w_lnode->n_nodes=w_lnode->conts!=NULL?w_lnode->conts->n_replicas:w_lnode->n_ilin; // If there is no species tree node (fixed locus tree), the locus tree doesn't change, and therefore we can use the number of input lineages in the same way as conts->n_replicas
                        break;
                }
                // *
                /// l_node::g_nodes are filled with g_node leaves by a loop with as much iterations as l_node::n_nodes, and each of these \ref g_node "g_nodes" is set with the same l_node::sp_index and with it own g_node::replica. </dd></dl></dd></dl>
                for (j=0;j<w_lnode->n_nodes;++j)
                {
                    w_gnode=gene_tree->m_node+next_leaf;
                    ++next_leaf;
                    *(w_lnode->g_nodes+j)=w_gnode;
                    w_gnode->sp_index=w_lnode->sp_index;
                    w_gnode->paralog=w_lnode->paralog;
                    w_gnode->replica=j;
                    w_gnode->contl=w_lnode;
                    w_gnode->conts=w_lnode->conts;
                    w_gnode->n_gen=w_lnode->n_gen;
                }
                break;
                
            default:
                w_lnode->n_nodes=0;
                break;
        }
    }
    
    gene_tree->species_tree=locus_tree->species_tree;
    gene_tree->locus_tree=locus_tree;
    gene_tree->gen_time=locus_tree->gen_time;
    locus_tree->gene_tree=gene_tree;
    
    return (NO_ERROR);
}

long int MatchTreesMLC(l_tree *locus_tree, g_tree *gene_tree, int reset_gtree)
{
    int i=0,j=0,next_leaf=0;
    l_node *w_lnode=NULL;
    g_node *w_gnode=NULL;
    
    // ****
    /// <dl><dt> Function structure </dt><dd>
    
    // ***
    ///  Errors control
    if (locus_tree==NULL || gene_tree==NULL || gene_tree->m_node==NULL)
    {
        return (MEM_ERROR);
    }
    else if (locus_tree->m_node==NULL)
    {
        if (locus_tree->root == NULL)
            return MEM_ERROR;
        
        CollapseLTree(locus_tree,1,0,1); //Post-order
    }
    
    // ***
    /// <dl><dt>Gene tree reset if it is required</dt><dd>
    
    // **
    /// <dl><dt>Gene tree loop</dt><dd>
    switch (reset_gtree)
    {
        case 1:
            for (i=0;i<gene_tree->n_nodes;++i)
            {
                // *
                /// Gene node info and pointers reset</dd></dl></dd></dl>
                w_gnode=gene_tree->m_node+i;
                w_gnode->index=i;
                w_gnode->sp_index=0;
                w_gnode->paralog=0;
                w_gnode->replica=0;
                w_gnode->conts=NULL;
                w_gnode->contl=NULL;
                w_gnode->n_child=0;
                w_gnode->anc_node=NULL;
                w_gnode->n_gen=0.0;
                w_gnode->bl=0.0;
                w_gnode->gen_length=0.0;
            }
            
            // **
            /// Gene tree root deletion</dd></dl>
            gene_tree->root=NULL;
            break;
    }
    
    // ***
    /// <dl><dt>Locus tree loop</dt><dd>
    for (i=0;i<locus_tree->n_nodes;++i)
    {
        w_lnode=locus_tree->m_node+i;
        
        //**
        /// Adding fixed number of output lineages for bounded coalescences
        if(w_lnode->anc_node!=NULL && ((w_lnode->anc_node->kind_node==DUP || w_lnode->anc_node->kind_node==TRFR || w_lnode->anc_node->kind_node==GC) && *(w_lnode->anc_node->children+1)==w_lnode) && w_lnode->fmax_nlin>0) //bounded and with lineages
            w_lnode->n_olin=1;
        else //=> nodes where we have to sample the number of lineages or it is 0.
            w_lnode->n_olin=0;
        
        // **
        /// <dl><dt>Leaf nodes association </dt><dd>
        switch (w_lnode->n_child)
        {
            case 0:
                // *
                /// l_node::g_nodes are filled with g_node leaves by a loop with as much iterations as l_node::n_nodes, and each of these \ref g_node "g_nodes" is set with the same l_node::sp_index and with it own g_node::replica. </dd></dl></dd></dl>
                for (j=0;j<w_lnode->n_ilin;++j)
                {
                    w_gnode=gene_tree->m_node+next_leaf;
                    ++next_leaf;
                    *(w_lnode->g_nodes+j)=w_gnode;
                    w_gnode->sp_index=w_lnode->sp_index;
                    w_gnode->paralog=w_lnode->paralog;
                    w_gnode->replica=j;
                    w_gnode->contl=w_lnode;
                    w_gnode->conts=w_lnode->conts;
                    w_gnode->n_gen=w_lnode->n_gen;
                }
                break;
                
            default:
                w_lnode->n_ilin=0;
                break;
        }
        w_lnode->n_nodes=w_lnode->n_ilin;
    }
    
    gene_tree->species_tree=locus_tree->species_tree;
    gene_tree->locus_tree=locus_tree;
    gene_tree->gen_time=locus_tree->gen_time;
    locus_tree->gene_tree=gene_tree;
    
    return (NO_ERROR);
}

// ** Tree conversion ** //

long int CollapseSTree (s_tree * in_tree, int post_order)
{
    s_node * m_node=NULL, * p_m_node=NULL, * w_node=NULL, * p_root=NULL;
    int index=0,i=0;
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Input tree root control
    if (in_tree->root == NULL)
        return MEM_ERROR;
    
    // *
    /// New s_node array allocation by \ref NewSTree
    m_node=NewSNodes(in_tree->n_nodes, in_tree->max_children);
    
    // *
    /// Unreference previous (in case it exists) input tree m_node
    if (in_tree->m_node!=NULL)
    {
        p_m_node=in_tree->m_node;
        in_tree->m_node=NULL;
    }
    else
    {
        p_root=in_tree->root;
    }
    
    // *
    /// Reordenation in a post or pre-order and node copy by \ref PostCollapseSNodes or \ref PreCollapseSNodes according to the post_order variable
    if (post_order==0)
    {
        PreReorderSNodes(in_tree->root, &index);
        PreCollapseSNodes(m_node,in_tree->root);
        
        // *
        /// Asigns the new root to the tree.
        in_tree->root=m_node; //The first node is the root.
        
    }
    else
    {
        PostReorderSNodes(in_tree->root, &index);
        PostCollapseSNodes(m_node,in_tree->root);
        
        // Asigns the new root to the tree.
        in_tree->root=(m_node+in_tree->n_nodes-1); //The last node is the root.
    }
    
    // *
    /// Asigns the new nodes to the tree
    in_tree->m_node=m_node;
    
    // *
    /// Frees the memory of the previous tree</dd></dl>
    if (p_m_node!=NULL)
    {
        for (i=0; i<in_tree->n_nodes;++i)
        {
            w_node=p_m_node+i;
            if (w_node->children !=NULL)
            {
                free(w_node->children);
                w_node->children=NULL;
            }
        }
        free(p_m_node);
        p_m_node=NULL;
    }
    else
    {
        FreeSNodes(p_root);
    }
    
    return NO_ERROR;
    
}

long int ReindexSTree (s_tree * in_tree, int post_order)
{
    int i=0;
    if (post_order==1)
        PostReorderSNodes(in_tree->root,&i);
    else
        PreReorderSNodes(in_tree->root, &i);
    
    return NO_ERROR;
}

long int CollapseLTree (l_tree * in_tree, int post_order, int relink, int probs)
{
    l_node * m_node=NULL, * p_m_node=NULL, *w_node=NULL, *p_root=NULL;
    int index=0, i=0;
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Input tree root control
    if (in_tree->root == NULL)
        return MEM_ERROR;
    
    // *
    /// New l_node array allocation by \ref NewLTree
    m_node=NewLNodes(in_tree->n_nodes, in_tree->n_gleaves, in_tree->max_children);
    
    // *
    /// Unreference previous (in case it exists) input tree m_node
    if (in_tree->m_node!=NULL)
    {
        p_m_node=in_tree->m_node;
        in_tree->m_node=NULL;
    }
    else
    {
        p_root=in_tree->root;
    }
    
    // *
    /// Reordenation in a post or pre-order and node copy by \ref PostCollapseLNodes or \ref PreCollapseLNodes according to the post_order variable
    if (post_order==0)
    {
        PreReorderLNodes(in_tree->root, &index);
        PreCollapseLNodes(m_node,in_tree->root, in_tree->n_gleaves, relink,probs);

        // *
        /// Asigns the new root to the tree.
        in_tree->root=m_node; //The first node is the root.
    }
    else
    {
        PostReorderLNodes(in_tree->root, &index);
        PostCollapseLNodes(m_node,in_tree->root, in_tree->n_gleaves, relink,probs);
        
        // Asigns the new root to the tree.
        in_tree->root=(m_node+in_tree->n_nodes-1); //The last node is the root.
    }
    
    // *
    /// Relinkage of \ref l_node::lat_node and the \ref s_node::l_nodes when it is necessary.
    if (relink==1)
    {
        for (i=0;i<in_tree->n_nodes;++i)
        {
            w_node=(m_node+i);
            if (w_node->conts!=NULL && w_node->conts->l_nodes->index == i)
                w_node->conts->l_nodes=w_node;
            if (w_node->lat_node!=NULL)
                w_node->lat_node=(m_node+w_node->lat_node->index);            
        }
    }

    // *
    /// Asigns the new nodes to the tree
    in_tree->m_node=m_node;
    
    // *
    /// Frees the memory of the previous tree</dd></dl>
    if (p_m_node!=NULL)
    {
        for (i=0; i<in_tree->n_nodes;++i)
        {
            w_node=p_m_node+i;
            if (w_node->children !=NULL)
            {
                free(w_node->children);
                w_node->children=NULL;
            }
            if (w_node->i_probs!=NULL)
            {
                free(w_node->i_probs);
                w_node->i_probs=NULL;
            }
            if (w_node->i_combprobs!=NULL)
            {
                free(w_node->i_combprobs);
                w_node->i_combprobs=NULL;
            }
            if (w_node->o_probs !=NULL)
            {
                free(w_node->o_probs);
                w_node->o_probs=NULL;
            }
            if (w_node->g_nodes !=NULL)
            {
                free(w_node->g_nodes);
                w_node->g_nodes=NULL;
            }
        }
        free(p_m_node);
        p_m_node=NULL;
    }
    else
    {
        FreeLNodes(p_root,1);
    }
    
    return NO_ERROR;
}

long int CollapseResizeLTree (l_tree * in_tree, int new_nnodes, int new_nleaves, int post_order, int relink, int probs)
{
    l_node * m_node=NULL, * p_m_node=NULL, *w_node=NULL, *p_root=NULL;
    int index=0, i=0;
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Input tree root control
    if (in_tree->root == NULL)
        return MEM_ERROR;
    
    // *
    /// New l_node array allocation by \ref NewLTree
    m_node=NewLNodes(new_nnodes, in_tree->n_gleaves, in_tree->max_children);
    
    // *
    /// Unreference previous (in case it exists) input tree m_node
    if (in_tree->m_node!=NULL)
    {
        p_m_node=in_tree->m_node;
        in_tree->m_node=NULL;
    }
    else
    {
        p_root=in_tree->root;
    }
    
    // *
    /// Reordenation in a post or pre-order and node copy by \ref PostCollapseLNodes or \ref PreCollapseLNodes according to the post_order variable
    if (post_order==0)
    {
        PreReorderLNodes(in_tree->root, &index);
        PreCollapseLNodes(m_node,in_tree->root, in_tree->n_gleaves, relink,probs);
        
        // *
        /// Asigns the new root to the tree.
        in_tree->root=m_node; //The first node is the root.
    }
    else
    {
        PostReorderLNodes(in_tree->root, &index);
        PostCollapseLNodes(m_node,in_tree->root, in_tree->n_gleaves, relink,probs);
        
        // Asigns the new root to the tree.
        in_tree->root=(m_node+new_nnodes-1); //The last node is the root.
    }
    
    // *
    /// Relinkage of \ref l_node::lat_node and the \ref s_node::l_nodes when it is necessary.
    if (relink==1)
    {
        for (i=0;i<new_nnodes;++i)
        {
            w_node=(m_node+i);
            if (w_node->conts!=NULL && w_node->conts->l_nodes->index == i)
                w_node->conts->l_nodes=w_node;
            if (w_node->lat_node!=NULL)
                w_node->lat_node=(m_node+w_node->lat_node->index);
        }
    }
    
    // *
    /// Asigns the new nodes to the tree
    in_tree->m_node=m_node;
    
    // *
    /// Frees the memory of the previous tree</dd></dl>
    if (p_m_node!=NULL)
    {
        for (i=0; i<in_tree->n_nodes;++i)
        {
            w_node=p_m_node+i;
            if (w_node->children !=NULL)
            {
                free(w_node->children);
                w_node->children=NULL;
            }
            if (w_node->i_probs!=NULL)
            {
                free(w_node->i_probs);
                w_node->i_probs=NULL;
            }
            if (w_node->i_combprobs!=NULL)
            {
                free(w_node->i_combprobs);
                w_node->i_combprobs=NULL;
            }
            if (w_node->o_probs !=NULL)
            {
                free(w_node->o_probs);
                w_node->o_probs=NULL;
            }
            if (w_node->g_nodes !=NULL)
            {
                free(w_node->g_nodes);
                w_node->g_nodes=NULL;
            }
        }
        free(p_m_node);
        p_m_node=NULL;
    }
    else
    {
        FreeLNodes(p_root,1);
    }
    
    in_tree->n_nodes=new_nnodes;
    in_tree->n_leaves=new_nleaves;
    
    return NO_ERROR;
}

long int CollapseGTree (g_tree * in_tree, int post_order, int relink)
{
    g_node * m_node=NULL, * p_m_node=NULL, *w_node=NULL, *p_root=NULL;
    int index=0, i=0;
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Input tree root control
    if (in_tree->root == NULL)
        return MEM_ERROR;
    
    // *
    /// New g_node array allocation by \ref NewGTree
    m_node=NewGNodes(in_tree->n_nodes,in_tree->max_children);
    
    // *
    /// Unreference previous (in case it exists) input tree m_node
    if (in_tree->m_node!=NULL)
    {
        p_m_node=in_tree->m_node;
        in_tree->m_node=NULL;
    }
    else
    {
        p_root=in_tree->root;
    }
    
    // *
    /// Reordenation in a post or pre-order and node copy by \ref PostCollapseGNodes or \ref PreCollapseGNodes according to the post_order variable
    if (post_order==0)
    {
        PreReorderGNodes(in_tree->root, &index);
        PreCollapseGNodes(m_node,in_tree->root);
        
        // *
        /// Asigns the new root to the tree.
        in_tree->root=m_node; //The first node is the root.
    }
    else
    {
        PostReorderGNodes(in_tree->root, &index);
        PostCollapseGNodes(m_node,in_tree->root);
        
        // Asigns the new root to the tree.
        in_tree->root=(m_node+in_tree->n_nodes-1); //The last node is the root.
    }
    
    if (relink==1)
    {
        RelinkLGTrees(in_tree->locus_tree,m_node);
    }
    
    // *
    /// Asigns the new nodes to the tree
    in_tree->m_node=m_node;
    
    // *
    /// Frees the memory of the previous tree</dd></dl>
    if (p_m_node!=NULL)
    {
        for (i=0; i<in_tree->n_nodes;++i)
        {
            w_node=p_m_node+i;
            if (w_node->children !=NULL)
            {
                free(w_node->children);
                w_node->children=NULL;
            }
        }
        free(p_m_node);
        p_m_node=NULL;
    }
    else
    {
        FreeGNodes(p_root,1);
    }
    
    return NO_ERROR;
}
// ** Gene tree branch length modification ** //

long int TemporalizeLTree(l_tree *tree)
{
    unsigned int i=0;
    
    tree->root->time=0;
    for (i=0; i<tree->root->n_child;++i)
    {
        TemporalizeLNodes(*(tree->root->children+i), tree->gen_time);
    }
    return NO_ERROR;
}

long int Rateheter_lineagespec(s_tree *sp_tree, double alpha, gsl_rng *seed, FILE * gammadump)
{
    s_node *w_node;
    int i=0;
    double sampled_gamma=0;
    
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Input tree control
    if (sp_tree->root==NULL || sp_tree->m_node == NULL)
        return MEM_ERROR;
    
    // *
    /// s_node::mu_multi sample using \ref RndGammaE1 </dd></dl>
    
    if (gammadump!=NULL)
    {
        for (i=0; i<sp_tree->n_nodes;++i)
        {
            w_node=sp_tree->m_node+i;
            sampled_gamma=gsl_ran_gamma(seed,alpha,1/alpha);
            fprintf(gammadump,"%.8lf\t",sampled_gamma);
            w_node->mu_mult*=sampled_gamma;
        }
        fprintf(gammadump,"\n");
    }
    else
    {
        for (i=0; i<sp_tree->n_nodes;++i)
        {
            w_node=sp_tree->m_node+i;
            w_node->mu_mult*=gsl_ran_gamma(seed,alpha,1/alpha);
        }
    }

    
    return NO_ERROR;
}

long int Rateheter_genespec(l_tree *locus_tree,double value)
{
    
    locus_tree->mu=value;
    return NO_ERROR;
}

long int Rateheter_GTbranchspec(g_tree *tree, double alpha, gsl_rng *seed, FILE * gammadump)
{
    g_node *w_node;
    int i=0;
    double sampled_gamma=0;
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Input tree control
    if (tree->root==NULL || tree->m_node == NULL)
        return MEM_ERROR;
    
    // *
    /// g_node::bl modification using sampled multis from \ref RndGammaE1 </dd></dl>
    if (gammadump!=NULL)
    {
        for (i=0; i<tree->n_nodes;++i)
        {
            w_node=tree->m_node+i;
            sampled_gamma=gsl_ran_gamma(seed,alpha,1/alpha);
            fprintf(gammadump,"%.8lf\t",sampled_gamma);
            w_node->bl*=sampled_gamma;
        }
        fprintf(gammadump,"\n");
    }
    else
    {
        for (i=0; i<tree->n_nodes;++i)
        {
            w_node=tree->m_node+i;
            w_node->bl*=gsl_ran_gamma(seed,alpha,1/alpha);
        }
    }

    
    return NO_ERROR;
}

// *** Tree features obtention *** //

long int Count_duplications(l_tree *tree, int *n_dup)
{
    int i=0;
    l_node *w_lnode=NULL;
    
    if (tree->root==NULL)
    {
        return MEM_ERROR;
    }
    
    *n_dup=0;
    
    if (tree->m_node==NULL)
    {
        *n_dup=Count_duplications_lnodes(tree->root);
    }
    else
    {
        for (i=0; i<tree->n_nodes; ++i)
        {
            w_lnode=tree->m_node+i;
            if (w_lnode->kind_node==DUP)
                ++(*n_dup);
        }
    }
    
    return NO_ERROR;
        
}

long int Count_losses(l_tree *tree, int *n_loss)
{
    int i=0;
    l_node *w_lnode=NULL;
    
    if (tree->root==NULL)
    {
        return MEM_ERROR;
    }
    
    *n_loss=0;
    
    if (tree->m_node==NULL)
    {
        *n_loss=Count_losses_lnodes(tree->root);
    }
    else
    {
        for (i=0; i<tree->n_nodes; ++i)
        {
            w_lnode=tree->m_node+i;
            if (w_lnode->kind_node==LOSS)
                ++(*n_loss);
        }
    }
    
    return NO_ERROR;
    
}

long int Count_transfers(l_tree *tree, int *n_trans)
{
    int i=0;
    l_node *w_lnode=NULL;
    
    if (tree->root==NULL)
    {
        return MEM_ERROR;
    }
    
    *n_trans=0;
    
    if (tree->m_node==NULL)
    {
        *n_trans=Count_transfers_lnodes(tree->root);
    }
    else
    {
        for (i=0; i<tree->n_nodes; ++i)
        {
            w_lnode=tree->m_node+i;
            if (w_lnode->kind_node==TRFR)
                ++(*n_trans);
        }
    }
    
    return NO_ERROR;
    
}

long int Count_gc(l_tree *tree, int *n_gc)
{
    int i=0;
    l_node *w_lnode=NULL;
    
    if (tree->root==NULL)
    {
        return MEM_ERROR;
    }
    *n_gc=0;
    if (tree->m_node==NULL)
    {
        *n_gc=Count_gc_lnodes(tree->root);
    }
    else
    {
        for (i=0; i<tree->n_nodes; ++i)
        {
            w_lnode=tree->m_node+i;
            if (w_lnode->kind_node==GC)
                ++(*n_gc);
        }
    }
    
    return NO_ERROR;
    
}

long int Measure_ST_height(s_tree *tree,long double *height, int unit)
{
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Input tree memory control
    if (tree->root == NULL)
        return MEM_ERROR;
    *height=0;
    
    // *
    /// Recursive function selection depending on the n_gen unit.</dd></dl>
    switch(unit)
    {
        case CU:
            *height=Measure_s_node_cu_height(tree->root,tree->Ne);
            return NO_ERROR;
            break;
        case GL:
            *height=Measure_s_node_gl_height(tree->root);
            return NO_ERROR;
            break;
        default:
            return UNEXPECTED_VALUE;
            break;
    }
}

long int Measure_ST_length(s_tree *tree,long double *length,int unit)
{
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Input tree memory control
    if (tree->root == NULL)
        return MEM_ERROR;
    *length=0;
    
    // *
    /// Recursive function selection depending on the length unit.</dd></dl>
    switch(unit)
    {
        case CU:
            *length=Measure_s_node_cu_length(tree->root,tree->Ne);
            return NO_ERROR;
            break;
        case GL:
            *length=Measure_s_node_gl_length(tree->root);
            return NO_ERROR;
            break;
        default:
            return UNEXPECTED_VALUE;
            break;
    }

}

long int Measure_GT_height(g_tree *tree,long double *height,int unit)
{
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Input tree memory control
    if (tree->root == NULL)
        return MEM_ERROR;
    *height=0;
    
    // *
    /// Recursive function selection depending on the n_gen unit.</dd></dl>
    switch(unit)
    {
        case CU:
            *height=Measure_g_node_cu_height(tree->root,tree->locus_tree->Ne);
            return NO_ERROR;
            break;
        case BL:
            *height=Measure_g_node_bl_height(tree->root);
            return NO_ERROR;
            break;
        default:
            return UNEXPECTED_VALUE;
            break;
    }
}

long int Measure_GT_length(g_tree *tree,long double *length,int unit)
{
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Input tree memory control
    if (tree->root == NULL)
        return MEM_ERROR;
    *length=0;
    
    // *
    /// Recursive function selection depending on the length unit.</dd></dl>
    switch(unit)
    {
        case CU:
            *length=Measure_g_node_cu_length(tree->root,tree->locus_tree->Ne);
            return NO_ERROR;
            break;
        case BL:
            *length=Measure_g_node_bl_length(tree->root);
            return NO_ERROR;
            break;
        default:
            return UNEXPECTED_VALUE;
            break;
    }
}

long int MeasureMRCAEVdistance(g_tree *wg_tree,int event,double **distances, int n_dup, int unit)
{
    l_node * w_lnode=NULL, *b_lnode=NULL, *p_blnode=NULL;
    g_node * w_gnode=NULL, *b_gnode=NULL, **wg_pointers=NULL;
    int i=0, j=0, it=0, pn_nodes=0,j_bgnode=0, next_b=0,d=0,done=0;
    *distances=realloc(*distances,n_dup*sizeof(double));
    if (*distances==NULL)
        return MEM_ERROR;
    
    // ****
    /// <dl><dt> Function structure </dt><dd>
    
    // ***
    /// Input tree memory control
    if (wg_tree->root == NULL)
        return MEM_ERROR;
    if (unit!=GL)
    {
        fprintf(stderr, "The only implemented distance for the MRCA/event error is number of generations (GL)\n");
        return SETTINGS_ERROR;
    }
    // ***
    /// <dl><dt>Locus tree node loop</dt><dd>
    
    for (i=0; i<wg_tree->locus_tree->n_nodes;++i)
    {
        w_lnode=wg_tree->locus_tree->m_node+i;

        // **
        /// Discarding nodes where the measure doesn't make any sense
        if (w_lnode->kind_node!=event) //Different events
            continue;
        else if (w_lnode->n_ilin<=1) //Less than two lineages
        {
            *(*distances+d)=-1;
            ++d;
            continue;
        }

        // **
        /// Memory realocation if necessary
        if (pn_nodes!=w_lnode->n_ilin)
        {
            wg_pointers=realloc(wg_pointers, sizeof(g_node*)*w_lnode->n_ilin);
            if (wg_pointers==NULL)
                return MEM_ERROR;
        }
        
        memmove(wg_pointers, w_lnode->g_nodes, sizeof(g_node*)*w_lnode->n_ilin);
        pn_nodes=w_lnode->n_ilin;
        
        // **
        /// Gene tree nodes loop, looking for one node coming alone from one locus tree lineage
        done=0;
        p_blnode=NULL;
        
        for (j=0; j<w_lnode->n_ilin; ++j)
        {
            it=0;
            w_gnode=*(wg_pointers+j);
            b_lnode=w_gnode->contl;
            while (b_lnode!=NULL && b_lnode->anc_node != w_lnode && it<MAX_IT)
            {
                b_lnode=b_lnode->anc_node;
                ++it;
            }
            if (it>=MAX_IT)
                return LOOP_ERROR;
            if (b_lnode==NULL)
                return MEM_ERROR;
            if (p_blnode!=b_lnode)
            {
                p_blnode=b_lnode;
                ++done;
            }
            if (b_lnode->n_olin==1) //This is a bounded l_node. If both are bounded, both would be valid(this is only to have one concrete lineage from one of the paralogs to compare to the lineage/s from the other, looking for the MRCA)
            {
                b_gnode=w_gnode;
                j_bgnode=j;
            }

        }
        
        // **
        /// It continues only if two locus tree lineages have output lineages
        if (done!=2)
        {
            *(*distances+d)=-1;
            ++d;
            continue;
        }
        
        done=0;
        
        // **
        /// Loop going back in time, looking for the first MRCA between the selected gene tree node and the rest 
        while (it<MAX_IT && done==0)
        {
            next_b=1;
            for (j=0; j<w_lnode->n_ilin;++j)
            {
                if (j==j_bgnode)
                    continue;
                
                w_gnode=*(wg_pointers+j);

                if (w_gnode->n_gen>b_gnode->n_gen) //We have to consider the ancestor of the candidate node, but the b_gnode can still be a good choice (next_b=0)
                {
                    *(wg_pointers+j)=w_gnode->anc_node;
                    next_b=0;
                }
                else if (w_gnode==b_gnode)
                {
                    *(*distances+d)=w_lnode->n_gen-w_gnode->n_gen; //GL
                    done=1;
                    break;
                }
            }
            if (next_b==1)
            {
                b_gnode=b_gnode->anc_node; //All the candidates below the t_node, so we need to consider its ancestor
            }
            ++it;
        }
        ++d;
    }

    free(wg_pointers);
    return NO_ERROR;
}

double ExpectedPrunedLtreeNleavesSNodes(s_tree *tree, double b_rate, double d_rate)
{
    return exp((b_rate-d_rate)*Measure_s_node_gl_height(tree->root))*tree->n_leaves;
}

long int CheckUltrametricitySTree(s_tree *tree)
{
    int i=0, is_set=0, n_leaves=0;
    double time=0;
    
    if (tree->m_node==NULL)
    {
        return CheckUltrametricitySNodes(tree->root);
    }
    
    for (i=0; i<tree->n_nodes;++i)
    {
        if (n_leaves==tree->n_leaves)
            break;
        switch ((tree->m_node+i)->n_child)
        {
            case 0:
                switch (is_set)
            {
                case 0:
                    is_set=1;
                    time=(tree->m_node+i)->time;
                    break;
                    
                default:
                    if (fabs((tree->m_node+i)->time-time)>FLOAT_PRECISION)
                        return UNEXPECTED_VALUE;
                    
                    break;
            }
                ++n_leaves;
                break;
        }
    }
    return NO_ERROR;
    
}

long int CheckUltrametricityLTree(l_tree *tree)
{
    int i=0, is_set=0, n_leaves=0;
    double time=0;
    
    if (tree->m_node==NULL)
    {
        return CheckUltrametricityLNodes(tree->root);
    }
    
    for (i=0; i<tree->n_nodes;++i)
    {
        if (n_leaves==tree->n_leaves)
            break;
        if ((tree->m_node+i)->n_child==0)
        {
            ++n_leaves;
            switch ((tree->m_node+i)->kind_node)
            {
                case SP:
                case DUP:
                case TRFR:
                case GC:
                    switch (is_set)
                {
                    case 0:
                        is_set=1;
                        time=(tree->m_node+i)->time;
                        break;
                        
                    default:
                        if (fabs((tree->m_node+i)->time-time)>FLOAT_PRECISION)
                            return UNEXPECTED_VALUE;
                        
                        break;
                }
                    
                    break;
            }
        }
    }
    return NO_ERROR;
    
}


// *** Tree I/O *** //

long int WriteSTree (s_tree *in_tree, name_c * names, int time, int int_labels)
{
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Writes the tree by \ref WriteSNodes (recursive)</dd></dl>
    if (in_tree->root==NULL)
        return (MEM_ERROR);
    else if(time>0)
    {
        if (in_tree->gen_time==0)
        {
            return(UNEXPECTED_VALUE);
        }
        switch (int_labels)
        {
            default:
                WriteSNodesTime(in_tree->root,names,in_tree->gen_time);
                break;
            case 1:
                WriteSNodesTimeIntlabel(in_tree->root,names,in_tree->gen_time);
                break;
        }
        
        return(NO_ERROR);
    }
    else
    {
        switch (int_labels)
        {
            default:
                WriteSNodesGen(in_tree->root,names);
                break;
            case 1:
                WriteSNodesGenIntlabel(in_tree->root,names);
                break;
        }
        return(NO_ERROR);
    }
    
}

long int WriteSTreeFile (FILE *file,s_tree *in_tree, name_c * names, int time, int int_labels)
{
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Writes the tree by \ref WriteSNodesFile (recursive)</dd></dl>
    if (in_tree->root==NULL)
        return (MEM_ERROR);
    else if(time>0)
    {
        if (in_tree->gen_time==0)
        {
            return(UNEXPECTED_VALUE);
        }
        switch (int_labels)
        {
            default:
                WriteSNodesFileTime(file,in_tree->root,names,in_tree->gen_time);
                break;
            case 1:
                WriteSNodesFileTimeIntlabel(file,in_tree->root,names,in_tree->gen_time);
                break;
        }
        
        return(NO_ERROR);
    }
    else
    {
        switch (int_labels)
        {
            default:
                WriteSNodesFileGen(file,in_tree->root,names);
                break;
            case 1:
                WriteSNodesFileGenIntlabel(file,in_tree->root,names);
                break;
        }

        return(NO_ERROR);
    }
    
}

long int WriteLTree (l_tree *in_tree, name_c * names, int time, int int_labels)
{
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Writes the tree by \ref WriteLNodes (recursive)</dd></dl>
    if (in_tree->root==NULL)
        return (MEM_ERROR);
    else if(time>0)
    {
        if (in_tree->gen_time==0)
        {
            return(UNEXPECTED_VALUE);
        }
        switch (int_labels)
        {
            default:
                WriteLNodesTime(in_tree->root,names,in_tree->gen_time);
                break;
            case 1:
                WriteLNodesTimeIntlabel(in_tree->root,names,in_tree->gen_time);
                break;
        }

        return(NO_ERROR);
    }
    else
    {
        switch (int_labels)
        {
            default:
                WriteLNodesGen(in_tree->root,names);
                break;
            case 1:
                WriteLNodesGenIntlabel(in_tree->root,names);
                break;
        }

        return(NO_ERROR);
    }
    
}

long int WriteLTreeFile (FILE *file,l_tree *in_tree, name_c * names, int time, int int_labels)
{
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Writes the tree by \ref WriteLNodesFile (recursive)</dd></dl>
    if (in_tree->root==NULL)
        return (MEM_ERROR);
    else if(time>0)
    {
        if (in_tree->gen_time==0)
        {
            return(UNEXPECTED_VALUE);
        }
        switch (int_labels)
        {
            default:
                WriteLNodesFileTime(file,in_tree->root,names,in_tree->gen_time);
                break;
            case 1:
                WriteLNodesFileTimeIntlabel(file,in_tree->root,names,in_tree->gen_time);
                break;
        }
        return(NO_ERROR);
    }
    else
    {
        switch (int_labels)
        {
            default:
                WriteLNodesFileGen(file,in_tree->root,names);
                break;
            case 1:
                WriteLNodesFileGenIntlabel(file,in_tree->root,names);
                break;
        }

        return(NO_ERROR);
    }
    
}

long int WriteDaughtersFile (FILE *file,l_tree *in_tree, name_c * names)
{
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Writes the daughters by \ref WriteDaughtersNodeFile (recursive)</dd></dl>
    if (in_tree->root==NULL)
        return (MEM_ERROR);
    else
    {
        WriteDaughtersNodesFile(file,in_tree->root,names);
        fseek(file, -1, SEEK_END);
        fprintf(file,"\n");
        return(NO_ERROR);
    }
    
}

long int WriteGTree (g_tree *in_tree, name_c * names, int int_labels)
{
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Writes the tree by \ref WriteGNodes (recursive)</dd></dl>
    if (in_tree->root==NULL)
        return (MEM_ERROR);
    else
    {
        switch (int_labels)
        {
            default:
                WriteGNodes(in_tree->root,names);
                break;
            case 1:
                WriteGNodesIntlabel(in_tree->root,names);
                break;
        }
        
        return(NO_ERROR);
    }
    
}

long int WriteGTreeStr (char * string, g_tree *in_tree, name_c * names, int int_labels)
{
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Writes the tree by \ref WriteGNodesStr (recursive)</dd></dl>
    if (in_tree->root==NULL)
        return (MEM_ERROR);
    else
    {
        switch (int_labels)
        {
            default:
                WriteGNodesStr(string,in_tree->root,names);
                break;
            case 1:
                WriteGNodesStrIntlabel(string,in_tree->root,names);
                break;
        }
        return(NO_ERROR);
    }
}

long int WriteGTreeFile (FILE * file, g_tree *in_tree, name_c * names, int int_labels)
{
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Writes the tree by \ref WriteGNodesFile (recursive)</dd></dl>
    if (in_tree->root==NULL)
        return (MEM_ERROR);
    else
    {
        switch (int_labels)
        {
            default:
                WriteGNodesFile(file,in_tree->root,names);
                break;
            case 1:
                WriteGNodesFileIntlabel(file,in_tree->root,names);
                break;
        }
        
        return(NO_ERROR);
    }
    
}

long int WriteMappingSL(s_tree *wsp_tree, l_tree *locus_tree, name_c *names, char *mapsl_outname)
{
    int i=0;
    l_node *wl_node=NULL;
    FILE *mapsl_outfile=NULL;
    
    if ((mapsl_outfile=fopen(mapsl_outname, "w"))==NULL)
    {
        return(IO_ERROR);
    }
    fprintf(mapsl_outfile,"Lt_node\tLt_paralog\tLt_event\tSt_node\n");
    
    PostReorderSNodes(wsp_tree->root,&i);
    
    for (i=0;i<locus_tree->n_nodes;++i)
    {
        wl_node=locus_tree->m_node+i;
        if (wl_node->conts!=NULL)
        {
            if (wl_node->n_child!=0) //Internal l_node
            {
                if (wl_node->conts->n_child!=0) //Everything internal
                {
                    switch (wl_node->kind_node)
                    {
                        case SP:
                            fprintf(mapsl_outfile,"%d\t%d\t%s\t%d\n",wl_node->index,wl_node->paralog,"Sp",wl_node->conts->index);
                            break;
                        case DUP:
                            fprintf(mapsl_outfile,"%d\t%d\t%s\t%d\n",wl_node->index,wl_node->paralog,"Dup",wl_node->conts->index);
                            break;
                        case TRFR:
                            fprintf(mapsl_outfile,"%d\t%d\t%s\t%d\n",wl_node->index,wl_node->paralog,"Transf",wl_node->conts->index);
                            break;
                        case GC:
                            fprintf(mapsl_outfile,"%d\t%d\t%s\t%d\n",wl_node->index,wl_node->paralog,"Gconv",wl_node->conts->index);
                            break;
                        default:
                            return UNEXPECTED_VALUE;
                            break;
                    }
                }
                else
                {
                    if (names==NULL)
                    {
                        switch (wl_node->kind_node)
                        {
                            case SP:
                                fprintf(mapsl_outfile,"%d\t%d\t%s\t\'%d\'\n",wl_node->index,wl_node->paralog,"Sp",wl_node->sp_index);
                                break;
                            case DUP:
                                fprintf(mapsl_outfile,"%d\t%d\t%s\t\'%d\'\n",wl_node->index,wl_node->paralog,"Dup",wl_node->sp_index);
                                break;
                            case TRFR:
                                fprintf(mapsl_outfile,"%d\t%d\t%s\t\'%d\'\n",wl_node->index,wl_node->paralog,"Transf",wl_node->sp_index);
                                break;
                            case GC:
                                fprintf(mapsl_outfile,"%d\t%d\t%s\t\'%d\'\n",wl_node->index,wl_node->paralog,"Gconv",wl_node->sp_index);
                                break;
                            default:
                                return UNEXPECTED_VALUE;
                                break;
                        }
                    }
                    else
                    {
                        switch (wl_node->kind_node)
                        {
                            case SP:
                                fprintf(mapsl_outfile,"%d\t%d\t%s\t\'%s\'\n",wl_node->index,wl_node->paralog,"Sp",(names->names+(wl_node->sp_index*names->max_lname)));
                                break;
                            case DUP:
                                fprintf(mapsl_outfile,"%d\t%d\t%s\t\'%s\'\n",wl_node->index,wl_node->paralog,"Dup",(names->names+(wl_node->sp_index*names->max_lname)));
                                break;
                            case TRFR:
                                fprintf(mapsl_outfile,"%d\t%d\t%s\t\'%s\'\n",wl_node->index,wl_node->paralog,"Transf",(names->names+(wl_node->sp_index*names->max_lname)));
                                break;
                            case GC:
                                fprintf(mapsl_outfile,"%d\t%d\t%s\t\'%s\'\n",wl_node->index,wl_node->paralog,"Gconv",(names->names+(wl_node->sp_index*names->max_lname)));
                                break;
                            default:
                                return UNEXPECTED_VALUE;
                                break;
                        }
                    }
                }
            }
            else if (wl_node->kind_node==SP)
            {
                if (names==NULL)
                    fprintf(mapsl_outfile,"\'%d_%d\'\t%d\tLeaf\t\'%d\'\n",wl_node->sp_index,wl_node->paralog,wl_node->paralog,wl_node->sp_index);
                else
                    fprintf(mapsl_outfile,"\'%s_%d\'\t%d\tLeaf\t\'%s\'\n",(names->names+(wl_node->sp_index*names->max_lname)),wl_node->paralog,wl_node->paralog,(names->names+(wl_node->sp_index*names->max_lname)));
            }
            else if (wl_node->kind_node==LOSS)
            {
                
                fprintf(mapsl_outfile,"\'Lost-%d_%d\'\t%d\tLoss\t%d\n",wl_node->conts->index,wl_node->paralog,wl_node->paralog,wl_node->conts->index);
            }
            else if (wl_node->kind_node==RTRFR)
            {
                
                fprintf(mapsl_outfile,"\'RTransf-%d_%d\'\t%d\tRTransf\t%d\n",wl_node->conts->index,wl_node->paralog,wl_node->paralog,wl_node->conts->index);
            }
            else
            {
                fprintf(mapsl_outfile,"\'Rgc-%d_%d\'\t%d\tRgc\t%d\n",wl_node->conts->index,wl_node->paralog,wl_node->paralog,wl_node->conts->index);
            }
            
        }
        else
            return UNEXPECTED_VALUE;
        
    }
    
    fclose(mapsl_outfile);
    
    return(NO_ERROR);
}

long int WriteMappingLG(g_tree *gene_tree, name_c *names, char *maplg_outname)
{
    FILE *maplg_outfile=NULL;
    g_node *wg_node=NULL;
    int i=0;
    
    if ((maplg_outfile=fopen(maplg_outname, "w"))==NULL)
    {
        return(IO_ERROR);
    }
    fprintf(maplg_outfile,"Gt_node\tLt_node\tLt_paralog\n");
    
    for (i=0;i<gene_tree->n_nodes;++i)
    {
        wg_node=gene_tree->m_node+i;
        if (wg_node->contl!=NULL)
        {
            if (wg_node->n_child!=0) //g_node internal
            {
                
                if (wg_node->contl->n_child!=0) // Everything internal
                {
                    fprintf(maplg_outfile,"%d\t%d\t%d\n",wg_node->index,wg_node->contl->index,wg_node->paralog);
                }
                else //External l_node
                {
                    switch (wg_node->contl->kind_node)
                    {
                        case LOSS:
                            fprintf(maplg_outfile,"%d\t\'Lost-%d_%d\'\t%d\n",wg_node->index,wg_node->conts->index,wg_node->paralog,wg_node->paralog);
                            break;
                        case RTRFR:
                            fprintf(maplg_outfile,"%d\t\'Rtransf-%d_%d\'\t%d\n",wg_node->index,wg_node->conts->index,wg_node->paralog,wg_node->paralog);
                            break;
                        case RGC:
                            fprintf(maplg_outfile,"%d\t\'Rgc-%d_%d\'\t%d\n",wg_node->index,wg_node->conts->index,wg_node->paralog,wg_node->paralog);
                            break;
                            
                        default:
                            if (names==NULL)
                                fprintf(maplg_outfile,"%d\t\'%d_%d\'\t%d\n",wg_node->index,wg_node->sp_index,wg_node->paralog,wg_node->paralog);
                            else
                                fprintf(maplg_outfile,"%d\t\'%s_%d\'\t%d\n",wg_node->index,(names->names+(wg_node->sp_index*names->max_lname)),wg_node->paralog,wg_node->paralog);
                            break;
                    }
                }
            }
            else //g_node and l_node external
            {
                switch (wg_node->contl->kind_node)
                {
                    case LOSS:
                        fprintf(maplg_outfile,"\'Lost-%d_%d_%d\'\t\'Lost-%d_%d\'\t%d\n",wg_node->conts->index,wg_node->paralog, wg_node->replica,wg_node->conts->index,wg_node->paralog,wg_node->paralog);
                        break;
                    case RTRFR:
                        fprintf(maplg_outfile,"\'Rtransf-%d_%d_%d\'\t\'Rtransf-%d_%d\'\t%d\n",wg_node->conts->index,wg_node->paralog, wg_node->replica,wg_node->conts->index,wg_node->paralog,wg_node->paralog);
                        break;
                    case RGC:
                        fprintf(maplg_outfile,"\'Rgc-%d_%d_%d\'\t\'Rgc-%d_%d\'\t%d\n",wg_node->conts->index,wg_node->paralog, wg_node->replica,wg_node->conts->index,wg_node->paralog,wg_node->paralog);
                        break;
                        
                    default:
                        if (names==NULL)
                            fprintf(maplg_outfile,"\'%d_%d_%d\'\t\'%d_%d\'\t%d\n",wg_node->sp_index,wg_node->paralog, wg_node->replica,wg_node->sp_index,wg_node->paralog,wg_node->paralog);
                        else
                            fprintf(maplg_outfile,"\'%s_%d_%d\'\t\'%s_%d\'\t%d\n",(names->names+(wg_node->sp_index*names->max_lname)),wg_node->paralog,wg_node->replica,(names->names+(wg_node->sp_index*names->max_lname)),wg_node->paralog,wg_node->paralog);
                        break;
                }
            }
        }
        else
            return (UNEXPECTED_VALUE);
        
    }
    
    fclose(maplg_outfile);
    return (NO_ERROR);
}

long int CheckNexusTree (char * tree)
{
	int i=0;
    int k=0,error=0,n_left=0,n_right=0,n_branches=0,n_nodes=0,node=0,f_leaf=0,leaf=0,comment=0;
    long int r_value=NO_ERROR;
    char current=0,previous='(',next=0;
    
	// **
    /// <dl><dt> Function structure </dt><dd>
    
    // **
    /// Initial tests
    if ((*tree)!='(')
    {
        error=1;
    }
    
    // **
    /// <dl><dt>Parsing loop, taking into account 3 characters (actual, pervious, next). It takes into account format problems, as well as structure problems:</dt><dd>
    
    
	for (i=0;i<strlen(tree);++i)
    {
        if (error == 1 || f_leaf==1)
        {
            // *
            /// Clade with only one leaf
            if (f_leaf==1)
            {
                fprintf(stderr,"\nERROR: This internal node does not have more than one leaf, Is there a missed comma?\n");
                f_leaf=0;
            }
            else
            {
                fprintf (stderr,"\nERROR: There is something wrong in the tree\n");
                error=0;
                
            }
            
			
			for (k=0; k<=i-1; k++)
				fprintf (stderr,"%c", tree[k]);
            
			fprintf (stderr," <- HERE");
            fflush(stderr);
            
            r_value=SETTINGS_ERROR;
            
        }
        
        switch (comment)
        {
            case 1:
                if (*(tree+i)==']')
                    comment=0;
                continue;
        }
        
        current=*(tree+i);
        previous=*(tree+(i-1>0?i-1:0));
        next=*(tree+i+1);
        
        switch (current)
        {
            case '(':
                if (!(previous == '(' || previous == ',')) /// "(" following something different to a "(" or ","
                    error=1;
                if (next!='(')
                {
                    ++n_nodes;
                    leaf=1;
                    if (isalnum(next)==0) /// "(" followed by alphanumeric.
                        error=1;
                }
                ++n_left;
                node=0;
                break;
                
            case ')':
                if (!(previous == ')' || previous == ']'|| isalnum(previous)!=0)) /// ")" following something different to ")", "]" or alphanumerics.
                    error=1;
                if (!(next == ')' || next== ',' || next== ':' || next== '[')) /// ")" followed by something different to a branch length, a comma or a comment
                    error=1;
                if (next==';' && i==strlen(tree)-2) /// Allowing the final ";"
                    error=0;
                ++n_right;
                ++n_nodes;
                if (node==0 && leaf==1)
                {
                    f_leaf=1;
                }
                leaf=0;
                node=0;
                break;
                
            case ',':
                if (!(previous == ')' || previous == ']' || isalnum(previous)!=0)) /// "," following something different to ")", "]" or alphanumerics.
                    error=1;
                if (next!='(') /// "," followed by something different than alphanumerics.
                {
                    ++n_nodes;
                    leaf=1;
                    if (isalnum(next)==0)
                        error=1;
                }
                node=1;
                break;
                
            case ':':
                ++n_branches;
                if (!(previous == ')' || previous == ']' || isalnum(previous)!=0))
                    error=1; /// ":" following something different than alphanumerics
                if (isdigit(next)==0)
                    error=1; /// ":" followed by something different to numbers
                break;
                
            case '[':
                if(!(previous == ')' || isalnum(previous)!=0))
                    error=1;
                if (next!='&')
                    error=1;
                comment=1;
                break;
            case '&':
                break;
            default:
                break;
        }
    }
    
    if (current!=';') ///Lake of final ";"
    {
		fprintf (stderr, "\n\tTree is not finished with ;");
        r_value=SETTINGS_ERROR;
    }
    if (n_branches!=n_nodes-1) ///Discordant number of branches and nodes
    {
        fprintf(stderr, "\n\tThere are %s of nodes related with the number of branches. Maybe there are a missed comma, or the root node have branch length (it must not have it)\n",n_branches>n_nodes-1?"a defect":"an excess");
        r_value=SETTINGS_ERROR;
    }
	if (n_left != n_right) ///Diferent number of "(" and ")" </dd></dl></dd></dl>
    {
		fprintf (stderr, "\n\tTree seems unbalanced (%d left and %d right parentheses)", n_left, n_right);
        r_value=SETTINGS_ERROR;
    }
    
    return r_value;
}

long int CheckNewickSTree (char * tree)
{
	int i=0;
    int k=0,error=0,n_left=0,n_right=0,n_branches=0,n_nodes=0,node=0,f_leaf=0,leaf=0;
    long int r_value=NO_ERROR;
    char current=0,previous='(',next=0;
    
	// **
    /// <dl><dt> Function structure </dt><dd>
    
    // **
    /// Initial tests
    if ((*tree)!='(')
    {
        error=1;
    }
    
    // **
    /// <dl><dt>Parsing loop, taking into account 3 characters (actual, pervious, next). It takes into account format problems, as well as structure problems:</dt><dd>
    
    
	for (i=0;i<strlen(tree);++i)
    {
        if (error == 1 || f_leaf==1)
        {
            // *
            /// Clade with only one leaf
            if (f_leaf==1)
            {
                fprintf(stderr,"\nERROR: This internal node does not have more than one leaf, Is there a missed comma?\n");
                f_leaf=0;
            }
            else
            {
                fprintf (stderr,"\nERROR: There is something wrong in the tree\n");
                error=0;
                
            }
            
			
			for (k=0; k<=i-1; k++)
				fprintf (stderr,"%c", tree[k]);
            
			fprintf (stderr," <- HERE");
            fflush(stderr);
            
            r_value=SETTINGS_ERROR;
            
        }
        
        current=*(tree+i);
        previous=*(tree+(i-1>0?i-1:0));
        next=*(tree+i+1);
        
        switch (current)
        {
            case '(':
                if (!(previous == '(' || previous == ',')) /// "(" following something different to a "(" or ","
                    error=1;
                if (next!='(')
                {
                    ++n_nodes;
                    leaf=1;
                    if (isalnum(next)==0) /// "(" followed by alphanumeric.
                        error=1;
                }
                ++n_left;
                node=0;
                break;
                
            case ')':
                if (!(previous == ')' || isalnum(previous)!=0)) /// ")" following something different to ")" or alphanumerics.
                    error=1;
                if (!(next == ')' || next== ',' || next== ':' || next=='#' || next=='/'|| next=='*' || next=='~')) /// ")" followed by something different to codes( ":","#","/", ")","*" or "," )
                    error=1;
                if (next==';' && i==strlen(tree)-2) /// Allowing the final ";"
                    error=0;
                ++n_right;
                ++n_nodes;
                if (node==0 && leaf==1)
                {
                    f_leaf=1;
                }
                leaf=0;
                node=0;
                break;
                
            case ',':
                if (!(previous == ')' || isalnum(previous)!=0)) /// "," following something different to ")" or alphanumerics.
                    error=1;
                if (next!='(') /// "," followed by something different than alphanumerics.
                {
                    ++n_nodes;
                    leaf=1;
                    if (isalnum(next)==0)
                        error=1;
                }
                node=1;
                break;
                
            case ':':
                ++n_branches;
            case '*':
            case '#':
                if (!(previous == ')' || isalnum(previous)!=0))
                    error=1; /// "#", "*" or ":" following something different than alphanumerics
                if (isdigit(next) ==0)
                    error=1; /// "#", "*" or ":" followed by something different to numbers
                break;
            case '/':
                if (!(previous == ')' || isalnum(previous)!=0))
                    error=1; /// "/" following something different than alphanumerics
                if (isdigit(next) ==0)
                    error=1; /// "/" followed by something different to numbers
                if (leaf==0)
                    error=1;
                break;
                
            default:
                break;
        }
    }
    
    if (current!=';') ///Lake of final ";"
    {
		fprintf (stderr, "\n\tTree is not finished with ;");
        r_value=SETTINGS_ERROR;
    }
    if (n_branches!=n_nodes-1) ///Discordant number of branches and nodes
    {
        fprintf(stderr, "\n\tThere are %s of nodes related with the number of branches. Maybe there are a missed comma, or the root node have branch length (it must not have it)\n",n_branches>n_nodes-1?"a defect":"an excess");
        r_value=SETTINGS_ERROR;
    }
	if (n_left != n_right) ///Diferent number of "(" and ")" </dd></dl></dd></dl>
    {
		fprintf (stderr, "\n\tTree seems unbalanced (%d left and %d right parentheses)", n_left, n_right);
        r_value=SETTINGS_ERROR;
    }
    
    return r_value;
}

long int CheckNewickLTree (char * tree)
{
	int i=0;
    int k=0,error=0,n_left=0,n_right=0,n_branches=0,n_nodes=0,node=0,f_leaf=0,leaf=0;
    long int r_value=NO_ERROR;
    char current=0,previous='(',next=0;
    
	// **
    /// <dl><dt> Function structure </dt><dd>
    
    // **
    /// Initial tests
    if ((*tree)!='(')
    {
        error=1;
    }
    
    // **
    /// <dl><dt>Parsing loop, taking into account 3 characters (actual, pervious, next). It takes into account format problems, as well as structure problems:</dt><dd>
    
    
	for (i=0;i<strlen(tree);++i)
    {
        if (error == 1 || f_leaf==1)
        {
            // *
            /// Clade with only one leaf
            if (f_leaf==1)
            {
                fprintf(stderr,"\nERROR: This internal node does not have more than one leaf, Is there a missed comma?\n");
                f_leaf=0;
            }
            else
            {
                fprintf (stderr,"\nERROR: There is something wrong in the tree\n");
                error=0;
                
            }
            
			
			for (k=0; k<=i-1; k++)
				fprintf (stderr,"%c", tree[k]);
            
			fprintf (stderr," <- HERE");
            fflush(stderr);
            
            r_value=SETTINGS_ERROR;
            
        }
        
        current=*(tree+i);
        previous=*(tree+(i-1>0?i-1:0));
        next=*(tree+i+1);
        
        switch (current)
        {
            case '(':
                if (!(previous == '(' || previous == ',')) /// "(" following something different to a "(" or ","
                    error=1;
                if (next!='(')
                {
                    ++n_nodes;
                    leaf=1;
                    if (isalnum(next)==0) /// "(" followed by alphanumeric.
                        error=1;
                }
                ++n_left;
                node=0;
                break;
                
            case ')':
                if (!(previous == ')' || isalnum(previous)!=0)) /// ")" following something different to ")" or alphanumerics.
                    error=1;
                if (!(next == ')' || next== ',' || next== ':' || next=='#' || next=='/' || next=='%' || next=='*' || next=='~')) /// ")" followed by something different to codes( ":","#","/", ")", "*", "%" or ",")
                    error=1;
                if (next==';' && i==strlen(tree)-2) /// Allowing the final ";"
                    error=0;
                ++n_right;
                ++n_nodes;
                if (node==0 && leaf==1)
                {
                    f_leaf=1;
                }
                leaf=0;
                node=0;
                break;
                
            case ',':
                if (!(previous == ')' || isalnum(previous)!=0)) /// "," following something different to ")" or alphanumerics.
                    error=1;
                if (next!='(') /// "," followed by something different than alphanumerics.
                {
                    ++n_nodes;
                    leaf=1;
                    if (isalnum(next)==0)
                        error=1;
                }
                node=1;
                break;
                
            case ':':
                ++n_branches;
            case '#':
            case '*':
            case '%':
                if (!(previous == ')' || isalnum(previous)!=0))
                    error=1; /// "#", ":", "*" or "%" following something different than alphanumerics
                if (isdigit(next) ==0)
                    error=1; /// "#", ":", "*" or "%" followed by something different to numbers
                break;
            case '/':
                if (!(previous == ')' || isalnum(previous)!=0))
                    error=1; /// "/" following something different than alphanumerics
                if (isdigit(next) ==0)
                    error=1; /// "/" followed by something different to numbers
                if (leaf==0)
                    error=1;
                break;
                
            default:
                break;
        }
    }
    
    if (current!=';') ///Lake of final ";"
    {
		fprintf (stderr, "\n\tTree is not finished with ;");
        r_value=SETTINGS_ERROR;
    }
    if (n_branches!=n_nodes-1) ///Discordant number of branches and nodes
    {
        fprintf(stderr, "\n\tThere are %s of nodes related with the number of branches. Maybe there are a missed comma, or the root node have branch length (it must not have it)\n",n_branches>n_nodes-1?"a defect":"an excess");
        r_value=SETTINGS_ERROR;
    }
	if (n_left != n_right) ///Diferent number of "(" and ")" </dd></dl></dd></dl>
    {
		fprintf (stderr, "\n\tTree seems unbalanced (%d left and %d right parentheses)", n_left, n_right);
        r_value=SETTINGS_ERROR;
    }
    
    return r_value;
}

void PrintXCharError(char *string, int x, char * errormsg1, char *errormsg)
{
    int i;
    printf("%s", errormsg1);
    for (i=1; i<x; ++i)
    {
        printf("%c", *(string+i));
    }
    printf(" %s \n", errormsg);
    fflush(stdout);
    
}

extern inline void ResetBuffer(char *buffer,size_t size)
{
    *buffer='\0';
    *(buffer+size-2)=TEST_CHAR;
}

extern inline void reallocBuffer(char **buffer,size_t *size, size_t newsize)
{
    *size=newsize;
    *buffer=(char *)realloc(*buffer, newsize*sizeof(char));
    ErrorReporter((long int) *buffer,NULL);
}

///@}

// **** Declarations of private functions **** //

// *** Recursive functions *** //

// ** Node deletion ** //

void FreeSNodes(s_node *node)
{
    int i=0;
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Post-order recursion (children,anc)
    for (i=0;i<node->n_child;++i)
    {
        FreeSNodes(*(node->children+i));
    }
    
    // *
    /// Frees memory (s_node::children, s_node)</dd></dl>
    if (node->children != NULL)
    {
        free(node->children);
        node->children=NULL;
    }
    free(node);
    node=NULL;
}

void FreeLNodes(l_node *node, int free_gnodes)
{
    int i=0;
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Post-order recursion (children,anc)
    for (i=0;i<node->n_child;++i)
    {
        FreeLNodes(*(node->children+i),free_gnodes);
    }
    
    // *
    /// Frees memory (l_node::g_nodes, l_node::children, l_node)</dd></dl>
    if (node->g_nodes != NULL && free_gnodes!=0)
    {
        free(node->g_nodes);
        node->g_nodes=NULL;
    }
    if (node->children != NULL)
    {
        free(node->children);
        node->children=NULL;
    }
    if (node->i_probs != NULL)
    {
        free(node->i_probs);
        node->i_probs=NULL;
    }
    if (node->i_combprobs != NULL)
    {
        free(node->i_combprobs);
        node->i_combprobs=NULL;
    }
    if (node->o_probs != NULL)
    {
        free(node->o_probs);
        node->o_probs=NULL;
    }
    free(node);
    node=NULL;
}

void FreeGNodes(g_node *node, int free_root)
{
    int i=0;
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Post-order recursion (children,anc)
    for (i=0;i<node->n_child;++i)
    {
        FreeGNodes(*(node->children+i),free_root);
    }
    
    // *
    /// Frees memory (s_node::children, s_node)</dd></dl>
    if (node->anc_node!=NULL || free_root==1)
    {
        if (node->children != NULL)
        {
            free(node->children);
            node->children=NULL;
        }
        free(node);
        node=NULL;
    }
    
}

// ** Node reorganitation ** //

void PostReorderSNodes(s_node * node, int * index)
{
    int i=0;
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Post-order recursion (children,anc)
    for (i=0;i<node->n_child;++i)
    {
        PostReorderSNodes(*(node->children+i),index);
    }
    // *
    /// Reorders the index, increasing the pointed variable in each call, and using it as s_node::index.</dd></dl>
    node->index= *index;
    ++*index;
}

void PostReorderLNodes(l_node * node, int * index)
{
    int i=0;
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Post-order recursion (children,anc)
    for (i=0;i<node->n_child;++i)
    {
        PostReorderLNodes(*(node->children+i),index);
    }
    // *
    /// Reorders the index, increasing the pointed variable in each call, and using it as l_node::index.</dd></dl>
    node->index= *index;
    ++*index;
}

void PostReorderGNodes(g_node * node, int * index)
{
    int i=0;
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Post-order recursion (children,anc)
    for (i=0;i<node->n_child;++i)
    {
        PostReorderGNodes(*(node->children+i),index);
    }
    // *
    /// Reorders the index, increasing the pointed variable in each call, and using it as g_node::index.</dd></dl>
    node->index= *index;
    ++*index;
}

void PreReorderSNodes (s_node * node, int * index)
{
    int i=0;
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Reorders the index, increasing the pointed variable in each call, and using it as s_node::index.
    node->index= *index;
    ++*index;
    
    // *
    /// Pre-order recursion (anc,children)</dd></dl>
    for (i=0;i<node->n_child;++i)
    {
        PreReorderSNodes(*(node->children+i),index);
    }
    
}


void PreReorderLNodes (l_node * node, int * index)
{
    int i=0;
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Reorders the index, increasing the pointed variable in each call, and using it as l_node::index.
    node->index= *index;
    ++*index;
    
    // *
    /// Pre-order recursion (anc,children)</dd></dl>
    for (i=0;i<node->n_child;++i)
    {
        PreReorderLNodes(*(node->children+i),index);
    }
    
}

void PreReorderGNodes (g_node * node, int * index)
{
    int i=0;
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Reorders the index, increasing the pointed variable in each call, and using it as g_node::index.
    node->index= *index;
    ++*index;
    
    // *
    /// Pre-order recursion (anc,children)</dd></dl>
    for (i=0;i<node->n_child;++i)
    {
        PreReorderGNodes(*(node->children+i),index);
    }
    
}

/** Node copy **/

void PostCollapseSNodes(s_node *output,s_node *input)
{
    s_node * w_output=NULL;
    s_node ** c_backup=NULL;
    int i=0;
    
    // ***
    /// <dl><dt> Function structure </dt><dd>
    
    // **
    /// Post-order recursion (children,anc)
    for (i=0;i<input->n_child;++i)
    {
        PostCollapseSNodes(output,*(input->children+i));
    }
    
    // **
    /// Saves the original s_tree::children of the output node
    w_output=(output+input->index);
    
    ErrorReporter((long int)w_output->children,NULL);
    c_backup=w_output->children;
    
    // **
    /// Copy of values
    memmove(w_output,input,sizeof(s_node));
    
    // **
    ///<dl><dt>Copy of saved pointers</dt><dd>
    
    // *
    /// children loop
    w_output->children=c_backup;
    for (i=0;i<input->n_child;++i)
    {
        *(w_output->children+i)=(output+(*(input->children+i))->index);
    }
    
    // **
    /// Copy of ancestor pointers </dd></dl></dd></dl>
    if (input->anc_node != NULL)
    {
        w_output->anc_node=(output+input->anc_node->index);
    }
    else
    {
        w_output->anc_node=NULL;
    }
    
    
}

static void CopyStoLNodes(l_node *output,s_node *input, int *l_id)
{
    int i=0;
    l_node *w_lnode=NULL;
    
    for (i=0; i<input->n_child;++i)
    {
        CopyStoLNodes(output, *(input->children+i),l_id);
    }
    w_lnode=output+*l_id;
    w_lnode->sp_index=input->sp_index;
    w_lnode->n_child=input->n_child;
    w_lnode->n_nodes=input->n_replicas;
    w_lnode->Ne=input->Ne;
    w_lnode->n_gen=input->n_gen;
    w_lnode->time=input->time;
    w_lnode->gen_length=input->gen_length;
    w_lnode->conts=input;
    w_lnode->mu_mult=input->mu_mult;
    w_lnode->gtime_mult=input->gtime_mult;
    input->n_lnodes=1;
    input->l_nodes=w_lnode;
    
    for (i=0;i<w_lnode->n_child; ++i)
    {
        *(w_lnode->children+i)=(output+((*(input->children+i))->index));
    }
    
    if (input->anc_node != NULL)
        w_lnode->anc_node=(output+(input->anc_node->index));
    else
        w_lnode->anc_node=NULL;
    
    ++*l_id;
}

void PostCollapseLNodes(l_node *output,l_node *input,int n_gleaves, int retain_lateral, int probs)
{
    l_node * w_output=NULL;
    g_node ** g_backup=NULL;
    l_node ** r_backup=NULL;
    
    int i=0;
    
    // ***
    /// <dl><dt> Function structure </dt><dd>
    
    // **
    /// Post-order recursion (children,anc)
    for (i=0;i<input->n_child;++i)
    {
        PostCollapseLNodes(output,*(input->children+i),n_gleaves,retain_lateral,probs);
    }
    
    // **
    /// Saves the original l_tree::g_nodes, l_tree::children, l_tree::i_probs and l_tree::o_probs of the output node, and allocates memory if necessary
    w_output=(output+input->index);
    
    ErrorReporter((long int)w_output->children,NULL);
    r_backup=w_output->children;
    if (w_output->g_nodes!=NULL)
        g_backup=w_output->g_nodes;
    
    // **
    /// Copy of values
    memmove(w_output,input,sizeof(l_node));
    if (retain_lateral!=1)
        w_output->lat_node=NULL;
    
    // **
    ///<dl><dt>Copy of saved pointers</dt><dd>
    
    //*
    /// N_lineages probabilities (if necessary)
    switch (probs)
    {
        case 1:
            w_output->i_probs=input->i_probs;
            input->i_probs=NULL;
            w_output->i_combprobs=input->i_combprobs;
            input->i_combprobs=NULL;
            w_output->o_probs=input->o_probs;
            input->o_probs=NULL;
            break;
        default:
            w_output->i_probs=NULL;
            w_output->i_combprobs=NULL;
            w_output->o_probs=NULL;
            break;
    }
    
    // *
    /// children loop
    w_output->children=r_backup;
    for (i=0;i<input->n_child;++i)
    {
        *(w_output->children+i)=(output+(*(input->children+i))->index);
    }
    
    // *
    /// g_nodes
    if (g_backup != NULL)
    {
        w_output->g_nodes=g_backup;
        
        if (input->g_nodes!=NULL)
        {
            memmove(w_output->g_nodes,input->g_nodes,sizeof(g_node*)*n_gleaves);
        }
        
    }
    
    // **
    /// Copy of ancestor pointers </dd></dl>
    if (input->anc_node != NULL)
    {
        w_output->anc_node=(output+input->anc_node->index);
    }
    else
    {
        w_output->anc_node=NULL;
    }
    
    
}

void PostCollapseGNodes(g_node *output,g_node *input)
{
    g_node * w_output=NULL;
    g_node ** c_backup=NULL;
    int i=0;
    
    // ***
    /// <dl><dt> Function structure </dt><dd>
    
    // **
    /// Post-order recursion (children,anc)
    for (i=0;i<input->n_child;++i)
    {
        PostCollapseGNodes(output,*(input->children+i));
    }
    
    // **
    /// Saves the original g_tree::children of the output node
    w_output=(output+input->index);
    
    ErrorReporter((long int)w_output->children,NULL);
    c_backup=w_output->children;
    
    // **
    /// Copy of values
    memmove(w_output,input,sizeof(g_node));
    
    // **
    ///<dl><dt>Copy of saved pointers</dt><dd>
    
    // *
    /// children loop
    w_output->children=c_backup;
    for (i=0;i<input->n_child;++i)
    {
        *(w_output->children+i)=(output+(*(input->children+i))->index);
    }
    
    // **
    /// Copy of ancestor pointers </dd></dl></dd></dl>
    if (input->anc_node != NULL)
    {
        w_output->anc_node=(output+input->anc_node->index);
    }
    else
    {
        w_output->anc_node=NULL;
    }
    
    
}

void PreCollapseSNodes(s_node *output,s_node *input)
{
    s_node * w_output=NULL;
    s_node ** c_backup=NULL;
    int i=0;
    
    // ***
    /// <dl><dt> Function structure </dt><dd>
    
    // **
    /// Saves the original s_tree::children of the output node
    w_output=(output+input->index);
    
    ErrorReporter((long int)w_output->children,NULL);
    c_backup=w_output->children;
    
    // **
    /// Copy of values
    memmove(w_output,input,sizeof(s_node));
    
    // **
    ///<dl><dt>Copy of saved pointers</dt><dd>
    
    // *
    /// children loop
    w_output->children=c_backup;
    for (i=0;i<input->n_child;++i)
    {
        *(w_output->children+i)=(output+(*(input->children+i))->index);
    }
    
    // **
    /// Copy of ancestor pointers
    if (input->anc_node != NULL)
    {
        w_output->anc_node=(output+input->anc_node->index);
    }
    else
    {
        w_output->anc_node=NULL;
    }
    
    // **
    /// Pre-order recursion (anc,children) </dd></dl></dd></dl>
    for (i=0;i<input->n_child;++i)
    {
        PreCollapseSNodes(output,*(input->children+i));
    }
    
    
}

void PreCollapseLNodes(l_node *output,l_node *input,int n_gleaves, int retain_lateral, int probs)
{
    l_node * w_output=NULL;
    g_node ** g_backup=NULL;
    l_node ** r_backup=NULL;
    
    int i=0;
    
    // ***
    /// <dl><dt> Function structure </dt><dd>
    
    // **
    /// Saves the original l_tree::g_nodes and l_tree::children of the output node
    w_output=(output+input->index);
    
    ErrorReporter((long int)w_output->children,NULL);
    r_backup=w_output->children;
    if (w_output->g_nodes!=NULL)
        g_backup=w_output->g_nodes;
    
    // **
    /// Copy of values
    memmove(w_output,input,sizeof(l_node));
    if (retain_lateral!=1)
    {
        w_output->lat_node=NULL;
    }
    
    // **
    ///<dl><dt>Copy of saved pointers</dt><dd>
    
    //*
    /// N_lineages probabilities (if necessary)
    switch (probs)
    {
        case 1:
            w_output->i_probs=input->i_probs;
            input->i_probs=NULL;
            w_output->i_combprobs=input->i_combprobs;
            input->i_combprobs=NULL;
            w_output->o_probs=input->o_probs;
            input->o_probs=NULL;
            break;
        default:
            w_output->i_probs=NULL;
            w_output->i_combprobs=NULL;
            w_output->o_probs=NULL;
            break;
    }
    // *
    /// children loop
    w_output->children=r_backup;
    for (i=0;i<input->n_child;++i)
    {
        *(w_output->children+i)=(output+(*(input->children+i))->index);
    }
    
    // *
    /// g_nodes </dd></dl>
    if (g_backup != NULL)
    {
        w_output->g_nodes=g_backup;
        if (input->g_nodes!=NULL)
        {
            memmove(w_output->g_nodes,input->g_nodes,sizeof(g_node*)*n_gleaves);
        }
        
    }
    
    // **
    /// Copy of ancestor pointers
    if (input->anc_node != NULL)
    {
        w_output->anc_node=(output+input->anc_node->index);
    }
    else
    {
        w_output->anc_node=NULL;
    }
    
    // **
    /// Pre-order recursion (anc,children)  </dd></dl>
    for (i=0;i<input->n_child;++i)
    {
        PreCollapseLNodes(output,*(input->children+i),n_gleaves,retain_lateral,probs);
    }
    
    
}

void PreCollapseGNodes(g_node *output,g_node *input)
{
    g_node * w_output=NULL;
    g_node ** c_backup=NULL;
    int i=0;
    
    // ***
    /// <dl><dt> Function structure </dt><dd>
    
    // **
    /// Saves the original g_tree::children of the output node
    w_output=(output+input->index);
    
    ErrorReporter((long int)w_output->children,NULL);
    c_backup=w_output->children;
    
    // **
    /// Copy of values
    memmove(w_output,input,sizeof(s_node));
    
    // **
    ///<dl><dt>Copy of saved pointers</dt><dd>
    
    // *
    /// children loop
    w_output->children=c_backup;
    for (i=0;i<input->n_child;++i)
    {
        *(w_output->children+i)=(output+(*(input->children+i))->index);
    }
    
    // **
    /// Copy of ancestor pointers
    if (input->anc_node != NULL)
    {
        w_output->anc_node=(output+input->anc_node->index);
    }
    else
    {
        w_output->anc_node=NULL;
    }
    
    // **
    /// Pre-order recursion (anc,children) </dd></dl></dd></dl>
    for (i=0;i<input->n_child;++i)
    {
        PreCollapseGNodes(output,*(input->children+i));
    }
    
}

// ** Node features obtention ** //
static int Count_duplications_lnodes(l_node *node)
{
    int i=0;
    int n_dup=0;
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Post-order recursion (children,anc)
    for (i=0;i<node->n_child;++i)
    {
        n_dup+=Count_duplications_lnodes(*(node->children+i));
    }
    
    // *
    /// Counting </dd></dl>
    if (node->kind_node==DUP)//Duplication
    {
        ++n_dup;
    }
    return n_dup;

}

static int Count_transfers_lnodes(l_node *node)
{
    int i=0;
    int n_trans=0;
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Post-order recursion (children,anc)
    for (i=0;i<node->n_child;++i)
    {
        n_trans+=Count_transfers_lnodes(*(node->children+i));
    }
    
    // *
    /// Counting </dd></dl>
    if (node->kind_node==TRFR)//Duplication
    {
        ++n_trans;
    }
    return n_trans;
    
}

static int Count_gc_lnodes(l_node *node)
{
    int i=0;
    int n_gc=0;
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Post-order recursion (children,anc)
    for (i=0;i<node->n_child;++i)
    {
        n_gc+=Count_gc_lnodes(*(node->children+i));
    }
    
    // *
    /// Counting </dd></dl>
    if (node->kind_node==GC)//Duplication
    {
        ++n_gc;
    }
    return n_gc;
    
}

static int Count_losses_lnodes(l_node *node)
{
    int i=0;
    int n_losses=0;
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Post-order recursion (children,anc)
    for (i=0;i<node->n_child;++i)
    {
        n_losses+=Count_losses_lnodes(*(node->children+i));
    }
    
    // *
    /// Counting </dd></dl>
    if (node->kind_node==LOSS)//Loss
    {
        ++n_losses;
    }
    
    return n_losses;
    
}


static long double Measure_g_node_cu_height(g_node *node, int g_Ne)
{
    int i=0;
    long double current=0, new=0;
    
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Post-order recursion selecting the longest accumulated branch length
    for (i=0;i<node->n_child;++i)
    {
        new=Measure_g_node_cu_height(*(node->children+i),g_Ne);
        
        if (new>current)
            current=new;
    }
    
    // **
    /// Returns the biggest branch lenght added to the own one.</dd></dl>
    return node->gen_length/(node->contl->Ne!=0?node->contl->Ne:g_Ne)+current;

}

static long double Measure_g_node_bl_height(g_node *node)
{
    int i=0;
    long double current=0, new=0;
    
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Post-order recursion selecting the longest accumulated branch length
    for (i=0;i<node->n_child;++i)
    {
        new=Measure_g_node_bl_height(*(node->children+i));
        
        if (new>current)
            current=new;
    }
    
    // **
    /// Returns the biggest branch lenght added to the own one.</dd></dl>
    return current+node->bl;
}

static long double Measure_g_node_cu_length(g_node *node, int g_Ne)
{
    int i=0;
    long double sum=0;
    
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Post-order recursion acumulating branch length
    for (i=0;i<node->n_child;++i)
    {
        sum+=Measure_g_node_cu_length(*(node->children+i),g_Ne);
    }
    
    // **
    /// Returns the biggest branch lenght added to the own one.</dd></dl>
    return node->gen_length/(node->contl->Ne!=0?node->contl->Ne:g_Ne)+sum;
    
}

static long double Measure_g_node_bl_length(g_node *node)
{
    int i=0;
    long double sum=0;
    
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Post-order recursion selecting the longest accumulated branch length
    for (i=0;i<node->n_child;++i)
    {
        sum+=Measure_g_node_bl_length(*(node->children+i));

    }
    
    // **
    /// Returns the biggest branch lenght added to the own one.</dd></dl>
    return sum+node->bl;
}

static long double Measure_s_node_cu_height(s_node *node, int g_Ne)
{
    int i=0;
    long double current=0, new=0;
    
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Post-order recursion selecting the longest accumulated branch length
    for (i=0;i<node->n_child;++i)
    {
        new=Measure_s_node_cu_height(*(node->children+i),g_Ne);
        
        if (new>current)
            current=new;
    }
    
    // **
    /// Returns the biggest branch lenght added to the own one.</dd></dl>
    return node->gen_length/(node->Ne!=0?node->Ne:g_Ne)+current;
    
}

static long double Measure_s_node_gl_height(s_node *node)
{
    int i=0;
    long double current=0, new=0;
    
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Post-order recursion selecting the longest accumulated branch length
    for (i=0;i<node->n_child;++i)
    {
        new=Measure_s_node_gl_height(*(node->children+i));
        
        if (new>current)
            current=new;
    }
    
    // **
    /// Returns the biggest branch lenght added to the own one.</dd></dl>
    return current+node->gen_length;
}

static long double Measure_s_node_cu_length(s_node *node, int g_Ne)
{
    int i=0;
    long double sum=0;
    
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Post-order recursion acumulating branch length
    for (i=0;i<node->n_child;++i)
    {
        sum+=Measure_s_node_cu_length(*(node->children+i),g_Ne);
    }
    
    // **
    /// Returns the biggest branch lenght added to the own one.</dd></dl>
    return node->gen_length/(node->Ne!=0?node->Ne:g_Ne)+sum;
    
}

static long double Measure_s_node_gl_length(s_node *node)
{
    int i=0;
    long double sum=0;
    
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    ///\attention This could be improved looking for the leaf with the biggest n_gen! (but we are not using loops here)
    
    // *
    /// Post-order recursion selecting the longest accumulated branch length
    for (i=0;i<node->n_child;++i)
    {
        sum+=Measure_s_node_gl_length(*(node->children+i));
        
    }
    
    // **
    /// Returns the biggest branch lenght added to the own one.</dd></dl>
    return sum+node->gen_length;
}



// ** Node modification ** //

static void RefineSNodes(s_node * node, int ind_persp, double gen_time)
{
    int i=0;
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Sets the correct s_node::n_nodes for leaf nodes or calculates and applies the s_node::n_gen for every node </dd></dl>
    
    if (node->anc_node!=NULL)
    {
        switch (node->n_child)
        {
            case 0:
                if (node->n_replicas==0)
                    node->n_replicas=ind_persp;
                break;
                
            default:
                if (node->n_replicas>0)
                    ErrorReporter(SETTINGS_ERROR,": number of replicas per taxa applied to internal branches are not allowed\n");
                break;
        }
        node->n_gen=node->anc_node->n_gen+node->gen_length;
        node->time=node->anc_node->time+node->gen_length*node->gtime_mult*gen_time;
    }
    /*else //Implicit
    {
        node->n_gen=0;
    }*/
    
    // *
    /// Pre-order recursion (anc, children)
    for (i=0;i<node->n_child;++i)
    {
        RefineSNodes(*(node->children+i),ind_persp,gen_time);
    }
    
}

static void RefineLNodes(l_node * node, int n_gleaves, int ind_persp, double gen_time, int *n_dup, int *n_loss, int *n_trans, int *n_gc)
{
    int i=0,j=0;
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Gets lnodes stats (number of duplications, losses, transferences and gc)
    switch (node->kind_node)
    {
        case DUP:
            (*n_dup)++;
            break;
        case LOSS:
            (*n_loss)++;
        case GC:
            (*n_gc)++;
            break;
        case TRFR:
            (*n_trans)++;
            break;
    }
    
    // *
    /// Sets the correct l_node::n_nodes for leaf nodes or calculates and applies l_node::n_gen for every node
    
    if (node->anc_node!=NULL) //No root
    {
        switch (node->n_child)
        {
            case 0:
                switch (node->kind_node)
                {
                    case LOSS:
                    case RTRFR:
                    case RGC:
                        node->n_nodes=0;
                        break;
                    default:
                        if (node->n_nodes==0)
                            node->n_nodes=ind_persp;
                        break;
                }
                node->n_ilin=node->n_nodes;
  
            default:
                node->n_gen=node->anc_node->n_gen+node->gen_length;
                node->time=node->anc_node->time+node->gen_length*node->gtime_mult*gen_time;
                break;
        }
    }

    // *
    /// Allocates l_node::g_nodes if it is necesary, and initializes them.</dd></dl>
    if (node->g_nodes==NULL)
    {
        node->g_nodes=calloc(n_gleaves,sizeof(g_node *));
        ErrorReporter((long int)node->g_nodes,NULL);
        for (j=0; j<n_gleaves; ++j)
        {
            *(node->g_nodes+j)=NULL;
        }
        
    }
    
    // *
    /// Pre-order recursion (anc,children)
    for (i=0;i<node->n_child;++i)
    {
        RefineLNodes(*(node->children+i),n_gleaves, ind_persp, gen_time, n_dup, n_loss, n_trans, n_gc);
    }
    
}
//\cond DOXYGEN_EXCLUDE
//static void CleanlossesLNodes(l_node * node, l_node ** root, int *n_deletions, int *n_leaves)
//{
//    int i=0,j=0,n_child=node->n_child;
//    l_node *w_lnode=NULL;
//    
//    // ***
//    /// <dl><dt> Function structure </dt><dd>
//    
//    // **
//    /// Post-order recursion (children,anc)
//    for (i=0;i<n_child;++i)
//    {
//        CleanlossesLNodes(*(node->children+i),root,n_deletions,n_leaves);
//    }
//    
//    // **
//    ///<dl><dt>No proper internal node</dt><dd>
//    if(node->n_child<2)
//    {
//        // *
//        // If real leave, count it.
//        if(node->kind_node==SP)
//        {
//            ++(*n_leaves);
//        }
//        else if (node->anc_node!=NULL)
//        {
//            w_lnode=node->anc_node;
//            // *
//            /// Deletion if the node has 0 children and is not a real leaf (SP).
//            if (node->n_child==0)
//            {
//                for (j=0;j<w_lnode->n_child;++j)
//                {
//                    if (*(w_lnode->children+j)==node)
//                        break;
//                }
//                *(w_lnode->children+j)=*(w_lnode->children+w_lnode->n_child-1);
//                --w_lnode->n_child;
//                FreeLNodes(node, 1);
//                ++(*n_deletions);
//            }
//            // *
//            /// Deletion and "by pass" if the node has 1 child
//            else if (node->n_child==1)
//            {
//                for (j=0;j<w_lnode->n_child;++j)
//                {
//                    if (*(w_lnode->children+j)==node)
//                        break;
//                }
//                *(w_lnode->children+j)=*node->children;
//                w_lnode=*node->children;
//                w_lnode->anc_node=node->anc_node;
//                w_lnode->gen_length=node->anc_node->n_gen-w_lnode->n_gen;
//                node->n_child=0; //To delete only this node
//                FreeLNodes(node, 1);
//                ++(*n_deletions);
//            }
//        }
//        // *
//        ///Change of the root of the tree if it should be deleted</dd></dl></dd></dl>
//        else
//        {
//            w_lnode=*node->children;
//            *root=w_lnode;
//            (*root)->anc_node=NULL;
//            node->n_child=0; //To delete only this node
//            FreeLNodes(node, 1);
//            ++(*n_deletions);
//        }
//        
//    }
//    
//}
//\endcond

static void AddConstantNgenSNodes(s_node *node, double constant, double gen_time)
{
    unsigned int i=0;
    
    node->n_gen+=constant;
    node->time=node->anc_node->time+node->gen_length*node->gtime_mult*gen_time;
    for (i=0; i<node->n_child;++i)
    {
        AddConstantNgenSNodes(*(node->children+i), constant, gen_time);
    }
}

static void TemporalizeLNodes(l_node *node, double gen_time)
{
    unsigned int i=0;
    
    node->time=node->anc_node->time+node->gen_length*node->gtime_mult*gen_time;
    for (i=0; i<node->n_child; ++i)
    {
        TemporalizeLNodes(*(node->children+i), gen_time);
    }
}

// *** Tree I/O *** //

void WriteSNodesGen (s_node * p, name_c * names)
{
    int i=0;
    
    // ***
    /// <dl><dt> Function structure </dt><dd>
    
    if (p != NULL)
    {
		if (p->n_child==0)
        {
            // **
            /// <dl><dt>If the node is a leaf:</dt><dd>
            // *
            /// Prints the node name and its branch length.</dd></dl>
            if (names==NULL)
                printf("%d:%.8lf",p->sp_index,p->gen_length);
            else
                printf("%s:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->gen_length);
        }
		else
        {
            // **
            /// <dl><dt>Else (internal node):</dt><dd>
            // *
            /// Prints "(" and starts the post-order recursion (child loop)
			printf("(");
            for (i=0;i<p->n_child-1;++i)
            {
                // *
                /// Calls itself using each different child. After each call (except last child) prints ","
                WriteSNodesGen(*(p->children+i),names);
                printf(",");
            }
			WriteSNodesGen(*(p->children+i), names);
            
            // **
            /// Prints the information of this internal node (closing it with ")") or finishes the tree if the node is the root</dd></dl></dd></dl>
			if (p->anc_node !=NULL)
				printf("):%.8lf",p->gen_length);
            else
                printf(");");
        }
    }
	
}

void WriteSNodesGenIntlabel (s_node * p, name_c * names)
{
    int i=0;
    
    // ***
    /// <dl><dt> Function structure </dt><dd>
    
    if (p != NULL)
    {
        if (p->n_child==0)
        {
            // **
            /// <dl><dt>If the node is a leaf:</dt><dd>
            // *
            /// Prints the node name and its branch length.</dd></dl>
            if (names==NULL)
                printf("%d:%.8lf",p->sp_index,p->gen_length);
            else
                printf("%s:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->gen_length);
        }
        else
        {
            // **
            /// <dl><dt>Else (internal node):</dt><dd>
            // *
            /// Prints "(" and starts the post-order recursion (child loop)
            printf("(");
            for (i=0;i<p->n_child-1;++i)
            {
                // *
                /// Calls itself using each different child. After each call (except last child) prints ","
                WriteSNodesGenIntlabel(*(p->children+i),names);
                printf(",");
            }
            WriteSNodesGenIntlabel(*(p->children+i), names);
            
            // **
            /// Prints the information of this internal node (closing it with ")") or finishes the tree if the node is the root</dd></dl></dd></dl>
            if (p->anc_node !=NULL)
                printf(")%d:%.8lf",p->index,p->gen_length);
            else
                printf(");");
        }
    }
    
}

void WriteSNodesTime (s_node * p, name_c * names, double gen_time)
{
    int i=0;
    
    // ***
    /// <dl><dt> Function structure </dt><dd>

    if (p != NULL)
    {
		if (p->n_child==0)
        {
            // **
            /// <dl><dt>If the node is a leaf:</dt><dd>
            // *
            /// Prints the node name and its branch length.</dd></dl>
            if (names==NULL)
                printf("%d:%.8lf",p->sp_index,p->gen_length*p->gtime_mult*gen_time);
            else
                printf("%s:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->gen_length*p->gtime_mult*gen_time);
        }
		else
        {
            // **
            /// <dl><dt>Else (internal node):</dt><dd>
            // *
            /// Prints "(" and starts the post-order recursion (child loop)
			printf("(");
            for (i=0;i<p->n_child-1;++i)
            {
                // *
                /// Calls itself using each different child. After each call (except last child) prints ","
                WriteSNodesTime(*(p->children+i),names,gen_time);
                printf(",");
            }
			WriteSNodesTime(*(p->children+i), names,gen_time);
            
            // **
            /// Prints the information of this internal node (closing it with ")") or finishes the tree if the node is the root</dd></dl></dd></dl>
			if (p->anc_node !=NULL)
				printf("):%.8lf",p->gen_length*p->gtime_mult*gen_time);
            else
                printf(");");
        }
    }
	
}

void WriteSNodesTimeIntlabel (s_node * p, name_c * names, double gen_time)
{
    int i=0;
    
    // ***
    /// <dl><dt> Function structure </dt><dd>
    
    if (p != NULL)
    {
        if (p->n_child==0)
        {
            // **
            /// <dl><dt>If the node is a leaf:</dt><dd>
            // *
            /// Prints the node name and its branch length.</dd></dl>
            if (names==NULL)
                printf("%d:%.8lf",p->sp_index,p->gen_length*p->gtime_mult*gen_time);
            else
                printf("%s:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->gen_length*p->gtime_mult*gen_time);
        }
        else
        {
            // **
            /// <dl><dt>Else (internal node):</dt><dd>
            // *
            /// Prints "(" and starts the post-order recursion (child loop)
            printf("(");
            for (i=0;i<p->n_child-1;++i)
            {
                // *
                /// Calls itself using each different child. After each call (except last child) prints ","
                WriteSNodesTimeIntlabel(*(p->children+i),names,gen_time);
                printf(",");
            }
            WriteSNodesTimeIntlabel(*(p->children+i), names,gen_time);
            
            // **
            /// Prints the information of this internal node (closing it with ")") or finishes the tree if the node is the root</dd></dl></dd></dl>
            if (p->anc_node !=NULL)
                printf(")%d:%.8lf",p->index,p->gen_length*p->gtime_mult*gen_time);
            else
                printf(");");
        }
    }
    
}

//void WriteSNodesNexus (s_node * p, name_c * names, double gen_time)
//{
//    int i=0;
//    
//    // ***
//    /// <dl><dt> Function structure </dt><dd>
//    
//    if (p != NULL)
//    {
//		if (p->n_child==0)
//        {
//            // **
//            /// <dl><dt>If the node is a leaf:</dt><dd>
//            // *
//            /// Prints the node name and its branch length.</dd></dl>
//            if (names==NULL)
//                printf("%d[&t_length=%.8lf]:%.8lf",p->sp_index,p->gen_length*p->gtime_mult*gen_time,p->gen_length);
//            else
//                printf("%s[&t_length=%.8lf]:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->gen_length*p->gtime_mult*gen_time,p->gen_length);
//        }
//		else
//        {
//            // **
//            /// <dl><dt>Else (internal node):</dt><dd>
//            // *
//            /// Prints "(" and starts the post-order recursion (child loop)
//			printf("(");
//            for (i=0;i<p->n_child-1;++i)
//            {
//                // *
//                /// Calls itself using each different child. After each call (except last child) prints ","
//                WriteSNodesNexus(*(p->children+i),names,gen_time);
//                printf(",");
//            }
//			WriteSNodesNexus(*(p->children+i), names,gen_time);
//            
//            // **
//            /// Prints the information of this internal node (closing it with ")") or finishes the tree if the node is the root</dd></dl></dd></dl>
//			if (p->anc_node !=NULL)
//				printf(")[&t_length=%.8lf]:%.8lf",p->gen_length*p->gtime_mult*gen_time,p->gen_length);
//            else
//                printf(");");
//        }
//    }
//	
//}

void WriteSNodesFileGen(FILE * file,s_node * p, name_c * names)
{
    int i=0;
    
    // ***
    /// <dl><dt> Function structure </dt><dd>
    
    if (p != NULL)
    {
		if (p->n_child==0)
        {
            // **
            /// <dl><dt>If the node is a leaf:</dt><dd>
            // *
            /// Prints the node name and its branch length.</dd></dl>
            if (names==NULL)
                fprintf(file,"%d:%.8lf",p->sp_index,p->gen_length);
            else
                fprintf(file,"%s:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->gen_length);
        }
		else
        {
            // **
            /// <dl><dt>Else (internal node):</dt><dd>
            // *
            /// Prints "(" and starts the post-order recursion (child loop)
			fprintf(file,"(");
            for (i=0;i<p->n_child-1;++i)
            {
                // *
                /// Calls itself using each different child. After each call (except last child) prints ","
                WriteSNodesFileGen(file,*(p->children+i),names);
                fprintf(file,",");
            }
			WriteSNodesFileGen(file,*(p->children+i), names);
            
            // **
            /// Prints the information of this internal node (closing it with ")") or finishes the tree if the node is the root</dd></dl></dd></dl>
			if (p->anc_node !=NULL)
				fprintf(file,"):%.8lf",p->gen_length);
            else
                fprintf(file,");\n");
        }
    }
	
}

void WriteSNodesFileGenIntlabel(FILE * file,s_node * p, name_c * names)
{
    int i=0;
    
    // ***
    /// <dl><dt> Function structure </dt><dd>
    
    if (p != NULL)
    {
        if (p->n_child==0)
        {
            // **
            /// <dl><dt>If the node is a leaf:</dt><dd>
            // *
            /// Prints the node name and its branch length.</dd></dl>
            if (names==NULL)
                fprintf(file,"%d:%.8lf",p->sp_index,p->gen_length);
            else
                fprintf(file,"%s:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->gen_length);
        }
        else
        {
            // **
            /// <dl><dt>Else (internal node):</dt><dd>
            // *
            /// Prints "(" and starts the post-order recursion (child loop)
            fprintf(file,"(");
            for (i=0;i<p->n_child-1;++i)
            {
                // *
                /// Calls itself using each different child. After each call (except last child) prints ","
                WriteSNodesFileGenIntlabel(file,*(p->children+i),names);
                fprintf(file,",");
            }
            WriteSNodesFileGenIntlabel(file,*(p->children+i), names);
            
            // **
            /// Prints the information of this internal node (closing it with ")") or finishes the tree if the node is the root</dd></dl></dd></dl>
            if (p->anc_node !=NULL)
                fprintf(file,")%d:%.8lf",p->index,p->gen_length);
            else
                fprintf(file,");\n");
        }
    }
    
}

void WriteSNodesFileTime (FILE * file,s_node * p, name_c * names, double gen_time)
{
    int i=0;
    
    // ***
    /// <dl><dt> Function structure </dt><dd>
    
    if (p != NULL)
    {
        if (p->n_child==0)
        {
            // **
            /// <dl><dt>If the node is a leaf:</dt><dd>
            // *
            /// Prints the node name and its branch length.</dd></dl>
            if (names==NULL)
                fprintf(file,"%d:%.8lf",p->sp_index,p->gen_length*p->gtime_mult*gen_time);
            else
                fprintf(file,"%s:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->gen_length*p->gtime_mult*gen_time);
        }
        else
        {
            // **
            /// <dl><dt>Else (internal node):</dt><dd>
            // *
            /// Prints "(" and starts the post-order recursion (child loop)
            fprintf(file,"(");
            for (i=0;i<p->n_child-1;++i)
            {
                // *
                /// Calls itself using each different child. After each call (except last child) prints ","
                WriteSNodesFileTime(file,*(p->children+i),names,gen_time);
                fprintf(file,",");
            }
            WriteSNodesFileTime(file,*(p->children+i), names,gen_time);
            
            // **
            /// Prints the information of this internal node (closing it with ")") or finishes the tree if the node is the root</dd></dl></dd></dl>
            if (p->anc_node !=NULL)
                fprintf(file,"):%.8lf",p->gen_length*p->gtime_mult*gen_time);
            else
                fprintf(file,");\n");
        }
    }
    
}

void WriteSNodesFileTimeIntlabel (FILE * file,s_node * p, name_c * names, double gen_time)
{
    int i=0;
    
    // ***
    /// <dl><dt> Function structure </dt><dd>
    
    if (p != NULL)
    {
		if (p->n_child==0)
        {
            // **
            /// <dl><dt>If the node is a leaf:</dt><dd>
            // *
            /// Prints the node name and its branch length.</dd></dl>
            if (names==NULL)
                fprintf(file,"%d:%.8lf",p->sp_index,p->gen_length*p->gtime_mult*gen_time);
            else
                fprintf(file,"%s:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->gen_length*p->gtime_mult*gen_time);
        }
		else
        {
            // **
            /// <dl><dt>Else (internal node):</dt><dd>
            // *
            /// Prints "(" and starts the post-order recursion (child loop)
			fprintf(file,"(");
            for (i=0;i<p->n_child-1;++i)
            {
                // *
                /// Calls itself using each different child. After each call (except last child) prints ","
                WriteSNodesFileTimeIntlabel(file,*(p->children+i),names,gen_time);
                fprintf(file,",");
            }
			WriteSNodesFileTimeIntlabel(file,*(p->children+i), names,gen_time);
            
            // **
            /// Prints the information of this internal node (closing it with ")") or finishes the tree if the node is the root</dd></dl></dd></dl>
			if (p->anc_node !=NULL)
				fprintf(file,")%d:%.8lf",p->index,p->gen_length*p->gtime_mult*gen_time);
            else
                fprintf(file,");\n");
        }
    }
	
}

void WriteLNodesGen (l_node * p, name_c * names)
{
    int i=0;
    
    // ***
    /// <dl><dt> Function structure </dt><dd>
    
    if (p != NULL)
    {
        switch (p->n_child)
        {
            case 0:
                // **
                /// <dl><dt>If the node is a leaf:</dt><dd>
                // *
                /// Prints the node name and its branch length.</dd></dl>
                switch (p->kind_node)
                {
                        
                    default:
                        if (names==NULL)
                            printf("%d_%d:%.8lf",p->sp_index,p->paralog,p->gen_length);
                        else
                            printf("%s_%d:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->gen_length);
                        break;
                        
                    case LOSS:
                        if (p->conts==NULL)
                        {
                            printf("Lost-%d_%d:%.8lf",p->index,p->paralog,p->gen_length);
                        }
                        else if(p->conts->n_child==0 && names!=NULL)
                        {
                            printf("Lost-%s_%d:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->gen_length);
                        }
                        else
                        {
                            printf("Lost-%d_%d:%.8lf",p->conts->index,p->paralog,p->gen_length);
                        }
                        break;
                    case RTRFR:
                        if (p->conts==NULL)
                        {
                            printf("Rtransf-%d_%d:%.8lf",p->index,p->paralog,p->gen_length);
                        }
                        else if(p->conts->n_child==0 && names!=NULL)
                        {
                            printf("Rtransf-%s_%d:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->gen_length);
                        }
                        else
                        {
                            printf("Rtransf-%d_%d:%.8lf",p->conts->index,p->paralog,p->gen_length);
                        }
                        break;
                    case RGC:
                        if (p->conts==NULL)
                        {
                            printf("Rgc-%d_%d:%.8lf",p->index,p->paralog,p->gen_length);
                        }
                        else if (p->conts->n_child==0 && names!=NULL)
                        {
                            printf("Rgc-%s_%d:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->gen_length);
                        }
                        else
                        {
                            printf("Rgc-%d_%d:%.8lf",p->conts->index,p->paralog,p->gen_length);
                        }
                        break;
                        
                }
                break;
                
            default:
                // **
                /// <dl><dt>Else (internal node):</dt><dd>
                // *
                /// Prints "(" and starts the post-order recursion (child loop)
                printf("(");
                for (i=0;i<p->n_child-1;++i)
                {
                    // *
                    /// Calls itself using each different child. After each call (except last child) prints ","
                    WriteLNodesGen(*(p->children+i),names);
                    printf(",");
                }
                WriteLNodesGen (*(p->children+i), names);
                
                // **
                /// Prints the information of this internal node (closing it with ")") or finishes the tree if the node is the root</dd></dl></dd></dl>
                if (p->anc_node !=NULL)
                    printf("):%.8lf",p->gen_length);
                else
                    printf(");");
                break;
        }
    }
}

void WriteLNodesGenIntlabel (l_node * p, name_c * names)
{
    int i=0;
    
    // ***
    /// <dl><dt> Function structure </dt><dd>
    
    if (p != NULL)
    {
        switch (p->n_child)
        {
            case 0:
                // **
                /// <dl><dt>If the node is a leaf:</dt><dd>
                // *
                /// Prints the node name and its branch length.</dd></dl>
                switch (p->kind_node)
            {
                    
                default:
                    if (names==NULL)
                        printf("%d_%d:%.8lf",p->sp_index,p->paralog,p->gen_length);
                    else
                        printf("%s_%d:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->gen_length);
                    break;
                    
                case LOSS:
                    if (p->conts==NULL)
                    {
                        printf("Lost-%d_%d:%.8lf",p->index,p->paralog,p->gen_length);
                    }
                    else if(p->conts->n_child==0 && names!=NULL)
                    {
                        printf("Lost-%s_%d:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->gen_length);
                    }
                    else
                    {
                        printf("Lost-%d_%d:%.8lf",p->conts->index,p->paralog,p->gen_length);
                    }
                    break;
                case RTRFR:
                    if (p->conts==NULL)
                    {
                        printf("Rtransf-%d_%d:%.8lf",p->index,p->paralog,p->gen_length);
                    }
                    else if(p->conts->n_child==0 && names!=NULL)
                    {
                        printf("Rtransf-%s_%d:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->gen_length);
                    }
                    else
                    {
                        printf("Rtransf-%d_%d:%.8lf",p->conts->index,p->paralog,p->gen_length);
                    }
                    break;
                case RGC:
                    if (p->conts==NULL)
                    {
                        printf("Rgc-%d_%d:%.8lf",p->index,p->paralog,p->gen_length);
                    }
                    else if (p->conts->n_child==0 && names!=NULL)
                    {
                        printf("Rgc-%s_%d:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->gen_length);
                    }
                    else
                    {
                        printf("Rgc-%d_%d:%.8lf",p->conts->index,p->paralog,p->gen_length);
                    }
                    break;
                    
            }
                break;
                
            default:
                // **
                /// <dl><dt>Else (internal node):</dt><dd>
                // *
                /// Prints "(" and starts the post-order recursion (child loop)
                printf("(");
                for (i=0;i<p->n_child-1;++i)
                {
                    // *
                    /// Calls itself using each different child. After each call (except last child) prints ","
                    WriteLNodesGenIntlabel(*(p->children+i),names);
                    printf(",");
                }
                WriteLNodesGenIntlabel(*(p->children+i), names);
                
                // **
                /// Prints the information of this internal node (closing it with ")") or finishes the tree if the node is the root</dd></dl></dd></dl>
                if (p->anc_node !=NULL)
                    printf(")%d:%.8lf",p->index,p->gen_length);
                else
                    printf(");");
                break;
        }
    }
}

void WriteLNodesTime (l_node * p, name_c * names, double gen_time)
{
    int i=0;
    
    // ***
    /// <dl><dt> Function structure </dt><dd>
    
    if (p != NULL)
    {
        switch (p->n_child)
        {
            case 0:
                // **
                /// <dl><dt>If the node is a leaf:</dt><dd>
                // *
                /// Prints the node name and its branch length.</dd></dl>
                switch (p->kind_node)
            {
                    
                default:
                    if (names==NULL)
                        printf("%d_%d:%.8lf",p->sp_index,p->paralog,p->gen_length*gen_time*p->gtime_mult);
                    else
                        printf("%s_%d:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->gen_length*gen_time*p->gtime_mult);
                    break;
                    
                case LOSS:
                    if (p->conts->n_child==0 && names!=NULL)
                    {
                        printf("Lost-%s_%d:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->gen_length*gen_time*p->gtime_mult);
                    }
                    else
                    {
                        printf("Lost-%d_%d:%.8lf",p->conts->index,p->paralog,p->gen_length*gen_time*p->gtime_mult);
                    }
                    break;
                case RTRFR:
                    if (p->conts->n_child==0 && names!=NULL)
                    {
                        printf("Rtransf-%s_%d:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->gen_length*gen_time*p->gtime_mult);
                    }
                    else
                    {
                        printf("Rtransf-%d_%d:%.8lf",p->conts->index,p->paralog,p->gen_length*gen_time*p->gtime_mult);
                    }
                    break;
                case RGC:
                    if (p->conts->n_child==0 && names!=NULL)
                    {
                        printf("Rgc-%s_%d:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->gen_length*gen_time*p->gtime_mult);
                    }
                    else
                    {
                        printf("Rgc-%d_%d:%.8lf",p->conts->index,p->paralog,p->gen_length*gen_time*p->gtime_mult);
                    }
                    break;
                    
            }
                break;
                
            default:
                // **
                /// <dl><dt>Else (internal node):</dt><dd>
                // *
                /// Prints "(" and starts the post-order recursion (child loop)
                printf("(");
                for (i=0;i<p->n_child-1;++i)
                {
                    // *
                    /// Calls itself using each different child. After each call (except last child) prints ","
                    WriteLNodesTime(*(p->children+i),names,gen_time);
                    printf(",");
                }
                WriteLNodesTime(*(p->children+i), names,gen_time);
                
                // **
                /// Prints the information of this internal node (closing it with ")") or finishes the tree if the node is the root</dd></dl></dd></dl>
                if (p->anc_node !=NULL)
                    printf("):%.8lf",p->gen_length*gen_time*p->gtime_mult);
                else
                    printf(");");
                break;
        }
    }
}

void WriteLNodesTimeIntlabel (l_node * p, name_c * names, double gen_time)
{
    int i=0;
    
    // ***
    /// <dl><dt> Function structure </dt><dd>
    
    if (p != NULL)
    {
        switch (p->n_child)
        {
            case 0:
                // **
                /// <dl><dt>If the node is a leaf:</dt><dd>
                // *
                /// Prints the node name and its branch length.</dd></dl>
                switch (p->kind_node)
            {
                    
                default:
                    if (names==NULL)
                        printf("%d_%d:%.8lf",p->sp_index,p->paralog,p->gen_length*gen_time*p->gtime_mult);
                    else
                        printf("%s_%d:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->gen_length*gen_time*p->gtime_mult);
                    break;
                    
                case LOSS:
                    if (p->conts->n_child==0 && names!=NULL)
                    {
                        printf("Lost-%s_%d:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->gen_length*gen_time*p->gtime_mult);
                    }
                    else
                    {
                        printf("Lost-%d_%d:%.8lf",p->conts->index,p->paralog,p->gen_length*gen_time*p->gtime_mult);
                    }
                    break;
                case RTRFR:
                    if (p->conts->n_child==0 && names!=NULL)
                    {
                        printf("Rtransf-%s_%d:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->gen_length*gen_time*p->gtime_mult);
                    }
                    else
                    {
                        printf("Rtransf-%d_%d:%.8lf",p->conts->index,p->paralog,p->gen_length*gen_time*p->gtime_mult);
                    }
                    break;
                case RGC:
                    if (p->conts->n_child==0 && names!=NULL)
                    {
                        printf("Rgc-%s_%d:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->gen_length*gen_time*p->gtime_mult);
                    }
                    else
                    {
                        printf("Rgc-%d_%d:%.8lf",p->conts->index,p->paralog,p->gen_length*gen_time*p->gtime_mult);
                    }
                    break;
                    
            }
                break;
                
            default:
                // **
                /// <dl><dt>Else (internal node):</dt><dd>
                // *
                /// Prints "(" and starts the post-order recursion (child loop)
                printf("(");
                for (i=0;i<p->n_child-1;++i)
                {
                    // *
                    /// Calls itself using each different child. After each call (except last child) prints ","
                    WriteLNodesTimeIntlabel(*(p->children+i),names,gen_time);
                    printf(",");
                }
                WriteLNodesTimeIntlabel(*(p->children+i), names,gen_time);
                
                // **
                /// Prints the information of this internal node (closing it with ")") or finishes the tree if the node is the root</dd></dl></dd></dl>
                if (p->anc_node !=NULL)
                    printf(")%d:%.8lf",p->index,p->gen_length*gen_time*p->gtime_mult);
                else
                    printf(");");
                break;
        }
    }
}

void WriteLNodesFileGen (FILE * file,l_node * p, name_c * names)
{
    int i=0;
    
    // ***
    /// <dl><dt> Function structure </dt><dd>
    
    if (p != NULL)
    {
        switch (p->n_child)
        {
            case 0:
                // **
                /// <dl><dt>If the node is a leaf:</dt><dd>
                // *
                /// Prints the node name and its branch length.</dd></dl>
                
                switch (p->kind_node)
                {
                    default:
                        if (names==NULL)
                            fprintf(file,"%d_%d:%.8lf",p->sp_index,p->paralog,p->gen_length);
                        else
                            fprintf(file,"%s_%d:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->gen_length);
                        break;
                    case LOSS:
                        if (p->conts==NULL)
                        {
                            fprintf(file,"Lost-%d_%d:%.8lf",p->index,p->paralog,p->gen_length);
                        }
                        else if (p->conts->n_child==0 && names!=NULL)
                        {
                            fprintf(file,"Lost-%s_%d:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->gen_length);
                        }
                        else
                        {
                            fprintf(file,"Lost-%d_%d:%.8lf",p->conts->index,p->paralog,p->gen_length);
                        }
                        break;
                    case RTRFR:
                        if (p->conts==NULL)
                        {
                            fprintf(file,"Rtransf-%d_%d:%.8lf",p->index,p->paralog,p->gen_length);
                        }
                        else if (p->conts->n_child==0 && names!=NULL)
                        {
                            fprintf(file,"Rtransf-%s_%d:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->gen_length);
                        }
                        else
                        {
                            fprintf(file,"Rtransf-%d_%d:%.8lf",p->conts->index,p->paralog,p->gen_length);
                        }
                        break;
                    case RGC:
                        if (p->conts==NULL)
                        {
                            fprintf(file,"Rgc-%d_%d:%.8lf",p->index,p->paralog,p->gen_length);
                        }
                        else if (p->conts->n_child==0 && names!=NULL)
                        {
                            fprintf(file,"Rgc-%s_%d:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->gen_length);
                        }
                        else
                        {
                            fprintf(file,"Rgc-%d_%d:%.8lf",p->conts->index,p->paralog,p->gen_length);
                        }
                        break;
                }
                break;
                
            default:
                // **
                /// <dl><dt>Else (internal node):</dt><dd>
                // *
                /// Prints "(" and starts the post-order recursion (child loop)
                fprintf(file,"(");
                for (i=0;i<p->n_child-1;++i)
                {
                    // *
                    /// Calls itself using each different child. After each call (except last child) prints ","
                    WriteLNodesFileGen(file,*(p->children+i),names);
                    fprintf(file,",");
                }
                WriteLNodesFileGen(file,*(p->children+i), names);
                
                // **
                /// Prints the information of this internal node (closing it with ")") or finishes the tree if the node is the root</dd></dl></dd></dl>
                if (p->anc_node !=NULL)
                    fprintf(file,"):%.8lf",p->gen_length);
                else
                    fprintf(file,");\n");
                break;
        }
    }
}

void WriteLNodesFileGenIntlabel (FILE * file,l_node * p, name_c * names)
{
    int i=0;
    
    // ***
    /// <dl><dt> Function structure </dt><dd>
    
    if (p != NULL)
    {
        switch (p->n_child)
        {
            case 0:
                // **
                /// <dl><dt>If the node is a leaf:</dt><dd>
                // *
                /// Prints the node name and its branch length.</dd></dl>
                
                switch (p->kind_node)
            {
                default:
                    if (names==NULL)
                        fprintf(file,"%d_%d:%.8lf",p->sp_index,p->paralog,p->gen_length);
                    else
                        fprintf(file,"%s_%d:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->gen_length);
                    break;
                case LOSS:
                    if (p->conts==NULL)
                    {
                        fprintf(file,"Lost-%d_%d:%.8lf",p->index,p->paralog,p->gen_length);
                    }
                    else if (p->conts->n_child==0 && names!=NULL)
                    {
                        fprintf(file,"Lost-%s_%d:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->gen_length);
                    }
                    else
                    {
                        fprintf(file,"Lost-%d_%d:%.8lf",p->conts->index,p->paralog,p->gen_length);
                    }
                    break;
                case RTRFR:
                    if (p->conts==NULL)
                    {
                        fprintf(file,"Rtransf-%d_%d:%.8lf",p->index,p->paralog,p->gen_length);
                    }
                    else if (p->conts->n_child==0 && names!=NULL)
                    {
                        fprintf(file,"Rtransf-%s_%d:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->gen_length);
                    }
                    else
                    {
                        fprintf(file,"Rtransf-%d_%d:%.8lf",p->conts->index,p->paralog,p->gen_length);
                    }
                    break;
                case RGC:
                    if (p->conts==NULL)
                    {
                        fprintf(file,"Rgc-%d_%d:%.8lf",p->index,p->paralog,p->gen_length);
                    }
                    else if (p->conts->n_child==0 && names!=NULL)
                    {
                        fprintf(file,"Rgc-%s_%d:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->gen_length);
                    }
                    else
                    {
                        fprintf(file,"Rgc-%d_%d:%.8lf",p->conts->index,p->paralog,p->gen_length);
                    }
                    break;
            }
                break;
                
            default:
                // **
                /// <dl><dt>Else (internal node):</dt><dd>
                // *
                /// Prints "(" and starts the post-order recursion (child loop)
                fprintf(file,"(");
                for (i=0;i<p->n_child-1;++i)
                {
                    // *
                    /// Calls itself using each different child. After each call (except last child) prints ","
                    WriteLNodesFileGenIntlabel(file,*(p->children+i),names);
                    fprintf(file,",");
                }
                WriteLNodesFileGenIntlabel(file,*(p->children+i), names);
                
                // **
                /// Prints the information of this internal node (closing it with ")") or finishes the tree if the node is the root</dd></dl></dd></dl>
                if (p->anc_node !=NULL)
                    fprintf(file,")%d:%.8lf",p->index,p->gen_length);
                else
                    fprintf(file,");\n");
                break;
        }
    }
}

void WriteLNodesFileTime (FILE * file,l_node * p, name_c * names, double gen_time)
{
    int i=0;
    
    // ***
    /// <dl><dt> Function structure </dt><dd>
    
    if (p != NULL)
    {
        switch (p->n_child)
        {
            case 0:
                // **
                /// <dl><dt>If the node is a leaf:</dt><dd>
                // *
                /// Prints the node name and its branch length.</dd></dl>
                
                switch (p->kind_node)
            {
                default:
                    if (names==NULL)
                        fprintf(file,"%d_%d:%.8lf",p->sp_index,p->paralog,p->gen_length*gen_time*p->gtime_mult);
                    else
                        fprintf(file,"%s_%d:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->gen_length*gen_time*p->gtime_mult);
                    break;
                case LOSS:
                    if (p->conts->n_child==0 && names!=NULL)
                    {
                        fprintf(file,"Lost-%s_%d:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->gen_length*gen_time*p->gtime_mult);
                    }
                    else
                    {
                        fprintf(file,"Lost-%d_%d:%.8lf",p->conts->index,p->paralog,p->gen_length*gen_time*p->gtime_mult);
                    }
                    break;
                case RTRFR:
                    if (p->conts->n_child==0 && names!=NULL)
                    {
                        fprintf(file,"Rtransf-%s_%d:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->gen_length*gen_time*p->gtime_mult);
                    }
                    else
                    {
                        fprintf(file,"Rtransf-%d_%d:%.8lf",p->conts->index,p->paralog,p->gen_length*gen_time*p->gtime_mult);
                    }
                    break;
                case RGC:
                    if (p->conts->n_child==0 && names!=NULL)
                    {
                        fprintf(file,"Rgc-%s_%d:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->gen_length*gen_time*p->gtime_mult);
                    }
                    else
                    {
                        fprintf(file,"Rgc-%d_%d:%.8lf",p->conts->index,p->paralog,p->gen_length*gen_time*p->gtime_mult);
                    }
                    break;
            }
                break;
                
            default:
                // **
                /// <dl><dt>Else (internal node):</dt><dd>
                // *
                /// Prints "(" and starts the post-order recursion (child loop)
                fprintf(file,"(");
                for (i=0;i<p->n_child-1;++i)
                {
                    // *
                    /// Calls itself using each different child. After each call (except last child) prints ","
                    WriteLNodesFileTime(file,*(p->children+i),names,gen_time);
                    fprintf(file,",");
                }
                WriteLNodesFileTime(file,*(p->children+i), names,gen_time);
                
                // **
                /// Prints the information of this internal node (closing it with ")") or finishes the tree if the node is the root</dd></dl></dd></dl>
                if (p->anc_node !=NULL)
                    fprintf(file,"):%.8lf",p->gen_length*gen_time*p->gtime_mult);
                else
                    fprintf(file,");\n");
                break;
        }
    }
}

void WriteLNodesFileTimeIntlabel (FILE * file,l_node * p, name_c * names, double gen_time)
{
    int i=0;
    
    // ***
    /// <dl><dt> Function structure </dt><dd>
    
    if (p != NULL)
    {
        switch (p->n_child)
        {
            case 0:
                // **
                /// <dl><dt>If the node is a leaf:</dt><dd>
                // *
                /// Prints the node name and its branch length.</dd></dl>
                
                switch (p->kind_node)
            {
                default:
                    if (names==NULL)
                        fprintf(file,"%d_%d:%.8lf",p->sp_index,p->paralog,p->gen_length*gen_time*p->gtime_mult);
                    else
                        fprintf(file,"%s_%d:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->gen_length*gen_time*p->gtime_mult);
                    break;
                case LOSS:
                    if (p->conts->n_child==0 && names!=NULL)
                    {
                        fprintf(file,"Lost-%s_%d:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->gen_length*gen_time*p->gtime_mult);
                    }
                    else
                    {
                        fprintf(file,"Lost-%d_%d:%.8lf",p->conts->index,p->paralog,p->gen_length*gen_time*p->gtime_mult);
                    }
                    break;
                case RTRFR:
                    if (p->conts->n_child==0 && names!=NULL)
                    {
                        fprintf(file,"Rtransf-%s_%d:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->gen_length*gen_time*p->gtime_mult);
                    }
                    else
                    {
                        fprintf(file,"Rtransf-%d_%d:%.8lf",p->conts->index,p->paralog,p->gen_length*gen_time*p->gtime_mult);
                    }
                    break;
                case RGC:
                    if (p->conts->n_child==0 && names!=NULL)
                    {
                        fprintf(file,"Rgc-%s_%d:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->gen_length*gen_time*p->gtime_mult);
                    }
                    else
                    {
                        fprintf(file,"Rgc-%d_%d:%.8lf",p->conts->index,p->paralog,p->gen_length*gen_time*p->gtime_mult);
                    }
                    break;
            }
                break;
                
            default:
                // **
                /// <dl><dt>Else (internal node):</dt><dd>
                // *
                /// Prints "(" and starts the post-order recursion (child loop)
                fprintf(file,"(");
                for (i=0;i<p->n_child-1;++i)
                {
                    // *
                    /// Calls itself using each different child. After each call (except last child) prints ","
                    WriteLNodesFileTimeIntlabel(file,*(p->children+i),names,gen_time);
                    fprintf(file,",");
                }
                WriteLNodesFileTimeIntlabel(file,*(p->children+i), names,gen_time);
                
                // **
                /// Prints the information of this internal node (closing it with ")") or finishes the tree if the node is the root</dd></dl></dd></dl>
                if (p->anc_node !=NULL)
                    fprintf(file,")%d:%.8lf",p->index,p->gen_length*gen_time*p->gtime_mult);
                else
                    fprintf(file,");\n");
                break;
        }
    }
}

void WriteDaughtersNodesFile(FILE * file,l_node * p, name_c *names)
{
    int i=0;
    l_node * daughter=NULL;
    
    // ***
    /// <dl><dt> Function structure </dt><dd>
    
    if (p != NULL)
    {
        switch (p->n_child)
        {
 
            default:
                // **
                /// <dl><dt>Internal node</dt><dd>
                // *
                /// Post-order recursion (child loop)
                for (i=0;i<p->n_child;++i)
                {
                    // *
                    /// Calls itself using each different child.
                    WriteDaughtersNodesFile(file,*(p->children+i),names);
                }
                // **
                /// Prints the identity of the daughter node</dd></dl></dd></dl>
                switch (p->kind_node)
                {
                    case DUP:
                    case TRFR:
                    case GC:
                        daughter=*(p->children+1);
                        if (daughter->n_child!=0)
                        {
                            fprintf(file,"%d,",daughter->index);
                        }
                        else if (names!=NULL)
                        {
                            fprintf(file,"\'%s_%d\',",(names->names+(daughter->sp_index*names->max_lname)),daughter->paralog);
                        }
                        else
                        {
                            fprintf(file,"\'%d_%d\',",daughter->conts->index,daughter->paralog);
                        }
                        break;
                }
                break;
            
            case 0:
                break;
               
        }
    }
}

void WriteGNodes (g_node * p, name_c * names)
{
    int i=0;
    
    // ***
    /// <dl><dt> Function structure </dt><dd>
    if (p != NULL)
    {
        switch (p->n_child)
        {
            case 0:
                // **
                /// <dl><dt>If the node is a leaf:</dt><dd>
                // *
                /// Prints the node name and its branch length.</dd></dl>
                switch (p->contl->kind_node)
                {
                    default:
                        if (names==NULL)
                            printf("%d_%d_%d:%.8lf",p->sp_index,p->paralog, p->replica,p->bl);
                        
                        else
                            printf("%s_%d_%d:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->replica,p->bl);
                        break;
                    case LOSS:
                        if (p->conts->n_child==0 && names!=NULL)
                        {
                            printf("Lost-%s_%d_0:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->bl);
                        }
                        else
                        {
                            printf("Lost-%d_%d_0:%.8lf",p->conts->index,p->paralog,p->bl);
                        }
                        break;
                    case RTRFR:
                        if (p->conts->n_child==0 && names!=NULL)
                        {
                            printf("Rtransf-%s_%d_0:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->bl);
                        }
                        else
                        {
                            printf("Rtransf-%d_%d_0:%.8lf",p->conts->index,p->paralog,p->bl);
                        }
                        break;
                    case RGC:
                        if (p->conts->n_child==0 && names!=NULL)
                        {
                            printf("Rgc-%s_%d_0:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->bl);
                        }
                        else
                        {
                            printf("Rgc-%d_%d_0:%.8lf",p->conts->index,p->paralog,p->bl);
                        }
                        break;
                }
                break;
                
            default:
                // **
                /// <dl><dt>Else (internal node):</dt><dd>
                // *
                /// Prints "(" and starts the post-order recursion (child loop)
                printf("(");
                for (i=0;i<p->n_child-1;++i)
                {
                    // *
                    /// Calls itself using each different child. After each call (except last child) prints ","
                    WriteGNodes(*(p->children+i),names);
                    printf(",");
                }
                WriteGNodes(*(p->children+i),names);
                
                // **
                /// Prints the information of this internal node (closing it with ")") or finishes the tree if the node is the root</dd></dl></dd></dl>
                if (p->anc_node !=NULL)
                    printf("):%.8lf",p->bl);
                
                else
                {
                    //printf("):%.8lf);",p->gen_length);//To print the length of the p
                    printf(");");
                }
                break;
        }
    }
}

void WriteGNodesIntlabel (g_node * p, name_c * names)
{
    int i=0;
    
    // ***
    /// <dl><dt> Function structure </dt><dd>
    if (p != NULL)
    {
        switch (p->n_child)
        {
            case 0:
                // **
                /// <dl><dt>If the node is a leaf:</dt><dd>
                // *
                /// Prints the node name and its branch length.</dd></dl>
                switch (p->contl->kind_node)
            {
                default:
                    if (names==NULL)
                        printf("%d_%d_%d:%.8lf",p->sp_index,p->paralog, p->replica,p->bl);
                    
                    else
                        printf("%s_%d_%d:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->replica,p->bl);
                    break;
                case LOSS:
                    if (p->conts->n_child==0 && names!=NULL)
                    {
                        printf("Lost-%s_%d_0:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->bl);
                    }
                    else
                    {
                        printf("Lost-%d_%d_0:%.8lf",p->conts->index,p->paralog,p->bl);
                    }
                    break;
                case RTRFR:
                    if (p->conts->n_child==0 && names!=NULL)
                    {
                        printf("Rtransf-%s_%d_0:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->bl);
                    }
                    else
                    {
                        printf("Rtransf-%d_%d_0:%.8lf",p->conts->index,p->paralog,p->bl);
                    }
                    break;
                case RGC:
                    if (p->conts->n_child==0 && names!=NULL)
                    {
                        printf("Rgc-%s_%d_0:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->bl);
                    }
                    else
                    {
                        printf("Rgc-%d_%d_0:%.8lf",p->conts->index,p->paralog,p->bl);
                    }
                    break;
            }
                break;
                
            default:
                // **
                /// <dl><dt>Else (internal node):</dt><dd>
                // *
                /// Prints "(" and starts the post-order recursion (child loop)
                printf("(");
                for (i=0;i<p->n_child-1;++i)
                {
                    // *
                    /// Calls itself using each different child. After each call (except last child) prints ","
                    WriteGNodesIntlabel(*(p->children+i),names);
                    printf(",");
                }
                WriteGNodesIntlabel(*(p->children+i),names);
                
                // **
                /// Prints the information of this internal node (closing it with ")") or finishes the tree if the node is the root</dd></dl></dd></dl>
                if (p->anc_node !=NULL)
                    printf(")%d:%.8lf",p->index,p->bl);
                
                else
                {
                    //printf("):%.8lf);",p->gen_length);//To print the length of the p
                    printf(");");
                }
                break;
        }
    }
}

void WriteGNodesStr (char * str, g_node * p, name_c * names)
{
    int i=0;
    
    // ***
    /// <dl><dt> Function structure </dt><dd>
    
    if (p != NULL)
    {
        switch (p->n_child)
        {
            case 0:
                // **
                /// <dl><dt>If the node is a leaf:</dt><dd>
                // *
                /// Prints the node name and its branch length.</dd></dl>
                switch (p->contl->kind_node)
                {
                    default:
                        if (names==NULL)
                            sprintf(str,"%s%d_%d_%d:%.8lf",str,p->sp_index,p->paralog,p->replica,p->bl);
                        
                        else
                            sprintf(str,"%s%s_%d_%d:%.8lf",str,(names->names+(p->sp_index*names->max_lname)),p->paralog,p->replica,p->bl);
                        break;
                        
                    case LOSS:
                        if (p->conts->n_child==0 && names!=NULL)
                        {
                            sprintf(str,"%sLost-%s_%d_0:%.8lf",str,(names->names+(p->sp_index*names->max_lname)),p->paralog,p->bl);
                        }
                        else
                        {
                            sprintf(str,"%sLost-%d_%d_0:%.8lf",str,p->conts->index,p->paralog,p->bl);
                        }
                        break;
                    case RTRFR:
                        if (p->conts->n_child==0 && names!=NULL)
                        {
                            sprintf(str,"%sRtransf-%s_%d_0:%.8lf",str,(names->names+(p->sp_index*names->max_lname)),p->paralog,p->bl);
                        }
                        else
                        {
                            sprintf(str,"%sRtransf-%d_%d_0:%.8lf",str,p->conts->index,p->paralog,p->bl);
                        }
                        break;
                    case RGC:
                        if (p->conts->n_child==0 && names!=NULL)
                        {
                            sprintf(str,"%sRgc-%s_%d_0:%.8lf",str,(names->names+(p->sp_index*names->max_lname)),p->paralog,p->bl);
                        }
                        else
                        {
                            sprintf(str,"%sRgc-%d_%d_0:%.8lf",str,p->conts->index,p->paralog,p->bl);
                        }
                        break;
                        

                }
                break;
                
            default:
                // **
                /// <dl><dt>Else (internal node):</dt><dd>
                // *
                /// Prints "(" and starts the post-order recursion (child loop)
                sprintf(str,"%s(",str);
                for (i=0;i<p->n_child-1;++i)
                {
                    // *
                    /// Calls itself using each different child. After each call (except last child) prints ","
                    WriteGNodesStr(str,*(p->children+i),names);
                    printf(",");
                    
                }
                WriteGNodesStr(str,*(p->children+i),names);
                
                // **
                /// Prints the information of this internal node (closing it with ")") or finishes the tree if the node is the root</dd></dl></dd></dl>
                if (p->anc_node !=NULL)
                    sprintf(str,"%s):%.8lf",str, p->bl);
                
                else
                {
                    //                sprintf(str,"(%s):%.8lf);",str, p->gen_length);//To print the length of the p
                    sprintf(str,"%s;",str);
                }
                break;
        }
    }
}

void WriteGNodesStrIntlabel (char * str, g_node * p, name_c * names)
{
    int i=0;
    
    // ***
    /// <dl><dt> Function structure </dt><dd>
    
    if (p != NULL)
    {
        switch (p->n_child)
        {
            case 0:
                // **
                /// <dl><dt>If the node is a leaf:</dt><dd>
                // *
                /// Prints the node name and its branch length.</dd></dl>
                switch (p->contl->kind_node)
            {
                default:
                    if (names==NULL)
                        sprintf(str,"%s%d_%d_%d:%.8lf",str,p->sp_index,p->paralog,p->replica,p->bl);
                    
                    else
                        sprintf(str,"%s%s_%d_%d:%.8lf",str,(names->names+(p->sp_index*names->max_lname)),p->paralog,p->replica,p->bl);
                    break;
                    
                case LOSS:
                    if (p->conts->n_child==0 && names!=NULL)
                    {
                        sprintf(str,"%sLost-%s_%d_0:%.8lf",str,(names->names+(p->sp_index*names->max_lname)),p->paralog,p->bl);
                    }
                    else
                    {
                        sprintf(str,"%sLost-%d_%d_0:%.8lf",str,p->conts->index,p->paralog,p->bl);
                    }
                    break;
                case RTRFR:
                    if (p->conts->n_child==0 && names!=NULL)
                    {
                        sprintf(str,"%sRtransf-%s_%d_0:%.8lf",str,(names->names+(p->sp_index*names->max_lname)),p->paralog,p->bl);
                    }
                    else
                    {
                        sprintf(str,"%sRtransf-%d_%d_0:%.8lf",str,p->conts->index,p->paralog,p->bl);
                    }
                    break;
                case RGC:
                    if (p->conts->n_child==0 && names!=NULL)
                    {
                        sprintf(str,"%sRgc-%s_%d_0:%.8lf",str,(names->names+(p->sp_index*names->max_lname)),p->paralog,p->bl);
                    }
                    else
                    {
                        sprintf(str,"%sRgc-%d_%d_0:%.8lf",str,p->conts->index,p->paralog,p->bl);
                    }
                    break;
                    
                    
            }
                break;
                
            default:
                // **
                /// <dl><dt>Else (internal node):</dt><dd>
                // *
                /// Prints "(" and starts the post-order recursion (child loop)
                sprintf(str,"%s(",str);
                for (i=0;i<p->n_child-1;++i)
                {
                    // *
                    /// Calls itself using each different child. After each call (except last child) prints ","
                    WriteGNodesStrIntlabel(str,*(p->children+i),names);
                    printf(",");
                    
                }
                WriteGNodesStrIntlabel(str,*(p->children+i),names);
                
                // **
                /// Prints the information of this internal node (closing it with ")") or finishes the tree if the node is the root</dd></dl></dd></dl>
                if (p->anc_node !=NULL)
                    sprintf(str,"%s)%d:%.8lf",str,p->index,p->bl);
                
                else
                {
                    //                sprintf(str,"(%s):%.8lf);",str, p->gen_length);//To print the length of the p
                    sprintf(str,"%s;",str);
                }
                break;
        }
    }
}

void WriteGNodesFile (FILE * file, g_node * p, name_c * names)
{
    int i=0;
    
    // ***
    /// <dl><dt> Function structure </dt><dd>
    
    if (p != NULL)
    {
        switch (p->n_child)
        {
            case 0:
                // **
                /// <dl><dt>If the node is a leaf:</dt><dd>
                // *
                /// Prints the node name and its branch length.</dd></dl>
                switch (p->contl->kind_node)
                {
                    default:
                        if (names==NULL)
                            fprintf(file,"%d_%d_%d:%.8lf",p->sp_index,p->paralog, p->replica,p->bl);
                        
                        else
                            fprintf(file,"%s_%d_%d:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->replica,p->bl);
                        break;
                    case LOSS:
                        if (p->conts->n_child==0 && names!=NULL)
                        {
                            fprintf(file,"Lost-%s_%d_0:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->bl);
                        }
                        else
                        {
                            fprintf(file,"Lost-%d_%d_0:%.8lf",p->conts->index,p->paralog,p->bl);
                        }
                        break;
                    case RTRFR:
                        if (p->conts->n_child==0 && names!=NULL)
                        {
                            fprintf(file,"Rtransf-%s_%d_0:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->bl);
                        }
                        else
                        {
                            fprintf(file,"Rtransf-%d_%d_0:%.8lf",p->conts->index,p->paralog,p->bl);
                        }
                        break;
                    case RGC:
                        if (p->conts->n_child==0 && names!=NULL)
                        {
                            fprintf(file,"Rgc-%s_%d_0:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->bl);
                        }
                        else
                        {
                            fprintf(file,"Rgc-%d_%d_0:%.8lf",p->conts->index,p->paralog,p->bl);
                        }
                        break;
                }
                break;
                
            default:
                // **
                /// <dl><dt>Else (internal node):</dt><dd>
                // *
                /// Prints "(" and starts the post-order recursion (child loop)
                fprintf(file,"(");
                for (i=0;i<p->n_child-1;++i)
                {
                    WriteGNodesFile(file, *(p->children+i), names);
                    fprintf(file,",");
                }
                WriteGNodesFile(file, *(p->children+i), names);
                
                // **
                /// Prints the information of this internal node (closing it with ")") or finishes the tree if the node is the root</dd></dl></dd></dl>
                if (p->anc_node !=NULL)
                    fprintf(file,"):%.8lf",p->bl);
                
                else
                {
                    //                fprintf(file,"):%.8lf);\n",p->gen_length); //To print the length of the p
                    fprintf(file,");\n");
                }
                break;
        }
    }
}

void WriteGNodesFileIntlabel (FILE * file, g_node * p, name_c * names)
{
    int i=0;
    
    // ***
    /// <dl><dt> Function structure </dt><dd>
    
    if (p != NULL)
    {
        switch (p->n_child)
        {
            case 0:
                // **
                /// <dl><dt>If the node is a leaf:</dt><dd>
                // *
                /// Prints the node name and its branch length.</dd></dl>
                switch (p->contl->kind_node)
            {
                default:
                    if (names==NULL)
                        fprintf(file,"%d_%d_%d:%.8lf",p->sp_index,p->paralog, p->replica,p->bl);
                    
                    else
                        fprintf(file,"%s_%d_%d:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->replica,p->bl);
                    break;
                case LOSS:
                    if (p->conts->n_child==0 && names!=NULL)
                    {
                        fprintf(file,"Lost-%s_%d_0:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->bl);
                    }
                    else
                    {
                        fprintf(file,"Lost-%d_%d_0:%.8lf",p->conts->index,p->paralog,p->bl);
                    }
                    break;
                case RTRFR:
                    if (p->conts->n_child==0 && names!=NULL)
                    {
                        fprintf(file,"Rtransf-%s_%d_0:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->bl);
                    }
                    else
                    {
                        fprintf(file,"Rtransf-%d_%d_0:%.8lf",p->conts->index,p->paralog,p->bl);
                    }
                    break;
                case RGC:
                    if (p->conts->n_child==0 && names!=NULL)
                    {
                        fprintf(file,"Rgc-%s_%d_0:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->bl);
                    }
                    else
                    {
                        fprintf(file,"Rgc-%d_%d_0:%.8lf",p->conts->index,p->paralog,p->bl);
                    }
                    break;
            }
                break;
                
            default:
                // **
                /// <dl><dt>Else (internal node):</dt><dd>
                // *
                /// Prints "(" and starts the post-order recursion (child loop)
                fprintf(file,"(");
                for (i=0;i<p->n_child-1;++i)
                {
                    WriteGNodesFileIntlabel(file, *(p->children+i), names);
                    fprintf(file,",");
                }
                WriteGNodesFileIntlabel(file, *(p->children+i), names);
                
                // **
                /// Prints the information of this internal node (closing it with ")") or finishes the tree if the node is the root</dd></dl></dd></dl>
                if (p->anc_node !=NULL)
                    fprintf(file,")%d:%.8lf",p->index,p->bl);
                
                else
                {
                    //                fprintf(file,"):%.8lf);\n",p->gen_length); //To print the length of the p
                    fprintf(file,");\n");
                }
                break;
        }
    }
}

long int RelinkLGTrees(l_tree * wl_tree, g_node * m_node)
{
    int i=0, j=0;
    l_node * w_lnode=NULL;
    
    // **
    /// <dl><dt>Locus tree memory check</dt><dd>
    if (wl_tree->m_node==NULL)
    {
        return MEM_ERROR;
    }
    
    //**
    /// <dl><dt>Locus tree nodes loop</dt><dd>
    for (i=0;i<wl_tree->n_nodes;++i)
    {
        w_lnode=(wl_tree->m_node+i);
        // *
        /// Locus node g_node pointer loop</dd></dl></dd></dl>
        for (j=0; j<w_lnode->n_ilin;++j)
        {
            *(w_lnode->g_nodes+j)=m_node+(*(w_lnode->g_nodes+j))->index;
        }
    }
    return NO_ERROR;
}

l_node * ChooseLNodePeriod(l_node **l_pointers, int n_nodes, l_node * t_node, double u_num, int verbosity)
{
    int i=0, j=0;
    int done=0, next_t=1;
    double *t_times=calloc(n_nodes,sizeof(double));
    double t_prob=0, c_prob=0;
    l_node *w_lnode=NULL,*w_tnode=t_node, **wl_pointers;
    
    wl_pointers=calloc(n_nodes,sizeof(l_node *));
    
    if (t_times==NULL || wl_pointers==NULL)
        ErrorReporter(MEM_ERROR,NULL);
    
    for (j=0; j<n_nodes; ++j) //Reset
    {
        *(t_times+j)=-1;
    }
    memmove(wl_pointers,l_pointers,sizeof(l_node *)*n_nodes);
    
    while (i<MAX_IT && done==0) //while it gets all the distances or the loop reaches MAX_IT
    {
        done=1;
        next_t=1;
        for (j=0; j<n_nodes;++j)
        {
            if (*(t_times+j)!=-1) //We already have the time for this node. If all of them enter here, done==1, and the program continue, leaving the while loop.
                continue;
            w_lnode=*(wl_pointers+j);
            done=0;
            if ((w_lnode->time>w_tnode->time) && (w_lnode != w_tnode)) //We have to consider the ancestor of the receptor, but the t_node can still be a good choice (next_t=0)
            {
                *(wl_pointers+j)=w_lnode->anc_node;
                next_t=0;
            }
            else if (w_lnode==w_tnode)
            {
                *(t_times+j)=w_lnode->time;
            }
        }
        if (next_t==1 && done==0)
        {
            w_tnode=w_tnode->anc_node; //All the candidates below the t_node, so we need to consider its ancestor
        }
        ++i;
    }
    
    if (verbosity>5)
    {
        printf("\n\t\t\t\t\t\tSampling receptor inversely proportional to MRCA time: (Node, inverse distance) ");
#ifdef DBG
        fflush(stdout);
#endif
    }
    for (j=0; j<n_nodes;++j) //Calculating the inverse of the distance and the total inverse distance
    {
        if (t_node->time==*(t_times+j))
            *(t_times+j)=0;
        else
        {
            *(t_times+j)=1/(t_node->time-*(t_times+j)); //Inverse distance
            t_prob+=*(t_times+j);
        }
        if (verbosity>5)
        {
            printf("%u %lf, ", (*(l_pointers+j))->index, *(t_times+j));
#ifdef DBG
            fflush(stdout);
#endif
        }
    }
    
    w_lnode=NULL;
    
    if (verbosity>5)
    {
        printf("\n\t\t\t\t\t\tSampled probability %lf, relative probabilities (Node, relative probability, cumulative probability): ",u_num);
#ifdef DBG
        fflush(stdout);
#endif
    }
    
    for (j=0; j<n_nodes;++j) //Samples the probabilities
    {
        if (u_num>c_prob+(*(t_times+j)/t_prob))
        {
            c_prob+=(*(t_times+j)/t_prob);
            if (verbosity>5)
            {
                printf("%u %lf %lf, ", (*(l_pointers+j))->index, (*(t_times+j)/t_prob), c_prob);
#ifdef DBG
                fflush(stdout);
#endif
            }
        }
        else
        {
            if (verbosity>5)
            {
                printf("Selected: %u %lf %lf, ", (*(l_pointers+j))->index, (*(t_times+j)/t_prob), c_prob+(*(t_times+j)/t_prob));
#ifdef DBG
                fflush(stdout);
#endif
            }
            w_lnode=*(l_pointers+j);
            break;
        }
    }
    
    free(t_times);
    free(wl_pointers);
    
    if (w_lnode==NULL)
        ErrorReporter(UNEXPECTED_VALUE,NULL);
    
    return w_lnode;
}

int firstnoblank(char *string)
{
    char i=1;
    int it=0;
    
    while (i!='\0' && it<MAX_IT)
    {
        i=*(string+it);
        if (isblank(i)==0)
            break;
        ++it;
    }
    if (it>=MAX_IT)
    {
        ErrorReporter(LOOP_ERROR,NULL);
        return EXIT_FAILURE;
    }
    else
    {
        return it;
    }
}

long int GetSnodeParamsFromNexusComments(char *string,s_node *node,int *n_char)
{
    int i=0,c_param,out_param=1, step=0;
    const int n_params=4;
    const char *params[4]={"pop_size","n_ind","u_mult","g_mult"};
    const int s_params[4]={8,5,6,6};
    char code, *buffer=NULL;
    size_t tbuffer=NUM_BUFFER;
    
    reallocBuffer(&buffer,&tbuffer,tbuffer);
    ResetBuffer(buffer, tbuffer);
    
    while (*(string+*n_char)!=']' && i<MAX_IT)
    {
        if (out_param==1)
        {
            c_param=-1;
            for (step=0; step<n_params;++step)
            {
                if (strncmp(params[step], string+*n_char, s_params[step])==0)
                {
                    c_param=step;
                }
            }
            if (c_param==-1)
            {
                PrintXCharError(string, *n_char+8, "\nNEXUS TREE PARSING ERROR\n", "|<- Unvalid param, please, revisit the manual to check the proper parameter spelling\n"); //+8, maximum parameter width
                return SETTINGS_ERROR;
            }
            else if (*(string+*n_char+s_params[c_param])=='=')
            {
                out_param=0;
                *n_char+=s_params[c_param]+1;
            }
        }
        else
        {
            ResetBuffer(buffer, tbuffer);
            step=0;
            code=*(string+*n_char);
            while (code<58 && (code>47 || code==46) && step<MAX_IT && step<tbuffer)
            {
                *(buffer+step)=code;
                ++*n_char;
                ++step;
                code=*(string+*n_char);
            }
            if (step==MAX_IT)
            {
                PrintXCharError(string, *n_char-step, "\nNEWICK PARSING ERROR\n", "|<- Unvalid value, please, check the NEXUS file\n");
                return LOOP_ERROR;
            }
            else if (step==tbuffer)
            {
                *n_char-=step;
                reallocBuffer(&buffer, &tbuffer, tbuffer*10);
            }
            else
            {
                *(buffer+step)=0;
                switch (c_param)
                {
                    case 0:
                        if(sscanf(buffer,"%d",&node->Ne)!=1)
                        {
                            PrintXCharError(string, *n_char, "\nNEWICK PARSING ERROR\n", "|<- Unvalid value, please, check the NEXUS file\n");
                            return SETTINGS_ERROR;
                        }
                        break;
                    case 1:
                        if(sscanf(buffer,"%d",&node->n_replicas)!=1)
                        {
                            PrintXCharError(string, *n_char, "\nNEWICK PARSING ERROR\n", "|<- Unvalid value, please, check the NEXUS file\n");
                            return SETTINGS_ERROR;
                        }
                        break;
                    case 2:
                        if(sscanf(buffer,"%lf",&node->mu_mult)!=1)
                        {
                            PrintXCharError(string, *n_char, "\nNEWICK PARSING ERROR\n", "|<- Unvalid value, please, check the NEXUS file\n");
                            return SETTINGS_ERROR;
                        }
                        break;
                    case 3:
                        if(sscanf(buffer,"%lf",&node->gtime_mult)!=1)
                        {
                            PrintXCharError(string, *n_char, "\nNEWICK PARSING ERROR\n", "|<- Unvalid value, please, check the NEXUS file\n");
                            return SETTINGS_ERROR;
                        }
                        break;
                        
                    default:
                        return SETTINGS_ERROR;
                        break;
                }
                if (*(string+*n_char)==',')
                    ++*n_char;
                out_param=1;
            }
        }
        
        ++i;
    }
    if(i==MAX_IT)
        return LOOP_ERROR;
    
    free(buffer);
    return NO_ERROR;
}

long int GetLnodeParamsFromNexusComments(char *string,l_node *node, int *n_char)
{
    int i=0,c_param,out_param=1, step=0;
    const int n_params=6;
    const char *params[6]={"pop_size","kind_n","paralog","n_ind","u_mult","g_mult"};
    const int s_params[6]={8,6,7,5,6,6};
    char code, *buffer=NULL;
    size_t tbuffer=NUM_BUFFER;
    
    reallocBuffer(&buffer,&tbuffer,tbuffer);
    ResetBuffer(buffer, tbuffer);
    
    while (*(string+*n_char)!=']' && i<=MAX_IT)
    {
        if (out_param==1)
        {
            c_param=-1;
            for (step=0; step<n_params;++step)
            {
                if (strncmp(params[step], string+*n_char, s_params[step])==0)
                {
                    c_param=step;
                }
            }
            if (c_param==-1)
            {
                PrintXCharError(string, *n_char+8, "\nNEWICK PARSING ERROR\n", "|<- Unvalid param, please, revisit the manual to check the proper parameter spelling\n"); //+8, maximum parameter width
                return SETTINGS_ERROR;
            }
            else if (*(string+*n_char+s_params[c_param])=='=')
            {
                out_param=0;
                *n_char+=s_params[c_param]+1;
            }
        }
        else
        {
            ResetBuffer(buffer, tbuffer);
            step=0;
            code=*(string+*n_char);
            while (code<58 && (code>47 || code==46) && step<MAX_IT && step<tbuffer)
            {
                *(buffer+step)=code;
                ++*n_char;
                ++step;
                code=*(string+*n_char);
            }
            if (step==MAX_IT)
            {
                PrintXCharError(string, *n_char-step, "\nNEWICK PARSING ERROR\n", "|<- Unvalid value, please, check the NEXUS file\n");
                return LOOP_ERROR;
            }
            else if (step==tbuffer)
            {
                *n_char-=step;
                reallocBuffer(&buffer, &tbuffer, tbuffer*10);
            }
            else
            {
                *(buffer+step)=0;
                switch (c_param)
                {
                    case 0:
                        if(sscanf(buffer,"%d",&node->Ne)!=1)
                        {
                            PrintXCharError(string, *n_char, "\nNEWICK PARSING ERROR\n", "|<- Unvalid value, please, check the NEXUS file\n");
                            return SETTINGS_ERROR;
                        }
                        break;
                    case 1:
                        if(sscanf(buffer,"%d",&node->kind_node)!=1)
                        {
                            PrintXCharError(string, *n_char, "\nNEWICK PARSING ERROR\n", "|<- Unvalid value, please, check the NEXUS file\n");
                            return SETTINGS_ERROR;
                        }
                        break;
                    case 2:
                        if(sscanf(buffer,"%d",&node->paralog)!=1)
                        {
                            PrintXCharError(string, *n_char, "\nNEWICK PARSING ERROR\n", "|<- Unvalid value, please, check the NEXUS file\n");
                            return SETTINGS_ERROR;
                        }
                        break;
                    case 3:
                        if(sscanf(buffer,"%d",&node->n_nodes)!=1)
                        {
                            PrintXCharError(string, *n_char, "\nNEWICK PARSING ERROR\n", "|<- Unvalid value, please, check the NEXUS file\n");
                            return SETTINGS_ERROR;
                        }
                        node->n_ilin=node->n_nodes;
                        break;
                    case 4:
                        if(sscanf(buffer,"%lf",&node->mu_mult)!=1)
                        {
                            PrintXCharError(string, *n_char, "\nNEWICK PARSING ERROR\n", "|<- Unvalid value, please, check the NEXUS file\n");
                            return SETTINGS_ERROR;
                        }
                        break;
                    case 5:
                        if(sscanf(buffer,"%lf",&node->gtime_mult)!=1)
                        {
                            PrintXCharError(string, *n_char, "\nNEWICK PARSING ERROR\n", "|<- Unvalid value, please, check the NEXUS file\n");
                            return SETTINGS_ERROR;
                        }
                        break;
                        
                    default:
                        return SETTINGS_ERROR;
                        break;
                }
                if (*(string+*n_char)==',')
                    ++*n_char;
                out_param=1;
            }
        }
        
        ++i;
    }
    if(i>MAX_IT)
        return LOOP_ERROR;
    
    free(buffer);
    return NO_ERROR;
}

double CheckUltrametricitySNodes(s_node *node)
{
    int i=0;
    double time=0, n_time=0;
    
    if (node->n_child!=0)
    {
        time=CheckUltrametricitySNodes(*(node->children));
        if (time==-1)
            return -1;
        else
            for (i=1; i<node->n_child; ++i)
            {
                n_time=CheckUltrametricitySNodes(*(node->children+i));
                if(n_time == -1 || fabs(time-n_time)>FLOAT_PRECISION)
                    return -1;
            }
        return time;
    }
    else
    {
        return node->time;
    }
}

double CheckUltrametricityLNodes(l_node *node)
{
    int i=0;
    double time=0, n_time=0;
    
    if (node->n_child!=0)
    {
        time=CheckUltrametricityLNodes(*(node->children));
        if (time==-1)
            return -1;
        else
            for (i=1; i<node->n_child; ++i)
            {
                if ((*(node->children+i))->kind_node!=LOSS && (*(node->children+i))->kind_node!=RTRFR && (*(node->children+i))->kind_node!=RGC)
                {
                    n_time=CheckUltrametricityLNodes(*(node->children+i));
                    if(n_time == -1 || fabs(time-n_time)>FLOAT_PRECISION)
                        return -1;
                }
            }
        return time;
    }
    else
    {
        return node->time;
    }
}

