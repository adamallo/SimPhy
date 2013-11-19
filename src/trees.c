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
 * including s_node::childs.
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
 * including l_node::childs and l_node::g_nodes (this could be avoid using the
 * free_gnodes tag).
 *
 * \param node
 *  l_node pointer (root of the group of l_nodes in the first call).
 * \param free_gnodes
 *  Logical flag. If =0 the g_node pointers (\ref l_node::g_nodes) are not freed.
 *
 *******************************************************************************/
static void FreeLNodes (l_node * node, unsigned int free_gnodes);

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
static void FreeGNodes (g_node * node, unsigned int complete);

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
void PostReorderSNodes (s_node * node, unsigned int * index); ///\todo Undo this provisional trick (removing the static status) to have this function accesible from the main function of the program.

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
static void PostReorderLNodes (l_node * node, unsigned int * index);

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
static void PostReorderGNodes (g_node * node, unsigned int * index); 

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
static void PreReorderSNodes (s_node * node, unsigned int * index);

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
static void PreReorderLNodes (l_node * node, unsigned int * index);

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
static void PreReorderGNodes (g_node * node, unsigned int * index);

// ** Node copy ** //

/**
 * Copies a group of s_nodes (with tree structure) in an array of s_nodes following
 * a post-order.
 *
 * It clones a tree/subtree (group of \ref s_node "s_nodes" with tree structure) in an array
 * of \ref s_node "s_nodes". This is an important part of \ref CollapseSTree public function.
 *
 * \param output
 *  Pre-allocated memory of s_nodes (with enougth memory for the number of
 *  \ref s_node "s_nodes" present in the input).
 * \param input
 *  Node of the original group of nodes, to be copied into the output (root in
 *  the first call).
 *
 * \attention It requires that the group of s_nodes have their index in a post-order.
 *******************************************************************************/
static void PostCollapseSNodes(s_node *output,s_node *input);

/**
 * Copies a group of l_nodes (with tree structure) in an array of l_nodes following
 * a post-order.
 *
 * It clones a tree/subtree (group of \ref l_node "l_nodes" with tree structure) in an array
 * of \ref l_node "l_nodes". This is an important part of \ref CollapseLTree public function.
 * The l_node::g_nodes are also copied if there is allocated memory for them.
 *
 * \param output
 *  Pre-allocated memory of l_nodes (with enougth memory for the number of
 *  \ref l_node "l_nodes" present in the input.
 * \param input
 *  Node of the original group of nodes, to be copied into the output (root in
 *  the first call).
 * \param n_gleaves
 *  Number of leaves of the g_tree pointed by these \ref l_node "l_nodes". It is used to copy
 *  l_node::g_nodes.
 * \param retain_lateral
 *  Logical flag. If retain_lateral = 1, the \ref l_node::lat_node will not be set as NULL, retained an outdated value (DANGEROUS).
 *
 * \attention It requires that the group of l_nodes have their index in a post-order.
 *******************************************************************************/
static void PostCollapseLNodes(l_node *output,l_node *input,unsigned int n_gleaves, unsigned int retain_lateral);

/**
 * Copies a group of g_nodes (with tree structure) in an array of g_nodes following
 * a post-order.
 *
 * It clones a tree/subtree (group of \ref g_node "g_nodes" with tree structure) in an array
 * of \ref g_node "g_nodes". This is an important part of \ref CollapseGTree public function.
 *
 * \param output
 *  Pre-allocated memory of g_nodes (with enougth memory for the number of
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
 *  Pre-allocated memory of s_nodes (with enougth memory for the number of
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
 *  Pre-allocated memory of l_nodes (with enougth memory for the number of
 *  \ref l_node "l_nodes" present in the input.
 * \param input
 *  Node of the original group of nodes, to be copied into the output (root in
 *  the first call).
 * \param n_gleaves
 *  Number of leaves of the g_tree pointed by these \ref l_node "l_nodes". It is used to copy
 *  l_node::g_nodes.
 * \param retain_lateral
 *  Logical flag. If retain_lateral = 1, the \ref l_node::lat_node will not be set as NULL, retained an outdated value (DANGEROUS).
 *
 * \attention It requires that the group of l_nodes have their index in a pre-order.
 *******************************************************************************/
static void PreCollapseLNodes(l_node *output,l_node *input,unsigned int n_gleaves, unsigned int retain_lateral);

/**
 * Copies a group of g_nodes (with tree structure) in an array of g_nodes following
 * a pre-order.
 *
 * It clones a tree/subtree (group of \ref g_node "g_nodes" with tree structure) in an array
 * of \ref g_node "g_nodes". This is an important part of \ref CollapseGTree public function.
 *
 * \param output
 *  Pre-allocated memory of g_nodes (with enougth memory for the number of
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
static unsigned int Count_duplications_lnodes(l_node *node);

/**
 * Gets the number of losses of a bunch of l_nodes following a post-order recursion.
 *
 * \param node
 *  l_node to analyze.
 * \return
 *  Number of losses.
 *******************************************************************************/
static unsigned int Count_losses_lnodes(l_node *node);

/**
 * Gets the height of the tree in coalescent units (n_gen/Ne).
 *
 * \param node
 *  s_node to analyze.
 * \param g_Ne
 *  Global effective population size.
 * \return Tree height in coalescent units.
 *******************************************************************************/
static long double Measure_s_node_cu_height(s_node *node, unsigned int g_Ne);

/**
 * Gets the height of the tree in generations.
 *
 * \param node
 *  g_node to analyze.
 * \return Tree height in number of generations.
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
static long double Measure_s_node_cu_length(s_node *node, unsigned int g_Ne);

/**
 * Gets the length of the tree in expected number of substitutions per site.
 *
 * \param node
 *  g_node to analyze.
 * \return Tree length in units of expected number of substitutions per site.
 *******************************************************************************/
static long double Measure_s_node_gl_length(s_node *node);

/**
 * Gets the height of the tree in coalescent units (n_gen/Ne).
 *
 * \param node
 *  g_node to analyze.
 * \param g_Ne
 *  Global effective population size.
 * \return Tree height in coalescent units.
 *******************************************************************************/
static long double Measure_g_node_cu_height(g_node *node, unsigned int g_Ne);

/**
 * Gets the height of the tree in expected number of substitutions per site.
 *
 * \param node
 *  g_node to analyze.
 * \return Tree height in units of expected number of substitutions per site.
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
static long double Measure_g_node_cu_length(g_node *node, unsigned int g_Ne);

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
 * \attention This function assumes ultrametricity (s_node::n_gen is calculated
 *  using the s_node::gen_length of the first offspring node)
 *******************************************************************************/
static void RefineSNodes(s_node * node, unsigned int ind_persp);

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
 * \attention This function assumes ultrametricity (l_node::n_gen is calculated
 *  using the l_node::gen_length of the first offspring node)
 *******************************************************************************/
static void RefineLNodes(l_node * node, unsigned int n_gleaves, unsigned int ind_persp);

/**
 * Deletes superfluous nodes due to losses.
 *
 * This function deletes superfluous nodes due to losses of a \ref l_tree after
 * it generation by a birth-death process.
 *
 * \param node
 *  Node to delete if it is a superfluous one.
 * \param root
 *  Point to the root of the tree, to change if it is deleted.
 * \param n_deletions
 *  Number of deleted nodes after the work of the function.
 * \param n_leaves
 *  Number of leaves after the work of the function.
 *******************************************************************************/
static void CleanlossesLNodes(l_node * node, l_node ** root,unsigned int * n_deletions, unsigned int * n_leaves);

/**
 * Transforms the branch lenght of a gene node from number of generations to time recursively following a post-order.
 *
 * \param node
 *  Input g_node (root in the first recursion).
 * \param gen_time
 *  Generation time.
 * \return \ref NO_ERROR on OK or an \ref ERRORS "error code" if any error
 *  ocurrs.
 * \attention Only applicable to ultrametric trees.
 *  *******************************************************************************/
static void Temporalize_GNodes(g_node * node,double gen_time);

// ** Node I/O ** //

/**
 * Writes a given group of s_nodes with tree structure in Newick format in
 * stdout.
 *
 * \param root
 *  s_node to print (root in the first call).
 * \param names
 *  Names (name_c *).
 * \param gen_time
 *  Generation time (1 if everything is working in generations)
 *******************************************************************************/
static void WriteSNodes (s_node * root, name_c * names, double gen_time);

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
 *  Generation time (1 if everything is working in generations)
 * \note The FILE * should be previously opened and checked to avoid errors.
 *******************************************************************************/
static void WriteSNodesFile (FILE * file, s_node * root, name_c * names, double gen_time);

/**
 * Writes a given group of l_nodes with tree structure in Newick format in
 * stdout.
 *
 * \param root
 *  l_node to print (root in the first call).
 * \param names
 *  Names (name_c *).
 * \param gen_time
 *  Generation time (1 if everything is working in generations)
 *******************************************************************************/
static void WriteLNodes (l_node * root, name_c * names, double gen_time);

/**
 * Writes a given group of l_nodes with tree structure in Newick format in a
 * file.
 *
 * \param file
 *  Output opened file.
 * \param root
 *  l_node to print (root in the first call).
 * \param names
 *  Names (name_c *).
 * \param gen_time
 *  Generation time (1 if everything is working in generations)
 * \note The FILE * should be previously opened and checked to avoid errors.
 *******************************************************************************/
static void WriteLNodesFile (FILE * file, l_node * root, name_c * names, double gen_time);

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
s_node * NewSNodes(unsigned int n_nodes, unsigned int max_childs)
{
    s_node * nodes=NULL, *w_node=NULL;
    unsigned int i=0,j=0;
    
    // ***
    /// <dl><dt> Function structure </dt><dd>
    
    // **
    /// \ref s_node "S_node" memory allocation
    nodes=calloc(n_nodes,sizeof(s_node));
    ErrorReporter((int long)nodes);
    
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
        w_node->mu_mult=1;
        
        // *
        /// s_node::childs memory allocation and initialization</dd></dl></dd></dl>
        w_node->childs=calloc(max_childs,sizeof(s_node *));
        ErrorReporter((long int) w_node->childs);
        
        for (j=0; j<max_childs; ++j)
        {
            *(w_node->childs+j)=NULL;
        }
        
    }
    
    return (nodes);
}

// ** L_Node creation ** //
l_node * NewLNodes(unsigned int n_nodes, unsigned int n_gleaves, unsigned int max_childs)
{
    l_node * nodes=NULL, *w_node=NULL;
    unsigned int i=0,j=0;
    
    // ***
    /// <dl><dt> Function structure </dt><dd>
    
    // **
    /// \ref l_node "L_node" memory allocation
    nodes=calloc(n_nodes,sizeof(l_node));
    ErrorReporter((int long)nodes);
    
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
         w_node->paralog=0;*/ // Implicit in calloc
        w_node->anc_node=NULL;
        w_node->lat_node=NULL;
        w_node->conts=NULL;
        w_node->gen_length=0.0;
        w_node->n_gen=0.0;
        w_node->mu_mult=1;
        
        // *
        /// l_node::childs memory allocation and initialization
        w_node->childs=calloc(max_childs,sizeof(l_node *));
        ErrorReporter((long int) w_node->childs);
        
        for (j=0; j<max_childs; ++j)
        {
            *(w_node->childs+j)=NULL;
        }
        
        // *
        /// l_node::g_nodes allocation and initialization (if it is required by n_gleaves)</dd></dl></dd></dl>
        if (n_gleaves!=0)
        {
            w_node->g_nodes=calloc(n_gleaves,sizeof(g_node *));
            
            for (j=0; j<n_gleaves; ++j)
            {
                *(w_node->g_nodes+j)=NULL;
            }
        }
        else
        {
            w_node->g_nodes=NULL;
        }
        
        
    }
    
    return (nodes);
}

// ** G_Node creation ** //

g_node * NewGNodes(unsigned int n_nodes, unsigned int max_childs)
{
    g_node * nodes=NULL, *w_node=NULL;
    unsigned int i=0,j=0;
    
    // ***
    /// <dl><dt> Function structure </dt><dd>
    
    // **
    /// \ref g_node "G_node" memory allocation
    if (n_nodes==0)
        return (NULL);
    
    nodes=calloc(n_nodes,sizeof(g_node));
    ErrorReporter((int long)nodes);
    
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
        w_node->childs=NULL;
        w_node->anc_node=NULL;
        w_node->contl=NULL;
        w_node->conts=NULL;
        w_node->height=0.0;
        w_node->bl=0.0;
        w_node->gen_length=0.0;
        
        // *
        /// g_node::childs memory allocation and initialization</dd></dl>
        w_node->childs=calloc(max_childs,sizeof(g_node *));
        ErrorReporter((long int) w_node->childs);
        
        for (j=0; j<max_childs; ++j)
        {
            *(w_node->childs+j)=NULL;
        }
    }
    
    return (nodes);
}

period * NewPeriods(unsigned int n_periods, unsigned int max_nodes)
{
    period *periods=NULL, *w_period;
    unsigned int i=0;
    
    periods=calloc(n_periods, sizeof(struct period));
    
    for (i=0; i<n_periods; ++i)
    {
        w_period=periods+i;
        w_period->r_bound=0;
        w_period->l_nodes=calloc(max_nodes,sizeof(l_node *));
        w_period->n_lnodes=0;
    }
    return periods;
}


// *** Tree memory manage *** //

// ** Tree creation ** //

s_tree * NewSTree (unsigned int n_nodes, unsigned int n_leaves, unsigned int n_gleaves, unsigned int max_childs, double gen_time, unsigned int Ne, double mu)
{
    s_tree * tree=NULL;
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Tree memory allocation
    tree=calloc(1,sizeof(s_tree));
    ErrorReporter((long int)tree);
    
    // *
    /// Tree initialization
    tree->n_nodes=n_nodes;
    tree->n_leaves=n_leaves;
    tree->n_gleaves=n_gleaves;
    tree->max_childs=max_childs;
    tree->gen_time=gen_time;
    tree->Ne=Ne;
    tree->mu=mu;
    tree->locus_tree=NULL;
    tree->gene_tree=NULL;
    
    // *
    /// Tree nodes allocation and initialization by \ref NewSNodes </dd></dl>
    
    if (n_nodes>1)
    {
        tree->m_node=NewSNodes(n_nodes,max_childs);
        tree->root=NULL;
    }
    else if (n_nodes==1)
    {
        tree->root=NewSNodes(n_nodes, max_childs);
        tree->m_node=NULL;
    }
    else
    {
        tree->root=NULL;
        tree->m_node=NULL;
    }
    
    return (tree);
}

s_tree * ReadNewickSTree (char * newick,name_c **names_ptr, unsigned int verbosity, double gen_time, unsigned int Ne, double mu, unsigned int ind_persp)
{
    s_node *current_node=NULL, *anc_node=NULL,*root=NULL;
    s_tree *tree=NULL;
    name_c * names=NULL;
    char code=' ';
    unsigned int step=0;
    unsigned int n_char=0, iteration=0, ffree_codename=1, n_leaves=0, n_gleaves=0, n_inodes=0, n_nodes=0, max_childs=0, max_lname=0, n_replica=0,index=0, n_priv_ngleaves=0;
    char bl_buffer[DBL_DIG+3]=""; //DBL_DIG= precission of double, + 3 positions (0 separator and \0)
    char * ui_buffer;
    char name_buffer[MAX_NAME]="";
    
    // ****
    /// <dl><dt> Function structure </dt><dd>
    
    // ***
    /// Test of the newick tree string by \ref CheckNewickSTree
    
    ErrorReporter(CheckNewickSTree(newick));
    
    // ***
    /// Error control
    
    if (gen_time==0)
        ErrorReporter(UNEXPECTED_VALUE);
    
    // ***
    /// Integer buffer allocation and initialization
    
    ui_buffer=calloc((int)log10(UINT_MAX)+2,sizeof(char)); // An string of the max number of digits of an unsigned int + 1 (\0)
    strcpy(ui_buffer, "");
    
    // ***
    /// First read of the Newick tree. Obtains the number of internal nodes (")"), s_tree::n_leaves ("(" or "," not followed by "(") and s_tree::n_gleaves (s_tree::n_leaves + (replicas "/" -1 for each node)).
    while (*(newick+step)!=';')
    {
        if((*(newick+step)=='(' ||*(newick+step)==',' )&&*(newick+step+1)!='(') ++n_leaves;
        if(*(newick+step)==')') ++n_inodes;
        if(*(newick+step)=='/')
        {
            // * Reseting variables * //
            strcpy(ui_buffer,"");
            n_replica=0;
            n_char=0;
            ++step;
            code=*(newick+step);
            // * Reading replicas digits * //
            while (code<58 && code>47) // Int numbers
            {
                *(ui_buffer+n_char)=code;
                ++n_char;
                ++step;
                code=*(newick+step);
            }
            --step;
            // * Translating * //
            sscanf(ui_buffer,"%ui",&n_replica);
            n_gleaves+=n_replica;
            ++n_priv_ngleaves;
        }
        ++step;
        if (step==MAX_IT) ErrorReporter(LOOP_ERROR); // Avoids ininite loops
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
                    current_node=NewSNodes(1,MAX_CHILDS);
                    // *
                    /// Points the pointers of this new node and its ancestor
                    *(anc_node->childs+anc_node->n_child)=current_node;
                    ++anc_node->n_child;
                    current_node->anc_node=anc_node;
                    // *
                    /// Searches for the maximum number of childs in the tree.</dd></dl>
                    if(max_childs<anc_node->n_child)
                    {
                        max_childs=anc_node->n_child;
                    }
                    
                }
                else //New root node
                {
                    current_node=NewSNodes(1,MAX_CHILDS);
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
                strcpy(bl_buffer,"");
                n_char=0;
                ++step;
                code=*(newick+step);
                while (code<58 && (code>47 || code==46))
                {
                    *(bl_buffer+n_char)=code;
                    ++n_char;
                    ++step;
                    code=*(newick+step);
                }
                --step;
                *(bl_buffer+n_char)=0; //Sets the end of the string to avoid problems in the sscanf
                // *
                /// Translates the buffer in a float number and asigns it as s_node::gen_length converted in number of generations of the current node</dd></dl>
                sscanf(bl_buffer,"%lf",&current_node->gen_length);
                current_node->gen_length/=gen_time;
                
                break;
            case '*':
                // **
                /// <dl><dt>Lineage specific substitution rate multi (code="*").</dt><dd>
                
                // *
                /// Reads all the following integers and . in a buffer
                strcpy(bl_buffer,"");
                n_char=0;
                ++step;
                code=*(newick+step);
                while (code<58 && (code>47 || code==46))
                {
                    *(bl_buffer+n_char)=code;
                    ++n_char;
                    ++step;
                    code=*(newick+step);
                }
                --step;
                *(bl_buffer+n_char)=0; //Sets the end of the string to avoid problems in the sscanf
                // *
                /// Translates the buffer in a float number and asigns it as s_node::mu_mult </dd></dl>
                sscanf(bl_buffer,"%lf",&current_node->mu_mult);
                
                break;
            case '#':
                // **
                /// <dl><dt>New effective population size (Ne) (code="#").</dt><dd>
                
                // *
                /// Reads all the following integers in a buffer
                strcpy(ui_buffer,"");
                n_char=0;
                ++step;
                code=*(newick+step);
                while (code<58 && code>47)
                {
                    *(ui_buffer+n_char)=code;
                    ++n_char;
                    ++step;
                    code=*(newick+step);
                }
                --step;
                *(ui_buffer+n_char)=0; //Sets the end of the string to avoid problems in the sscanf
                // *
                /// Translates the buffer in a integer number and asigns it as s_node::Ne of the current node</dd></dl>
                sscanf(ui_buffer,"%ui",&current_node->Ne);
                if (current_node->Ne==0)
                    ErrorReporter(SETTINGS_ERROR);
                
                break;
            case '/':
                // **
                /// <dl><dt>New number of replicas (code="/").</dt><dd>
                
                // *
                /// Reads all the following integers in a buffer
                strcpy(ui_buffer,"");
                n_replica=0;
                n_char=0;
                ++step;
                code=*(newick+step);
                while (code<58 && code>47)
                {
                    *(ui_buffer+n_char)=code;
                    ++n_char;
                    ++step;
                    code=*(newick+step);
                }
                --step;
                *(ui_buffer+n_char)=0; //Sets the end of the string to avoid problems in the sscanf
                // *
                /// Translates the buffer in a integer number and asigns it as s_node::n_nodes of the current node</dd></dl>
                sscanf(ui_buffer,"%ui",&n_replica);
                if (n_replica==0)
                    ErrorReporter(SETTINGS_ERROR);
                current_node->n_replicas=n_replica;
                break;
                
            default:
                // **
                /// <dl><dt>New leaf(code!= former ones).</dt><dd>
                
                // *
                /// Reads the name of the leaf in a buffer
                
                strcpy(name_buffer,"");
                n_char=0;
                while (code!='(' && code!=')' && code!= ',' && code!= ';' && code!= ':')
                {
                    *(name_buffer+n_char)=code;
                    ++n_char;
                    ++step;
                    code=*(newick+step);
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
                current_node=NewSNodes(1,MAX_CHILDS);
                // *
                /// Points the pointers of this new node and its ancestor
                *(anc_node->childs+anc_node->n_child)=current_node;
                ++anc_node->n_child;
                current_node->anc_node=anc_node;
                // *
                /// Asigns the new name code (position) to the current node (s_node::sp_index)
                current_node->sp_index=ffree_codename;
                ++ffree_codename;
                // *
                /// Searches for the maximum number of childs in the tree.</dd></dl></dd></dl>
                if(max_childs<anc_node->n_child)
                {
                    max_childs=anc_node->n_child;
                }
                
                break;
        }
        
        ++step;
        code=*(newick+step);
        ++iteration;
        if (iteration==MAX_IT) ErrorReporter(LOOP_ERROR); // Avoids ininite loops
        
    }
    
    // ***
    /// Refines the readed nodes by \ref RefineSNodes
    
    RefineSNodes(root,ind_persp);
    max_lname++; //One extra character for \0.
    
    // ***
    /// Reallocates the name_c memory by \ref ReallocNames
    ReallocNames(names,max_lname);
    
    // ***
    /// Allocates memory for the readed tree and completes it.
    
    tree=NewSTree(0, n_leaves, n_gleaves, max_childs,gen_time,Ne,mu); //Nodes have already been allocated.
    tree->root=root;
    tree->n_nodes=n_nodes;
    tree->locus_tree=NULL;
    tree->gene_tree=NULL;
    
    // ***
    /// Fills the node indexes following a post_order.
    PostReorderSNodes(tree->root,&index);
    
    // ***
    /// Frees dynamic memory (buffers)</dd></dl>
    free(ui_buffer);
    ui_buffer=NULL;
    
    if (verbosity>2)
    {
        printf("\n\t\t %d-node species tree correctly built",(n_leaves*2)-1);
        if (verbosity>3)
        {
            printf(": ");
            WriteSNodes(root,names,gen_time);
        }
        printf("\n");
        
    }
    
    return (tree);
}

long int NewBDSTree (s_tree ** out_tree, unsigned int leaves, double time, double b_rate, double d_rate, double gen_time, unsigned int Ne, double mu, unsigned int ind_per_sp, double outgroup, unsigned int complete, unsigned int two_lineages, gsl_rng *seed, int verbosity)
{
    // ***** Declaration and initialization of variables ***** //
    
    unsigned int eq_rates=0, j=0, node_index=0, n_leaves=0, n_nodes=0, avail_leaves=0, extra_nodes=0, iter=0;
    //unsigned long int stats_leaves=0;
    //double stats_time=0;
    double w_prob=d_rate+b_rate, b_prob=(b_rate/(d_rate+b_rate)), random=0, inv_leaves=0, res1=0, res2=0, *i_nodes=NULL, s_time=0, t_time=0;
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
    /// Method election, BDSA or SSA algorithms.
    
    if (leaves>0)
    {
        // *****
        /// <dl><dt>BDSA (fixed number of leaves). Faster than GSA, but without information about extinct lineages.</dt><dd>
        if (verbosity>3)
        {
            printf("\n\t\tUsing BDSA algorithm (fixed number of leaves) to simulate the trees...");
            if (verbosity>4 && two_lineages==0)
                printf("\n\t\tTwo_lineages parameter disabled has no sense in BDSA context, and it will be ignored.");
#ifdef DEBUG
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
        ErrorReporter((long int)i_nodes);
        node_ptrs=calloc((*out_tree)->n_nodes,sizeof(s_node *));
        ErrorReporter((long int)node_ptrs);
        for (j=0;j<(*out_tree)->n_nodes;++j)
        {
            *(node_ptrs+j)=(*out_tree)->m_node+j;
        }
        
        n_leaves=leaves;
        
        // ****
        /// Calculates the start of the tree (From user options or sampled from the inverse of the pdf of the time of the origin of the tree (flat prior), conditional on having n species at the present, Hartmann et al., 2010; Gernhard, 2008.)
        if (time==0)
        {
            random=gsl_rng_uniform_pos(seed); //If it is 0 we get negative times, although in the paper they say [0,1]
            inv_leaves=1/(double)leaves;
            
            if (eq_rates>0)
                time=1/(b_rate*(pow(random,-inv_leaves)-1));
            else
            {
                res1=pow(random,inv_leaves);
                time=log((1-(res1*d_rate/b_rate))/(1-res1))/(b_rate-d_rate);
            }
            if (verbosity>4)
            {
                printf("\n\t\t\tRoot time (sampled using the inverse pdf of the origin of the tree conditional on having %d leaves (flat prior) = %.8lf",leaves,time);
#ifdef DEBUG
                fflush(stdout);
#endif
            }
            
        }
        else if (verbosity>4)
        {
            printf("\n\t\t\tRoot time (user defined) = %.8lf",time);
#ifdef DEBUG
            fflush(stdout);
#endif
        }
        //        if (verbosity>1)
        //        {
        //            stats_time+=(time/gen_time);
        //        }
        
        // ****
        /// Calculates the internal node branches, or sampled from the inverse of the pdf of the speciation events, Hartmann et al., 2010; Gernhard, 2008.
        
        for (j=0;j<leaves-2;++j)
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
        *(i_nodes+leaves-2)=time;
        
        // ****
        /// Orders the times in ascending order using a quicksort algorithm
        qsort(i_nodes, leaves-1, sizeof(*i_nodes), Compare_DBL);
        
        if (verbosity>4)
        {
            printf("\n\t\t\tReconstructing backwards the tree from sampled duplication times...");
            if (verbosity == 6)
            {
                printf("\n\t\t\t\tDuplication times:");
                for (j=0; j<leaves-1; ++j)
                {
                    printf(" %.8lf,",*(i_nodes+j));
                }
            }
#ifdef DEBUG
            fflush(stdout);
#endif
        }
        
        // ****
        /// Reconstructs the tree from the birth-death times
        for (j=0; j<leaves-1;++j)
        {
            // * First offspring node * //
            node_index= gsl_rng_uniform(seed)*n_leaves;//Samples nodes
            w_node1=*(node_ptrs+node_index);
            --n_leaves;
            *(node_ptrs + node_index) = *(node_ptrs + n_leaves);/* The pointer to this node now points to the last node_ptr, inaccesible due to --n_leaves*/
            
            // * Second offspring node * //
#ifndef NO_VAR
            if (n_leaves>1)
                node_index= gsl_rng_uniform(seed)*n_leaves;//Samples nodes
            else
                node_index=0;
#endif
#ifdef NO_VAR
            node_index=gsl_rng_uniform(seed)*n_leaves;
#endif
            w_node2= *(node_ptrs + node_index);
            *(node_ptrs + node_index)= (*out_tree)->m_node + j + leaves; //Now the node_ptr that pointed to w_node2 points to the new ancestral node, for using it in the next loop iteration
            
            // ** Ancestral node ** //
            anc_node= (*out_tree)->m_node + j + leaves; //Uses as anc node the next avaliable memory block for internal nodes
            
            // ** Filling working nodes ** //
            // * Pointers * //
            *(anc_node->childs)=w_node1;
            *(anc_node->childs+1)=w_node2;
            w_node1->anc_node=anc_node;
            w_node2->anc_node=anc_node;
            // * Branch lenghts and times * //
            anc_node->n_gen=*(i_nodes+j)/gen_time; //Time stored as number of generations.
            anc_node->n_child=2;
            w_node1->gen_length=anc_node->n_gen-w_node1->n_gen;
            w_node2->gen_length=anc_node->n_gen-w_node2->n_gen;
            
        }
        (*out_tree)->root=anc_node;
        
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
            w_node1=(*out_tree)->m_node+leaves*2;
            w_node2=(*out_tree)->m_node+leaves*2-1;
            w_node1->anc_node=w_node2;
            *w_node2->childs=(*out_tree)->root;
            *(w_node2->childs+1)=w_node1;
            w_node2->n_gen=(*out_tree)->root->n_gen+(outgroup*time/2);
            w_node2->n_child=2;
            (*out_tree)->root->gen_length=(outgroup*time/2);
            (*out_tree)->root->anc_node=w_node2;
            w_node1->gen_length=w_node2->n_gen;
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

        if (verbosity>4)
        {
            printf("\n\t\t\tDone");
#ifdef DEBUG
            fflush(stdout);
#endif
        }
        //        if (verbosity>4)
        //        {
        //            if (verbosity>5)
        //            {
        //                printf("\n\t\tSTATS: Time length of the tree %.8lf\n\t", stats_time);
        //            }
        //            printf("Done \n");
        //            #ifdef DEBUG
        //fflush(stdout);
//#endif
        //        }
        
        // ****
        /// Frees allocated memory</dd></dl>
        free (i_nodes);
        i_nodes=NULL;
        free (node_ptrs);
        node_ptrs=NULL;
    }
    else
    {
        
        // *****
        /// <dl><dt>SSA algorithm (GSA is not required with fixed time)</dt><dd>
        if (verbosity>3)
        {
            printf("\n\t\tUsing SSA algorithm (fixed generations) to simulate the species tree...");
#ifdef DEBUG
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
            ErrorReporter((long int)node_ptrs);
            (*out_tree)=NewSTree(1, 1, 0, 2,gen_time,Ne, mu); // Recursion oriented tree (root).
            
            
            // ****
            /// Iteration-dependent variable initialization
            
            // ****
            /// Start of the tree (one or two lineages)
            w_node1=NewSNodes(1, 2);
            *((*out_tree)->root->childs)=w_node1;//This root node is only a way to set the start of the tree. It is deleted at the end.
            w_node1->anc_node=(*out_tree)->root;
            anc_node=w_node1;
            extra_nodes=0;
            s_time=-log(gsl_rng_uniform_pos(seed))/(w_prob);
            t_time=time;
            
            if (two_lineages==0)
            {
                *(node_ptrs)=w_node1; // The actual root of the tree
                n_leaves=1;
                avail_leaves=1;
                n_nodes=1;
                (*out_tree)->root->n_gen=t_time/gen_time;
            }
            else
            {
                w_node1=NewSNodes(1,2);
                w_node2=NewSNodes(1,2);
                *(anc_node->childs)=w_node1;
                *(anc_node->childs+1)=w_node2;
                w_node1->anc_node=anc_node;
                w_node2->anc_node=anc_node;
                anc_node->n_gen=time/gen_time;
                anc_node->n_child=2;
                n_leaves=2;
                avail_leaves=2;
                n_nodes=3;
                (*out_tree)->root->n_gen=(time+s_time)/gen_time;
                s_time=-log(gsl_rng_uniform_pos(seed))/(w_prob*2);
                *(node_ptrs)=w_node1;
                *(node_ptrs+1)=w_node2;
                
            }
            
            // ****
            /// <dl><dt>Birth-death process loop (while final time is not exceeded)</dt><dd>
            while (t_time-s_time>0 && avail_leaves!=0)
            {
                // ***
                /// Memory error test (maximum number of leaves)
                if (n_leaves>=MAX_LEAVES)
                {
                    fprintf(stderr,"\nMaximum number of leaves reached. This limit is hard-coded as a constant to avoid memory problems, so if you have memory enough, you can recompile the software with a bigger limit\n");
                    return(MEM_ERROR);
                }
                
                // ***
                /// Chooses the leave for the next event
#ifndef NO_VAR
                if(avail_leaves>1)
                    node_index=gsl_rng_uniform(seed)*avail_leaves;
                else
                    node_index=0;
#endif
#ifdef NO_VAR
                node_index=gsl_rng_uniform(seed)*avail_leaves;
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
                    *(anc_node->childs)=w_node1;
                    *(anc_node->childs+1)=w_node2;
                    w_node1->anc_node=anc_node;
                    w_node2->anc_node=anc_node;
                    
                    anc_node->n_gen=(t_time-s_time)/gen_time;
                    anc_node->gen_length=anc_node->anc_node->n_gen-anc_node->n_gen;
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
#ifdef DEBUG
                        fflush(stdout);
#endif
                    }
                    else if (verbosity==6)
                    {
                        printf("\n\t\t\tNew duplication, time %.8lf",anc_node->n_gen);
#ifdef DEBUG
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
                        if (*(anc_node->childs) == w_node1)
                            w_node2=*(anc_node->childs+1);
                        else
                            w_node2=*(anc_node->childs);
                        
                        if (complete==0)
                        {
                            // * Reconfiguring pointers * //
                            w_node2->anc_node=anc_node->anc_node;
                            if (*(w_node2->anc_node->childs)==anc_node)
                                *(w_node2->anc_node->childs)=w_node2;
                            else
                                *(w_node2->anc_node->childs+1)=w_node2;
                            
                            // * Reconfiguring branches * //
                            if (w_node2->n_gen!=0)
                                w_node2->gen_length=w_node2->anc_node->n_gen-w_node2->n_gen;
                            
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
                        else
                        {
                            // * Reconfiguring branches * //
                            w_node1->n_gen=(t_time-s_time)/gen_time;
                            w_node1->gen_length=anc_node->n_gen-w_node1->n_gen;
                            
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
#ifdef DEBUG
                        fflush(stdout);
#endif
                    }
                    else if (verbosity==6)
                    {
                        printf("\n\t\t\tNew loss, time %.8lf",(t_time-s_time)/gen_time);
#ifdef DEBUG
                        fflush(stdout);
#endif
                    }
                    
                }
                
                // ***
                /// New step time </dl>
                t_time-=s_time;
                s_time=-log(gsl_rng_uniform_pos(seed))/(w_prob*avail_leaves);
                
            }
            // ****
            /// If the tree ends without leaves, it is discarted and the algorithm retries it construction.
            if (avail_leaves<2)
            {
                FreeSTree(out_tree);
                if (verbosity>2)
                {
                    printf("\n\t\t\tSpecies tree with less than 2 leaves, restart of the simulation of this tree. Try %d of %d",iter,MAX_IT);
                    if (verbosity>4)
                        printf("\n");
#ifdef DEBUG
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
            w_node1->n_gen=0;
            w_node1->gen_length=w_node1->anc_node->n_gen;
            w_node1->n_replicas=ind_per_sp;
            w_node1->sp_index=j+n_leaves+1-avail_leaves;
        }
        
        // ****
        /// Deletion of the dummie root of the tree if there is no outgroup addition
        if (outgroup>0)
        {
            outgroup=(*(*out_tree)->root->childs)->n_gen*outgroup/2; //From a deviation to the half of the tree heigth to a real internal branch length
            
            w_node1=NewSNodes(1,2);
            *((*out_tree)->root->childs+1)=w_node1;
            w_node1->anc_node=(*out_tree)->root;
            (*out_tree)->root->n_gen=(*(*out_tree)->root->childs)->n_gen+outgroup;
            (*out_tree)->root->n_child=2;
            (*(*out_tree)->root->childs)->gen_length=outgroup;
            w_node1->gen_length=(*(*out_tree)->root->childs)->n_gen+outgroup;
            w_node1->n_replicas=1;
            
  
            (*out_tree)->n_nodes=n_nodes+extra_nodes+2;
            (*out_tree)->n_leaves=avail_leaves+(extra_nodes/2)+1;
            (*out_tree)->n_gleaves=(avail_leaves+(extra_nodes/2))*ind_per_sp+1;
            (*out_tree)->locus_tree=NULL;
            (*out_tree)->gene_tree=NULL;

        }
        else
        {
            w_node1=*((*out_tree)->root->childs);
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
#ifdef DEBUG
                fflush(stdout);
#endif
            }
            return (LOOP_ERROR);
        }
    }
    if (verbosity>2)
    {
        printf("\n\t\t%d-node species tree correctly simulated",(*out_tree)->n_nodes);
        if (verbosity>3)
        {
            printf(": ");
            WriteSTree((*out_tree),(name_c *)NULL,(*out_tree)->gen_time!=1?1:0);
        }
        printf("\n");
#ifdef DEBUG
        fflush(stdout);
#endif
    }
    // ******
    /// Return
    return (NO_ERROR);
}

long int SimBDLTree(s_tree *wsp_tree,l_tree **wlocus_tree, l_node **node_ptrs, double b_rate,double d_rate,gsl_rng *seed, unsigned int min_lleaves, unsigned int min_lsleaves, double gen_time, unsigned int verbosity, unsigned int *st_losses, unsigned int *st_dups, unsigned int *st_leaves, unsigned int *st_gleaves)
{
    
    // *******
    /// <dl><dt>Declarations</dt><dd>
    
    // ******
    ///Structures
    l_node dummy_node,*aux_lnode=NULL,*g_lnode=NULL,*anc_lnode=NULL,*w_lnode=NULL,*w_lnode2=NULL;
    s_node *w_snode=NULL;
    
    // ******
    /// B-D process related variables
    unsigned int n_leaves=0,lt_true_leaves=0,lt_diffs_true_leaves=0,diffs_true_leaves=0,tn_nodes=0,extra_nodes=0,n_losses=0;
    unsigned int next_paralog=0,node_index=0,avail_leaves=0;
    int n_nodes=0;
    double w_prob=d_rate+b_rate,b_prob=(b_rate/(d_rate+b_rate)),current_ngen=0,max_ngen=0,sampled_ngen=0;
    
    // ******
    /// Loop related variables</dd></dl>
    unsigned int i=0,j=0,k=0,ltree_iter=0,bd_iter=0,l_tree_retries=0, maxnleaves_reached=0,done=0;
    
    while (done!=1 && ltree_iter<=MAX_IT)
    {
        // ******
        /// Locus tree re/initialization
        if (*wlocus_tree!=NULL)
            FreeLTree(wlocus_tree);
        *wlocus_tree=NewLTree(1, 0, 0, wsp_tree->max_childs,wsp_tree->gen_time,wsp_tree->Ne,wsp_tree->mu); //Tree with only a root node and without g_node pointers
        
        // ******
        /// Species tree reset by \ref ResetSTreeSimL
        
        if(ResetSTreeSimL(wsp_tree)!= NO_ERROR)
            return MEM_ERROR;
        
        // ******
        /// Reseting of B-D process variables and linking the root of both trees
        wsp_tree->root->l_nodes=(*wlocus_tree)->root;
        wsp_tree->root->n_lnodes=1;
        (*wlocus_tree)->root->n_gen=wsp_tree->root->n_gen;
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
#ifdef DEBUG
                fflush(stdout);
#endif
            }
            
            w_snode=wsp_tree->m_node+i;
            
            if (w_snode->anc_node->n_lnodes==0) //There is no l_nodes in the ancestor of this s_node (lost lineage)
                continue;
            
            dummy_node.lat_node=w_snode->anc_node->l_nodes; //dummy node to allow the loop to start in the first iteration, due to the use of different structures (first iteration: s_node->l_nodes, rest: l_node->lat_node)
            g_lnode=&dummy_node;
            diffs_true_leaves=0;
            
            // *****
            /// <dl><dt>Loop of associated \ref s_node::l_nodes "l_nodes" to the working s_node, i.e paralogs</dt><dd>
            
            
            for (j=0;j<w_snode->anc_node->n_lnodes;++j)
            {
                if (verbosity>4)
                {
                    printf("\n\t\t\t Locus tree node %u of %u...",j+1,w_snode->anc_node->n_lnodes);
#ifdef DEBUG
                    fflush(stdout);
#endif
                }
                
                // ****
                /// Waiting-time sampling
                g_lnode=g_lnode->lat_node;
                anc_lnode=g_lnode;
                
                current_ngen=w_snode->anc_node->n_gen;
                max_ngen=w_snode->n_gen;
                sampled_ngen=-log(gsl_rng_uniform_pos(seed))/(w_prob);
                
                if (sampled_ngen>=current_ngen-max_ngen) //No events
                {
                    if (verbosity>4)
                    {
                        printf("\n\t\t\tThere is neither birth nor death event");
#ifdef DEBUG
                        fflush(stdout);
#endif
                    }
                    
                    // ****
                    /// New l-node allocation and configuration if there is neither birth nor death events along this branch
                    
                    w_lnode=NewLNodes(1, 0, wsp_tree->max_childs); //New node
                    //Locus tree pointers
                    w_lnode->anc_node=anc_lnode;
                    *(anc_lnode->childs+anc_lnode->n_child)=w_lnode;
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
                        printf("\n\t\t\tBirth-death process");
#ifdef DEBUG
                        fflush(stdout);
#endif
                    }
                    
                    // ****
                    /// <dl><dt>Else, Birth-death process along this branch</dt><dd>
                    
                    // ****
                    /// Iteration-dependent variable initialization
                    
                    w_lnode=NewLNodes(1,0,wsp_tree->max_childs); //New node
                    
                    //Pointers
                    w_lnode->anc_node=anc_lnode;
                    *(anc_lnode->childs+anc_lnode->n_child)=w_lnode;
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
                    while(current_ngen-sampled_ngen>max_ngen && avail_leaves!=0 && bd_iter<=MAX_IT && bd_iter<MAX_LEAVES-1)
                    {
                        ++bd_iter;
                        // ***
                        /// Chooses the leave for the next event
                        if(avail_leaves>1)
                            node_index=gsl_rng_uniform(seed)*avail_leaves;
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
                            w_lnode=NewLNodes(1,0,wsp_tree->max_childs);
                            w_lnode2=NewLNodes(1,0,wsp_tree->max_childs);
                            
                            //Tree pointers and branch lengths reconfiguration
                            anc_lnode=*(node_ptrs+node_index);
                            *(anc_lnode->childs)=w_lnode;
                            *(anc_lnode->childs+1)=w_lnode2;
                            w_lnode->anc_node=anc_lnode;
                            w_lnode2->anc_node=anc_lnode;
                            
                            anc_lnode->n_gen=current_ngen-sampled_ngen;
                            anc_lnode->gen_length=anc_lnode->anc_node->n_gen-anc_lnode->n_gen;
                            
                            //Ancestor node info
                            anc_lnode->kind_node=DUP; //Duplication
                            anc_lnode->n_child=2;
                            anc_lnode->sp_index=w_snode->sp_index;
                            anc_lnode->Ne=w_snode->Ne;
                            anc_lnode->mu_mult=w_snode->mu_mult;
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
#ifdef DEBUG
                                fflush(stdout);
#endif
                            }
                            else if (verbosity==6)
                            {
                                printf("\n\t\t\tNew duplication, time %lf",anc_lnode->n_gen*gen_time);
#ifdef DEBUG
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
                            w_lnode->n_gen=current_ngen-sampled_ngen;
                            w_lnode->gen_length=anc_lnode->n_gen-w_lnode->n_gen;
                            
                            //Setting the leaf as non avaliable
                            *(node_ptrs+node_index)=*(node_ptrs+avail_leaves-1);
                            --avail_leaves;
                            --n_leaves;
                            
                            w_lnode->sp_index=w_snode->sp_index;
                            w_lnode->kind_node=LOSS;
                            w_lnode->Ne=w_snode->Ne;
                            w_lnode->mu_mult=w_snode->mu_mult;
                            w_lnode->n_nodes=0;
                            w_lnode->conts=w_snode;
                            n_nodes-=2; //2 true nodes lost (w_lnode an its ancestor) becoming extra_nodes (they have allocated memory, and are part of the l_tree, but they will not be represented in the g_tree)
                            extra_nodes+=2;
                            ++n_losses;
                            ++(*st_losses);
                            
                            if (verbosity==5)
                            {
                                printf("\n\t\t\tNew loss");
#ifdef DEBUG
                                fflush(stdout);
#endif
                            }
                            else if (verbosity==6)
                            {
                                printf("\n\t\t\tNew loss, time %lf",w_lnode->n_gen*gen_time);
#ifdef DEBUG
                                fflush(stdout);
#endif
                            }
                            
                            
                        }
                        
                        // ***
                        /// New step time </dd></dl>
                        current_ngen-=sampled_ngen;
                        sampled_ngen=-log(gsl_rng_uniform_pos(seed))/(w_prob*avail_leaves);
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
                            w_lnode->gen_length=w_lnode->anc_node->n_gen-w_lnode->n_gen;
                            w_lnode->sp_index=w_snode->sp_index;
                            w_lnode->Ne=w_snode->Ne;
                            w_lnode->mu_mult=w_snode->mu_mult;
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
                        w_lnode->gen_length=w_lnode->anc_node->n_gen-w_lnode->n_gen;
                        w_lnode->sp_index=w_snode->sp_index;
                        w_lnode->Ne=w_snode->Ne;
                        w_lnode->mu_mult=w_snode->mu_mult;
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
#ifdef DEBUG
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
                    printf("\n\tThe locus tree simulation reached the maximum number of locus tree lineages (%d) inside a species tree branch. Try %d of %d",MAX_LEAVES,l_tree_retries,MAX_IT);
#ifdef DEBUG
                    fflush(stdout);
#endif
                }
                else if (lt_true_leaves<min_lleaves)
                {
                    printf("\n\tLocus tree with %u leaves, less than the minimum %u , restart of the simulation of this tree. Try %d of %d",lt_true_leaves,min_lleaves,l_tree_retries,MAX_IT);
#ifdef DEBUG
                    fflush(stdout);
#endif
                }
                else
                {
                    printf("\n\tLocus tree with %u leaves from different species, less than the minimum %u , restart of the simulation of this tree. Try %d of %d",lt_diffs_true_leaves,min_lsleaves,l_tree_retries,MAX_IT);
#ifdef DEBUG
                    fflush(stdout);
#endif
                }
            }
            
            
            ++l_tree_retries;
        }
        else
        {
            done=1;
        }
    }
    
    
    if (ltree_iter>MAX_IT)
    {
        ErrorReporter(LOOP_ERROR);
        return (NO_ERROR);
    }
    else
    {
        // ******
        /// Filling information of the just simulated l_tree
        (*wlocus_tree)->n_leaves=*st_leaves;
        (*wlocus_tree)->n_nodes=tn_nodes;
        (*wlocus_tree)->n_gleaves=*st_gleaves;
        (*wlocus_tree)->species_tree=wsp_tree;
        wsp_tree->locus_tree=*wlocus_tree;
        return (NO_ERROR);
    }
    
}

long int SimBDLHTree(s_tree *wsp_tree,l_tree **wlocus_tree, l_node **node_ptrs, double b_rate,double d_rate, double h_rate, gsl_rng *seed, unsigned int min_lleaves, unsigned int min_lsleaves, double gen_time, unsigned int verbosity, unsigned int *st_losses, unsigned int *st_dups, unsigned int *st_transfr, unsigned int *st_leaves, unsigned int *st_gleaves)
{
    
    //////////////////////////////DEBUG/////////////////////////////////////////
    
    unsigned int dist_dependent=1;
    ////////////////////////////////////////////////////////////////////////////
    
    
    // *******
    /// <dl><dt>Declarations</dt><dd>
    
    // ******
    ///Structures
    l_node dummy_node,*aux_lnode=NULL,*g_lnode=NULL,*anc_lnode=NULL,*w_lnode=NULL,*w_lnode2=NULL;
    s_node *w_snode=NULL;
    period *periods=NULL, *w_period=NULL, *w_period2=NULL;
    
    // ******
    /// B-D process related variables
    unsigned int n_leaves=0,lt_true_leaves=0,lt_diffs_true_leaves=0,diffs_true_leaves=0,tn_nodes=0,extra_nodes=0,n_losses=0, n_periods=0;
    unsigned int next_paralog=0,node_index=0,avail_leaves=0,n_transfer=0;
    int n_nodes=0;
    double w_prob=d_rate+b_rate+h_rate,b_prob=(b_rate/(d_rate+b_rate+h_rate)), bd_prob=(b_rate+d_rate/(b_rate+d_rate+h_rate)),current_ngen=0,max_ngen=0,sampled_ngen=0, rnumber=0;
    
    // ******
    /// Loop related variables</dd></dl>
    unsigned int ltree_iter=0,bd_iter=0,l_tree_retries=0, maxnleaves_reached=0,done=0;
    int i=0,j=0,k=0;
    
    while (done!=1 && ltree_iter<=MAX_IT)
    {
        // ******
        /// Locus tree re/initialization
        if (*wlocus_tree!=NULL)
            FreeLTree(wlocus_tree);
        *wlocus_tree=NewLTree(1, 0, 0, wsp_tree->max_childs,wsp_tree->gen_time,wsp_tree->Ne,wsp_tree->mu); //Tree with only a root node and without g_node pointers
        
        // ******
        /// Species tree reset by \ref ResetSTreeSimL
        
        if(ResetSTreeSimL(wsp_tree)!= NO_ERROR)
            return MEM_ERROR;
        
        // ******
        /// Reseting of B-D process variables and linking the root of both trees
        wsp_tree->root->l_nodes=(*wlocus_tree)->root;
        wsp_tree->root->n_lnodes=1;
        (*wlocus_tree)->root->n_gen=wsp_tree->root->n_gen;
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
        *st_transfr=0;
        ++ltree_iter;
        maxnleaves_reached=0;
        
        // ******
        /// <dl><dt>Pre-order s_node loop, to do a SSA (Simple sampling approach) algorithm along each branch of the species tree (GSA (General Sampling approach) is not required with fixed time)</dt><dd>
        for (i=1;i<wsp_tree->n_nodes;++i) //The root is not iterated.
        {
            
            if (verbosity>4)
            {
                printf("\n\t\tBranch %u of %u...",i,wsp_tree->n_nodes-1);
#ifdef DEBUG
                fflush(stdout);
#endif
            }
            
            w_snode=wsp_tree->m_node+i;
            
            if (w_snode->anc_node->n_lnodes==0) //There is no l_nodes in the ancestor of this s_node (lost lineage)
                continue;
            
            dummy_node.lat_node=w_snode->anc_node->l_nodes; //dummy node to allow the loop to start in the first iteration, due to the use of different structures (first iteration: s_node->l_nodes, rest: l_node->lat_node)
            g_lnode=&dummy_node;
            diffs_true_leaves=0;
            
            // *****
            /// <dl><dt>Loop of associated \ref s_node::l_nodes "l_nodes" to the working s_node, i.e paralogs</dt><dd>
            
            
            for (j=0;j<w_snode->anc_node->n_lnodes;++j)
            {
                if (verbosity>4)
                {
                    printf("\n\t\t\t Locus tree node %u of %u...",j+1,w_snode->anc_node->n_lnodes);
#ifdef DEBUG
                    fflush(stdout);
#endif
                }
                
                // ****
                /// Waiting-time sampling
                g_lnode=g_lnode->lat_node;
                anc_lnode=g_lnode;
                
                current_ngen=w_snode->anc_node->n_gen;
                max_ngen=w_snode->n_gen;
                sampled_ngen=-log(gsl_rng_uniform_pos(seed))/(w_prob);
                
                if (sampled_ngen>=current_ngen-max_ngen) //No events
                {
                    if (verbosity>4)
                    {
                        printf("\n\t\t\tThere is neither birth nor death event");
#ifdef DEBUG
                        fflush(stdout);
#endif
                    }
                    
                    // ****
                    /// New l-node allocation and configuration if there is neither birth nor death events along this branch
                    
                    w_lnode=NewLNodes(1, 0, wsp_tree->max_childs); //New node
                    //Locus tree pointers
                    w_lnode->anc_node=anc_lnode;
                    *(anc_lnode->childs+anc_lnode->n_child)=w_lnode;
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
                        printf("\n\t\t\tBirth-death process");
#ifdef DEBUG
                        fflush(stdout);
#endif
                    }
                    
                    // ****
                    /// <dl><dt>Else, Birth-death process along this branch</dt><dd>
                    
                    // ****
                    /// Iteration-dependent variable initialization
                    
                    w_lnode=NewLNodes(1,0,wsp_tree->max_childs); //New node
                    
                    //Pointers
                    w_lnode->anc_node=anc_lnode;
                    *(anc_lnode->childs+anc_lnode->n_child)=w_lnode;
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
                    while(current_ngen-sampled_ngen>max_ngen && avail_leaves!=0 && bd_iter<=MAX_IT && bd_iter<MAX_LEAVES-1)
                    {
                        ++bd_iter;
                        // ***
                        /// Chooses the leave for the next event
                        if(avail_leaves>1)
                            node_index=gsl_rng_uniform(seed)*avail_leaves;
                        else
                            node_index=0;
                        
                        // ***
                        /// <dl><dt>Chooses the process</dt><dd>
                        rnumber=gsl_rng_uniform(seed);
                        if (rnumber>bd_prob)
                        {
                            // **
                            /// <dl><dt>Transfer (donnor)</dt><dd>
                            
                            // *
                            /// New nodes allocation and reconfiguration of tree pointers/info</dd></dl>
                            w_lnode=NewLNodes(1,0,wsp_tree->max_childs);
                            
                            //Tree pointers and branch lengths reconfiguration
                            anc_lnode=*(node_ptrs+node_index);
                            *(anc_lnode->childs)=w_lnode;
                            w_lnode->anc_node=anc_lnode;
                            
                            anc_lnode->n_gen=current_ngen-sampled_ngen;
                            anc_lnode->gen_length=anc_lnode->anc_node->n_gen-anc_lnode->n_gen;
                            
                            //Ancestor node info
                            anc_lnode->kind_node=TRFR; //Transfer (donnor)
                            anc_lnode->n_child=1;
                            anc_lnode->sp_index=w_snode->sp_index;
                            anc_lnode->Ne=w_snode->Ne;
                            anc_lnode->mu_mult=w_snode->mu_mult;
                            anc_lnode->conts=w_snode;
                            
                            //New nodes info
                            w_lnode->paralog=anc_lnode->paralog;
                            
                            //New leaves addition
                            *(node_ptrs+node_index)=w_lnode;
                            ++n_nodes;
                            ++(*st_transfr);
                            
                            if (verbosity==5)
                            {
                                printf("\n\t\t\tNew transfer donnor");
#ifdef DEBUG
                                fflush(stdout);
#endif
                            }
                            else if (verbosity==6)
                            {
                                printf("\n\t\t\tNew transfer donnor, time %lf",anc_lnode->n_gen*gen_time);
#ifdef DEBUG
                                fflush(stdout);
#endif
                                
                            }
                        }
                        else if (rnumber<b_prob)
                        {
                            // **
                            /// <dl><dt>Birth</dt><dd>
                            
                            // *
                            /// New nodes allocation and reconfiguration of tree pointers/info</dd></dl>
                            w_lnode=NewLNodes(1,0,wsp_tree->max_childs);
                            w_lnode2=NewLNodes(1,0,wsp_tree->max_childs);
                            
                            //Tree pointers and branch lengths reconfiguration
                            anc_lnode=*(node_ptrs+node_index);
                            *(anc_lnode->childs)=w_lnode;
                            *(anc_lnode->childs+1)=w_lnode2;
                            w_lnode->anc_node=anc_lnode;
                            w_lnode2->anc_node=anc_lnode;
                            
                            anc_lnode->n_gen=current_ngen-sampled_ngen;
                            anc_lnode->gen_length=anc_lnode->anc_node->n_gen-anc_lnode->n_gen;
                            
                            //Ancestor node info
                            anc_lnode->kind_node=DUP; //Duplication
                            anc_lnode->n_child=2;
                            anc_lnode->sp_index=w_snode->sp_index;
                            anc_lnode->Ne=w_snode->Ne;
                            anc_lnode->mu_mult=w_snode->mu_mult;
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
#ifdef DEBUG
                                fflush(stdout);
#endif
                            }
                            else if (verbosity==6)
                            {
                                printf("\n\t\t\tNew duplication, time %lf",anc_lnode->n_gen*gen_time);
#ifdef DEBUG
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
                            w_lnode->n_gen=current_ngen-sampled_ngen;
                            w_lnode->gen_length=anc_lnode->n_gen-w_lnode->n_gen;
                            
                            //Setting the leaf as non avaliable
                            *(node_ptrs+node_index)=*(node_ptrs+avail_leaves-1);
                            --avail_leaves;
                            --n_leaves;
                            
                            w_lnode->sp_index=w_snode->sp_index;
                            w_lnode->kind_node=LOSS;
                            w_lnode->Ne=w_snode->Ne;
                            w_lnode->mu_mult=w_snode->mu_mult;
                            w_lnode->n_nodes=0;
                            w_lnode->conts=w_snode;
                            n_nodes-=2; //2 true nodes lost (w_lnode an its ancestor) becoming extra_nodes (they have allocated memory, and are part of the l_tree, but they will not be represented in the g_tree)
                            extra_nodes+=2;
                            ++n_losses;
                            ++(*st_losses);
                            
                            if (verbosity==5)
                            {
                                printf("\n\t\t\tNew loss");
#ifdef DEBUG
                                fflush(stdout);
#endif
                            }
                            else if (verbosity==6)
                            {
                                printf("\n\t\t\tNew loss, time %lf",w_lnode->n_gen*gen_time);
#ifdef DEBUG
                                fflush(stdout);
#endif
                            }
                            
                            
                        }
                        
                        // ***
                        /// New step time </dd></dl>
                        current_ngen-=sampled_ngen;
                        sampled_ngen=-log(gsl_rng_uniform_pos(seed))/(w_prob*avail_leaves);
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
                            w_lnode->gen_length=w_lnode->anc_node->n_gen-w_lnode->n_gen;
                            w_lnode->sp_index=w_snode->sp_index;
                            w_lnode->Ne=w_snode->Ne;
                            w_lnode->mu_mult=w_snode->mu_mult;
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
                        w_lnode->gen_length=w_lnode->anc_node->n_gen-w_lnode->n_gen;
                        w_lnode->sp_index=w_snode->sp_index;
                        w_lnode->Ne=w_snode->Ne;
                        w_lnode->mu_mult=w_snode->mu_mult;
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
#ifdef DEBUG
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
                    printf("\n\tThe locus tree simulation reached the maximum number of locus tree lineages (%d) inside a species tree branch. Try %d of %d",MAX_LEAVES,l_tree_retries,MAX_IT);
#ifdef DEBUG
                    fflush(stdout);
#endif
                }
                else if (lt_true_leaves<min_lleaves)
                {
                    printf("\n\tLocus tree with %u leaves, less than the minimum %u , restart of the simulation of this tree. Try %d of %d",lt_true_leaves,min_lleaves,l_tree_retries,MAX_IT);
#ifdef DEBUG
                    fflush(stdout);
#endif
                }
                else
                {
                    printf("\n\tLocus tree with %u leaves from different species, less than the minimum %u , restart of the simulation of this tree. Try %d of %d",lt_diffs_true_leaves,min_lsleaves,l_tree_retries,MAX_IT);
#ifdef DEBUG
                    fflush(stdout);
#endif
                }
            }
            
            
            ++l_tree_retries;
        }
        else
        {
            
            done=1;
        }
    }
    
    
    if (ltree_iter>MAX_IT)
    {
        ErrorReporter(LOOP_ERROR);
        return (NO_ERROR);
    }
    else
    {
        if (*st_transfr>0)
        {
            if (verbosity>4)
            {
                printf("\n\n\t\tCompleting transfers... ");
#ifdef DEBUG
                fflush(stdout);
#endif
            }
            (*wlocus_tree)->n_leaves=*st_leaves;
            (*wlocus_tree)->n_nodes=tn_nodes+*st_transfr;//One extra node will be needed for each transfer
            (*wlocus_tree)->n_gleaves=*st_gleaves;
            (*wlocus_tree)->species_tree=wsp_tree;
            ErrorReporter(CollapseLTree(*wlocus_tree,1,0));
            
            n_periods=tn_nodes-*st_leaves+*st_losses+1; //Maximum number of periods of an ultrametric tree (without the root) =Internal nodes + losses (tip_dates). I add a dummy one, with r_bound==0 to avoid some pointer problems
            
            //Here I'm not using NewPeriods to perform the initialization in a more efficient way, saving one extra for loop
            periods=calloc(n_periods, sizeof(struct period));
            periods->r_bound=0;
            periods->l_nodes=NULL;
            periods->n_lnodes=0;
            j=1;
            
            if (verbosity>4)
            {
                printf("\n\t\t\tInitializing periods... ");
#ifdef DEBUG
                fflush(stdout);
#endif
            }
            for (i=0; i<tn_nodes; ++i)
            {
                w_lnode=((*wlocus_tree)->m_node+i);
                if (w_lnode->n_gen==0)
                    continue;
                w_period=periods+j;
                w_period->l_nodes=calloc(*st_leaves, sizeof(l_node *)); //The maximum number of lineages in a period is the number of leaves (star tree)
                w_period->n_lnodes=w_lnode->n_child;
                for (k=0; k<w_lnode->n_child; ++k)
                {
                    *(w_period->l_nodes+k)=*(w_lnode->childs+k);
                }
                w_period->r_bound=w_lnode->n_gen;
                ++j;
            }
            if (verbosity>4)
            {
                printf("Done");
#ifdef DEBUG
                fflush(stdout);
#endif
            }
            
            if (verbosity>4)
            {
                printf("\n\t\t\tOrdering periods... ");
#ifdef DEBUG
                fflush(stdout);
#endif
            }
            
            qsort(periods, n_periods, sizeof(struct period), Compare_periods);
            
            if (verbosity>4)
            {
                printf("Done");
#ifdef DEBUG
                fflush(stdout);
#endif
            }
            
            if (verbosity>4)
            {
                printf("\n\t\t\tTraversing periods and performing transferences... ");
#ifdef DEBUG
                fflush(stdout);
#endif
            }
            for (i=n_periods-1; i>=0; --i)
            {
                w_period=periods+i;
                w_period2=periods+i-1;
                if (verbosity>4)
                {
                    printf("\n\t\t\t\tPeriod %d, lower bound %lf, n_lineages %d",i, w_period->r_bound, w_period->n_lnodes);
#ifdef DEBUG
                    fflush(stdout);
#endif
                }
                for (j=0; j<w_period->n_lnodes; ++j)
                {
                    w_lnode=*(w_period->l_nodes+j);
                    if (w_lnode->n_gen<w_period2->r_bound)
                    {
                        if (verbosity>5)
                        {
                            printf("\n\t\t\t\t\tCopying node %d to the previous period",w_lnode->index);
#ifdef DEBUG
                            fflush(stdout);
#endif
                        }
                        *(w_period2->l_nodes+w_period2->n_lnodes)=w_lnode;
                        w_period2->n_lnodes+=1;
                    }
                    else if (w_lnode->kind_node==TRFR && w_lnode->n_gen==w_period2->r_bound)
                    {
                        //Perform the transfer
                        if (w_period->n_lnodes>2)
                        {
                            if (dist_dependent!=0)
                                node_index=ChooseLNodePeriod(w_period,w_lnode,gsl_rng_uniform_pos(seed));
                            else
                            {
                                node_index=gsl_rng_uniform_pos(seed)*(w_period->n_lnodes-2);
                                if (node_index>=j)
                                    node_index++;
                            }
                        }
                        else if (j==0)
                            node_index=1;
                        else
                            node_index=0;


                        
                        //w_lnode -> donnor
                        //w_lnode2 -> receptor after reception
                        //g_lnode -> receptor (lost)
                        
                        w_lnode2=*(w_period->l_nodes+node_index);
                        
                        g_lnode=(*wlocus_tree)->m_node+tn_nodes+n_transfer; //This will be the receptor
                        anc_lnode=w_lnode2->anc_node;
                        g_lnode->anc_node=anc_lnode;
                        g_lnode->n_gen=w_lnode->n_gen;
                        g_lnode->gen_length=g_lnode->anc_node->n_gen-g_lnode->n_gen;
                        g_lnode->kind_node=RTRFR;
                        g_lnode->paralog=w_lnode2->paralog;
                        g_lnode->conts=w_lnode2->conts;
                        g_lnode->mu_mult=w_lnode2->mu_mult;
                        g_lnode->Ne=w_lnode2->Ne;
                        
                        for (k=0; k<anc_lnode->n_child; ++k)
                        {
                            if (*(anc_lnode->childs+k)==w_lnode2)
                                *(anc_lnode->childs+k)=g_lnode;
                        }
                        
                        w_lnode2->gen_length=w_lnode->n_gen-w_lnode2->n_gen;
                        w_lnode2->anc_node=w_lnode;
                        
                        
                        *(w_lnode->childs+1)=w_lnode2; //The second node in a TRFR node is allways the transfered lineage (bound!!).
                        w_lnode->n_child=2;
                        
                        ++n_transfer;
                        if (verbosity>5)
                        {
                            printf("\n\t\t\t\t\tTransference from %d to %d, time %lf",w_lnode->index, w_lnode2->index, w_period2->r_bound);
#ifdef DEBUG
                            fflush(stdout);
#endif
                        }
                    }
                }
                
            }
            if (verbosity>4)
            {
                printf("\n\t\t\tDone");
#ifdef DEBUG
                fflush(stdout);
#endif
            }
            FreePeriods(periods, n_periods);
            (*wlocus_tree)->root=((*wlocus_tree)->m_node+tn_nodes-1);
        }
        // ******
        /// Filling information of the just simulated l_tree
        (*wlocus_tree)->n_leaves=*st_leaves+*st_transfr;
        (*wlocus_tree)->n_nodes=tn_nodes+*st_transfr;
        (*wlocus_tree)->n_gleaves=*st_gleaves;
        (*wlocus_tree)->species_tree=wsp_tree;
        wsp_tree->locus_tree=*wlocus_tree;
        
        if (verbosity>4)
        {
            printf("\n\t\tDone\n");
#ifdef DEBUG
            fflush(stdout);
#endif
        }
        return (NO_ERROR);
    }
    
}

long int SimGTree(l_tree *wlocus_tree, g_tree **gene_tree, name_c * names, float epsilon_brent,float min_cu_bc, int precision, gsl_rng *seed, unsigned int *n_lcoals,unsigned int simlosses, unsigned int verbosity, double gen_time)
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
    unsigned int next_avail_inode=0,max_coals=0,n_anc_nodes=0,p_Ne=0,k=0;
    unsigned int bounded=0, politomy=0;
    unsigned int is_dup=0;
    double p_mu=0;
    double sampled_ngen=0, max_ngen=0, current_ngen=0;
    unsigned int avail_leaves=0,node_index=0, extra_nodes=0;
    
    // ******
    /// Loop related variables</dd></dl>
    unsigned int i=0,j=0,l=0;
    
    // ******
    /// Initialization of other general variables</dd></dl>
    if (verbosity>4)
        iobuffer=calloc(count_intdigits((long)UINT_MAX, 0), sizeof(char));
    
    // ****
    /// Association between locus tree and gene tree and gene tree reset
    MatchTrees(wlocus_tree,*gene_tree,1,simlosses);
    
    w_gnodes=(*gene_tree)->m_node;
    next_avail_inode=wlocus_tree->n_gleaves;
    extra_nodes=0;
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
        bounded=0;
        is_dup=0;
        politomy=0;
        
        // ** Taking into account a new locus tree node ** //
        w_lnode=(wlocus_tree->m_node+i);
        anc_lnode=w_lnode->anc_node;
        k=w_lnode->n_nodes;
        avail_leaves=k;
        w_gnodes_ptr=w_lnode->g_nodes;
        
        if(w_lnode->Ne!=0) //Node with a private Ne.
            p_Ne=w_lnode->Ne;
        else
            p_Ne=wlocus_tree->Ne; //Global Ne.
        
        p_mu=wlocus_tree->mu*w_lnode->mu_mult;
        
        if (k==0)
            continue;
        
        // ** Constraining coalescent time and number of events ** //
        current_ngen=w_lnode->n_gen;
        max_ngen=w_lnode->gen_length+current_ngen;
        max_coals=w_lnode->n_nodes-1;//Max num of coalescent events = number of nodes to coalesce -1
        
        // ***
        /// Selection of Bounded coalescent or regular coalescent if the ancestor is a duplication, depending on the paralog id \ref l_node::paralog
        if(anc_lnode!=NULL)
        {
            is_dup=anc_lnode->kind_node;
            
            if((anc_lnode->kind_node==DUP || anc_lnode->kind_node==TRFR) && *(anc_lnode->childs+1) == w_lnode)
                bounded=1;
        }
        
        if (verbosity==5)
        {
            if (w_lnode->conts==NULL) //There is no species tree
            {
                if(names!=NULL)
                    printf("\n\t\tLocus tree node %s_%d... ",w_lnode->n_gen==0?(names->names+(w_lnode->sp_index*names->max_lname)):"Internal node",w_lnode->paralog);
                else
                {
                    sprintf(iobuffer, "%d",w_lnode->sp_index);
                    printf("\n\t\tLocus tree node %s_%d... ",w_lnode->n_gen==0?iobuffer:"Internal node",w_lnode->paralog);
                }
            }
            else
            {
                if(names!=NULL)
                    printf("\n\t\tSpecies tree node %s, locus tree node %s_%d... ",(names->names+(w_lnode->conts->sp_index*names->max_lname)),w_lnode->n_gen==0?(names->names+(w_lnode->sp_index*names->max_lname)):"Internal node",w_lnode->paralog);
                else
                {
                    sprintf(iobuffer, "%d",w_lnode->conts->sp_index);
                    printf("\n\t\tSpecies tree node %s, locus tree node %s_%d... ",w_lnode->conts->sp_index==0?"Internal node":iobuffer,w_lnode->n_gen==0?iobuffer:"Internal node",w_lnode->paralog);
                }
            }
            
            
#ifdef DEBUG
            fflush(stdout);
#endif
        }
        else if (verbosity==6)
        {
            if (w_lnode->conts==NULL)
            {
                if(names!=NULL)
                    printf("\n\t\tLocus tree node %s_%d (index %d, sp_index %d, n_nodes %d, kind_ancestor %d) ",w_lnode->n_gen==0?(names->names+(w_lnode->sp_index*names->max_lname)):"Internal node",w_lnode->paralog,w_lnode->index, w_lnode->sp_index, k,is_dup);
                else
                {
                    sprintf(iobuffer, "%d",w_lnode->sp_index);
                    printf("\n\t\tLocus tree node %s_%d (index %d, sp_index %d, n_nodes %d, kind_ancestor %d) ",w_lnode->n_gen==0?iobuffer:"Internal node",w_lnode->paralog,w_lnode->index, w_lnode->sp_index, k,is_dup);
                }
            }
            else
            {
                if(names!=NULL)
                    printf("\n\t\tSpecies tree node %s, locus tree node %s_%d (index %d, sp_index %d, n_nodes %d, kind_ancestor %d) ",(names->names+(w_lnode->conts->sp_index*names->max_lname)),w_lnode->n_gen==0?(names->names+(w_lnode->sp_index*names->max_lname)):"Internal node",w_lnode->paralog,w_lnode->index, w_lnode->sp_index, k,is_dup);
                else
                {
                    sprintf(iobuffer, "%d",w_lnode->conts->sp_index);
                    printf("\n\t\tSpecies tree node %s, locus tree node %s_%d (index %d, sp_index %d, n_nodes %d, kind_ancestor %d) ",w_lnode->conts->sp_index==0?"Internal node":iobuffer,w_lnode->n_gen==0?iobuffer:"Internal node",w_lnode->paralog,w_lnode->index, w_lnode->sp_index, k,is_dup);
                }
            }
            
#ifdef DEBUG
            fflush(stdout);
#endif
        }
        
        // ***
        /// <dl><dt> Coalescent simulation loop. Each iteration means one coalescence. Ends if the end of the branch has reached,there is no more posible coalescences (only 1 active g_node leaves) or a politomy is required</dt><dd>
        
        for (j=0; j<max_coals;++j)
        {
            // **
            /// Waiting time (coalescent time) sample. Regular coalescent process:  \f$ Time\sim exp({k\choose 2}) \f$. Bounded coalescent process: Inverse transform sampling over CDF of waiting times of BMC. \f$ CDF= P(x|a,b=1,t,N)= \frac{P(b=1|a-1,t-x,N)P(x|a,N)}{P(b=1|a,t,N)} \f$. Root finding algorithm \ref brent_root "Brent method" over the \ref sample_bounded_coalescent "BMC waiting times sampling function"
            if (bounded==1)
            {
                sampled_ngen=brent_root(*(sample_bounded_coalescent),0.0,max_ngen-current_ngen,gsl_rng_uniform_pos(seed), avail_leaves, p_Ne, max_ngen-current_ngen,precision,epsilon_brent,verbosity);
                if(sampled_ngen>=max_ngen-current_ngen)
                    return UNEXPECTED_VALUE;
                
            }
            else
            {
                sampled_ngen=((-log(gsl_rng_uniform_pos(seed)))*(2*p_Ne))/(k*(k-1));// - (1/lambda) * ln uniform random variable. Lambda = parameter (k over 2). K, number of nodes avaliable to coalesce.
            }
            
            if(sampled_ngen<=0)
                return UNEXPECTED_VALUE;
            else if (sampled_ngen/p_Ne<min_cu_bc)
            {
                politomy=1;
            }
            
            current_ngen+=sampled_ngen;
            --k;
            
            //End of the branch reached. There is no restriction if anc_resnode==NULL (root)
            if (current_ngen>max_ngen && anc_lnode!=NULL)
            {
                
                (*n_lcoals)+=(avail_leaves-1);
                
                if (verbosity==5)
                {
                    printf("\n\t\t\t %d extra lineages going deeper using a %s    ",avail_leaves-1,bounded==1?"bounded coalescent process":"regular coalescent process");
#ifdef DEBUG
                    fflush(stdout);
#endif
                }
                if (verbosity==6)
                {
                    printf("\n\t\t\t %d extra lineages going deeper using a %s, max coalescent time %f, coalescent time %f",avail_leaves-1,bounded==1?"bounded coalescent process":"regular coalescent process",max_ngen*gen_time,current_ngen*gen_time);
#ifdef DEBUG
                    fflush(stdout);
#endif
                }
                break;
                
            }
            else if (politomy==1)
            {
                if (verbosity==5)
                {
                    printf("\n\t\t\t New politomy of %d gene tree nodes, coalescent time %f\n",avail_leaves,current_ngen*gen_time);
#ifdef DEBUG
                    fflush(stdout);
#endif
                }
                if (verbosity==6)
                {
                    printf("\n\t\t\t New politomy of %d gene tree nodes, coalescent time %f, due to very low cu to the bound (biologicaly-senseless short branches).",avail_leaves,current_ngen*gen_time);
#ifdef DEBUG
                    fflush(stdout);
#endif
                }
                
                //New polytomic ancestor
                anc_gnode= w_gnodes + next_avail_inode;
                
                //Ancestor configuration
                anc_gnode->height=max_ngen;
                anc_gnode->contl=w_lnode;
                anc_gnode->conts=w_lnode->conts;
                anc_gnode->paralog=w_lnode->paralog;
                
                //Ancestor and tree reconfiguration due to the polytomy
                if ((*gene_tree)->max_childs<avail_leaves)
                    (*gene_tree)->max_childs=avail_leaves;
                anc_gnode->childs=realloc(anc_gnode->childs, sizeof(g_node *)*avail_leaves); //\attention The use of politomyes in the gene tree is implemented creating "improper" gene trees, as the size of g_node::childs is diferent in each node. This is not a problem in this context, but It would deal memory errors if these trees were used in other context. In this case, a reallocation of the full tree would be compulsory.
                
                //Child pointers and info configuration
                for (l=0;l<avail_leaves;++l)
                {
                    off1=*(w_gnodes_ptr+l);
                    *(anc_gnode->childs+l)=off1;
                    off1->anc_node=anc_gnode;
                    off1->gen_length= anc_gnode->height - off1->height;
                    off1->bl+=(anc_gnode->height-(off1->height<w_lnode->n_gen?w_lnode->n_gen:off1->height))*p_mu;
                }
                anc_gnode->n_child=avail_leaves;
                
                //Current working l_node reconfiguration
                w_lnode->n_nodes=1;
                *(w_gnodes_ptr)= anc_gnode;
                
                //Common variables reconfiguration
                extra_nodes+=avail_leaves-2;
                avail_leaves=1; //Only one remaining
                ++next_avail_inode; //One node used
                
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
            *(anc_gnode->childs+anc_gnode->n_child)=off1;
            *(anc_gnode->childs+anc_gnode->n_child+1)=off2;
            off1->anc_node=anc_gnode;
            off2->anc_node=anc_gnode;
            anc_gnode->contl=w_lnode;
            anc_gnode->conts=w_lnode->conts;
            anc_gnode->n_child+=2;
            
            //Branch lengths and times
            anc_gnode->height=current_ngen;
            off1->gen_length= anc_gnode->height - off1->height;
            off1->bl+=(anc_gnode->height-(off1->height<w_lnode->n_gen?w_lnode->n_gen:off1->height))*p_mu;
            off2->gen_length= anc_gnode->height - off2->height;
            off2->bl+=(anc_gnode->height-(off2->height<w_lnode->n_gen?w_lnode->n_gen:off2->height))*p_mu;
            
            //Info
            anc_gnode->paralog=w_lnode->paralog;
            
            if (verbosity==5)
            {
                if (bounded==0)
                {
                    printf("\n\t\t\tNew regular coalescence\t\t    ");
#ifdef DEBUG
                    fflush(stdout);
#endif
                }
                else
                {
                    printf("\n\t\t\tNew bounded coalescence\t\t    ");
#ifdef DEBUG
                    fflush(stdout);
#endif
                }
                
                
            }
            else if (verbosity==6)
            {
                if (w_lnode->anc_node==NULL)
                {
                    printf("\n\t\t\tNew regular coalescence, nodes %u and %u, coalescent time %f",off1->index,off2->index,current_ngen*gen_time);
#ifdef DEBUG
                    fflush(stdout);
#endif
                }
                else if (bounded==0)
                {
                    printf("\n\t\t\tNew regular coalescence, nodes %u and %u, max coalescent time %f, coalescent time %f",off1->index,off2->index,max_ngen*gen_time,current_ngen*gen_time);
#ifdef DEBUG
                    fflush(stdout);
#endif
                }
                else
                {
                    printf("\n\t\t\tNew bounded coalescence, nodes %u and %u, max coalescent time %f, coalescent time %f",off1->index, off2->index,max_ngen*gen_time,current_ngen*gen_time);
#ifdef DEBUG
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
                off1->bl+=(off1->height<w_lnode->n_gen?w_lnode->gen_length:max_ngen-off1->height)*p_mu; //Addition of the branch length of this period, which can't be added at the end of the g_node branch due to the lineage especific substitution rates.
                ++anc_lnode->n_nodes;
            }
            
            if (verbosity==6)
            {
                printf("\n\t\t\tConfiguring ancestor restrictions: n_nodes=%d",anc_lnode->n_nodes);
#ifdef DEBUG
                fflush(stdout);
#endif
            }
        }
        else //Last
        {
            (*gene_tree)->root=anc_gnode; //Saves the root of the gene tree
            (*gene_tree)->n_nodes-=extra_nodes; //Unused nodes due to the presence of polytomies.
        }
        
        if (verbosity>4)
        {
            printf("\n\t\tDone");
#ifdef DEBUG
            fflush(stdout);
#endif
        }
        
    }
    
    if (verbosity>4)
        free(iobuffer);
    
    return(NO_ERROR);
}

l_tree * NewLTree (unsigned int n_nodes, unsigned int n_leaves, unsigned int n_gleaves, unsigned int max_childs, double gen_time, unsigned int Ne, double mu)
{
    l_tree * tree=NULL;
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Tree memory allocation
    tree=calloc(1,sizeof(l_tree));
    ErrorReporter((long int)tree);
    
    // *
    /// Tree initialization
    tree->n_nodes=n_nodes;
    tree->n_leaves=n_leaves;
    tree->n_gleaves=n_gleaves;
    tree->max_childs=max_childs;
    tree->gen_time=gen_time;
    tree->Ne=Ne;
    tree->mu=mu;
    tree->species_tree=NULL;
    tree->gene_tree=NULL;
    
    // *
    /// Tree nodes allocation and initialization by \ref NewLNodes </dd></dl>
    
    if (n_nodes>1)
    {
        tree->m_node=NewLNodes(n_nodes,n_gleaves,max_childs);
        tree->root=NULL;
    }
    else if (n_nodes==1)
    {
        tree->root=NewLNodes(n_nodes, n_gleaves, max_childs);
        tree->m_node=NULL;
    }
    else
    {
        tree->root=NULL;
        tree->m_node=NULL;
    }
    
    return (tree);
}

l_tree * ReadNewickLTree (char * newick,name_c **names_ptr, unsigned int verbosity, double gen_time, unsigned int Ne, double mu, unsigned int ind_persp)
{
    l_node *current_node=NULL, *anc_node=NULL,*root=NULL;
    l_tree *tree=NULL;
    name_c * names=NULL;
    char code=' ';
    unsigned int register step=0;
    unsigned int n_char=0, iteration=0, ffree_codename=1, n_leaves=0, n_gleaves=0, n_inodes=0, n_nodes=0, max_childs=0, max_lname=0, n_replica=0,index=0, is_dup=0, n_paralog=0, n_priv_ngleaves=0;
    char bl_buffer[DBL_DIG+3]=""; //DBL_DIG= precission of double, + 3 positions (0 separator and \0)
    char * ui_buffer;
    char name_buffer[MAX_NAME]="";
    
    // ****
    /// <dl><dt> Function structure </dt><dd>
    
    // ***
    /// Test of the newick tree string by \ref CheckNewickLTree
    
    ErrorReporter(CheckNewickLTree(newick));
    
    // ***
    /// Error control
    
    if (gen_time==0)
        ErrorReporter(UNEXPECTED_VALUE);
    
    // ***
    /// Integer buffer allocation and initialization
    
    ui_buffer=calloc((int)log10(UINT_MAX)+2,sizeof(char)); // An string of the max number of digits of an unsigned int + 1 (\0)
    strcpy(ui_buffer, "");
    
    // ***
    /// First read of the Newick tree. Obtains the number of internal nodes (")"), s_tree::n_leaves ("(" or "," not followed by "(") and s_tree::n_gleaves (s_tree::n_leaves + (replicas "/" -1 for each node)).
    while (*(newick+step)!=';')
    {
        if((*(newick+step)=='(' ||*(newick+step)==',' )&&*(newick+step+1)!='(') ++n_leaves;
        if(*(newick+step)==')') ++n_inodes;
        if(*(newick+step)=='/')
        {
            // * Reseting variables * //
            strcpy(ui_buffer,"");
            n_replica=0;
            n_char=0;
            ++step;
            code=*(newick+step);
            // * Reading replicas digits * //
            while (code<58 && code>47) // Int numbers
            {
                *(ui_buffer+n_char)=code;
                ++n_char;
                ++step;
                code=*(newick+step);
            }
            --step;
            // * Translating * //
            sscanf(ui_buffer,"%ui",&n_replica);
            n_gleaves+=n_replica;
            ++n_priv_ngleaves;
        }
        ++step;
        if (step==MAX_IT) ErrorReporter(LOOP_ERROR); // Avoids ininite loops
    }
    n_nodes=n_leaves+n_inodes;
    n_gleaves+=(n_leaves-n_priv_ngleaves)*ind_persp; //Addition of the ngleaves of the nodes using the common number of individuals (ind_persp)
    
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
                    current_node=NewLNodes(1, 0, MAX_CHILDS);
                    // *
                    /// Points the pointers of this new node and its ancestor
                    *(anc_node->childs+anc_node->n_child)=current_node;
                    ++anc_node->n_child;
                    current_node->anc_node=anc_node;
                    // *
                    /// Searches for the maximum number of childs in the tree.</dd></dl>
                    if(max_childs<anc_node->n_child)
                    {
                        max_childs=anc_node->n_child;
                    }
                    
                }
                else //New root node
                {
                    current_node=NewLNodes(1,0,MAX_CHILDS);
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
                strcpy(bl_buffer,"");
                n_char=0;
                ++step;
                code=*(newick+step);
                while (code<58 && (code>47 || code==46))
                {
                    *(bl_buffer+n_char)=code;
                    ++n_char;
                    ++step;
                    code=*(newick+step);
                }
                --step;
                *(bl_buffer+n_char)=0; //Sets the end of the string to avoid problems in the sscanf
                // *
                /// Translates the buffer in a float number and asigns it as s_node::gen_length of the current node</dd></dl>
                sscanf(bl_buffer,"%lf",&current_node->gen_length);
                current_node->gen_length/=gen_time;
                
                break;
            case '*':
                // **
                /// <dl><dt>Lineage specific substitution rate multi(code="*").</dt><dd>
                
                // *
                /// Reads all the following integers and . in a buffer
                strcpy(bl_buffer,"");
                n_char=0;
                ++step;
                code=*(newick+step);
                while (code<58 && (code>47 || code==46))
                {
                    *(bl_buffer+n_char)=code;
                    ++n_char;
                    ++step;
                    code=*(newick+step);
                }
                --step;
                *(bl_buffer+n_char)=0; //Sets the end of the string to avoid problems in the sscanf
                // *
                /// Translates the buffer in a float number and asigns it as l_node::mu_mult </dd></dl>
                sscanf(bl_buffer,"%lf",&current_node->mu_mult);
                
                break;
            case '#':
                // **
                /// <dl><dt>New effective population size (Ne) (code="#").</dt><dd>
                
                // *
                /// Reads all the following integers in a buffer
                strcpy(ui_buffer,"");
                n_char=0;
                ++step;
                code=*(newick+step);
                while (code<58 && code>47)
                {
                    *(ui_buffer+n_char)=code;
                    ++n_char;
                    ++step;
                    code=*(newick+step);
                }
                --step;
                *(ui_buffer+n_char)=0; //Sets the end of the string to avoid problems in the sscanf
                // *
                /// Translates the buffer in a integer number and asigns it as l_node::Ne of the current node</dd></dl>
                sscanf(ui_buffer,"%ui",&current_node->Ne);
                if (current_node->Ne==0)
                    ErrorReporter(SETTINGS_ERROR);
                
                break;
            case '/':
                // **
                /// <dl><dt>New number of replicas (code="/").</dt><dd>
                
                // *
                /// Reads all the following integers in a buffer
                strcpy(ui_buffer,"");
                n_replica=0;
                n_char=0;
                ++step;
                code=*(newick+step);
                while (code<58 && code>47)
                {
                    *(ui_buffer+n_char)=code;
                    ++n_char;
                    ++step;
                    code=*(newick+step);
                }
                --step;
                *(ui_buffer+n_char)=0; //Sets the end of the string to avoid problems in the sscanf
                // *
                /// Translates the buffer in a integer number and asigns it as l_node::n_nodes of the current node</dd></dl>
                sscanf(ui_buffer,"%ui",&n_replica);
                if (n_replica==0)
                    ErrorReporter(SETTINGS_ERROR);
                current_node->n_nodes=n_replica;
                break;
            case '%':
                // **
                /// <dl><dt>New duplication info (code="%").</dt><dd>
                
                // *
                /// Reads all the following integers in a buffer
                strcpy(ui_buffer,"");
                n_replica=0;
                n_char=0;
                ++step;
                code=*(newick+step);
                while (code<58 && code>47)
                {
                    *(ui_buffer+n_char)=code;
                    ++n_char;
                    ++step;
                    code=*(newick+step);
                }
                --step;
                *(ui_buffer+n_char)=0; //Sets the end of the string to avoid problems in the sscanf
                // *
                /// Translates the buffer in a integer number and asigns it as l_node::kind_node of the current node</dd></dl>
                sscanf(ui_buffer,"%ui",&is_dup);
                if (is_dup!=DUP)
                {
                    //DEBUG
                    fprintf(stderr,"WARNING!!!: The user is setting private values to the l_node event info\n");
                    //ErrorReporter(SETTINGS_ERROR);
                }
                current_node->kind_node=is_dup;
                break;
            case '_':
                
                // **
                /// <dl><dt>New paralog info (code="_").</dt><dd>
                
                // *
                /// Reads all the following integers in a buffer
                strcpy(ui_buffer,"");
                n_paralog=0;
                n_char=0;
                ++step;
                code=*(newick+step);
                while (code<58 && code>47)
                {
                    *(ui_buffer+n_char)=code;
                    ++n_char;
                    ++step;
                    code=*(newick+step);
                }
                --step;
                *(ui_buffer+n_char)=0; //Sets the end of the string to avoid problems in the sscanf
                // *
                /// Translates the buffer in a integer number and asigns it as l_node::paralog of the current node</dd></dl>
                sscanf(ui_buffer,"%ui",&n_paralog);
                current_node->paralog=n_paralog;
                break;
                
            default:
                // **
                /// <dl><dt>New leaf(code!= former ones).</dt><dd>
                
                // *
                /// Reads the name of the leaf in a buffer
                
                strcpy(name_buffer,"");
                n_char=0;
                while (code!='(' && code!=')' && code!= ',' && code!= ';' && code!= ':' && code!= '_')
                {
                    *(name_buffer+n_char)=code;
                    ++n_char;
                    ++step;
                    code=*(newick+step);
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
                current_node=NewLNodes(1,0,MAX_CHILDS);
                // *
                /// Points the pointers of this new node and its ancestor
                *(anc_node->childs+anc_node->n_child)=current_node;
                ++anc_node->n_child;
                current_node->anc_node=anc_node;
                // *
                /// Asigns the new name code (position) to the current node (l_node::lp_index)
                current_node->sp_index=ffree_codename;
                ++ffree_codename;
                // *
                /// Searches for the maximum number of childs in the tree.</dd></dl></dd></dl>
                if(max_childs<anc_node->n_child)
                {
                    max_childs=anc_node->n_child;
                }
                
                break;
        }
        
        ++step;
        code=*(newick+step);
        ++iteration;
        if (iteration==MAX_IT) ErrorReporter(LOOP_ERROR); // Avoids ininite loops
        
    }
    
    // ***
    /// Refines the readed nodes by \ref RefineLNodes
    
    RefineLNodes(root,n_gleaves,ind_persp);
    max_lname++; //One extra character for \0.
    
    // ***
    /// Reallocates the name_c memory by \ref ReallocNames
    ReallocNames(names,max_lname);
    
    // ***
    /// Allocates memory for the readed tree and completes it.
    
    tree=NewLTree(0, n_leaves, n_gleaves, max_childs, gen_time, Ne, mu);//Nodes have already been allocated.
    tree->root=root;
    tree->n_nodes=n_nodes;
    tree->species_tree=NULL;
    tree->gene_tree=NULL;
    
    // ***
    /// Fills the node indexes following a post_order.
    PostReorderLNodes(tree->root,&index);
    
    // ***
    /// Frees dynamic memory (buffers)</dd></dl>
    free(ui_buffer);
    ui_buffer=NULL;
    
    if (verbosity>2)
    {
        printf("\n\t\t %d-node locus tree correctly built",(n_leaves*2)-1);
        if (verbosity>3)
        {
            printf(": ");
            WriteLNodes(root,names,gen_time);
        }
        printf("\n");
        
    }
    
    return (tree);
}

g_tree * NewGTree (unsigned int n_nodes, unsigned int max_childs)
{
    g_tree * tree=NULL;
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Tree memory allocation
    tree=calloc(1,sizeof(g_tree));
    ErrorReporter((long int) tree);
    
    // *
    /// Tree initialization
    tree->n_nodes=n_nodes;
    tree->max_childs=max_childs;
    tree->species_tree=NULL;
    tree->locus_tree=NULL;
    
    // *
    /// Tree nodes allocation and initialization by \ref NewGNodes </dd></dl>
    if (n_nodes<2)
    {
        tree->root=NewGNodes(n_nodes,max_childs);
        tree->m_node=NULL;
    }
    else
    {
        tree->m_node=NewGNodes(n_nodes,max_childs);
        tree->root=NULL;
    }
    
    
    return (tree);
}

// ** Tree copy ** //

long int CopySTree (s_tree ** out_tree_ptr,s_tree * in_tree,unsigned int tree_struct, unsigned int l_nodes_ptr)
{
    s_tree * out_tree=NULL;
    s_node * w_output=NULL, * w_input=NULL;
    unsigned int i=0, j=0,post_order=0;
    
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
        *out_tree_ptr=NewSTree(in_tree->n_nodes, in_tree->n_leaves,in_tree->n_gleaves, in_tree->max_childs,in_tree->gen_time, in_tree->Ne, in_tree->mu);
    }
    // **
    /// Reallocation (by \ref FreeSTree and \ref NewSTree) if the output tree pointer has a sparse tree, and s_tree::n_nodes or s_tree::max_childs are different. </dd></dl>
    else if ((*out_tree_ptr)->m_node==NULL || (*out_tree_ptr)->n_nodes!=in_tree->n_nodes || (*out_tree_ptr)->max_childs!=in_tree->max_childs)
    {
        FreeSTree(out_tree_ptr);
        *out_tree_ptr=NewSTree(in_tree->n_nodes, in_tree->n_leaves,in_tree->n_gleaves, in_tree->max_childs,in_tree->gen_time, in_tree->Ne, in_tree->mu);
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
        w_output->gen_length=w_input->gen_length;
        w_output->mu_mult=w_input->mu_mult;
        
        
        // **
        /// <dl><dt>Copy of tree struct pointers</dt><dd>
        if (tree_struct>0)
        {
            if (w_output->childs==NULL)
                return MEM_ERROR;
            // *
            /// Child loop
            for (j=0;j<w_input->n_child;++j)
            {
                *(w_output->childs+j)=(out_tree->m_node+(*(w_input->childs+j))->index);
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

long int CopyLTree (l_tree **out_tree_ptr, l_tree *in_tree, unsigned int tree_struct, unsigned int l_nodes_ptr, unsigned int g_nodes_ptr, unsigned int relink)
{
    l_tree * out_tree=NULL;
    l_node * w_output=NULL, * w_input=NULL;
    unsigned int i=0, j=0, post_order=0;
    
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
        CollapseLTree(in_tree,post_order,1);
    }
    
    // ***
    /// <dl><dt>Output tree memory control</dt><dd>
    
    // **
    /// Allocation if *out_tree_ptr==NULL
    if (*out_tree_ptr==NULL)
    {
        *out_tree_ptr=NewLTree(in_tree->n_nodes, in_tree->n_leaves,in_tree->n_gleaves, in_tree->max_childs, in_tree->gen_time, in_tree->Ne, in_tree->mu);
    }
    // **
    /// Reallocation (by \ref FreeSTree and \ref NewLTree) if the output tree pointer has a sparse tree, and l_tree::n_nodes or l_tree::max_childs are different. </dd></dl>
    else if ((*out_tree_ptr)->m_node==NULL || (*out_tree_ptr)->n_nodes!=in_tree->n_nodes || (*out_tree_ptr)->max_childs!=in_tree->max_childs)
    {
        FreeLTree(out_tree_ptr);
        *out_tree_ptr=NewLTree(in_tree->n_nodes, in_tree->n_leaves,in_tree->n_gleaves, in_tree->max_childs, in_tree->gen_time, in_tree->Ne, in_tree->mu);
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
        w_output->gen_length=w_input->gen_length;
        w_output->conts=w_input->conts;
        w_output->mu_mult=w_input->mu_mult;
        
        if (tree_struct>0)
        {
            if (w_output->childs==NULL)
                return MEM_ERROR;
            
            for (j=0;j<w_input->n_child;++j)
            {
                *(w_output->childs+j)=(out_tree->m_node+(*(w_input->childs+j))->index);
            }
            w_output->anc_node=NULL;
        }
        
        if (g_nodes_ptr>0)
        {
            for (j=0; j>w_input->n_nodes; ++j)
            {
                *(w_output->childs+j)=out_tree->m_node+((*(w_input->childs+j))->index);
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
        w_output->gen_length=w_input->gen_length;
        w_output->conts=w_input->conts;
        w_output->mu_mult=w_input->mu_mult;
        
        // **
        /// <dl><dt>Copy of tree struct pointers</dt><dd>
        if (tree_struct>0)
        {
            if (w_output->childs==NULL)
                return MEM_ERROR;
            
            // *
            /// Child loop
            for (j=0;j<w_input->n_child;++j)
            {
                *(w_output->childs+j)=(out_tree->m_node+(*(w_input->childs+j))->index);
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
                *(w_output->childs+j)=out_tree->m_node+((*(w_input->childs+j))->index);
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

inline long int CopyStoLTree(s_tree *sp_tree, l_tree *locus_tree)
{
    unsigned int i=0,j=0;
    s_node *w_snode=NULL;
    l_node *w_lnode=NULL;
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Input tree memory control
    
    // *
    /// Tree equivalence 
    if (sp_tree->m_node==NULL||sp_tree->max_childs!=locus_tree->max_childs||sp_tree->n_gleaves!=locus_tree->n_gleaves||sp_tree->n_leaves!=locus_tree->n_leaves||sp_tree->n_nodes!=locus_tree->n_nodes)
        return MEM_ERROR;
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
        w_lnode->gen_length=w_snode->gen_length;
        w_lnode->conts=w_snode;
        w_lnode->mu_mult=w_snode->mu_mult;
        
        for (j=0;j<w_lnode->n_child; ++j)
        {
            *(w_lnode->childs+j)=(locus_tree->m_node+((*(w_snode->childs+j))->index));
        }
        
        if (w_snode->anc_node != NULL)
            w_lnode->anc_node=(locus_tree->m_node+(w_snode->anc_node->index));
            
    }

    // *
    /// Tree info copy </dd></dl>
    locus_tree->root=locus_tree->m_node+locus_tree->n_nodes-1;
    locus_tree->gen_time=sp_tree->gen_time;
    locus_tree->Ne=sp_tree->Ne;
    locus_tree->mu=sp_tree->mu;
    locus_tree->species_tree=sp_tree;
    sp_tree->locus_tree=locus_tree;
    
    return NO_ERROR;
}



// ** Tree edition ** //

long int CleanlossesLTree(l_tree *locus_tree)
{
    unsigned int n_deletions=0;
    unsigned int n_leaves=0;
    if (locus_tree->m_node!=NULL || locus_tree->root == NULL)
        return MEM_ERROR;
    
    CleanlossesLNodes(locus_tree->root,&locus_tree->root,&n_deletions,&n_leaves);
    
    locus_tree->n_nodes-=n_deletions;
    locus_tree->n_leaves=n_leaves;
    
    return NO_ERROR;
}

// ** Tree reset ** //

inline long int ResetSTreeSimL (s_tree *tree)
{
    unsigned int i=0;
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

inline long int ResetGTree (g_tree *tree)
{
    unsigned int i=0,j=0;
    g_node * w_node=NULL;
    
    if (tree->m_node==NULL)
        return MEM_ERROR;
    
    for (i=0;i<tree->n_nodes;++i)
    {
        w_node=tree->m_node+i;
        
        w_node->sp_index=0;
        w_node->replica=0;
        w_node->paralog=0;
        w_node->height=0;
        w_node->bl=0;
        w_node->gen_length=0.0;
        for (j=0;j<tree->max_childs;++j)
        {
            *(w_node->childs+j)=NULL;
        }
        w_node->anc_node=NULL;
        
    }
    
    return NO_ERROR;
}

// ** Tree deletion ** //

void FreeSTree(s_tree ** tree)
{
    register unsigned int i;
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
            /// Frees s_node::childs </dd></dl>
            if (w_node->childs !=NULL)
            {
                free(w_node->childs);
                w_node->childs=NULL;
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
    register unsigned int i;
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
            /// Frees l_node::g_nodes and l_node::childs </dd></dl>
            if (w_node->g_nodes !=NULL)
            {
                free(w_node->g_nodes);
                w_node->g_nodes=NULL;
            }
            if (w_node->childs !=NULL)
            {
                free(w_node->childs);
                w_node->childs=NULL;
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

void FreeGTree(g_tree ** tree, unsigned int complete)
{
    unsigned int i;
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
            /// Frees g_node::childs</dd></dl>
            if (w_node->childs !=NULL)
            {
                free(w_node->childs);
                w_node->childs=NULL;
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
        (*tree)->root->height=0;
        (*tree)->root->n_child=0;
        (*tree)->root->bl=0;
        (*tree)->root->gen_length=0.0;
        (*tree)->root->contl=NULL;
        (*tree)->root->conts=NULL;
        (*tree)->root->anc_node=NULL;
    }
}

void FreePeriods(period * periods, unsigned int n_periods)
{
    unsigned int i=0;
    
    for (i=0; i<n_periods; ++i)
    {
        if ((periods+i)->l_nodes!=NULL)
            free((periods+i)->l_nodes);
    }
    
    free(periods);
}

// *** Names memory manage *** //

name_c * NewNames (unsigned int n_names, unsigned int max_lname)
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
    ErrorReporter((long int) names);
    
    // *
    /// Inizialitation of name_c
    names->n_names=n_names;
    names->max_lname=max_lname;
    
    // *
    /// Allocation of the string of names (name_c::names)
    names->names=calloc((n_names+1)*max_lname,sizeof(char));
    ErrorReporter((long int)names->names);
    
    // *
    /// Inizialitation of the string </dd></dl>
    strncpy(names->names,"Internal node",max_lname);
    
    return (names);
    
}

void ReallocNames (name_c * names,unsigned int max_lname)
{
    char * new_names=NULL;
    unsigned int register i=0;
    unsigned int min_max_lname=0;
    
    // *
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Control of max_lname variable
    min_max_lname=(unsigned int)strlen(names->names)+1; // strlen is measuring only the first name, usually "Internal node"". +1 due to \0.
    
    if (max_lname<min_max_lname)
        max_lname=min_max_lname;
    
    // *
    /// Allocation of the new string of names
    new_names=calloc((names->n_names+1)*max_lname,sizeof(char));
    ErrorReporter((long int)new_names);
    
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

long int MatchTrees(l_tree *locus_tree, g_tree *gene_tree, unsigned int reset_gtree, unsigned int includelosses)
{
    unsigned int i=0;
    unsigned int j=0,next_leaf=0;
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
        
        CollapseLTree(locus_tree,1,0); //Post-order
    }
    if (includelosses>1)
        return UNEXPECTED_VALUE;
    
    // ***
    /// <dl><dt>Gene tree reset if it is required</dt><dd>
    
    // **
    /// <dl><dt>Gene tree loop</dt><dd>
    if (reset_gtree==1)
    {
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
            w_gnode->height=0.0;
            w_gnode->bl=0.0;
            w_gnode->gen_length=0.0;
        }
        
        // **
        /// Gene tree root deletion</dd></dl>
        gene_tree->root=NULL;
    }
    
    // ***
    /// <dl><dt>Locus tree loop</dt><dd>
    for (i=0;i<locus_tree->n_nodes;++i)
    {
        w_lnode=locus_tree->m_node+i;
        
        // **
        /// <dl><dt>Leaf nodes association </dt><dd>
        if (w_lnode->n_child==0)
        {
            // *
            /// The l_node::n_nodes is reset using l_node::conts s_node::n_replicas
            if (w_lnode->conts!=NULL)
                w_lnode->n_nodes=w_lnode->conts->n_replicas;
            if(w_lnode->kind_node==LOSS || w_lnode->kind_node==RTRFR)
                w_lnode->n_nodes=includelosses;
            
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
            }
        }
        else
        {
            w_lnode->n_nodes=0;
        }
    }
    
    gene_tree->species_tree=locus_tree->species_tree;
    gene_tree->locus_tree=locus_tree;
    locus_tree->gene_tree=gene_tree;
    
    return (NO_ERROR);
}

// ** Tree conversion ** //

long int CollapseSTree (s_tree * in_tree, unsigned int post_order)
{
    s_node * m_node=NULL, * p_m_node=NULL, * w_node=NULL, * p_root=NULL;
    unsigned int index=0,i=0;
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Input tree root control
    if (in_tree->root == NULL)
        return MEM_ERROR;
    
    // *
    /// New s_node array allocation by \ref NewSTree
    m_node=NewSNodes(in_tree->n_nodes, in_tree->max_childs);
    
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
            if (w_node->childs !=NULL)
            {
                free(w_node->childs);
                w_node->childs=NULL;
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

long int CollapseLTree (l_tree * in_tree, unsigned int post_order, unsigned int relink)
{
    l_node * m_node=NULL, * p_m_node=NULL, *w_node=NULL, *p_root=NULL;
    unsigned int index=0, i=0;
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Input tree root control
    if (in_tree->root == NULL)
        return MEM_ERROR;
    
    // *
    /// New l_node array allocation by \ref NewLTree
    m_node=NewLNodes(in_tree->n_nodes, in_tree->n_gleaves, in_tree->max_childs);
    
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
        PreCollapseLNodes(m_node,in_tree->root, in_tree->n_gleaves, relink);

        // *
        /// Asigns the new root to the tree.
        in_tree->root=m_node; //The first node is the root.
    }
    else
    {
        PostReorderLNodes(in_tree->root, &index);
        PostCollapseLNodes(m_node,in_tree->root, in_tree->n_gleaves, relink);
        
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
            if (w_node->childs !=NULL)
            {
                free(w_node->childs);
                w_node->childs=NULL;
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

long int CollapseGTree (g_tree * in_tree, unsigned int post_order, unsigned int relink)
{
    g_node * m_node=NULL, * p_m_node=NULL, *w_node=NULL, *p_root=NULL;
    unsigned int index=0, i=0;
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Input tree root control
    if (in_tree->root == NULL)
        return MEM_ERROR;
    
    // *
    /// New g_node array allocation by \ref NewGTree
    m_node=NewGNodes(in_tree->n_nodes,in_tree->max_childs);
    
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
            if (w_node->childs !=NULL)
            {
                free(w_node->childs);
                w_node->childs=NULL;
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

long int Temporalize_GTree(g_tree *gene_tree,double gen_time)
{
    g_node *w_node;
    unsigned int i=0;
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Input tree root control
    if (gene_tree->root == NULL && gene_tree->m_node == NULL)
        return MEM_ERROR;
    // *
    /// Branch length conversion inside a loop (collapsed \ref g_tree ) or using a recursive function (\ref Temporalize_GNodes) </dd></dl>
    else if (gene_tree->m_node!=NULL)
    {
        for (i=0; i<gene_tree->n_nodes;++i)
        {
            w_node=gene_tree->m_node+i;
            w_node->bl*=gen_time;
            w_node->height*=gen_time;
            w_node->gen_length*=gen_time;
        }
    }
    else
    {
        Temporalize_GNodes(gene_tree->root, gen_time);
    }

    return NO_ERROR;
}

long int Rateheter_lineagespec(s_tree *sp_tree, double alpha, gsl_rng *seed, FILE * gammadump)
{
    s_node *w_node;
    unsigned int i=0;
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
    unsigned int i=0;
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


//long int Mutate_GTree(g_tree *gene_tree,double mu)
//{
//    g_node *w_node;
//    unsigned int i=0;
//    
//    // **
//    /// <dl><dt> Function structure </dt><dd>
//    
//    // *
//    /// Input tree control
//    if (gene_tree->root==NULL || gene_tree->m_node == NULL)
//        return MEM_ERROR;
//    if (gene_tree->m_node->contl==NULL || gene_tree->m_node==gene_tree->root) //Pre-order
//        return MEM_ERROR;
//   
//    // *
//    /// Branch length conversion inside a post-order loop
//
//    for (i=0; i<gene_tree->n_nodes;++i)
//    {
//        w_node=gene_tree->m_node+i;
//        w_node->bl*=w_node->contl->mu!=0?w_node->contl->mu:mu;
//        w_node->height=0;
//    }
//    
//    // * Height conversion inside an inverse post-order loop</dd></dl>
//    ///
//    
//    (gene_tree->m_node+gene_tree->n_nodes-1)->height=0;
//    
//    for (i=gene_tree->n_nodes-1; i>0;--i)
//    {
//        w_node=gene_tree->m_node+i-1;
//        w_node->height=w_node->anc_node->height-w_node->bl;
//    }
//
//    return NO_ERROR;
//}
//
//long int Evolve_GTree(g_tree *gene_tree,double gen_time, double mu)
//{
//    g_node *w_node;
//    unsigned int i=0;
//    
//    // **
//    /// <dl><dt> Function structure </dt><dd>
//    
//    // *
//    /// Input tree control
//    if (gene_tree->root==NULL || gene_tree->m_node == NULL)
//        return MEM_ERROR;
//    if (gene_tree->m_node->contl==NULL || gene_tree->m_node==gene_tree->root) //Pre-order
//        return MEM_ERROR;
//    
//    // *
//    /// Branch length conversion inside a post-order loop
//    
//    for (i=0; i<gene_tree->n_nodes;++i)
//    {
//        w_node=gene_tree->m_node+i;
//        w_node->bl*=gen_time*(w_node->contl->mu!=0?w_node->contl->mu:mu);
//        w_node->height=0;
//    }
//    
//    // * Height conversion inside an inverse post-order loop</dd></dl>
//    ///
//    
//    (gene_tree->m_node+gene_tree->n_nodes-1)->height=0;
//    
//    for (i=gene_tree->n_nodes-1; i>0;--i)
//    {
//        w_node=gene_tree->m_node+i-1;
//        w_node->height=w_node->anc_node->height-w_node->bl;
//    }
//    
//    return NO_ERROR;
//}

// *** Tree features obtention *** //

long int Count_duplications(l_tree *tree, unsigned int *n_dup)
{
    unsigned int i=0;
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

long int Count_losses(l_tree *tree, unsigned int *n_losses)
{
    unsigned int i=0;
    l_node *w_lnode=NULL;
    
    if (tree->root==NULL)
    {
        return MEM_ERROR;
    }
    *n_losses=0;
    if (tree->m_node==NULL)
    {
        *n_losses=Count_losses_lnodes(tree->root);
    }
    else
    {
        for (i=0; i<tree->n_nodes; ++i)
        {
            w_lnode=tree->m_node+i;
            if (w_lnode->kind_node==LOSS)
                ++(*n_losses);
        }
    }
    
    return NO_ERROR;
    
}

long int Measure_ST_height(s_tree *tree,long double *height, unsigned int unit)
{
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Input tree memory control
    if (tree->root == NULL)
        return MEM_ERROR;
    *height=0;
    
    // *
    /// Recursive function selection depending on the height unit.</dd></dl>
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

long int Measure_ST_length(s_tree *tree,long double *length,unsigned int unit)
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

long int Measure_GT_height(g_tree *tree,long double *height,unsigned int unit)
{
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Input tree memory control
    if (tree->root == NULL)
        return MEM_ERROR;
    *height=0;
    
    // *
    /// Recursive function selection depending on the height unit.</dd></dl>
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

long int Measure_GT_length(g_tree *tree,long double *length,unsigned int unit)
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


// *** Tree I/O *** //

long int WriteSTree (s_tree *in_tree, name_c * names, unsigned int time)
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
        WriteSNodes(in_tree->root,names,in_tree->gen_time);
        return(NO_ERROR);
    }
    else
    {
        WriteSNodes(in_tree->root,names,1);
        return(NO_ERROR);
    }
    
}
long int WriteSTreeFile (FILE *file,s_tree *in_tree, name_c * names, unsigned int time)
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
        WriteSNodesFile(file,in_tree->root,names,in_tree->gen_time);
        return(NO_ERROR);
    }
    else
    {
        WriteSNodesFile(file,in_tree->root,names,1);
        return(NO_ERROR);
    }
    
}

long int WriteLTree (l_tree *in_tree, name_c * names, unsigned int time)
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
        WriteLNodes(in_tree->root,names,in_tree->gen_time);
        return(NO_ERROR);
    }
    else
    {
        WriteLNodes(in_tree->root,names,1);
        return(NO_ERROR);
    }
    
}

long int WriteLTreeFile (FILE *file,l_tree *in_tree, name_c * names, unsigned int time)
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
        WriteLNodesFile(file,in_tree->root,names,in_tree->gen_time);
        return(NO_ERROR);
    }
    else
    {
        WriteLNodesFile(file,in_tree->root,names,1);
        return(NO_ERROR);
    }
    
}

long int WriteGTree (g_tree *in_tree, name_c * names)
{
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Writes the tree by \ref WriteGNodes (recursive)</dd></dl>
    if (in_tree->root==NULL)
        return (MEM_ERROR);
    else
    {
        WriteGNodes(in_tree->root,names);
        return(NO_ERROR);
    }
    
}

long int WriteGTreeStr (char * string, g_tree *in_tree, name_c * names)
{
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Writes the tree by \ref WriteGNodesStr (recursive)</dd></dl>
    if (in_tree->root==NULL)
        return (MEM_ERROR);
    else
    {
        WriteGNodesStr(string,in_tree->root,names);
        return(NO_ERROR);
    }
}

long int WriteGTreeFile (FILE * file, g_tree *in_tree, name_c * names)
{
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Writes the tree by \ref WriteGNodesFile (recursive)</dd></dl>
    if (in_tree->root==NULL)
        return (MEM_ERROR);
    else
    {
        WriteGNodesFile(file,in_tree->root,names);
        return(NO_ERROR);
    }
    
}

long int WriteReconSL(s_tree *wsp_tree, l_tree *locus_tree, name_c *names, char *reconsl_outname)
{
    unsigned int i=0;
    l_node *wl_node=NULL;
    FILE *reconsl_outfile=NULL;
    
    if ((reconsl_outfile=fopen(reconsl_outname, "w"))==NULL)
    {
        return(IO_ERROR);
    }
    fprintf(reconsl_outfile,"Lt_node\tLt_paralog\tLt_event\tSt_node\n");
    
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
                            fprintf(reconsl_outfile,"%d\t%d\t%s\t%d\n",wl_node->index,wl_node->paralog,"Sp",wl_node->conts->index);
                            break;
                        case DUP:
                            fprintf(reconsl_outfile,"%d\t%d\t%s\t%d\n",wl_node->index,wl_node->paralog,"Dup",wl_node->conts->index);
                            break;
                        case TRFR:
                            fprintf(reconsl_outfile,"%d\t%d\t%s\t%d\n",wl_node->index,wl_node->paralog,"Transf",wl_node->conts->index);
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
                                fprintf(reconsl_outfile,"%d\t%d\t%s\t\'%d\'\n",wl_node->index,wl_node->paralog,"Sp",wl_node->sp_index);
                                break;
                            case DUP:
                                fprintf(reconsl_outfile,"%d\t%d\t%s\t\'%d\'\n",wl_node->index,wl_node->paralog,"Dup",wl_node->sp_index);
                                break;
                            case TRFR:
                                fprintf(reconsl_outfile,"%d\t%d\t%s\t\'%d\'\n",wl_node->index,wl_node->paralog,"Transf",wl_node->sp_index);
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
                                fprintf(reconsl_outfile,"%d\t%d\t%s\t\'%s\'\n",wl_node->index,wl_node->paralog,"Sp",(names->names+(wl_node->sp_index*names->max_lname)));
                                break;
                            case DUP:
                                fprintf(reconsl_outfile,"%d\t%d\t%s\t\'%s\'\n",wl_node->index,wl_node->paralog,"Dup",(names->names+(wl_node->sp_index*names->max_lname)));
                                break;
                            case TRFR:
                                fprintf(reconsl_outfile,"%d\t%d\t%s\t\'%s\'\n",wl_node->index,wl_node->paralog,"Transf",(names->names+(wl_node->sp_index*names->max_lname)));
                                break;
                            default:
                                return UNEXPECTED_VALUE;
                                break;
                        }
                    }
                }
            }
            else if (wl_node->kind_node==LOSS)
            {
                
                fprintf(reconsl_outfile,"\'Lost–%d_%d\'\t%d\tLoss\t%d\n",wl_node->conts->index,wl_node->paralog,wl_node->paralog,wl_node->conts->index);
            }
            else if (wl_node->kind_node==RTRFR)
            {
                
                fprintf(reconsl_outfile,"\'RTransf–%d_%d\'\t%d\tRTransf\t%d\n",wl_node->conts->index,wl_node->paralog,wl_node->paralog,wl_node->conts->index);
            }
            else
            {
                
                if (names==NULL)
                    fprintf(reconsl_outfile,"\'%d_%d\'\t%d\tLeaf\t\'%d\'\n",wl_node->sp_index,wl_node->paralog,wl_node->paralog,wl_node->sp_index);
                else
                    fprintf(reconsl_outfile,"\'%s_%d\'\t%d\tLeaf\t\'%s\'\n",(names->names+(wl_node->sp_index*names->max_lname)),wl_node->paralog,wl_node->paralog,(names->names+(wl_node->sp_index*names->max_lname)));
            }
            
        }
        else
            return UNEXPECTED_VALUE;
        
    }
    
    fclose(reconsl_outfile);
    
    return(NO_ERROR);
}

long int WriteReconLG(g_tree *gene_tree, name_c *names, char *reconlg_outname)
{
    FILE *reconlg_outfile=NULL;
    g_node *wg_node=NULL, **gn_pointers=NULL;
    unsigned int i=0;
    
    
    gn_pointers=calloc(gene_tree->max_childs,sizeof(g_node *));
    
    CollapseGTree(gene_tree,1,1);
    
    if ((reconlg_outfile=fopen(reconlg_outname, "w"))==NULL)
    {
        return(IO_ERROR);
    }
    fprintf(reconlg_outfile,"Gt_node\tLt_node\tLt_paralog\n");
    
    for (i=0;i<gene_tree->n_nodes;++i)
    {
        wg_node=gene_tree->m_node+i;
        if (wg_node->contl!=NULL)
        {
            if (wg_node->n_child!=0) //g_node internal
            {
                
                if (wg_node->contl->n_child!=0) // Everything internal
                {
                    fprintf(reconlg_outfile,"%d\t%d\t%d\n",wg_node->index,wg_node->contl->index,wg_node->paralog);
                }
                else //External l_node
                {
                    if (wg_node->contl->kind_node==LOSS) //Lost l_node
                        fprintf(reconlg_outfile,"%d\t\'Lost–%d_%d\'\t%d\n",wg_node->index,wg_node->conts->index,wg_node->paralog,wg_node->paralog);
                    else if(wg_node->contl->kind_node==RTRFR) //Rtransfer l_node
                        fprintf(reconlg_outfile,"%d\t\'Rtransf–%d_%d\'\t%d\n",wg_node->index,wg_node->conts->index,wg_node->paralog,wg_node->paralog);
                    else
                    {
                        if (names==NULL)
                            fprintf(reconlg_outfile,"%d\t\'%d_%d\'\t%d\n",wg_node->index,wg_node->sp_index,wg_node->paralog,wg_node->paralog);
                        else
                            fprintf(reconlg_outfile,"%d\t\'%s_%d\'\t%d\n",wg_node->index,(names->names+(wg_node->sp_index*names->max_lname)),wg_node->paralog,wg_node->paralog);
                    }
                    
                }
                
            }
            else //g_node and l_node external
            {
                if (wg_node->contl->kind_node==LOSS)
                    fprintf(reconlg_outfile,"\'Lost–%d_%d_%d\'\t\'Lost–%d_%d\'\t%d\n",wg_node->conts->index,wg_node->paralog, wg_node->replica,wg_node->conts->index,wg_node->paralog,wg_node->paralog);
                else if (wg_node->contl->kind_node==RTRFR) //Rtransfer l_node
                    fprintf(reconlg_outfile,"\'Rtransf–%d_%d_%d\'\t\'Rtransf–%d_%d\'\t%d\n",wg_node->conts->index,wg_node->paralog, wg_node->replica,wg_node->conts->index,wg_node->paralog,wg_node->paralog);
                else //Everything external
                {
                    if (names==NULL)
                        fprintf(reconlg_outfile,"\'%d_%d_%d\'\t\'%d_%d\'\t%d\n",wg_node->sp_index,wg_node->paralog, wg_node->replica,wg_node->sp_index,wg_node->paralog,wg_node->paralog);
                    else
                        fprintf(reconlg_outfile,"\'%s_%d_%d\'\t\'%s_%d\'\t%d\n",(names->names+(wg_node->sp_index*names->max_lname)),wg_node->paralog,wg_node->replica,(names->names+(wg_node->sp_index*names->max_lname)),wg_node->paralog,wg_node->paralog);
                }
                
            }
            
        }
        else
            return (UNEXPECTED_VALUE);
        
    }
    
    fclose(reconlg_outfile);
    return (NO_ERROR);
}

//long int WriteReconLosses(s_tree *wsp_tree, g_tree *gene_tree, name_c *names, char *recon_outname)
//{
//    unsigned int i=0, j=0, is_dup=0, locus=0;
//    g_node *wg_node=NULL;
//    FILE *recon_outfile=NULL;
//    
//    if ((recon_outfile=fopen(recon_outname, "w"))==NULL)
//    {
//        return(IO_ERROR);
//    }
//    fprintf(recon_outfile,"Gt_node\tLocus\tSt_node\tEvent\n");
//    
//    
//    PostReorderSNodes(wsp_tree->root,&i);
//    CollapseGTree(gene_tree,1,1);
//    
//    for (i=0;i<gene_tree->n_nodes;++i)
//    {
//        wg_node=gene_tree->m_node+i;
//        if (wg_node->conts!=NULL)
//        {
//            if (wg_node->n_child!=0) //Internal g_node
//            {
//                is_dup=0;
//                locus=(*wg_node->childs)->paralog;
//                
//                for (j=1; j<wg_node->n_child; ++j)
//                {
//                    if (locus!=(*(wg_node->childs+j))->paralog)
//                    {
//                        is_dup=1;
//                        break;
//                    }
//                }
//                if (wg_node->conts->n_child==0) //External s_node
//                {
//                    if (names==NULL)
//                        fprintf(recon_outfile,"%d\t%d\t\'%d\'\t%s\n",wg_node->index,wg_node->paralog,wg_node->sp_index,is_dup==1?"Dup":"Coal");
//                    else
//                        fprintf(recon_outfile,"%d\t%d\t\'%s\'\t%s\n",wg_node->index,wg_node->paralog,(names->names+(wg_node->sp_index*names->max_lname)),is_dup==1?"Dup":"Coal");
//                }
//                else
//                    fprintf(recon_outfile,"%d\t%d\t%d\t%s\n",wg_node->index,wg_node->paralog,wg_node->conts->index,is_dup==1?"Dup":"Coal");
//            }
//            else //External g_node
//            {
//                if (wg_node->contl->kind_node==LOSS)
//                {
//                    if (wg_node->conts->n_child==0) //Everything external
//                    {
//                        if (names==NULL)
//                            fprintf(recon_outfile,"\'Lost–%d_%d_%d\'\t%d\t\'%d\'\tLoss\n",wg_node->conts->index,wg_node->paralog,wg_node->replica,wg_node->paralog,wg_node->sp_index);
//                        else
//                            fprintf(recon_outfile,"\'Lost–%s_%d_%d\'\t%d\t\'%s\'\tLoss\n",(names->names+(wg_node->sp_index*names->max_lname)),wg_node->paralog,wg_node->replica,wg_node->paralog,(names->names+(wg_node->sp_index*names->max_lname)));
//                    }
//                    else //Loss inside an internal species branch
//                    {
//                        fprintf(recon_outfile,"\'Lost–%d_%d_%d\'\t%d\t\%d\tLoss\n",wg_node->conts->index,wg_node->paralog,wg_node->replica,wg_node->paralog,wg_node->conts->index);
//                    }
//                    
//                }
//                else //Everything external
//                {
//                    if (names==NULL)
//                        fprintf(recon_outfile,"\'%d_%d_%d\'\t%d\t\'%d\'\tLeaf\n",wg_node->sp_index,wg_node->paralog, wg_node->replica,wg_node->paralog,wg_node->sp_index);
//                    else
//                        fprintf(recon_outfile,"\'%s_%d_%d\'\t%d\t\'%s\'\tLeaf\n",(names->names+(wg_node->sp_index*names->max_lname)),wg_node->paralog,wg_node->replica,wg_node->paralog,(names->names+(wg_node->sp_index*names->max_lname)));
//                }
//                
//            }
//            
//        }
//        else
//            return UNEXPECTED_VALUE;
//        
//    }
//    
//    fclose(recon_outfile);
//    
//    return(NO_ERROR);
//}

long int CheckNewickSTree (char * tree)
{
	int i=0;
    unsigned int k=0,error=0,n_left=0,n_right=0,n_branches=0,n_nodes=0,node=0,f_leaf=0,leaf=0;
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
                if (!(next == ')' || next== ',' || next== ':' || next=='#' || next=='/'|| next=='*')) /// ")" followed by something different to codes( ":","#","/", ")","*" or "," )
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
    unsigned int k=0,error=0,n_left=0,n_right=0,n_branches=0,n_nodes=0,node=0,f_leaf=0,leaf=0;
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
                if (!(next == ')' || next== ',' || next== ':' || next=='#' || next=='/' || next=='%' || next=='*')) /// ")" followed by something different to codes( ":","#","/", ")", "*", "%" or ",")
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

///@}

// **** Declarations of private functions **** //

// *** Recursive functions *** //

// ** Node deletion ** //

void FreeSNodes(s_node *node)
{
    unsigned int i=0;
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Post-order recursion (childs,anc)
    for (i=0;i<node->n_child;++i)
    {
        FreeSNodes(*(node->childs+i));
    }
    
    // *
    /// Frees memory (s_node::childs, s_node)</dd></dl>
    if (node->childs != NULL)
    {
        free(node->childs);
        node->childs=NULL;
    }
    free(node);
    node=NULL;
}

void FreeLNodes(l_node *node, unsigned int free_gnodes)
{
    unsigned int i=0;
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Post-order recursion (childs,anc)
    for (i=0;i<node->n_child;++i)
    {
        FreeLNodes(*(node->childs+i),free_gnodes);
    }
    
    // *
    /// Frees memory (l_node::g_nodes, l_node::childs, l_node)</dd></dl>
    if (node->g_nodes != NULL && free_gnodes!=0)
    {
        free(node->g_nodes);
        node->g_nodes=NULL;
    }
    if (node->childs != NULL)
    {
        free(node->childs);
        node->childs=NULL;
    }
    free(node);
    node=NULL;
}

void FreeGNodes(g_node *node, unsigned int free_root)
{
    unsigned int i=0;
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Post-order recursion (childs,anc)
    for (i=0;i<node->n_child;++i)
    {
        FreeGNodes(*(node->childs+i),free_root);
    }
    
    // *
    /// Frees memory (s_node::childs, s_node)</dd></dl>
    if (node->anc_node!=NULL || free_root==1)
    {
        if (node->childs != NULL)
        {
            free(node->childs);
            node->childs=NULL;
        }
        free(node);
        node=NULL;
    }
    
}

// ** Node reorganitation ** //

void PostReorderSNodes(s_node * node, unsigned int * index)
{
    unsigned int i=0;
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Post-order recursion (childs,anc)
    for (i=0;i<node->n_child;++i)
    {
        PostReorderSNodes(*(node->childs+i),index);
    }
    // *
    /// Reorders the index, increasing the pointed variable in each call, and using it as s_node::index.</dd></dl>
    node->index= *index;
    ++*index;
}

void PostReorderLNodes(l_node * node, unsigned int * index)
{
    unsigned int i=0;
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Post-order recursion (childs,anc)
    for (i=0;i<node->n_child;++i)
    {
        PostReorderLNodes(*(node->childs+i),index);
    }
    // *
    /// Reorders the index, increasing the pointed variable in each call, and using it as l_node::index.</dd></dl>
    node->index= *index;
    ++*index;
}

void PostReorderGNodes(g_node * node, unsigned int * index)
{
    unsigned int i=0;
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Post-order recursion (childs,anc)
    for (i=0;i<node->n_child;++i)
    {
        PostReorderGNodes(*(node->childs+i),index);
    }
    // *
    /// Reorders the index, increasing the pointed variable in each call, and using it as g_node::index.</dd></dl>
    node->index= *index;
    ++*index;
}

void PreReorderSNodes (s_node * node, unsigned int * index)
{
    unsigned int i=0;
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Reorders the index, increasing the pointed variable in each call, and using it as s_node::index.
    node->index= *index;
    ++*index;
    
    // *
    /// Pre-order recursion (anc,childs)</dd></dl>
    for (i=0;i<node->n_child;++i)
    {
        PreReorderSNodes(*(node->childs+i),index);
    }
    
}


void PreReorderLNodes (l_node * node, unsigned int * index)
{
    unsigned int i=0;
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Reorders the index, increasing the pointed variable in each call, and using it as l_node::index.
    node->index= *index;
    ++*index;
    
    // *
    /// Pre-order recursion (anc,childs)</dd></dl>
    for (i=0;i<node->n_child;++i)
    {
        PreReorderLNodes(*(node->childs+i),index);
    }
    
}

void PreReorderGNodes (g_node * node, unsigned int * index)
{
    unsigned int i=0;
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Reorders the index, increasing the pointed variable in each call, and using it as g_node::index.
    node->index= *index;
    ++*index;
    
    // *
    /// Pre-order recursion (anc,childs)</dd></dl>
    for (i=0;i<node->n_child;++i)
    {
        PreReorderGNodes(*(node->childs+i),index);
    }
    
}

/** Node copy **/

void PostCollapseSNodes(s_node *output,s_node *input)
{
    s_node * w_output=NULL;
    s_node ** c_backup=NULL;
    unsigned int i=0;
    
    // ***
    /// <dl><dt> Function structure </dt><dd>
    
    // **
    /// Post-order recursion (childs,anc)
    for (i=0;i<input->n_child;++i)
    {
        PostCollapseSNodes(output,*(input->childs+i));
    }
    
    // **
    /// Saves the original s_tree::childs of the output node
    w_output=(output+input->index);
    
    ErrorReporter((long int)w_output->childs);
    c_backup=w_output->childs;
    
    // **
    /// Copy of values
    memcpy(w_output,input,sizeof(s_node));
    
    // **
    ///<dl><dt>Copy of saved pointers</dt><dd>
    
    // *
    /// Childs loop
    w_output->childs=c_backup;
    for (i=0;i<input->n_child;++i)
    {
        *(w_output->childs+i)=(output+(*(input->childs+i))->index);
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

void PostCollapseLNodes(l_node *output,l_node *input,unsigned int n_gleaves, unsigned int retain_lateral)
{
    l_node * w_output=NULL;
    g_node ** g_backup=NULL;
    l_node ** r_backup=NULL;
    unsigned int i=0;
    
    // ***
    /// <dl><dt> Function structure </dt><dd>
    
    // **
    /// Post-order recursion (childs,anc)
    for (i=0;i<input->n_child;++i)
    {
        PostCollapseLNodes(output,*(input->childs+i),n_gleaves,retain_lateral);
    }
    
    // **
    /// Saves the original l_tree::g_nodes and l_tree::childs of the output node
    w_output=(output+input->index);
    
    ErrorReporter((long int)w_output->childs);
    r_backup=w_output->childs;
    if (w_output->g_nodes!=NULL)
    {
        g_backup=w_output->g_nodes;
    }
    
    // **
    /// Copy of values
    memcpy(w_output,input,sizeof(l_node));
    if (retain_lateral!=1)
    {
        w_output->lat_node=NULL;
    }
    
    // **
    ///<dl><dt>Copy of saved pointers</dt><dd>
    
    // *
    /// Childs loop
    w_output->childs=r_backup;
    for (i=0;i<input->n_child;++i)
    {
        *(w_output->childs+i)=(output+(*(input->childs+i))->index);
    }
    
    // *
    /// g_nodes </dd></dl>
    if (g_backup != NULL)
    {
        w_output->g_nodes=g_backup;
        
        if (input->g_nodes!=NULL)
        {
            memcpy(w_output->g_nodes,input->g_nodes,sizeof(g_node*)*n_gleaves);
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
    unsigned int i=0;
    
    // ***
    /// <dl><dt> Function structure </dt><dd>
    
    // **
    /// Post-order recursion (childs,anc)
    for (i=0;i<input->n_child;++i)
    {
        PostCollapseGNodes(output,*(input->childs+i));
    }
    
    // **
    /// Saves the original g_tree::childs of the output node
    w_output=(output+input->index);
    
    ErrorReporter((long int)w_output->childs);
    c_backup=w_output->childs;
    
    // **
    /// Copy of values
    memcpy(w_output,input,sizeof(g_node));
    
    // **
    ///<dl><dt>Copy of saved pointers</dt><dd>
    
    // *
    /// Childs loop
    w_output->childs=c_backup;
    for (i=0;i<input->n_child;++i)
    {
        *(w_output->childs+i)=(output+(*(input->childs+i))->index);
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
    unsigned int i=0;
    
    // ***
    /// <dl><dt> Function structure </dt><dd>
    
    // **
    /// Saves the original s_tree::childs of the output node
    w_output=(output+input->index);
    
    ErrorReporter((long int)w_output->childs);
    c_backup=w_output->childs;
    
    // **
    /// Copy of values
    memcpy(w_output,input,sizeof(s_node));
    
    // **
    ///<dl><dt>Copy of saved pointers</dt><dd>
    
    // *
    /// Childs loop
    w_output->childs=c_backup;
    for (i=0;i<input->n_child;++i)
    {
        *(w_output->childs+i)=(output+(*(input->childs+i))->index);
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
    /// Pre-order recursion (anc,childs) </dd></dl></dd></dl>
    for (i=0;i<input->n_child;++i)
    {
        PreCollapseSNodes(output,*(input->childs+i));
    }
    
    
}

void PreCollapseLNodes(l_node *output,l_node *input,unsigned int n_gleaves, unsigned int retain_lateral)
{
    l_node * w_output=NULL;
    g_node ** g_backup=NULL;
    l_node ** r_backup=NULL;
    unsigned int i=0;
    
    // ***
    /// <dl><dt> Function structure </dt><dd>
    
    // **
    /// Saves the original l_tree::g_nodes and l_tree::childs of the output node
    w_output=(output+input->index);
    
    ErrorReporter((long int)w_output->childs);
    r_backup=w_output->childs;
    if (w_output->g_nodes!=NULL)
    {
        g_backup=w_output->g_nodes;
    }
    
    // **
    /// Copy of values
    memcpy(w_output,input,sizeof(l_node));
    if (retain_lateral!=1)
    {
        w_output->lat_node=NULL;
    }
    // **
    ///<dl><dt>Copy of saved pointers</dt><dd>
    
    // *
    /// Childs loop
    w_output->childs=r_backup;
    for (i=0;i<input->n_child;++i)
    {
        *(w_output->childs+i)=(output+(*(input->childs+i))->index);
    }
    
    // *
    /// g_nodes </dd></dl>
    if (g_backup != NULL)
    {
        w_output->g_nodes=g_backup;
        if (input->g_nodes!=NULL)
        {
            memcpy(w_output->g_nodes,input->g_nodes,sizeof(g_node*)*n_gleaves);
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
    /// Pre-order recursion (anc,childs)  </dd></dl>
    for (i=0;i<input->n_child;++i)
    {
        PreCollapseLNodes(output,*(input->childs+i),n_gleaves,retain_lateral);
    }
    
    
}

void PreCollapseGNodes(g_node *output,g_node *input)
{
    g_node * w_output=NULL;
    g_node ** c_backup=NULL;
    unsigned int i=0;
    
    // ***
    /// <dl><dt> Function structure </dt><dd>
    
    // **
    /// Saves the original g_tree::childs of the output node
    w_output=(output+input->index);
    
    ErrorReporter((long int)w_output->childs);
    c_backup=w_output->childs;
    
    // **
    /// Copy of values
    memcpy(w_output,input,sizeof(s_node));
    
    // **
    ///<dl><dt>Copy of saved pointers</dt><dd>
    
    // *
    /// Childs loop
    w_output->childs=c_backup;
    for (i=0;i<input->n_child;++i)
    {
        *(w_output->childs+i)=(output+(*(input->childs+i))->index);
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
    /// Pre-order recursion (anc,childs) </dd></dl></dd></dl>
    for (i=0;i<input->n_child;++i)
    {
        PreCollapseGNodes(output,*(input->childs+i));
    }
    
}

// ** Node features obtention ** //
static unsigned int Count_duplications_lnodes(l_node *node)
{
    unsigned int i=0;
    unsigned int n_dup=0;
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Post-order recursion (childs,anc)
    for (i=0;i<node->n_child;++i)
    {
        n_dup+=Count_duplications_lnodes(*(node->childs+i));
    }
    
    // *
    /// Counting </dd></dl>
    if (node->paralog==DUP)//Duplication
    {
        ++n_dup;
    }
    return n_dup;

}

static unsigned int Count_losses_lnodes(l_node *node)
{
    unsigned int i=0;
    unsigned int n_losses=0;
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Post-order recursion (childs,anc)
    for (i=0;i<node->n_child;++i)
    {
        n_losses+=Count_losses_lnodes(*(node->childs+i));
    }
    
    // *
    /// Counting </dd></dl>
    if (node->paralog==LOSS)//Loss
    {
        ++n_losses;
    }
    
    return n_losses;
    
}


static long double Measure_g_node_cu_height(g_node *node, unsigned int g_Ne)
{
    unsigned int i=0;
    long double current=0, new=0;
    
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Post-order recursion selecting the longest accumulated branch length
    for (i=0;i<node->n_child;++i)
    {
        new=Measure_g_node_cu_height(*(node->childs+i),g_Ne);
        
        if (new>current)
            current=new;
    }
    
    // **
    /// Returns the biggest branch lenght added to the own one.</dd></dl>
    return node->gen_length/(node->contl->Ne!=0?node->contl->Ne:g_Ne)+current;

}

static long double Measure_g_node_bl_height(g_node *node)
{
    unsigned int i=0;
    long double current=0, new=0;
    
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Post-order recursion selecting the longest accumulated branch length
    for (i=0;i<node->n_child;++i)
    {
        new=Measure_g_node_bl_height(*(node->childs+i));
        
        if (new>current)
            current=new;
    }
    
    // **
    /// Returns the biggest branch lenght added to the own one.</dd></dl>
    return current+node->bl;
}

static long double Measure_g_node_cu_length(g_node *node, unsigned int g_Ne)
{
    unsigned int i=0;
    long double sum=0;
    
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Post-order recursion acumulating branch length
    for (i=0;i<node->n_child;++i)
    {
        sum+=Measure_g_node_cu_height(*(node->childs+i),g_Ne);
    }
    
    // **
    /// Returns the biggest branch lenght added to the own one.</dd></dl>
    return node->gen_length/(node->contl->Ne!=0?node->contl->Ne:g_Ne)+sum;
    
}

static long double Measure_g_node_bl_length(g_node *node)
{
    unsigned int i=0;
    long double sum=0;
    
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Post-order recursion selecting the longest accumulated branch length
    for (i=0;i<node->n_child;++i)
    {
        sum+=Measure_g_node_bl_height(*(node->childs+i));

    }
    
    // **
    /// Returns the biggest branch lenght added to the own one.</dd></dl>
    return sum+node->bl;
}

static long double Measure_s_node_cu_height(s_node *node, unsigned int g_Ne)
{
    unsigned int i=0;
    long double current=0, new=0;
    
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Post-order recursion selecting the longest accumulated branch length
    for (i=0;i<node->n_child;++i)
    {
        new=Measure_s_node_cu_height(*(node->childs+i),g_Ne);
        
        if (new>current)
            current=new;
    }
    
    // **
    /// Returns the biggest branch lenght added to the own one.</dd></dl>
    return node->gen_length/(node->Ne!=0?node->Ne:g_Ne)+current;
    
}

static long double Measure_s_node_gl_height(s_node *node)
{
    unsigned int i=0;
    long double current=0, new=0;
    
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Post-order recursion selecting the longest accumulated branch length
    for (i=0;i<node->n_child;++i)
    {
        new=Measure_s_node_gl_height(*(node->childs+i));
        
        if (new>current)
            current=new;
    }
    
    // **
    /// Returns the biggest branch lenght added to the own one.</dd></dl>
    return current+node->gen_length;
}

static long double Measure_s_node_cu_length(s_node *node, unsigned int g_Ne)
{
    unsigned int i=0;
    long double sum=0;
    
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Post-order recursion acumulating branch length
    for (i=0;i<node->n_child;++i)
    {
        sum+=Measure_s_node_cu_height(*(node->childs+i),g_Ne);
    }
    
    // **
    /// Returns the biggest branch lenght added to the own one.</dd></dl>
    return node->gen_length/(node->Ne!=0?node->Ne:g_Ne)+sum;
    
}

static long double Measure_s_node_gl_length(s_node *node)
{
    unsigned int i=0;
    long double sum=0;
    
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Post-order recursion selecting the longest accumulated branch length
    for (i=0;i<node->n_child;++i)
    {
        sum+=Measure_s_node_gl_height(*(node->childs+i));
        
    }
    
    // **
    /// Returns the biggest branch lenght added to the own one.</dd></dl>
    return sum+node->gen_length;
}



// ** Node modification ** //

static void RefineSNodes(s_node * node, unsigned int ind_persp)
{
    unsigned int i=0;
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Post-order recursion (childs,anc)
    for (i=0;i<node->n_child;++i)
    {
        RefineSNodes(*(node->childs+i),ind_persp);
    }
    
    // *
    /// Sets the correct s_node::n_nodes for leaf nodes or calculates and applies the s_node::n_gen for internal nodes </dd></dl>
    if (node->n_child==0)//Leaf
    {
        if (node->n_replicas==0)
            node->n_replicas=ind_persp;
    }
    else
    {
        node->n_gen=(*(node->childs))->gen_length+(*(node->childs))->n_gen;
        if (node->n_replicas>0)
        {
            fprintf(stderr,"Number of replicas per taxa applied to internal branches are not allowed\n");
            ErrorReporter(SETTINGS_ERROR);
        }
    }
    
}

static void RefineLNodes(l_node * node, unsigned int n_gleaves, unsigned int ind_persp)
{
    unsigned int i=0,j=0;
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Post-order recursion (childs,anc)
    for (i=0;i<node->n_child;++i)
    {
        RefineLNodes(*(node->childs+i),n_gleaves, ind_persp);
    }
    
    // *
    /// Sets the correct l_node::n_nodes for leaf nodes or calculates and applies the l_node::n_gen for internal nodes
    if (node->n_child==0)//Leaf
    {
        if (node->n_nodes==0)
            node->n_nodes=ind_persp;
    }
    else
    {
        node->n_gen=(*(node->childs))->gen_length+(*(node->childs))->n_gen;
    }
    
    // *
    /// Allocates l_node::g_nodes if it is necesary, and initializes them.</dd></dl>
    if (node->g_nodes==NULL)
    {
        node->g_nodes=calloc(n_gleaves,sizeof(g_node *));
        
        for (j=0; j<n_gleaves; ++j)
        {
            *(node->g_nodes+j)=NULL;
        }
        
    }
    
}

static void CleanlossesLNodes(l_node * node, l_node ** root, unsigned int *n_deletions, unsigned int *n_leaves)
{
    unsigned int i=0,j=0,n_child=node->n_child;
    l_node *w_lnode=NULL;
    
    // ***
    /// <dl><dt> Function structure </dt><dd>
    
    // **
    /// Post-order recursion (childs,anc)
    for (i=0;i<n_child;++i)
    {
        CleanlossesLNodes(*(node->childs+i),root,n_deletions,n_leaves);
    }
    
    // **
    ///<dl><dt>No proper internal node</dt><dd>
    if(node->n_child<2)
    {
        // *
        // If real leave, count it.
        if(node->n_gen==0)
        {
            ++(*n_leaves);
        }
        else if (node->anc_node!=NULL)
        {
            w_lnode=node->anc_node;
            // *
            /// Deletion if the node has 0 childs and is not a real leaf (n_gen=0).
            if (node->n_child==0)
            {
                for (j=0;j<w_lnode->n_child;++j)
                {
                    if (*(w_lnode->childs+j)==node)
                        break;
                }
                *(w_lnode->childs+j)=*(w_lnode->childs+w_lnode->n_child-1);
                --w_lnode->n_child;
                FreeLNodes(node, 1);
                ++(*n_deletions);
            }
            // *
            /// Deletion and "by pass" if the node has 1 child
            else if (node->n_child==1)
            {
                for (j=0;j<w_lnode->n_child;++j)
                {
                    if (*(w_lnode->childs+j)==node)
                        break;
                }
                *(w_lnode->childs+j)=*node->childs;
                w_lnode=*node->childs;
                w_lnode->anc_node=node->anc_node;
                w_lnode->gen_length=node->anc_node->n_gen-w_lnode->n_gen;
                node->n_child=0; //To delete only this node
                FreeLNodes(node, 1);
                ++(*n_deletions);
            }
        }
        // *
        ///Change of the root of the tree if it should be deleted</dd></dl></dd></dl>
        else
        {
            w_lnode=*node->childs;
            *root=w_lnode;
            (*root)->anc_node=NULL;
            node->n_child=0; //To delete only this node
            FreeLNodes(node, 1);
            ++(*n_deletions);
        }
        
    }
    
}

static void Temporalize_GNodes(g_node * w_node,double gen_time)
{
    unsigned int i=0;
    
    w_node->bl*=gen_time;
    w_node->gen_length*=gen_time;
    w_node->height*=gen_time;
    for(i=0; i<w_node->n_child; ++i)
    {
        Temporalize_GNodes(*(w_node->childs+i), gen_time);
    }
}

// *** Tree I/O *** //

void WriteSNodes (s_node * p, name_c * names, double gen_time)
{
    unsigned int i=0;
    
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
                printf("%d:%.8lf",p->sp_index,p->gen_length*gen_time);
            else
                printf("%s:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->gen_length*gen_time);
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
                WriteSNodes(*(p->childs+i),names,gen_time);
                printf(",");
            }
			WriteSNodes (*(p->childs+i), names,gen_time);
            
            // **
            /// Prints the information of this internal node (closing it with ")") or finishes the tree if the node is the root</dd></dl></dd></dl>
			if (p->anc_node !=NULL)
				printf("):%.8lf",p->gen_length*gen_time);
            else
                printf(");");
        }
    }
	
}

void WriteSNodesFile (FILE * file,s_node * p, name_c * names, double gen_time)
{
    unsigned int i=0;
    
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
                fprintf(file,"%d:%.8lf",p->sp_index,p->gen_length*gen_time);
            else
                fprintf(file,"%s:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->gen_length*gen_time);
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
                WriteSNodesFile(file,*(p->childs+i),names,gen_time);
                fprintf(file,",");
            }
			WriteSNodesFile (file,*(p->childs+i), names,gen_time);
            
            // **
            /// Prints the information of this internal node (closing it with ")") or finishes the tree if the node is the root</dd></dl></dd></dl>
			if (p->anc_node !=NULL)
				fprintf(file,"):%.8lf",p->gen_length*gen_time);
            else
                fprintf(file,");\n");
        }
    }
	
}

void WriteLNodes (l_node * p, name_c * names, double gen_time)
{
    unsigned int i=0;
    
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
            if(p->kind_node!=LOSS && p->kind_node!=RTRFR)
            {
                if (names==NULL)
                    printf("%d_%d:%.8lf",p->sp_index,p->paralog,p->gen_length*gen_time);
                else
                    printf("%s_%d:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->gen_length*gen_time);
            }
            else if (p->kind_node==RTRFR)
            {
                if (p->conts->n_child==0 && names!=NULL)
                {
                    printf("Rtransf–%s_%d:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->gen_length*gen_time);
                }
                else
                {
                    printf("Rtransf–%d_%d:%.8lf",p->conts->index,p->paralog,p->gen_length*gen_time);
                }
            }
            else
            {
                if (p->conts->n_child==0 && names!=NULL)
                {
                    printf("Lost–%s_%d:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->gen_length*gen_time);
                }
                else
                {
                    printf("Lost–%d_%d:%.8lf",p->conts->index,p->paralog,p->gen_length*gen_time);
                }
            }

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
                WriteLNodes(*(p->childs+i),names,gen_time);
                printf(",");
            }
			WriteLNodes (*(p->childs+i), names,gen_time);
            
            // **
            /// Prints the information of this internal node (closing it with ")") or finishes the tree if the node is the root</dd></dl></dd></dl>
			if (p->anc_node !=NULL)
				printf("):%.8lf",p->gen_length*gen_time);
            else
                printf(");");
        }
    }
	
}

void WriteLNodesFile (FILE * file,l_node * p, name_c * names, double gen_time)
{
    unsigned int i=0;
    
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
            if(p->kind_node!=LOSS && p->kind_node!=RTRFR)
            {
                if (names==NULL)
                    fprintf(file,"%d_%d:%.8lf",p->sp_index,p->paralog,p->gen_length*gen_time);
                else
                    fprintf(file,"%s_%d:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->gen_length*gen_time);
            }
            else if (p->kind_node==RTRFR)
            {
                if (p->conts->n_child==0 && names!=NULL)
                {
                    fprintf(file,"Rtransf–%s_%d:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->gen_length*gen_time);
                }
                else
                {
                    fprintf(file,"Rtransf–%d_%d:%.8lf",p->conts->index,p->paralog,p->gen_length*gen_time);
                }
            }
            else
            {
                if (p->conts->n_child==0 && names!=NULL)
                {
                    fprintf(file,"Lost–%s_%d:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->gen_length*gen_time);
                }
                else
                {
                    fprintf(file,"Lost–%d_%d:%.8lf",p->conts->index,p->paralog,p->gen_length*gen_time);
                }
            }

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
                WriteLNodesFile(file,*(p->childs+i),names,gen_time);
                fprintf(file,",");
            }
			WriteLNodesFile (file,*(p->childs+i), names,gen_time);
            
            // **
            /// Prints the information of this internal node (closing it with ")") or finishes the tree if the node is the root</dd></dl></dd></dl>
			if (p->anc_node !=NULL)
				fprintf(file,"):%.8lf",p->gen_length*gen_time);
            else
                fprintf(file,");\n");
        }
    }
	
}

void WriteGNodes (g_node * p, name_c * names)
{
    unsigned int i=0;
    
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
            if (p->contl->kind_node!=LOSS && p->contl->kind_node!=RTRFR)
            {
                if (names==NULL)
                    printf("%d_%d_%d:%.8lf",p->sp_index,p->paralog, p->replica,p->bl);
                
                else
                    printf("%s_%d_%d:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->replica,p->bl);
            }
            else if (p->contl->kind_node==RTRFR)
            {
                if (p->conts->n_child==0 && names!=NULL)
                {
                    printf("Rtransf–%s_%d_0:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->bl);
                }
                else
                {
                    printf("Rtransf–%d_%d_0:%.8lf",p->conts->index,p->paralog,p->bl);
                }
            }
            else
            {
                if (p->conts->n_child==0 && names!=NULL)
                {
                    printf("Lost–%s_%d_0:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->bl);
                }
                else
                {
                    printf("Lost–%d_%d_0:%.8lf",p->conts->index,p->paralog,p->bl);
                }
            }

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
                WriteGNodes(*(p->childs+i),names);
                printf(",");
            }
            WriteGNodes(*(p->childs+i),names);
            
            // **
            /// Prints the information of this internal node (closing it with ")") or finishes the tree if the node is the root</dd></dl></dd></dl>
			if (p->anc_node !=NULL)
				printf("):%.8lf",p->bl);
            
            else
            {
                //printf("):%.8lf);",p->gen_length);//To print the length of the p
                printf(");");
            }
        }
    }
    
}

void WriteGNodesStr (char * str, g_node * p, name_c * names)
{
    unsigned int i=0;
    
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
            if (p->contl->kind_node!=LOSS && p->contl->kind_node==RTRFR)
            {
                if (names==NULL)
                    sprintf(str,"%s%d_%d_%d:%.8lf",str,p->sp_index,p->paralog,p->replica,p->bl);
                
                else
                    sprintf(str,"%s%s_%d_%d:%.8lf",str,(names->names+(p->sp_index*names->max_lname)),p->paralog,p->replica,p->bl);
            }
            else if (p->contl->kind_node==RTRFR)
            {
                if (p->conts->n_child==0 && names!=NULL)
                {
                    sprintf(str,"%sRtransf–%s_%d_0:%.8lf",str,(names->names+(p->sp_index*names->max_lname)),p->paralog,p->bl);
                }
                else
                {
                    sprintf(str,"%sRtransf–%d_%d_0:%.8lf",str,p->conts->index,p->paralog,p->bl);
                }
            }
            else
            {
                if (p->conts->n_child==0 && names!=NULL)
                {
                    sprintf(str,"%sLost–%s_%d_0:%.8lf",str,(names->names+(p->sp_index*names->max_lname)),p->paralog,p->bl);
                }
                else
                {
                    sprintf(str,"%sLost–%d_%d_0:%.8lf",str,p->conts->index,p->paralog,p->bl);
                }
            }
        }
		else
        {
            // **
            /// <dl><dt>Else (internal node):</dt><dd>
            // *
            /// Prints "(" and starts the post-order recursion (child loop)
            sprintf(str,"%s(",str);
            for (i=0;i<p->n_child-1;++i)
            {
                // *
                /// Calls itself using each different child. After each call (except last child) prints ","
                WriteGNodesStr(str,*(p->childs+i),names);
                printf(",");
                
            }
			WriteGNodesStr(str,*(p->childs+i),names);
            
            // **
            /// Prints the information of this internal node (closing it with ")") or finishes the tree if the node is the root</dd></dl></dd></dl>
			if (p->anc_node !=NULL)
				sprintf(str,"%s):%.8lf",str, p->bl);
            
            else
            {
                //                sprintf(str,"(%s):%.8lf);",str, p->gen_length);//To print the length of the p
                sprintf(str,"%s;",str);
            }
        }
    }
	
}

void WriteGNodesFile (FILE * file, g_node * p, name_c * names)
{
    unsigned int i=0;
    
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
            if (p->contl->kind_node!=LOSS && p->contl->kind_node!=RTRFR)
            {
                if (names==NULL)
                    fprintf(file,"%d_%d_%d:%.8lf",p->sp_index,p->paralog, p->replica,p->bl);
                
                else
                    fprintf(file,"%s_%d_%d:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->replica,p->bl);
            }
            else if (p->contl->kind_node==RTRFR)
            {
                if (p->conts->n_child==0 && names!=NULL)
                {
                    fprintf(file,"Rtransf–%s_%d_0:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->bl);
                }
                else
                {
                    fprintf(file,"Rtransf–%d_%d_0:%.8lf",p->conts->index,p->paralog,p->bl);
                }
            }
            else
            {
                if (p->conts->n_child==0 && names!=NULL)
                {
                    fprintf(file,"Lost–%s_%d_0:%.8lf",(names->names+(p->sp_index*names->max_lname)),p->paralog,p->bl);
                }
                else
                {
                    fprintf(file,"Lost–%d_%d_0:%.8lf",p->conts->index,p->paralog,p->bl);
                }
            }
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
                WriteGNodesFile(file, *(p->childs+i), names);
                fprintf(file,",");
            }
            WriteGNodesFile(file, *(p->childs+i), names);
            
            // **
            /// Prints the information of this internal node (closing it with ")") or finishes the tree if the node is the root</dd></dl></dd></dl>
			if (p->anc_node !=NULL)
				fprintf(file,"):%.8lf",p->bl);
            
            else
            {
                //                fprintf(file,"):%.8lf);\n",p->gen_length); //To print the length of the p
                fprintf(file,");\n");
            }
        }
    }
	
}

long int RelinkLGTrees(l_tree * wl_tree, g_node * m_node)
{
    unsigned int i=0, j=0;
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
        for (j=0; j<w_lnode->n_nodes;++j)
        {
            *(w_lnode->g_nodes+j)=(m_node+(*(w_lnode->g_nodes+j))->index);
        }
    }
    return NO_ERROR;
}

unsigned int ChooseLNodePeriod(period *w_period, l_node * t_node, double u_num)
{
    unsigned int i=0, j=0;
    l_node **l_pointers=calloc(w_period->n_lnodes,sizeof(l_node*));
    double *t_times=calloc(w_period->n_lnodes,sizeof(double));
    double t_prob=0, c_prob=0;
    l_node *w_lnode=NULL;
    l_node *w_tnode=t_node;
    unsigned int done=0, next_t=1;
    
    for (j=0; j<w_period->n_lnodes;++j)
    {
        *(l_pointers+j)=*(w_period->l_nodes+j);
    }
    while (i<MAX_IT && done==0)
    {
        done=1;
        next_t=1;
        for (j=0; j<w_period->n_lnodes;++j)
        {
            if (*(t_times+j)!=0)
                continue;
            w_lnode=*(l_pointers+j);
            if (w_lnode!=NULL)
            {
                done=0;
                if ((w_lnode->n_gen<w_tnode->n_gen) && (w_lnode != w_tnode))
                {
                    *(l_pointers+j)=w_lnode->anc_node;
                    next_t=0;
                }
                else if (w_lnode==w_tnode)
                {
                    *(t_times+j)=w_lnode->n_gen;
                }
            }
        }
        if (next_t==1 && done==0)
        {
            w_tnode=w_tnode->anc_node;
        }
        
    }
    
    for (j=0; j<w_period->n_lnodes;++j)
    {
        if (*(t_times+j)-t_node->n_gen<DBL_EPSILON)
            *(t_times+j)=0;
        else
        {
            *(t_times+j)=1/(*(t_times+j)-t_node->n_gen);
            t_prob+=*(t_times+j);
        }
    }
    
    w_lnode=NULL;
    
    done=0;
    for (j=0; j<w_period->n_lnodes;++j)
    {
        if (u_num>c_prob+(*(t_times+j)/t_prob))
        {
            c_prob+=(*(t_times+j)/t_prob);
        }
        else
        {
            w_lnode=*(w_period->l_nodes+j);
            done=1;
            break;
        }
    }
    
    free(l_pointers);
    free(t_times);
    
    if (done==0)
        ErrorReporter(UNEXPECTED_VALUE);
    
    return j;
}