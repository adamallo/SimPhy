
/**
 * 
 * \file trees.h
 * Header of the library of tree-manipulation functions (and tree types definitions).
 *******************************************************************************/ 

/**
 * \defgroup trees_pub Trees managing module
 *
 * Module with variables and functions to work with phylogenetic trees in a gene
 * tree / species tree context, mainly in simulation programs.
 *******************************************************************************/ 
 ///@{

#ifndef trees_h
#define trees_h
#include <stdio.h>
#ifndef stdlib_h
#include <stdlib.h>
#endif
#include <memory.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <ctype.h>
#ifndef num_methods_h
#include "num_methods.h"
#endif
#ifndef __GSL_RNG_H
#include <gsl/gsl_rng.h>
#endif
#ifndef __GSL_RANDIST_H
#include <gsl/gsl_randist.h>
#endif
#ifndef _STDARG_H
#include <stdarg.h>
#endif
#ifndef sampling_h
#include "sampling.h"
#endif
#ifndef string_h
#include <string.h>
#endif

//#define DBG
//\todo remove for final release

/**
 * Overestimation of the memory required to store one gene tree in Newick format.
 * Based on this formula \f[ ((max\_taxa\_name\_length \ast number\_of\_leaves+((number\_of\_leaves \ast 2)-1) \ast (max\_double\_true\_digits+4)+1) \f]
 *******************************************************************************/ 

#define STREE_MAXMEM ((max_lname*n_gleaves+((n_gleaves*2)-1)*(DBL_DIG+4))+1)  //Str memory overstimation. Max length of the name * n_names + n_nodes* max_blength_ndigits+2(separators) + ; * number of replicas (+ \0).

extern int MAX_IT;
extern int MAX_NAME;
extern int MAX_CHILDREN;
extern int MAX_LEAVES;
extern int NUM_BUFFER;
extern char TEST_CHAR;
extern int IO_BUFFER;
extern double FLOAT_PRECISION;

/**
 * \enum ERRORS
 * Error control constants
 *
 * This constants are used for \ref ErrorReporter. They are related to the
 * value of NULL in each implementation.
 *
 * \var NO_ERROR
 *  There is no error.
 * \var MEM_ERROR
 *  NULL value. It is a dynamic memory error (calloc returning NULL pointer)
 * \var IO_ERROR
 *  Problem with file I/O.
 * \var SETTINGS_ERROR
 *  Error with the settings parsing (values and newick tree).
 * \var LOOP_ERROR
 *  Infinite LOOP
 * \var UNEXPECTED_VALUE
 *  There is something wrong in any algorithm and/or math operation.
 *******************************************************************************/ 
enum ERRORS {NO_ERROR=(long int)NULL+1,MEM_ERROR=(long int) NULL,IO_ERROR=(long int)NULL-1,SETTINGS_ERROR=(long int)NULL+2, LOOP_ERROR=(long int)NULL+3, UNEXPECTED_VALUE=(long int)NULL+4, DB_ERROR=(long int)NULL+5, TERMINATE_NOERROR=(long int) NULL+6};



/**
 * \enum LIMITS
 * Limit constants
 *
 * This constants are used as limits of diferent processes.
 *
 * \var MAX_IT
 *  Limit to avoid infinite loops.
 * \var MAX_NAME
 *  Length limit of a taxa name. It is used to read a s_tree from Newick
 *  in \ref ReadNewickSTree. After this, the memory is reallocated in a more 
 *  efficient way in \ref ReallocNames.
 * \var MAX_CHILDREN
 *  Max number of children of a node. It is used to read a s_tree from 
 *  Newick in \ref ReadNewickSTree. After this, the memory is reallocated in a 
 *  more efficient way in \ref CollapseSTree, \ref CollapseLTree and \todo complete.
 * \var MAX_LEAVES
 *  Max number of leaves in a birth-death process. It is used to construct the
 *  locus tree from the species tree under a birth-death process.
 *******************************************************************************/
//enum LIMITS{MAX_IT=10000, MAX_NAME=100, MAX_CHILDREN=100, MAX_LEAVES=1000};

/**
 * \enum NODE_VALUES
 * Constants
 *
 * This constants are a kind of text codification of interesting int values.
 *
 * \var DUP
 *  Value of l_node::kind_node which identifies the node as a duplication.
 * \var LOSS
 *  Value of l_node::kind_node which identifies the node as a loss.
 * \var SP
 *  Value of l_node::kind_node which identifies the node as a speciation.
 * \var TRFR
 *  Value of l_node::kind_node which identifies the node as a transfer.
 * \var SP
 *  Value of l_node::kind_node which identifies the node as a transfer reception (loss).
 *******************************************************************************/
enum NODE_VALUES{SP=0,DUP=1,LOSS=2,TRFR=3,RTRFR=4,GC=5,RGC=6};

/**
 * \enum UNITS
 * Constants
 *
 * This constants are a kind of text codification of interesting units
 *
 * \var CU
 *  Coalescent units.
 * \var BL
 *  Branch length (substitutions per site)
 * \var GL
 *  Number of generations.
 *******************************************************************************/
enum UNITS{CU=0,BL=1,GL=2,TL=3};

/**
 * \enum OFORMAT
 * Output formats
 *
 * \var NEWICK
 *  Regular newick format.
 * \var NEXUS
 *  Nexus formats with branch lengths in generation units and extra information coded as comments.
 *******************************************************************************/
enum OFORMAT{NEWICK=0,NEXUS=1};

// *********************** Declaration of custom types ************************** //

typedef struct g_node g_node;
typedef struct g_tree g_tree;
typedef struct l_node l_node;
typedef struct l_tree l_tree;
typedef struct s_node s_node;
typedef struct s_tree s_tree;
typedef struct name_c name_c;
typedef struct period period;


/**
 * Gene tree node.
 *
 * These nodes are allocated as an array to form a g_tree, and they also conserve the identity of
 * his containers (l_node and/or g_node). They allow politomy.
 * Taxa names are stored in a name_c variable.
 *******************************************************************************/ 
struct g_node
{
    /** \name Codes **/
    /// @{
    int index;///< Memory position in the array.
    int sp_index; ///< Position of the name of this node in name_c
    int replica; ///< Code for the replica number of the node (Applicable in species trees with more than one replica in any taxa)
    int paralog; ///< Number to identify one "paralog"
    /// @}
    
    /** \name Node info **/
    /// @{
    int n_child; ///< Number of children (offspring) of this node.
    /// @}
    
    /** \name Data **/
    /// @{
    double n_gen; ///< Number of generations of this node (from the root).
    double gen_length;///< Number of generations from this node to the ancestor.
    double bl; ///< Branch length (substitutions per site per generation)
    /// @}
    
    /** \name Node pointers
     * Pointers to the offspring, ancestor, s_node and l_node.
     **/
    /// @{
    g_node **children; ///< Pointer to the array of children (pointers).
    g_node *anc_node; ///< Ancestor
    l_node *contl; ///< l_node that contains this g_node along his branch.
    s_node *conts; ///< s_node that contains this g_node along his branch.
    /// @}
    
};

/**
 * Gene tree.
 *
 * This structure holds a gene tree (g_node array) and some other relevant info of it.
 *******************************************************************************/
struct g_tree
{
    int n_nodes;///< Number of nodes
    int max_children; ///< Max number of children in each node. This variable is used to know (and/or modify) the length of the g_node::children allocated memory.
    double gen_time;///< Global generation time
    g_node * m_node; ///< Memory block of nodes
    g_node * root; ///< Pointer to the root of the tree
    l_tree *locus_tree; ///< Pointer to the locus tree.
    s_tree *species_tree; ///< Pointer to the species tree.

};

/**
 * Locus tree node.
 *
 * These nodes are allocated as an array or a simple tree to form a l_tree.
 * These nodes hold time restrictions (\ref l_node::gen_length) and combinatory restrictions (pointers to avaliable
 * gene nodes, \ref l_node::g_nodes "**g_nodes"). 
 * Polytomies are allowed, since nodes have an \ref children/ 1-\ref l_node::anc_node "ancestor" structure.
 * Taxa names are stored in a name_c variable.
 *******************************************************************************/ 
struct l_node 
{
    /** \name Codes **/
    /// @{
    int index; ///< Memory position in the array.
    int sp_index; ///< Position of the name of this node in name_c
    int kind_node; ///< Event tag of the node (SP, DUP, LOSS, TRFR, RTRFR). See \ref VALUES
    int paralog; ///< Number to identify one "paralog".
    /// @}
    
    /** \name Node info **/
    /// @{
    int n_child; ///< Number of children (offspring) of this node.
    int n_nodes; ///< Number of g_nodes wich are being pointed by this restriction tree (avaliable nodes).
    int n_ilin; ///< Number of input (backwards in time) lineages.
    int n_olin; ///< Number of output (backwards in time) lineages.
    int fmax_nlin; ///< Number of maximum gene lineages coming in and out from this locus tree branch (only used for bounded subtrees).
    /// @}
    
    /** \name Data **/
    /// @{
    int Ne; ///< Branch specific effective population size (0=> the global one).
    double n_gen; ///< Number of generations of this node (from the root).
    double time; ///< Time of this node (from the root)
    double gen_length;///< Number of generations from this node to the ancestor.
    double mu_mult; ///< Branch specific substitution rate multi.
    double gtime_mult; ///< Branch specific generation time multi. This parameter is species_tree branch dependent, but I'm replicating it at locus_tree level in order to allow fixed user-specified locus trees as input.
    /// @}
    
    /** \name Pointers **/
    ///@{
    l_node **children; ///< Pointer to the array of children (pointers).
    l_node *anc_node; ///< Ancestor
    l_node *lat_node; ///< Lateral node to not use a big array of memory to save the avaliable locus nodes in the birth-death process of the locus tree simulation.
    g_node **g_nodes; ///< Array of pointers to g_nodes (avaliable gene nodes).
    s_node *conts; ///< s_node that contains this g_node along his branch.
    double *i_probs; ///< Array of probabilities for each number of input lineages (0,1,2...)
    double *i_combprobs; ///< Array of probabilities for each combination of input lineages coming from it children (1:1,1:2,2:1...)
    double *o_probs; ///< Array of probabilities for each number of output lineages (0,1,2...)
    ///@}
    
};

/**
 * Locus tree. 
 * 
 * This structure holds a locus tree, built from independently allocated l_nodes
 * (root) and/or from an array of l_nodes (m_node).
 * It can be associated with a gene tree (using \ref MatchTrees) to use it as a guide tree
 * in a constrained coalescent simulation of gene trees.
 *******************************************************************************/
struct l_tree
{
    /** \name Tree info **/
    /// @{
    int n_nodes; ///< Number of nodes of the tree.
    int max_children; ///< Max number of children in each node. This variable is used to know (and/or modify) the length of the l_node::children allocated memory.
    int n_leaves; ///< Number of leaves of the l_tree. It is necessary due to allowing polytomy.
    int n_gleaves; ///< Number of leaves of the asociated gene tree. It is necessary due to allowing more than one replicate (g_node leaf) per taxa (l_node leaf).
    double gen_time;///< Global generation time
    int Ne; ///< Global effective population size.
    double mu; ///< Global substitution rate
    /// @}
    
    /** \name Pointers **/
    ///@{
    l_node * m_node; ///< s_node array.
    l_node * root; ///< root of a tree of spread s_nodes.
    s_tree * species_tree; ///< Pointer to the species tree.
    g_tree * gene_tree; ///< Pointer to the gene tree.
    /// @}
};

/**
 * Species tree node.
 *
 * These nodes are allocated as an array or a simple tree to form a s_tree.
 * These nodes hold time restrictions (\ref gen_length) and combinatory restrictions (pointers to avaliable 
 * locus nodes, \ref s_node::l_nodes "**l_nodes"). 
 * Polytomies are allowed, since nodes have an \ref children/ 1-\ref s_node::anc_node "ancestor" structure.
 * Taxa names are stored in a name_c variable.
 *******************************************************************************/ 
struct s_node 
{
    /** \name Codes **/
    /// @{
    int index; ///< Memory position in the array.
    int sp_index; ///< Position of the name of this node in name_c
    /// @}
    
    /** \name Node info **/
    /// @{
    int n_child; ///< Number of children (offspring) of this node.
    int n_lnodes; ///< Number of l_nodes wich are being pointed by this species tree (avaliable nodes).
    int n_replicas; ///< Number of g nodes that will evolve within this leaf (only for leaves). 

    /// @}
    
    /** \name Data **/
    /// @{
    int Ne; ///< Branch specific effective population size (0=> the global one).
    double n_gen; ///< Number of generations of this node (from the root).
    double time; ///< Time of this node (from the root)
    double gen_length;///< Number of generations from this node to the ancestor.
    double mu_mult; ///< Branch specific effective substitution rate multi.
    double gtime_mult; ///< Branch specific generation time multi.
    /// @}
    
    /** \name Pointers **/
    ///@{
    s_node **children; ///< Pointer to the array of children (pointers).
    s_node *anc_node; ///< Ancestor
    l_node *l_nodes; ///< Pointers to l_node (avaliable locus node, to walk by *lat_node).
    ///@}
    
};

/**
 * Species tree.
 *
 * This structure holds a species tree , built from independently
 * allocated s_nodes (root) and/or from an array of s_nodes (m_node). 
 *******************************************************************************/

struct s_tree
{
    /** \name Tree info **/
    /// @{
    int n_nodes; ///< Number of nodes of the tree.
    int max_children; ///< Max number of children in each node. This variable is used to know (and/or modify) the length of the s_node::children allocated memory.
    int n_leaves; ///< Number of leaves of the s_tree. It is necessary due to allowing polytomy.
    int n_gleaves; ///< Number of leaves of the asociated gene tree. It is necessary due to allowing more than one replicate (g_node leaf) per taxa (l_node leaf).
    double gen_time;///< Global generation time
    int Ne; ///< Global effective population size.
    double mu; ///< Global substitution rate
    /// @}
    
    /** \name Pointers **/
    ///@{
    s_node * m_node; ///< s_node array.
    s_node * root; ///< root of a tree of spread s_nodes.
    l_tree * locus_tree; ///< Pointer to the locus tree.
    g_tree * gene_tree; ///< Pointer to the gene tree.
    /// @}
};

/**
 * Name container.
 * 
 * This structure contains a string with taxa names. For each name there are
 * \ref max_lname positions in the string, and there are
 * \ref n_names+1 in the string. The first name is the name of a no-named node
 * (usually "Internal node" or something similar).
 *******************************************************************************/
 
struct name_c
{
    char * names; ///< String of names
    int max_lname; ///< Max name-length
    int n_names; ///< Number of taxa names
    
};

/**
 * Period.
 *
 * This structure contains the information of one period of the tree.
 *******************************************************************************/

struct period
{
    double r_bound;
    l_node **l_nodes;
    int n_lnodes;
    
};

// *********************** Prototype of functions ******************************* //

/** \name Node memory manage **/
///@{

// ** Node creation ** //

/**
 * Creates a memory block of n s_nodes, and initializes it.
 *
 * \param n_nodes
 *  Number of nodes.
 * \param max_children
 *  Number of maximum children in each node. It is used to allocate s_node::children.
 * \return s_node pointer to the allocated memory.
 * \note If an error ocurrs, it exits by \ref ErrorReporter.
 *******************************************************************************/
s_node * NewSNodes (int n_nodes,int max_children);

/**
 * Creates a memory block of n l_nodes, and initializes it.
 *
 * \param n_nodes
 *  Number of nodes.
 * \param n_treeleaves
 *  Number of leaves of the final gene tree (equivalent to l_tree::n_gleaves) 
 *  where these new nodes will point. It is used to allocate l_node::g_node.
 * \param max_children
 *  Number of maximum children in each node. It is used to allocate l_node::children.
 * \return l_node pointer to the allocated memory.
 * \note If an error ocurrs, it exits by \ref ErrorReporter.
 *******************************************************************************/
l_node * NewLNodes (int n_nodes,int n_treeleaves,int max_children);

/**
 * Creates a memory block of n g_nodes, and initializes it. 
 *
 * \param n_nodes
 *  Number of nodes.
 * \param max_children
 *  Number of maximum children in each node. It is used to allocate g_node::children.
 * \return g_node pointer to the allocated memory of NULL if n_nodes=0.
 * \note If an error ocurrs, it exits by \ref ErrorReporter.
 *******************************************************************************/
g_node * NewGNodes(int n_nodes, int max_children);

/**
 * Creates a memory block of n periods, and initializes it.
 *
 * \param n_periods
 *  Number of periods.
 * \param max_nodes
 *  Number of maximum l_nodes associated to each period (usually the number of leaves). It is used to allocate period::l_nodes.
 * \return period pointer to the allocated memory or NULL if n_periods=0.
 * \note If an error ocurrs, it exits by \ref ErrorReporter.
 *******************************************************************************/
period * NewPeriods(int n_periods, int max_nodes);

///@}

/** \name NEXUS trees IO functions **/

///@{

long double NNexusTrees(FILE *ifile, int *n_trees);
long double InitNexusParser(FILE *ifile);
long double NextNexusTree(FILE *ifile,char **species_tree_str);

s_tree *ParseNexusSTree(char * nexus,name_c **names_ptr, int verbosity, double gen_time, int Ne, double mu,int ind_persp);

/**
 * Tree parser, from a NEXIS tree to l_tree with spread nodes.
 *
 * Reads a string with a NEXUS tree and converts it in a l_tree, with
 * a l_tree::max_children given by \ref LIMITS::MAX_CHILDREN, and their names stored in
 * a name_c variable.
 * The resulting l_tree has the correct values of l_tree::max_children.
 * Thus, when the tree is converted in a l_tree with clustered nodes, this new tree becomes more memory
 * efficient. The name_c variable is reallocated for memory efficience
 * within this function (at the start, a name_c with \ref LIMITS::MAX_NAME is used). Moreover, count probabilities
 * are also calculated.
 *
 * \param nexus
 *  String with a nexus tree.
 * \param names
 *  name_c ** where save the taxa names (it allocates the memory and fills it).
 * \param verbosity
 *  Code of the amount of communication with the user.
 * \param gen_time
 *  Gen time is going to be used to convert time in generations (1 if the newick tree
 *  has been already given in number of generations).
 * \param Ne
 *   Global effective population size.
 * \param mu
 *   Global substitution rate.
 * \param ind_persp
 *  Common number of desired individuals per species (only applicable for external nodes without a private value specified by the Newick tree)
 * \param n_dup
 * Pointer to get the number of duplications.
 * \param n_loss
 * Pointer to get the number of losses.
 * \param n_trans
 * Pointer to get the number of transferences.
 * \param n_gc
 * Pointer to get the number of gene conversions.
 * \return New readed l_tree.
 * \note If an error ocurrs, it exits by \ref ErrorReporter.
 * \attention NEXUS FORMAT:
 * This function requires rooted NEXUS trees, with branch lengths in all but the
 * root node. <B>It parses Ne, generation times and number of replicates. This has
 * only biological sense if all the nodes of the same branch of the species tree
 * have the same info of Ne & number of replicas, and this is not checked!!!</B>
 * This biological info should be specified as nexus comments. 
 * See the manual and/or program usage for more info.
 *******************************************************************************/
l_tree *ParseNexusLTree(char * nexus,name_c **names_ptr, int verbosity, double gen_time, int Ne, double mu, int ind_persp, int *n_dup, int *n_loss, int *n_trans, int *n_gc);

///@}

/** \name Former Newick IO functions **/
///@{

/**
 * Tree parser, from modified Newick tree to s_tree with spread nodes.
 *
 * Reads a string with a Newick tree and converts it in a s_tree, with
 * a s_tree::max_children given by \ref LIMITS::MAX_CHILDREN, and their names stored in
 * a name_c variable.
 * The resulting s_tree has the correct values of s_tree::max_children.
 * Thus, when the tree is converted in a s_tree with clustered nodes, this new tree becomes more memory
 * efficient. The name_c variable is reallocated for memory efficience
 * within this function (at the start, a name_c with \ref LIMITS::MAX_NAME is used).
 *
 * \param newick
 *  String with a Newick tree.
 * \param names
 *  name_c ** where save the taxa names (it allocates the memory and fills it).
 * \param verbosity
 *  Code of the amount of communication with the user.
 * \param gen_time
 *  Generation time.
 * \param Ne
 *   Global effective population size.
 * \param mu
 *   Global substitution rate.
 * \param ind_persp
 *  Common number of desired individuals per species (only applicable for external nodes without a private value specified by the Newick tree)
 * \return New readed s_tree.
 * \note If an error ocurrs, it exits by \ref ErrorReporter.
 * \attention NEWICK FORMAT:
 * This function requires rooted Newick trees, with branch lengths in all but
 * root node. It allows Ne info in each node, and various replicas of each
 * taxa (only in the leaves). This info should be written in the same position
 * of branch lengths, just after or before it. Ne is indicated with a # and
 * the number of replicas with a /. Example: A:200#2000/3 its the same as
 * A/3#2000:200.
 *******************************************************************************/
s_tree * ParseNewickSTree(char* newick,name_c ** names,int verbosity, double gen_time, int Ne, double mu, int ind_persp);

/**
 * Tree parser, from modified Newick tree to l_tree with spread nodes.
 *
 * Reads a string with a Newick tree and converts it in a l_tree, with
 * a l_tree::max_children given by \ref LIMITS::MAX_CHILDREN, and their names stored in
 * a name_c variable.
 * The resulting l_tree has the correct values of l_tree::max_children.
 * Thus, when the tree is converted in a l_tree with clustered nodes, this new tree becomes more memory
 * efficient. The name_c variable is reallocated for memory efficience
 * within this function (at the start, a name_c with \ref LIMITS::MAX_NAME is used). Moreover, count probabilities
 * are also calculated.
 *
 * \param newick
 *  String with a Newick tree.
 * \param names
 *  name_c ** where save the taxa names (it allocates the memory and fills it).
 * \param verbosity
 *  Code of the amount of communication with the user.
 * \param gen_time
 *  Gen time is going to be used to convert time in generations (1 if the newick tree
 *  has been already given in number of generations).
 * \param Ne
 *   Global effective population size.
 * \param mu
 *   Global substitution rate.
 * \param ind_persp
 *  Common number of desired individuals per species (only applicable for external nodes without a private value specified by the Newick tree)
 * \param n_dup
 * Pointer to get the number of duplications.
 * \param n_loss
 * Pointer to get the number of losses.
 * \param n_trans
 * Pointer to get the number of transferences.
 * \param n_gc
 * Pointer to get the number of gene conversions.
 * \return New readed l_tree.
 * \note If an error ocurrs, it exits by \ref ErrorReporter.
 * \attention NEWICK FORMAT:
 * This function requires rooted Newick trees, with branch lengths in all but the
 * root node. <B>It parses Ne, generation times and number of replicates. This has
 * only biological sense if all the nodes of the same branch of the species tree
 * have the same info of Ne & number of replicas, and this is not checked!!!</B>
 * This biological info should be written in the same position as branch lengths,
 * just after or before it. Example: A:200#2000/3. See the manual and/or program usage for more info.
 *******************************************************************************/
l_tree * ParseNewickLTree(char* newick,name_c ** names,int verbosity, double gen_time, int Ne, double mu, int ind_persp, int *n_dup, int *n_loss, int *n_trans, int *n_gc);

///@}

/** \name Tree memory manage **/
///@{

// ** Tree creation ** //

/**
 * Creates a new s_tree, and initializes it.
 *
 * The nodes of the new tree are allocated as a hole (m_node) unless n_leaves
 * is 1. In this case, this node is the root.
 * \param n_nodes
 *   Number of nodes.
 * \param n_leaves
 *   Number of leaves.
 * \param n_gleaves
 *   Number of leaves of the final gene tree (equivalent to s_tree::n_gleaves).
 * \param max_children
 *   Number of max children per node.
 * \param gen_time
 *  Generation time.
 * \param Ne
 *   Global effective population size.
 * \param mu
 *   Global substitution rate.
 * \return
 *   s_tree pointer to the allocated memory.
 *  \note If an error ocurrs, it exits by \ref ErrorReporter.
 *******************************************************************************/
s_tree * NewSTree (int n_nodes, int n_leaves, int n_gleaves, int max_children, double gen_time, int Ne, double mu);

///@}
/**
 *  Makes a new non-collapsed \ref s_tree "species tree", using a birth-death process.
 *  The birth-death process can be stopped using a maximum time (SSA algorithm) and
 *  a desired number of leaves (conditioned birth-death using a BDSA algorithm).
 *  \note The conditioned birth-death is only possible with birth rates bigger or equal than the death rate. See Hartmann, 2010 and Gernhard, 2008.
 *
 * \param out_tree
 *   Pointer for the new s_tree.
 * \param leaves
 *   Number of leaves of the desired tree. If this option is used, a BSDA algorithm
 *      will be used. The time could be also fixed setting the proper parameter, or 
 *      it will be sampled from  the CDF.
 * \param time
 *   Length of the tree. If only this is fixed, the function will use a SSA algorithm,
 *      so the number of leaves is variable.
 * \param b_rate
 *   Birth rate (birth per generation).
 * \param d_rate
 *   Death rate (death per generation).
 * \param gen_time
 *  Generation time.
 * \param Ne
 *   Global effective population size.
 * \param mu
 *   Global substitution rate.
 * \param ind_per_sp
 *   Number of individuals per species.
 * \param outgroup
 *   Ratio between the ingroup height and the root height (internal branch generated by the addition of an outgroup)
 * \param complete
 *   Logical flag. If complete==1, the tree maintains it death lineages (complete tree). In the other case, the resulting tree is the reconstructed tree shape.
 * \param mrca_time
 *  Logical flag. If mrca_time==1, the trees start in a root with time==time,
 *   and two leaves (two birth-death processes). This will be the MRCA of the
 *   leaves of the tree.
 * \param seed
 *   Seed for the RandomNumber function.
 * \param out_time
 *   Whether the species tree will be written in number of generations (0) or time units (1)
 * \param verbosity
 *   Config about verbosity.
 *
 * \return NO_ERROR on OK or an ErrorCode if any error ocurrs.
 * \attention The resulting tree has to be collapsed or reindexed to be a proper tree (with proper indices and memory structure)
 *******************************************************************************/
long int NewBDSTree (s_tree ** out_tree,int leaves, double time, double b_rate, double d_rate, double gen_time, int Ne, double mu, int ind_per_sp, double outgroup, int complete, int mrca_time, gsl_rng *seed, int out_time,int labels, int verbosity);

/**
 *  Simulates a new locus tree under a birth death model (duplication-loss).
 * \attention The resulting locus tree requires a \ref CollapseLTree to work with it.
 *
 * \param wsp_tree
 *   Guide tree (s_tree).
 * \param wlocus_tree
 *   Resulting locus tree pointer.
 * \param node_ptrs
 *   Big bunch of l_node_ptrs to use in the simulation. This may be allocated and deallocated inside the function, but it is a lot of memory, so this way I avoid some overload in big simulation scenarios.
 * \param b_rate
 *   Birth rate (birth per generation).
 * \param d_rate
 *   Death rate (death per generation).
 * \param seed
 *   Seed for the random number generator.
 * \param min_lleaves
 *   Minimum number of locus tree leaves to consider the simulated tree as valid.
 * \param min_lsleaves
 *   Minimum number of gene tree leaves to consider the simulated locus tree as valid.
 * \param Ne
 *   Global effective population size.
 * \param verbosity
 *   Config about verbosity.
 * \paran st_losses
 *   Pointer to return the number of observed losses in the new tree.
 * \paran st_losses
 *   Pointer to return the number of observed losses in the new tree.
 * \paran st_dups
 *   Pointer to return the number of observed duplications in the new tree.
 * \paran st_leaves
 *   Pointer to return the number of observed number of leaves.
 * \paran st_gleaves
 *   Pointer to return the number of leaves of the future gene trees simulated using this new simulated locus tree.
 *
 * \return NO_ERROR on OK or an ErrorCode if any error ocurrs.
 * \attention The resulting tree has to be collapsed or reindexed to be a proper tree (with proper indices and memory structure)
 *******************************************************************************/
long int SimBDLTree(s_tree *wsp_tree,l_tree **wlocus_tree, l_node **node_ptrs, double b_rate,double d_rate,gsl_rng *seed, int min_lleaves, int min_lsleaves, int verbosity, int *st_losses, int *st_dups, int *st_leaves, int *st_gleaves);

/**
 *  Simulates a new locus tree taking into account GDL and HGT.
 * \attention The resulting locus tree requires a \ref CollapseLTree to work with it.
 *
 * \param wsp_tree
 *   Guide tree (s_tree).
 * \param wlocus_tree
 *   Resulting locus tree pointer.
 * \param node_ptrs
 *   Big bunch of l_node_ptrs to use in the simulation. This may be allocated and deallocated inside the function, but it is a lot of memory, so this way I avoid some overload in big simulation scenarios.
 * \param b_rate
 *   Birth rate (birth per generation).
 * \param d_rate
 *   Death rate (death per generation).
 * \param h_rate
 *   Transfer rate (transfers per generation).
 * \param gc_rate
 *   Gene conversion rate (transfers per generation).
 * \param t_kind
 *   Logical flag. 0=> RTRFR randomly sampled. 1=> RTRFR sampled with probability inversely related to distance (generations).
 * \param seed
 *   Seed for the random number generator.
 * \param min_lleaves
 *   Minimum number of locus tree leaves to consider the simulated tree as valid.
 * \param min_lsleaves
 *   Minimum number of gene tree leaves to consider the simulated locus tree as valid.
 * \param Ne
 *   Global effective population size.
 * \param verbosity
 *   Config about verbosity.
 * \param st_losses
 *   Pointer to return the number of observed losses in the new tree.
 * \param st_dups
 *   Pointer to return the number of observed duplications in the new tree.
 * \param st_transfr
 *   Pointer to return the number of observed transferences in the new tree.
 * \param st_gc
 *   Pointer to return the number of observed gene conversions in the new tree.
 * \param st_leaves
 *   Pointer to return the number of observed number of leaves.
 * \param st_gleaves
 *   Pointer to return the number of leaves of the future gene trees simulated using this new simulated locus tree.
 * \param st_ntrials
 *   Pointer to return the number of trials needed to generate the tree.
 *
 * \return NO_ERROR on OK or an ErrorCode if any error ocurrs.
 * \attention The resulting tree has to be collapsed or reindexed to be a proper tree (with proper indices and memory structure)
 *******************************************************************************/
long int SimBDLHTree(s_tree *wsp_tree,l_tree **wlocus_tree, l_node **node_ptrs, double b_rate,double d_rate, double h_rate, double gc_rate,int t_kind, gsl_rng *seed, int min_lleaves, int min_lsleaves, int verbosity, int *st_losses, int *st_dups, int *st_transfr, int *st_gc, int *st_leaves, int *st_gleaves, int *st_ntrials);

/**
 *  Simulates a new gene tree under the multispecies coalescent process along a locus tree.
 *
 * \param wlocus_tree
 *   Guide tree pointer (locus tree).
 * \param gene_tree
 *   Pointer to a gene tree ponter, where the program returns the simulated gene tree.
 * \param names
 *  Pointer to the name container.
 * \param epsilon_brent
 *  Epsilon of the convergence of the brent method for sampling Bounded multispecies coalescent
 * \param seed
 *   Seed for the random number generator.
 * \paran tn_lcoals
 *   Pointer to return the number of observed extra lineages.
 * \param simlosses
 *   Logical flag. If ==1, then the lost locus tree lineages are simulated at the gene tree level, using only one leave. This is mainly for DBGging and it is not recomended for general users.
 * \param verbosity
 *   Config about verbosity.
 * \param gen_time
 *   Generation time.
 * \param collapse
 *   Logical flag. If ==1, then the gene tree is collapsed and post-ordered.
 * \return NO_ERROR on OK or an ErrorCode if any error ocurrs.
 * \attention If collapse==0 the resulting tree has to be collapsed to be a proper tree (with proper memory structure)
 *******************************************************************************/
long int SimMSCGTree(l_tree *wlocus_tree, g_tree **gene_tree, name_c * names, float epsilon_brent, gsl_rng *seed, int *tn_lcoals, int simlosses,int verbosity, double gen_time, int collapse);

/**
 *  Simulates a new gene tree under the multilocus coalescent process along a locus tree.
 *
 * \param wlocus_tree
 *   Guide tree pointer (locus tree).
 * \param gene_tree
 *   Pointer to a gene tree ponter, where the program returns the simulated gene tree.
 * \param names
 *  Pointer to the name container.
 * \param epsilon_brent
 *  Epsilon of the convergence of the brent method for sampling Bounded multispecies coalescent
 * \param seed
 *   Seed for the random number generator.
 * \paran tn_lcoals
 *   Pointer to return the number of observed extra lineages.
 * \param verbosity
 *   Config about verbosity.
 * \param gen_time
 *   Generation time.
 * \param collapse
 *   Logical flag. If ==1, then the gene tree is collapsed and post-ordered.
 * \return NO_ERROR on OK or an ErrorCode if any error ocurrs.
 * \attention If collapse==0 the resulting tree has to be collapsed to be a proper tree (with proper memory structure)
 *******************************************************************************/
long int SimMLCGTree(l_tree *wlocus_tree, g_tree **gene_tree, name_c * names, float epsilon_brent, gsl_rng *seed, int *tn_lcoals,int verbosity, double gen_time, int collapse);

/**
 * Creates a new l_tree, and initializes it.
 *
 * The nodes of the new tree are allocated as a hole (m_node) unless n_leaves
 * is 1. In this case, this node is the root.
 * \param n_nodes
 *   Number of nodes.
 * \param n_leaves
 *   Number of leaves.
 * \param n_gleaves
 *   Number of leaves of the final gene tree (equivalent to l_tree::n_gleaves).
 * \param max_children
 *   Number of max children per node.
 * \param gen_time
 *   Generation time.
 * \param Ne
 *   Global effective population size.
 * \param mu
 *   Global substitution rate.
 * \return
 *   l_tree pointer to the allocated memory.
 *  \note If an error ocurrs, it exits by \ref ErrorReporter.
 *******************************************************************************/
l_tree * NewLTree (int n_nodes, int n_leaves, int n_gleaves, int max_children, double gen_time, int Ne, double mu);

/**
 * Creates a new g_tree and initializes it.
 *
 * This function creates a iterative tree (n_nodes>1, uses g_tree::m_node) or a recursive oriented
 * tree (n_nodes==1, uses g_tree::root_node)
 *
 * \param n_nodes
 *  Number of desired nodes.
 * \param max_children
 *  Number of maximum children in each node. It is used to allocate s_node::children.
 * \param gen_time
 *  Generation time.
 * \return Pointer to the new allocated g_tree.
 * \note If an error ocurrs, it exits by \ref ErrorReporter.
 *******************************************************************************/
g_tree * NewGTree (int n_nodes, int max_children, double gen_time);

// ** Tree copy ** //

/**
 * Copies or reallocates a s_tree.
 *
 * This function clones a s_tree in a new or old s_tree pointer (in the first case
 * allocating new memory, in the second case rewriting tree values). It can
 * copy only the main values, and/or the pointers of the tree structure and the
 * pointer of l_tree association. If the input tree has sparse nodes, it is 
 * collapsed by \ref CollapseSTree .
 *
 * \param out_tree_ptr
 *  Output s_tree **. At the first, it can has allocated memory or not.
 * \param in_tree
 *  Input s_tree.
 * \param tree_struct
 *  Logical flag. If =0 the s_tree pointers (\ref s_node::children) are not copied.
 * \param l_nodes_ptr
 *  Logical flag. If =0 the l_node pointers (\ref s_node::l_nodes) are not copied.
 * \return \ref NO_ERROR on OK or an \ref ERRORS "error code" if any error
 *  ocurrs.
 *******************************************************************************/
long int CopySTree (s_tree **out_tree_ptr, s_tree *in_tree, int tree_struct, int l_nodes_ptr);

/**
 * Copies or reallocates a l_tree.
 *
 * This function clones a l_tree in a new or old l_tree pointer (in the first case
 * allocating new memory, in the second case rewriting tree values). It can
 * copy just the main values, or also the pointers of the tree structure and the
 * pointer of g_tree association. If the input tree has sparse nodes, it is 
 * collapsed by \ref CollapseLTree .
 *
 * \param out_tree_ptr
 *  Output l_tree **. At the first, it can has allocated memory or not.
 * \param in_tree
 *  Input l_tree.
 * \param tree_struct
 *  Logical flag. If =0 the l_tree pointers (\ref l_node::children) are not copied.
 * \param l_nodes_ptr
 *  Logical flag. If =0 the l_node pointers (\ref l_node::lat_node) are not copied.
 * \param g_nodes_ptr
 *  Logical flag. If =0 the g_node pointers (\ref l_node::g_nodes) are not copied.
 * \param probs
 *  Logical flag. If =0 the l_node::i_probs, l_node::i_combprobs and l_node::o_probs will not be copied.
 * \param relink
 *  Logical flag. If relink = 1, the \ref s_node::l_nodes pointer of the \ref l_node::conts nodes will be updated linking the new locus tree to the species tree.
 * \param preserve
 *  Logical flag. If preserve =1, the original tree (in_tree) will not be modified. It it is 0 certain memory will be dereferenced in the in_tree and referenced by the out_tree in order to save time and memory.
 * \return \ref NO_ERROR on OK or an \ref ERRORS "error code" if any error
 *  ocurrs.
 *******************************************************************************/
long int CopyLTree (l_tree **out_tree_ptr, l_tree *in_tree, int tree_struct, int l_nodes_ptr, int g_nodes_ptr, int probs, int relink, int preserve);

/**
 * Direct copy of a species tree (pre or post-ordered) into a locus tree (post-ordered)
 *
 *
 * \param species_tree
 *  Collapsed species tree.
 * \param locus_tree
 *  Preallocated proper(same size) locus tree.
 * \param no_reindex
 * Logical flag. If reindex == 1 the resulting species tree will be reindexed in post-order (but not reordered in the array)
 * \return \ref NO_ERROR on OK or an \ref ERRORS "error code" if any error
 *  ocurrs.
 *******************************************************************************/
long int CopyStoLTree(s_tree *species_tree, l_tree *locus_tree, int reindex);
//\cond DOXYGEN_EXCLUDE
//// ** Tree edition ** //
//
///**
// * Deletes superfluous nodes due to losses.
// *
// * This function deletes superfluous nodes due to losses of a \ref l_tree after
// * their generation by a birth-death process. It needs an spread \ref l_tree .
// * \param locus_tree
// *  Tree to clean.
// * \return \ref NO_ERROR on OK or an \ref ERRORS "error code" if any error
// *  ocurrs.
// *******************************************************************************/
//long int CleanlossesLTree(l_tree *locus_tree);
//\endcon
// ** Tree reset ** //

/**
 * Resets locus tree simulation related s_tree and s_node values.
 *
 * Useful to use the same species tree to retry l_tree simulations.
 *
 * \param tree
 *  s_tree to restart.
 * \return \ref NO_ERROR on OK or an \ref ERRORS "error code" if any error
 *  ocurrs.
 *******************************************************************************/
long int ResetSTreeSimL (s_tree *tree);

/**
 * Resets a g_tree.
 *
 * This function restarts the values of each g_node of a g_tree into the default
 * value. Useful to use the same memory to simulate all the replicas of a g_tree.
 *
 * \param tree
 *  g_tree to restart.
 * \return \ref NO_ERROR on OK or an \ref ERRORS "error code" if any error
 *  ocurrs.
 *******************************************************************************/
long int ResetGTree (g_tree *tree);

// ** Tree deletion ** //

/**
 * Frees the allocated memory of a s_tree.
 *
 * \param tree
 *  s_tree to free.
 *******************************************************************************/
void FreeSTree (s_tree ** tree);

/**
 * Frees the allocated memory of a l_tree.
 *
 * \param tree
 *  l_tree to free.
 *******************************************************************************/
void FreeLTree (l_tree ** tree);

/**
 * Frees the allocated memory of a g_tree.
 * 
 * \param tree
 *  g_tree to free.
 * \param complete
 *  Logical flag. If complete == 1 all the tree memory is deleted, else the g_tree
 * memory and the root node alive, and only is deleted g_tree::m_nodes.
 *******************************************************************************/
void FreeGTree (g_tree ** tree, int complete);

/**
 * Frees a block of periods.
 *
 * \param periods
 *  Block of periods to free.
 * \param n_periods
 *  Number of periods that constitues this block.
 *******************************************************************************/
void FreePeriods(period * periods, int n_periods);

///@}

/** \name Names memory manage **/
///@{

/**
 * Creates a new name_c memory space where store the taxa names.
 *
 * \param n_names
 *  Number of names to store.
 * \param max_lname
 *  Maximum length of each name (if =0 it uses \ref LIMITS::MAX_NAME).
 * \return name_c pointer pointing to the allocated memory.
 * \note If an error ocurrs, it exits by \ref ErrorReporter.
 *******************************************************************************/
name_c * NewNames (int n_names, int max_lname);

/**
 * Reallocs the name string of a name_c struct.
 *
 * It is used to improve the memory efficience if a default max_lname has been
 * used to allocate the name_c memory. If max_lname is less than the necessary
 * space to save the first name (generic name, usually "Internal node"), this 
 * minimum space is set as max_lname.
 * 
 * \param names
 *  Pointer to the name container.
 * \param max_lname
 *  Maximum length of each name.
 * \note If an error ocurrs, it exits by \ref ErrorReporter.
 *******************************************************************************/
void ReallocNames (name_c * names, int max_lname);

/**
 * Frees the allocated memory of names_c.
 * 
 * \param names -- Names structure to be freed.
 *******************************************************************************/
void FreeNames (name_c ** names);
///@}

/** \name Tree data modification **/
///@{

// ** G and L tree reunification ** //

/**
 * Associates a l_tree and a g_tree to use \ref SimMSCGTree
 *
 * Resets the g_tree simulation dependent values of a l_tree, and associates
 * the leaves with g_nodes to start a new g_tree simulation. After so, each
 * restriction tree leaf points the corresponding g_tree leaf/leaves
 * (if there are more than one replicate) and the internal nodes of the l_tree
 * does not point any g_node.
 * Indeed, it can reset the g_tree (pointers and data).
 *
 * \param locus_tree
 *  l_tree (locus tree).
 * \param gene_tree
 *  g_tree (gene tree).
 * \param reset_gtree
 *  Logical flag. If =1 the g_tree is reset by the function.
 * \param includelosses
 *  Number of g_nodes associated to each lost l_node (0= normal behaviour, 1= simulate lost lineages).
 * \return \ref NO_ERROR on OK or an \ref ERRORS "error code" if any error
 *  ocurrs.
 * \note If the l_tree is spread (l_tree::root instead of l_tree::m_node) the
 *  tree is collapsed by  \ref CollapseLTree in a post-order.
 *******************************************************************************/
long int MatchTreesMSC(l_tree * locus_tree, g_tree * gene_tree, int reset_gtree, int includelosses);

/**
 * Associates a l_tree and a g_tree to use \ref SimMLCGTre
 *
 * Resets the g_tree simulation dependent values of a l_tree, and associates
 * the leaves with g_nodes to start a new g_tree simulation. After so, each
 * restriction tree leaf points the corresponding g_tree leaf/leaves
 * (if there are more than one replicate) and the internal nodes of the l_tree
 * does not point to any g_nodes. It also get the bounded subtrees, and initializes the lnode::n_ilin and lnode:n_olin for leaves and bounds.
 * Indeed, it can reset the g_tree (pointers and data).
 *
 * \param locus_tree
 *  l_tree (locus tree).
 * \param gene_tree
 *  g_tree (gene tree).
 * \param reset_gtree
 *  Number of g_nodes associated to each lost l_node (0= normal behaviour, 1= simulate lost lineages).
 * \return \ref NO_ERROR on OK or an \ref ERRORS "error code" if any error
 *  ocurrs.
 * \note If the l_tree is spread (l_tree::root instead of l_tree::m_node) the
 *  tree is collapsed by  \ref CollapseLTree in a post-order.
 *******************************************************************************/
long int MatchTreesMLC(l_tree *locus_tree, g_tree *gene_tree, int reset_gtree);

// ** Tree conversion ** //

/**
 * Collapse a sparse s_tree (root) in an s_tree with nodes in an array (m_node)
 * using a post-order or a pre-order.
 *
 * \param in_tree
 *  Input s_tree.
 * \param post_order
 *  Logical flag. If post_order = 1, the new tree will be in a post-order. Else,
 *  it will be in a pre-order.
 * \return \ref NO_ERROR on OK or an \ref ERRORS "error code" if any error
 *  ocurrs.
 *******************************************************************************/
long int CollapseSTree (s_tree * in_tree, int post_order);

/**
 * Reindex a s_tree using either a post-order or a pre-order.
 *
 * \param in_tree
 *  Input s_tree.
 * \param post_order
 *  Logical flag. If post_order = 1, the new tree will be in a post-order. Else,
 *  it will be in a pre-order.
 * \return \ref NO_ERROR on OK or an \ref ERRORS "error code" if any error
 *  ocurrs.
 * \attention THIS FUNCTION GENERATES AN UNPROPER TREE WITH NODE INDEXES THAT DO NOT CORESPOND TO THE POSITION IN THE MEMORY ARRAY
 *******************************************************************************/
long int ReindexSTree (s_tree * in_tree, int post_order);

/**
 * Collapse a sparse l_tree (root) in an l_tree with nodes in an array (m_node)
 *
 * \param in_tree
 *  Input l_tree.
 * \param post_order
 *  Logical flag. If post_order = 1, the new tree will be in a post-order. Else,
 *  it will be in a pre-order.
 * \param relink
 *  Logical flag. If relink = 1, the \ref l_node::lat_node and \ref s_node::l_nodes pointer of the \ref l_node::conts nodes will be updated. If = 0, the \ref l_node::lat_node will be erased to avoid malfunctioning due to the reference of freed memory.
 * \param probs
 *  Logical flag. If probs = 1, the new tree will have updated the l_tree::i_probs, l_tree::o_probs and ltree::i_combprobs.
 * \return \ref NO_ERROR on OK or an \ref ERRORS "error code" if any error
 *  ocurrs.
 *******************************************************************************/
long int CollapseLTree (l_tree * in_tree, int post_order, int relink, int probs);

/**
 * Collapse a sparse l_tree (root) in an l_tree with nodes in an array (m_node) allowing to increase or decrease the size of the tree
 * This avoids memory leaks when freezing unused memory (discarded transfers/gcs)
 *
 * \param in_tree
 *  Input l_tree.
 * \param new_nnodes
 *  New number of nodes.
 * \param new_nleaves
 *  New number of leaves.
 * \param post_order
 *  Logical flag. If post_order = 1, the new tree will be in a post-order. Else,
 *  it will be in a pre-order.
 * \param relink
 *  Logical flag. If relink = 1, the \ref l_node::lat_node and \ref s_node::l_nodes pointer of the \ref l_node::conts nodes will be updated. If = 0, the \ref l_node::lat_node will be erased to avoid malfunctioning due to the reference of freed memory.
 * \param probs
 *  Logical flag. If probs = 1, the new tree will have updated the l_tree::i_probs, l_tree::o_probs and ltree::i_combprobs.
 * \return \ref NO_ERROR on OK or an \ref ERRORS "error code" if any error
 *  ocurrs.
 *******************************************************************************/
long int CollapseResizeLTree (l_tree * in_tree,int new_nnodes, int new_nleaves, int post_order, int relink, int probs);

/**
 * Collapse a sparse g_tree (root) in an g_tree with nodes in an array (m_node)
 *
 * \param in_tree
 *  Input g_tree.
 * \param post_order
 *  Logical flag. If post_order = 1, the new tree will be in a post-order. Else,
 *  it will be in a pre-order.
 * \param relink
 *  Logical flag. If relink = 1, the \ref l_node::g_nodes pointers of the \ref g_tree::locus_tree nodes will be updated. If = 0, they will be retained (dangerous, references to freed memory).
 * \return \ref NO_ERROR on OK or an \ref ERRORS "error code" if any error
 *  ocurrs.
 *******************************************************************************/
long int CollapseGTree (g_tree * in_tree, int post_order, int relink);

// ** Branch length modification ** //
/**
 * Recursively (pre-order) updates the l_node::time of a bunch/tree of l_nodes.
 *
 * \param tree
 *  Tree to work with.
 * \return \ref NO_ERROR on OK or an \ref ERRORS "error code" if any error
 *  ocurrs.
 *******************************************************************************/
long int TemporalizeLTree(l_tree *tree);

/**
 * Modifies the branch specific substitution rate multiplier of species tree branches, creating lineage specific rate heterogeneity.
 * \attention This function has to be used before the locus tree simulation/copy.
 *
 * \param sp_tree
 *  Input s_tree.
 * \param alpha
 *  Alpha parameter of the gamma.
 * \param seed
 *  Pointer to the common seed for random number generators.
 * \param gammadump
 *  File to write raw values for each sampled gamma value. Just in case we need it. NULL avoids this feature.
 * \return \ref NO_ERROR on OK or an \ref ERRORS "error code" if any error
 *  ocurrs.
 *******************************************************************************/
long int Rateheter_lineagespec(s_tree *sp_tree,double alpha, gsl_rng *seed, FILE * gammadump);

/**
 * Modifies the general substitution rate of the locus tree, creating gene specific rate heterogeneity.
 * \attention This function has to be used before the gene tree simulation.
 *
 * \param locus_tree
 *  Input l_tree.
 * \param value
 *  VALUE (not multiplier) of the substitution rate.
 * \return \ref NO_ERROR on OK or an \ref ERRORS "error code" if any error
 *  ocurrs.
 *******************************************************************************/
long int Rateheter_genespec(l_tree *locus_tree,double value);

/**
 * Modifies the branch length of the gene tree, generating gene tree branch specific rate heterogeneity.
 *
 * \param tree
 *  Input g_tree.
 * \param alpha
 *  Alpha parameter of the gamma.
 * \param seed
 *  Pointer to the common seed for random number generators.
 * \param gammadump
 *  File to write raw values for each sampled gamma value. Just in case we need it. NULL avoids this feature.
 * \return \ref NO_ERROR on OK or an \ref ERRORS "error code" if any error
 *  ocurrs.
 *******************************************************************************/
long int Rateheter_GTbranchspec(g_tree *tree,double alpha, gsl_rng *seed,FILE * gammadump);

// \cond DOXYGEN_EXCLUDE

///**
// * Transforms the branch lenghts of a gene tree from time to expected number of substitutions per site, using a common substitution rate or the lineage specific ones if they have been given.
// *
// * \param gene_tree
// *  Input g_tree.
// * \param mu
// *  General substitution rate.
// * \return \ref NO_ERROR on OK or an \ref ERRORS "error code" if any error
// *  ocurrs.
// * \note This function only works with post-order collapsed gene trees.
// * \attention This function generates potential non-ultrametric trees. The height meaning changes, being 0 at the root, and becoming negative towards the tips.
// *******************************************************************************/
//long int Mutate_GTree(g_tree *gene_tree,double mu);
//
///**
// * Transforms the branch lenghts of a gene tree from generation time to expected number of substitutions per site, using the generation time and the common substitution rate or the lineage specific ones if they have been given.
// *
// * \param gene_tree
// *  Input g_tree.
// * \param gen_time
// *  Generation time.
// * \param mu
// *  General substitution rate.
// * \return \ref NO_ERROR on OK or an \ref ERRORS "error code" if any error
// *  ocurrs.
// * \note This function only works with post-order collapsed gene trees.
// * \attention This function generates potential non-ultrametric trees. The height meaning changes, being 0 at the root, and becoming negative towards the tips.
// *******************************************************************************/
//long int Evolve_GTree(g_tree *gene_tree,double gen_time, double mu);

//  \endcond

///@}

/** \name Tree features obtention **/
///@{

/**
 * Counts the number of duplications in a locus tree.
 *
 * \param tree
 *  locus tree to be analyzed.
 * \param n_dup
 *  number of duplications (result).
 * \return \ref NO_ERROR on OK or an \ref ERRORS "error code" if any error
 *  ocurrs.
 *******************************************************************************/
long int Count_duplications(l_tree *tree, int *n_dup);

/**
 * Counts the number of transfers in a locus tree.
 *
 * \param tree
 *  locus tree to be analyzed.
 * \param n_dup
 *  number of transfers (result).
 * \return \ref NO_ERROR on OK or an \ref ERRORS "error code" if any error
 *  ocurrs.
 *******************************************************************************/
long int Count_transfers(l_tree *tree, int *n_trans);

/**
 * Counts the number of gene conversions in a locus tree.
 *
 * \param tree
 *  locus tree to be analyzed.
 * \param n_dup
 *  number of gene conversions (result).
 * \return \ref NO_ERROR on OK or an \ref ERRORS "error code" if any error
 *  ocurrs.
 *******************************************************************************/
long int Count_gc(l_tree *tree, int *n_gc);

/**
 * Counts the number of losses in a locus tree.
 *
 * \param tree
 *  locus tree to be analyzed.
 * \param n_losses
 *  number of losses (result).
 * \return \ref NO_ERROR on OK or an \ref ERRORS "error code" if any error
 *  ocurrs.
 *******************************************************************************/
long int Count_losses(l_tree *tree, int *n_losses);

/**
 * Measures g_tree height under different units using a post-order recursion.
 *
 * \param tree
 *  Pointer to the g_tree to be analyzed.
 * \param height
 *  tree height pointer (result).
 * \param unit
 *  Length unit \ref UNITS
 * \return \ref NO_ERROR on OK or an \ref ERRORS "error code" if any error
 *  ocurrs.
 *******************************************************************************/
long int Measure_GT_height(g_tree *tree,long double *height,int unit);

/**
 * Measures g_tree length under different units using a post-order recursion.
 *
 * \param tree
 *  Pointer to the g_tree to be analyzed.
 * \param height
 *  tree height pointer (result).
 * \param unit
 *  Length unit \ref UNITS
 * \return \ref NO_ERROR on OK or an \ref ERRORS "error code" if any error
 *  ocurrs.
 *******************************************************************************/
long int Measure_GT_length(g_tree *tree, long double *length, int unit);

/**
 * Measures s_tree height under different units using a post-order recursion.
 *
 * \param tree
 *  Pointer to the s_tree to be analyzed.
 * \param height
 *  tree height pointer (result).
 * \param unit
 *  Length unit \ref UNITS
 * \return \ref NO_ERROR on OK or an \ref ERRORS "error code" if any error
 *  ocurrs.
 *******************************************************************************/
long int Measure_ST_height(s_tree *tree,long double *height, int unit);

/**
 * Measures s_tree length under different units using a post-order recursion.
 *
 * \param tree
 *  Pointer to the s_tree to be analyzed.
 * \param height
 *  tree height pointer (result).
 * \param unit
 *  Length unit \ref UNITS
 * \return \ref NO_ERROR on OK or an \ref ERRORS "error code" if any error
 *  ocurrs.
 *******************************************************************************/
long int Measure_ST_length(s_tree *tree,long double *length,int unit);

/**
 * Measures the distance between locus tree duplication, transfers or gene conversion nodes and the MRCA of the associated gene tree lineages (one from each locus tree), equivalent to the estimated event time using classical most parsimonious reconciliation methods.
 *
 * \param wg_tree
 *  Pointer to the gene tree.
 * \param event
 *  Event code, as in \ref NODE_VALUES .
 * \param distances
 *  Pointer to an array of doubles where the distance for each tree height pointer (result) will be written.
 * \param n_events
 *  Number of events to analyze. (Used to realloc distances array)
 * \param unit
 *  Unit in which the distances are to be measured (only GL has been implemented so far) \ref UNITS
 * \return \ref NO_ERROR on OK or an \ref ERRORS "error code" if any error
 *  ocurrs.
 * \attention -1 is given as the distance when there is no way of getting the MRCA (one loss or reception of transfer/gene conversion after the duplication, for example)
 *******************************************************************************/
long int MeasureMRCAEVdistance(g_tree *wg_tree,int event,double **distances, int n_dup, int unit);

/**
 * Calculates the expected locus tree number of leaves
 *
 * \param tree
 *  Species tree
 * \param b_rate
 *  Birth (duplication) rate.
 * \param d_rate
 *  Death (loss) rate.
 * \return Expected number of pruned leaves (leaves reaching to the present)
 * \attention -1 is given as the distance when there is no way of getting the MRCA (one loss or reception of transfer/gene conversion after the duplication, for example)
 *******************************************************************************/
double ExpectedPrunedLtreeNleavesSNodes(s_tree *tree, double b_rate, double d_rate);

long int CheckUltrametricitySTree(s_tree *tree);

long int CheckUltrametricityLTree(l_tree *tree);

///@}

/** \name Trees I/O **/
///@{

/**
 * Writes a given s_tree in Newick format in stdout.
 *
 * \param names
 *  name_c pointer with taxa names.
 * \param in_tree
 *  Tree to print.
 * \param time
 *  Logical flag. If ==1, the tree will be printed in time units instead of
 *  generations.
 * \param int_labels
 *  Logical flag. If ==1, the tree will be printed with internal nodes labeled with its index number (postorder starting at 0).
 * \return \ref NO_ERROR on OK or an \ref ERRORS "error code" if any error
 *  ocurrs.
 *  \note It requires a known root node into the r_tree structure.
 *******************************************************************************/
long int WriteSTree (s_tree *in_tree, name_c * names, int time, int int_labels);

/**
 * Writes a given s_tree in Newick format in a file.
 *
 * \param file
 *  File where the tree will be printed.
 * \param in_tree
 *  Tree to print.
 * \param names
 *  name_c pointer with taxa names.
 * \param time
 *  Logical flag. If ==1, the tree will be printed in time units instead of
 *  generations.
 * \return \ref NO_ERROR on OK or an \ref ERRORS "error code" if any error
 *  ocurrs.
 * \param int_labels
 *  Logical flag. If ==1, the tree will be printed with internal nodes labeled with its index number (postorder starting at 0).
 * \note The FILE * should be previously opened and checked to avoid errors, and
 *  closed after this function.
 * \note It requires a known root node into the s_tree structure.
 *******************************************************************************/
long int WriteSTreeFile(FILE * file,s_tree *in_tree, name_c * names, int time, int int_labels);

/**
 * Writes a given l_tree in Newick format in stdout.
 *
 * \param names
 *  name_c pointer with taxa names.
 * \param in_tree
 *  Tree to print.
 * \param time
 *  Logical flag. If ==1, the tree will be printed in time units instead of
 *  generations.
 * \param int_labels
 *  Logical flag. If ==1, the tree will be printed with internal nodes labeled with its index number (postorder starting at 0).
 * \return \ref NO_ERROR on OK or an \ref ERRORS "error code" if any error
 *  ocurrs.
 * \note It requires a known root node into the r_tree structure.
 *******************************************************************************/
long int WriteLTree (l_tree *in_tree, name_c * names, int time, int int_labels);

/**
 * Writes a given l_tree in Newick format in a file.
 *
 * \param file
 *  File where the tree will be printed.
 * \param in_tree
 *  Tree to print.
 * \param names
 *  name_c pointer with taxa names.
 * \param time
 *  Logical flag. If ==1, the tree will be printed in time units instead of
 *  generations.
 * \param int_labels
 *  Logical flag. If ==1, the tree will be printed with internal nodes labeled with its index number (postorder starting at 0).
 * \return \ref NO_ERROR on OK or an \ref ERRORS "error code" if any error
 *  ocurrs.
 * \note The FILE * should be previously opened and checked to avoid errors, and
 *  closed after this function.
 *  \note It requires a known root node into the l_tree structure.
 *******************************************************************************/
long int WriteLTreeFile(FILE * file,l_tree *in_tree, name_c * names, int time, int int_labels);


/**
 * Writes a comma-separated list of daughters locus tree lineages in a file.
 *
 * \param file
 *  File where the list will be printed.
 * \param in_tree
 *  Tree to analize.
 * \param names
 *  name_c pointer with taxa names.
 * \return \ref NO_ERROR on OK or an \ref ERRORS "error code" if any error
 *  ocurrs.
 * \note The FILE * should be previously opened and checked to avoid errors, and
 *  closed after this function.
 *  \note It requires a known root node into the l_tree structure.
 *******************************************************************************/
long int WriteDaughtersFile (FILE *file,l_tree *in_tree, name_c * names);

/**
 * Writes a given g_tree in Newick format in stdout.
 *
 * \param names
 *  name_c pointer with taxa names.
 * \param in_tree
 *  Tree to print.
 * \param int_labels
 *  Logical flag. If ==1, the tree will be printed with internal nodes labeled with its index number (postorder starting at 0).
 * \return \ref NO_ERROR on OK or an \ref ERRORS "error code" if any error
 *  ocurrs.
 * \note The FILE * should be previously opened and checked to avoid errors, and
 *  closed after this function.
 *  \note It requires a known root node into the g_tree structure.
 *******************************************************************************/
long int WriteGTree (g_tree *in_tree, name_c * names, int int_labels);

/**
 * Writes a given g_tree in Newick format in a string.
 *
 * \param string
 *  Output string.
 * \param in_tree
 *  Input tree pointer.
 * \param names
 *  names_c containing the taxa name of the tree.
 * \param int_labels
 *  Logical flag. If ==1, the tree will be printed with internal nodes labeled with its index number (postorder starting at 0).
 * \return \ref NO_ERROR on OK or an \ref ERRORS "error code" if any error
 *  ocurrs.
 * \note The FILE * should be previously opened and checked to avoid errors, and
 *  closed after this function.
 *  \note It requires a known root node into the g_tree structure.
 *******************************************************************************/
long int WriteGTreeStr (char * string, g_tree *in_tree, name_c * names, int int_labels);

/**
 * Writes a given g_tree in Newick format in a file.
 *
 * \param file
 *  Output opened file.
 * \param in_tree
 *  Input tree pointer.
 * \param names
 *  names_c containing the taxa name of the tree.
 * \param int_labels
 *  Logical flag. If ==1, the tree will be printed with internal nodes labeled with its index number (postorder starting at 0).
 * \return \ref NO_ERROR on OK or an \ref ERRORS "error code" if any error
 *  ocurrs.
 * \note The FILE * should be previously opened and checked to avoid errors, and
 *  closed after this function.
 *  \note It requires a known root node into the l_tree structure.
 *  This function is more efficient than using WriteGTreeStr followed by fprintf.
 *******************************************************************************/
long int WriteGTreeFile (FILE * file, g_tree *in_tree, name_c * names, int int_labels);

/**
 * Calculates and writes the mapping of the species and locus trees in a file (those trees have to be linked in the way that the program does along the locus tree simulation).
 *
 * \param wsp_tree
 *  Species tree.
 * \param locus_tree
 *  Locus tree.
 * \param names
 *  names_c containing the taxa names of the tree.
 * \param mapsl_outname
 *  Name of the file where the function is going to write the mapping.
 * \return \ref NO_ERROR on OK or an \ref ERRORS "error code" if any error
 *  ocurrs.
 *******************************************************************************/
long int WriteMappingSL(s_tree *wsp_tree, l_tree *locus_tree, name_c *names, char *mapsl_outname);

/**
 * Calculates and writes the mapping of the locus and gene trees in a file (those trees have to be linked in the way that the program does along the gene tree simulation).
 *
 * \param gene_tree
 *  Gene tree.
 * \param names
 *  names_c containing the taxa names of the tree.
 * \param maplg_outname
 *  Name of the file where the function is going to write the mapping.
 * \return \ref NO_ERROR on OK or an \ref ERRORS "error code" if any error
 *  ocurrs.
 * \attention The gene tree must have been previously collapsed in postorder.
 *******************************************************************************/
long int WriteMappingLG(g_tree *gene_tree, name_c *names, char *maplg_outname);

// \cond DOXYGEN_EXCLUDE

//DEPRECATED AND OUTDATED (It doesn't work with transfers)

/**
 * Calculates and writes the mapping of the gene and species trees in a file (those trees have to be linked in the way that the program does along the locus and gene tree simulation). This function if mainly for DBGging and it should not be used by a common user. The resulting mapping does not have to be a true representation of the reality (since you cannot skip the locus tree).
 *
 * \param wsp_tree
 *  Species tree.
 * \param gene_tree
 *  Gene tree.
 * \param names
 *  names_c containing the taxa names of the tree.
 * \param map_outname
 *  Name of the file where the function is going to write the mapping.
 * \return \ref NO_ERROR on OK or an \ref ERRORS "error code" if any error
 *  ocurrs.
 *******************************************************************************/
//long int WriteMappingLosses(s_tree *wsp_tree, g_tree *gene_tree, name_c *names, char *map_outname);

// \endcond

/**
 * Checks a Nexus tree string.
 *
 * Checks a Nexus tree string and warns if it finds any problem concerning to the structure of the tree. It does not test the fixed custom parameters given using Nexus comments [&xxx]. They are tested while being parsed using GetXnodeParamsFromNexusComments (S or L depending on the nature of the tree).
 * \param tree
 *  Nexus tree string.
 * \return \ref NO_ERROR on OK or an \ref ERRORS "error code" if any error
 *  ocurrs.
 * \note Function based on CheckTree from Mosaic 1.0 (Posada D.) 
 *******************************************************************************/
long int CheckNexusTree (char * tree);

/**
 * Checks a newick tree string taking it as a species tree.
 *
 * Checks a newick tree string and warns if it finds any problem.
 * \param tree
 *  Newick tree string.
 * \return \ref NO_ERROR on OK or an \ref ERRORS "error code" if any error
 *  ocurrs.
 * \note Function based on CheckTree from Mosaic 1.0 (Posada D.)
 *******************************************************************************/
long int CheckNewickSTree (char * tree);

/**
 * Checks a newick tree string taking it as a locus tree.
 *
 * Checks a newick tree string and warns if it finds any problem.
 * \param tree
 *  Newick tree string.
 * \return \ref NO_ERROR on OK or an \ref ERRORS "error code" if any error
 *  ocurrs.
 * \note Function based on CheckTree from Mosaic 1.0 (Posada D.)
 *******************************************************************************/
long int CheckNewickLTree (char * tree);


///@}

/** \name Miscellanea **/
///@{
extern void PrintUsage(void);
void PrintXCharError(char *string, int x, char *errormsg1, char *errormsg);
inline void ResetBuffer(char *buffer, size_t size);
inline void reallocBuffer(char **buffer,size_t *size, size_t newsize);

///@}


/**
 * Detector of errors. 
 * This function writes info of errors in stderr and closes
 * the program if it is necessary. The error_codes are declared in \ref ERRORS. 
 *
 * \param code
 *  Error code
 * \param MSN
 *  Mensage to write in stderr before exiting the program.
 ******************************************************************************/

static inline void ErrorReporter(long int code, char * MSN)
{
    switch (code) 
    {
        case NO_ERROR:
            break;
        case MEM_ERROR:
            fprintf(stderr,"\nDynamic memory error %s\n",MSN);
            fflush(stdout);
            fflush(stderr);
            exit (EXIT_FAILURE);
            break;
        case IO_ERROR:
            fprintf(stderr,"\nIO_error %s\n", MSN);
            fflush(stdout);
            fflush(stderr);
            exit (EXIT_FAILURE);
            break;
        case LOOP_ERROR:
            fprintf(stderr,"\nError: Infinite loop %s\n", MSN);
            fflush(stdout);
            fflush(stderr);
            exit (EXIT_FAILURE);
            break;
        case SETTINGS_ERROR:
            fprintf(stderr,"\nSettings error: %s\nExecute SimPhy with the option -h in order to print the usage information\n", MSN);
            fflush(stdout);
            fflush(stderr);
            exit (EXIT_FAILURE);
            break;
        case UNEXPECTED_VALUE:
            fprintf(stderr,"\nError: Unexpected value has been reached %s\n", MSN);
            fflush(stdout);
            fflush(stderr);
            exit (EXIT_FAILURE);
            break;
        case DB_ERROR:
            fprintf(stderr,"\nError with the DB manipulation %s\n", MSN);
            fflush(stdout);
            fflush(stderr);
            exit (EXIT_FAILURE);
            break;
        case TERMINATE_NOERROR:
            exit(EXIT_SUCCESS);
            break;
        default: //Pointers with no error.
            break;
    }
}
#endif

///@}
