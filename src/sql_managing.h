/**
 *
 * \file sql_managing.h
 * Library to manage the SQLite database
 *******************************************************************************/

/**
 * \defgroup sql_managing_pub SQL managing module
 *
 * Module with types and functions to drive SQLite databases managing
 *******************************************************************************/
///@{


#ifndef sql_managing_h
#define sql_managing_h

#ifndef _SQLITE3_H_
#include <sqlite3.h>
#endif

#ifndef _sampling_h
#include "sampling.h"
#endif

#ifndef trees_h
#include "trees.h"
#endif

// *********************** Prototype of functions ******************************* //

/** \name Sampling variables managing **/

///@{

// *** Sampling variables managing *** //

// *** Database Initialization functions ***//

long int InitDB(sqlite3 **database, char * db_filename);


// *** Database write values *** //

//long int WriteReplicateDB(sqlite3 **database,double gen_time, int n_iltrees, double sb_rate, double sd_rate, int bds_leaves, double bds_length, double b_rate, double d_rate, double t_rate, double gc_rate, int t_kind,double outgroup, int min_lleaves,int min_lsleaves,int ind_per_sp,int nl_trees,int ng_trees, int Ne, double mu,double alpha_s,double alpha_l,double alpha_g, int curr_stree);
long int WriteSTreeDB(sqlite3 **database, int n_leaves, double sb_rate, double sd_rate, int bds_leaves, double bds_length, double height, double length, double outgroup, int ind_per_sp, int n_loci, double b_rate,double d_rate, double t_rate, double gc_rate, double alpha_s, double alpha_l, double salpha_g, int Ne, double mu, double gen_time);
long int WriteLTreeDB(sqlite3 **database,int n_ltree, int SID, double b_rate, double d_rate, double t_rate, double gc_rate, int n_leaves, int n_dup, int n_loss, int n_transf, int n_gc, int n_ltrials,double gamma,double lalpha_g);
long int WriteGTreeDB(sqlite3 **database, int n_gtree, unsigned long long LID, int n_ltree, int SID, double alpha_g,int n_leaves, int extra_lineages, double tree_h_cu, double tree_l_ss);

long int CloseDB(sqlite3 **database);



///@}


#endif
///@}
