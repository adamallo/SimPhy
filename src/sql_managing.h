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

extern int MAX_IT;
extern int MAX_NAME;
extern int MAX_CHILDS;
extern int MAX_LEAVES;

// *********************** Prototype of functions ******************************* //

/** \name Sampling variables managing **/

///@{

// *** Sampling variables managing *** //

// *** Database Initialization functions ***//

long int InitDB(sqlite3 **database, char * db_filename);


// *** Database write values *** //

long int WriteSTreeDB(sqlite3 **database, int n_leaves, double height, double length, double outgroup, int ind_per_sp, int n_loci, double alpha_s, double alpha_l, double alpha_g, int Ne, double mu, double gen_time);
long int WriteLTreeDB(sqlite3 **database,int n_ltree, int SID, double b_rate, double d_rate, double t_rate, double gc_rate, int n_leaves, int n_dup, int n_loss, int n_transf, int n_gc, double gamma);
long int WriteGTreeDB(sqlite3 **database, int n_gtree, unsigned long long LID, int n_ltree, int SID, int n_leaves, int extra_lineages, double tree_h_cu, double tree_l_ss);

long int CloseDB(sqlite3 **database);



///@}


#endif
///@}
