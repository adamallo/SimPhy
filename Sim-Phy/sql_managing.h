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
#include "sqlite3.h"
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

long int WriteSTreeDB(sqlite3 **database, unsigned int n_leaves, double height, double length, double outgroup, unsigned int ind_per_sp, unsigned int n_loci, double alpha_s, double alpha_l, double alpha_g, unsigned int Ne, double mu, double gen_time);
long int WriteLTreeDB(sqlite3 **database,unsigned int n_ltree, unsigned int SID, double b_rate, double d_rate, double t_rate, unsigned int n_leaves, unsigned int n_dup, unsigned int n_loss, unsigned int n_transf, double gamma);
long int WriteGTreeDB(sqlite3 **database, unsigned int n_gtree, unsigned long long LID, unsigned int n_ltree, unsigned int SID, unsigned int n_leaves, unsigned int extra_lineages, double tree_h_cu, double tree_l_ss);

long int CloseDB(sqlite3 **database);



///@}


#endif
///@}
