/**
 *
 * \file sampling.h
 * Library of bayesian-like sampling to drive simulations
 *******************************************************************************/

/**
 * \defgroup sampling_pub Bayesian-like sampling module
 *
 * Module with types and functions to drive bayesian-like simulations.
 *******************************************************************************/
///@{


#ifndef sampling_h
#define sampling_h

#ifndef trees_h
#include "trees.h"
#endif

///GSL Library
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#ifndef _STDARG_H
#include <stdarg.h>
#endif

/**
 * Function-like macro to return the value with the proper type.
 ********************/
#define get_sampling(a) (a.vtype==UI?a.value.i:a.value.d)

/**
 * Function-like macro to set the value with the proper type (it has to had been initialized before).
 ********************/
#define set_propsampling(a,b) (a->vtype==UI?(a->value=(ui_d)(int)b):(a->value=(ui_d)(double)b))

/**
 * Macro to initialize SUs
 ********************/
#define INIT_SU {FIXED,{0,0,0,0,0},{UI,UI,UI,UI,UI},0,UI,0}

extern int MAX_IT;

/**
 * \enum DISTRIBUTIONS
 * Codes for distribution sampling using GSL.
 *
 *
 * \var FIXED
 *  No sampling (1 parameter, value).
 * \var UNIFORM
 *  Uniform distribution (2 parameters)
 * \var NORMAL
 *  Normal distribution (2 parameters).
 *******************************************************************************/

enum DISTRIBUTIONS {FIXED=0, UNIFORM=1, NORMAL=2, EXPONENTIAL=3, GAMMA=4, LOGNORMAL=5, LOGNORMAL_MULT=7};

/**
 * \enum TYPES
 * Possible types of the ui_d union.
 *
 * This constants are a kind of text codification of interesting int values.
 *
 * \var UI
 *  int
 * \var D
 *  Double
 *******************************************************************************/
enum TYPES{UI=0,D=1,SU=2};

// *********************** Declaration of custom types ************************** //

typedef union ui_d ui_d;
typedef struct sampling_unit sampling_unit;
typedef union u_idsu u_idsu;
typedef struct sampling_duple sampling_duple;
typedef struct sampling_table sampling_table;

/**
 * Union of double and int
 *******************************************************************************/

union ui_d
{
    int i;
    double d;
};

/**
 * Union of int, double and sampling unit
 *******************************************************************************/

union u_idsu
{
    int i;
    double d;
    sampling_unit *p;
};

/**
 * Data struct to guide parameter sampling
 *******************************************************************************/

struct sampling_unit
{
    int distribution_code; ///< Distribution id
    u_idsu params[5]; ///< parameters
    int params_type[5];
    ui_d value; ///< sampled value (current value)
    int vtype;
    int dependent_index;
    
};

/**
 * Data structure to allow sampling unit interdependence
 *******************************************************************************/

struct sampling_duple
{
    char name[3];
    sampling_unit *p;
};

/**
 * Data structure to allow sampling unit interdependence
 *******************************************************************************/
struct sampling_table
{
    sampling_duple *table;
    int n_duples;
};

// *********************** Prototype of functions ******************************* //

/** \name Sampling variables managing **/

///@{

void set_sampling_uint(sampling_unit *variable, int value);
void set_sampling_double(sampling_unit *variable, double value);
void set_sampling_pointeruint(sampling_unit *variable, sampling_unit *dependent);
void set_sampling_pointerdouble(sampling_unit *variable, sampling_unit *dependent);
long int sample_distr(gsl_rng *r, int n_arg,...);
long int ParseSampling(char * p, sampling_unit * sample, const sampling_table sampling_vars);
void Print_Sampling(sampling_unit *sample, char * buffer,const sampling_table stable);
int is_variable(sampling_unit value);
int is_sampling_set(sampling_unit value);

///@}

#endif
///@}
