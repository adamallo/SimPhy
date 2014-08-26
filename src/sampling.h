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
enum TYPES{UI=0,D=1};

// *********************** Declaration of custom types ************************** //

typedef union ui_d ui_d;
typedef struct sampling_unit sampling_unit;

/**
 * Union of double and int
 *******************************************************************************/

union ui_d
{
    int i;
    double d;
};

/**
 * Data struct to guide parameter sampling
 *******************************************************************************/

struct sampling_unit
{
    int distribution_code; ///< Distribution id
    double params[5]; ///< parameters
    ui_d value; ///< sampled value (current value)
    int vtype;
    
};

// *********************** Prototype of functions ******************************* //

/** \name Sampling variables managing **/

///@{

void set_sampling_uint(sampling_unit *variable, int value);
void set_sampling_double(sampling_unit *variable, double value);
long int sample_distr(gsl_rng *r, int n_arg,...);
long int ParseSampling(char * p, sampling_unit * sample);
void Print_Sampling(sampling_unit sample, char * buffer);
int is_sampling_set(sampling_unit value);

///@}


#endif
///@}
