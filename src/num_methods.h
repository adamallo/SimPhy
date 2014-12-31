/**
 * 
 * \file num_methods.h
 * Header of the library of numerical methods.
 *******************************************************************************/

/**
 * \defgroup mum_methods Numerical methods
 *
 * Module with variables and functions of numerical methods.
 *******************************************************************************/ 
///@{

#ifndef num_methods_h
#define num_methods_h
#endif
#include <math.h>
#ifndef trees_h
#include "trees.h"
#endif
#ifndef __GSL_SF_GAMMA_H
#include <gsl/gsl_sf_gamma.h>
#endif
#ifndef __GSL_RNG_H
#include <gsl/gsl_rng.h>
#endif
#ifndef __GSL_SF_LOG_H
#include <gsl/gsl_sys.h>
#endif

////MPFR Library
//#ifndef  _STDARG_H
//#include <stdarg.h>
//#endif
#ifndef __MPFR_H
#include <mpfr.h>
#endif
//#include <float.h>

extern int MAX_IT;

/** \name Root finding methods **/
///@{

/**
 * Finds the root of a function using the Brent's method.
 *
 * Function based in the pseudo-C code found in 
 * http://en.wikipedia.org/wiki/Brent's_method modified both attending to the implementation and the addition of 2 extra convergence conditions to avoid the return of the initial points (a and b), as it is required in this context.
 * \param *function 
 *   Pointer to the function.
 * \param a
 *   Point a.
 * \param b
 *   Point b.
 * \param epsilon
 *   Convergence value.
 * \param result
 *   Pointer to return the result.
 * \param verbosity
 *   How many info to display by stdout.
 * \return
 *   Root.
 *  \note If an error ocurrs, it exits by \ref ErrorReporter.
 *******************************************************************************/

long int brent_root (double (*f)(double value, int n, va_list ap), double a, double b, float epsilon, double *result, int verbosity, int n_arg, ...);

///@}

/** 
 * \name Coalescent theory related methods 
 *******************************************************************************/
///@{

/**
 * CDF function of a bounded coalescent process, derived by Leo Martins, implemented by me.
 *
 * \param w_time 
 *   Waiting time.
 * \param n_leaves
 *   Original number of leaves to coalesce.
 * \param pop_size
 *   Ne of the branch.
 * \param bound_time
 *   Bound of the process (time).
 * \return
 *   Cumulative density.
 *  \note If an error ocurrs, it exits by \ref ErrorReporter.
 *******************************************************************************/

double CdfBoundedCoalescent(double w_time, int n_leaves,int pop_size, double bound_time);

#ifdef __MPFR_H
long double cdf_bounded_coalescent_mpfr(long double time, int n_leaves,int pop_size, long double bound_time, int precision);
#endif

/**
 * Sampler of a bounded coalescent process.
 *
 * \param w_time 
 *   Waiting time.
 * \param n_arg
 *   Number of arguments
 * \param ap
 *   Rest of arguments, they must be:
 * \param density
 *   Density (sampled uniform random number) of the w_time to find.
 * \param n_leaves
 *   Original number of leaves to coalesce.
 * \param pop_size
 *   Ne of the branch.
 * \param bound_time
 *   Bound of the process (time).
 * \return
 *   Cumulative density.
 *  \note If an error ocurrs, it exits by \ref ErrorReporter.
 *******************************************************************************/

double SampleBoundedCoalescent(double w_time, int n_arg, va_list ap);

/**
 * Calculates de probability of going from x to y lineages along a locus tree branch, defined by a time and a population size.
 *
 * \param i_lineages
 *   Input lineages.
 * \param o_lineages
 *   Output lineages.
 * \param bound_time
 *   Length of the tree branch measured in number of generations.
 * \param pop_size
 *   Ne of the branch.
 * \return
 *   Probability.
 *  \attention This function is not usable for big number of lineages (around 150) due to big number issues.
 *  \note If an error ocurrs, it exits by \ref ErrorReporter.
 *******************************************************************************/

double ProbCoalFromXtoYLineages(int i_lineages,int o_lineages,double bound_time,int pop_size);

/**
 * Calculates de probability of going from x to y lineages along a locus tree branch, defined by a time and a population size. Internal calculations in logspace in order to allow bigger problems.
 *
 * \param i_lineages
 *   Input lineages.
 * \param o_lineages
 *   Output lineages.
 * \param bound_time
 *   Length of the tree branch measured in number of generations.
 * \param pop_size
 *   Ne of the branch.
 * \return
 *   Probability.
 *  \note If an error ocurrs, it exits by \ref ErrorReporter.
 *******************************************************************************/
double LogscaleProbCoalFromXtoYLineages(int i_lin, int o_lin, double bound_time, int Ne);

/**
 * Calculates de probability of going from x to y lineages along a locus tree branch, defined by a time and a population size. Internal calculations in logspace in order to allow bigger problems, coupled with compensated summations (Kahan) in order to reduce precision problems and with certain optimizations.
 *
 * \param i_lineages
 *   Input lineages.
 * \param o_lineages
 *   Output lineages.
 * \param bound_time
 *   Length of the tree branch measured in number of generations.
 * \param pop_size
 *   Ne of the branch.
 * \return
 *   Probability.
 *  \note If an error ocurrs, it exits by \ref ErrorReporter.
 *******************************************************************************/
double KahanLogscaleProbCoalFromXtoYLineages(int i_lin, int o_lin, double bound_time, int Ne);

/**
 * Calculates de logprobability of going from x to y lineages along a locus tree branch, defined by a time and a population size. Internal calculations in logspace in order to allow bigger problems, coupled with compensated summations (Kahan) in order to reduce precision problems and with certain optimizations.
 *
 * \param i_lineages
 *   Input lineages.
 * \param o_lineages
 *   Output lineages.
 * \param bound_time
 *   Length of the tree branch measured in number of generations.
 * \param pop_size
 *   Ne of the branch.
 * \return
 *   Probability.
 *  \note If an error ocurrs, it exits by \ref ErrorReporter.
 *******************************************************************************/
double KahanLogscaleLogProbCoalFromXtoYLineages(int i_lin, int o_lin, double bound_time, int Ne);


/**
 * Calculates de probability of going from x to y lineages along a locus tree branch, defined by a time and a population size. Original implementation, much slower.
 *
 * \param i_lineages
 *   Input lineages.
 * \param o_lineages
 *   Output lineages.
 * \param bound_time
 *   Length of the tree branch measured in number of generations.
 * \param pop_size
 *   Ne of the branch.
 * \return
 *   Probability.
 *  \attention This function is not usable for big number of lineages (around 150) due to big number issues.
 *  \note If an error ocurrs, it exits by \ref ErrorReporter.
 *******************************************************************************/
double OriginalProbCoalFromXtoYLineages(int i_lin, int o_lin, double bound_time, int Ne);

//double LogscaleOriginalProbCoalFromXtoYLineages(int i_lin, int o_lin, double bound_time, int Ne); //Unfinished

#ifdef __MPFR_H
/**
 * Calculates de probability of going from x to y lineages along a locus tree branch, defined by a time and a population size. Arbitrary precision version over the original implementation. (Slowest)
 *
 * \param i_lineages
 *   Input lineages.
 * \param o_lineages
 *   Output lineages.
 * \param bound_time
 *   Length of the tree branch measured in number of generations.
 * \param pop_size
 *   Ne of the branch.
 * \param precision
 *   Precision used for MPFR calculations.
 * \return
 *   Probability.
 *  \note If an error ocurrs, it exits by \ref ErrorReporter.
 *******************************************************************************/

double MPFROriginalProbcoalFromXtoYLineages (int i_lin, int o_lin, double bound_time, int Ne, int precision);
double MPFRSampleCoalTimeMLCFromXtoYLineages(int i_lin, int o_lin, double bound_time, int pop_size, double brent_epsilon, double density, int verbosity, int precision);
double MPFRSampleCDFCoalTimeMLCFromXtoYLineages(double w_time, int n_arg, va_list ap);
#endif

double SampleCoalTimeMLCFromXtoYLineages(int i_lin, int o_lin, double bound_time, int pop_size, double brent_epsilon, double density, int verbosity);
double LogscaleSampleCoalTimeMLCFromXtoYLineages(int i_lin, int o_lin, double bound_time, int pop_size, double brent_epsilon, double density, int verbosity);
double SampleCDFCoalTimeMLCFromXtoYLineages(double w_time, int n_arg, va_list ap);
double LogscaleSampleCDFCoalTimeMLCFromXtoYLineages(double w_time, int n_arg, va_list ap);

///@}
/**
 * \name Gamma related functions
 *******************************************************************************/
///@{

/**
 * Random sampling of a Gamma distribution with scale = 1
 *
 * Edited from PAML
 *
 * \param s
 *   Shape.
 * \param seed
 *   Seed for the pseudorandom number generator.
 * \return
 *   Sample
 *  \note If an error ocurrs, it exits by \ref ErrorReporter.
 *******************************************************************************/

double RndGamma (double s, long int *seed);

/**
 * Random sampling of a Gamma distribution with expected mean 1
 *
 * Edited from PAML
 *
 * \param s
 *   Shape.
 * \param seed
 *   Seed for the pseudorandom number generator.
 * \return
 *   Sample
 *  \note If an error ocurrs, it exits by \ref ErrorReporter.
 *******************************************************************************/

double RndGammaE1 (double s, long int *seed);

///@}
/**
 *\name Logscale functions
 *******************************************************************************/
///@{
/**
 * Logscale addition log(x+y) given logx and logy
 *
 * \param logx
 *   First term (log).
 * \param logy
 *   Second term (log).
 * \return
 *   Result of the addition.
 *******************************************************************************/
double LogscaleAdd(double logx, double logy);

/**
 * Logscale addition log(x+y) given logx and logy using Kahan summation to reduce precision problems.
 *
 * \param logx
 *   First term (log).
 * \param logy
 *   Second term (log).
 * \param comp_sum
 *   Pointer to a double to perform the compensated summation (previously lost, updated for a next iteration)
 * \return
 *   Result of the addition.
 *******************************************************************************/
double KahanLogscaleAdd(double logx, double logy, double * comp_sum);

/**
 * Logscale substraction log(x-y) given logx and logy
 *
 * \param logx
 *   First term (log).
 * \param logy
 *   Second term (log).
 * \return
 *   Result of the substraction.
 *******************************************************************************/
double LogscaleSub(double logx, double logy);

/**
 * Logscale factorial log(log(x)!) given logx using tabulated values form 1 to 256 and then Stirling's approximation
 *
 * Adapted from from http://www.johndcook.com/blog/2010/08/16/how-to-compute-log-factorial/
 *
 * \param logx
 *   Value.
 * \return
 *   Result of the factorial.
 *******************************************************************************/
double LogscaleFact(int n);

///@}

/**
 * \name Miscelanea
 *******************************************************************************/
///@{

/**
 * MINSTD Pseudorandom number generator.
 *
 * Simple linear congruence algorithm with certain a and m values.
 *
 * Described in:
 * Park, S. K. and K. W. Miller.  1988.  Random number generators: good
 * ones are hard to find.  Communications of the ACM, 31(10):1192-1201.
 *
 * \todo Change it to use a more sophisticated method, like SFMT.
 *
 * \param seed
 *  Seed to generate numbers.
 * \return Pseudorandom number from 0 to 1.
 *******************************************************************************/
double RandomNumber (long int *seed);

/**
 * Function to obtain the number of digits of a int, taking or not into account the sign.
 *
 * \param val
 *   Value.
 * \param sign
 *   Logical value, 1 = take the sign into account. 
 * \return
 *   Number of digits.
 * \note The int should be cast to long, to avoid the exception if val=int_max
 *******************************************************************************/

long count_intdigits(long val, int sign);

/**
 * Function to compare doubles.
 *
 * This function is used in the qsort (quicksort algorithm) to compare doubles
 *  and to order it in ascending order.
 *
 * \param n1
 *  Pointer to the first double to compare.
 * \param n2
 *  Pointer to the second double to compare.
 * \return 0 if the numbers are equal (difference less than DBL_EPSILON), 1 if
 *  n1>n2 and -1 if n1<n2.
 *******************************************************************************/
int Compare_DBL (const void * n1, const void *n2);

/**
 * Function to compare period::r_bound.
 *
 * This function is used in the qsort (quicksort algorithm) to order periods.
 *
 * \param n1
 *  Pointer to the first period to compare.
 * \param n2
 *  Pointer to the second period to compare.
 * \return 0 if the numbers are equal (difference less than DBL_EPSILON), 1 if
 *  n1>n2 and -1 if n1<n2.
 *******************************************************************************/
int Compare_periods (const void * n1, const void *n2);

/**
 * Function to sample values from an array, proportional to their values.
 *
 * \param n
 *  Number of values.
 * \param array
 *  Pointer to the array of values to sample from.
 * \param seed
 *  Random number generator seed.
 * \return Array index of the sampled value.
 *******************************************************************************/
size_t SampleNDoubles(size_t n, double * array, gsl_rng *seed);

/**
 * Performs a compensated sumation using Kahan's algorithm across an array of doubles.
 *
 * \param array
 *  Pointer to the array with the values.
 * \param n_elements
 *  Number of elements to add.
 * \return Result.
 *******************************************************************************/
double VKahanSum(double * array, int n_elements);

/**
 * Performs a compensated sumation using Kahan's algorithm externally, to use in loops.
 *
 * \param sum
 *  Pointer to the provisional result.
 * \param input
 *  Term to add.
 * \param compensations
 *  Pointer to a variable where store the precision lost. It gets updated in each iteration and must be initializated to 0 in the first.
 *******************************************************************************/
void SKahanSum(double *sum, double input, double *compensation);
///@}
///@}

