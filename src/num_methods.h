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
#ifndef trees_h
#include "trees.h"
#endif

//MPFR Library
#ifndef  _STDARG_H
#include <stdarg.h>
#endif
#include <mpfr.h>
#include <float.h>


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
 * \param c1
 *   Constant 1.
 * \param c2
 *   Constant 2.
 * \param c3
 *   Constant 3.
 * \param c4
 *   Constant 4.
 * \param c5
 *   Constant 5.
 * \param epsilon
 *   Convergence value.
 * \param verbosity
 *   How many info to display by stdout.
 * \return
 *   Root.
 *  \note If an error ocurrs, it exits by \ref ErrorReporter.
 *******************************************************************************/

extern double brent_root (double (*function) (double value, double c1, unsigned int c2,unsigned int c3, double c4, int c5),double a, double b,double c1, unsigned int c2, unsigned int c3, double c4, int c5, float epsilon,unsigned int verbosity);
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

extern double cdf_bounded_coalescent(double w_time, unsigned int n_leaves,unsigned int pop_size, double bound_time);

/**
 * Sampler of a bounded coalescent process.
 *
 * \param w_time 
 *   Waiting time.
 * \param density
 *   Density (sampled uniform random number) of the w_time to find.
 * \param n_leaves
 *   Original number of leaves to coalesce.
 * \param pop_size
 *   Ne of the branch.
 * \param bound_time
 *   Bound of the process (time).
 * \param precision
 *   Number of bits of the MPFR mantissa.
 * \return
 *   Cumulative density.
 *  \note If an error ocurrs, it exits by \ref ErrorReporter.
 *******************************************************************************/

extern double sample_bounded_coalescent(double w_time, double density, unsigned int n_leaves, unsigned int pop_size, double bound_time, int precision);

/**
 * Arbitrary-precision CDF function of a bounded coalescent process.
 *
 * \param time
 *   Waiting time.
 * \param n_leaves
 *   Original number of leaves to coalesce.
 * \param pop_size
 *   Ne of the branch.
 * \param bound_time
 *   Bound of the process (time).
 * \param precision
 *   Desired precision (manatissa's bits).
 * \return
 *   Cumulative density.
 *  \note If an error ocurrs, it exits by \ref ErrorReporter.
 *******************************************************************************/

extern long double cdf_bounded_coalescent_mpfr(long double time, unsigned int n_leaves,unsigned int pop_size, long double bound_time, int precision);

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

extern double RndGamma (double s, long int *seed);

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

extern double RndGammaE1 (double s, long int *seed);

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
extern double RandomNumber (long int *seed);

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

extern long count_intdigits(long val, int sign);

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
extern int Compare_DBL (const void * n1, const void *n2);

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
extern int Compare_periods (const void * n1, const void *n2);

///@}
///@}

