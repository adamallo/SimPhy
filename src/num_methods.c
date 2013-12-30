/**
* 
* \file num_methods.c
* Source of the library of numerical methods.
*******************************************************************************/

#include "num_methods.h"

// **** Prototipes of private functions **** //

/**
 * Random sampling of a Gamma distribution with scale = 1 and small (<1) shapes.
 *
 * Taken from PAML
 *
 * \param s
 *   Shape.
 * \param seed
 *   Seed for the pseudorandom number generator.
 * \return
 *   Sample
 *  \note If an error ocurrs, it exits by \ref ErrorReporter.
 *******************************************************************************/

double RndGamma1 (double s, long int *seed);

/**
 * Random sampling of a Gamma distribution with scale = 1 and big (>1) normal shapes.
 *
 * Taken from PAML
 *
 * \param s
 *   Shape.
 * \param seed
 *   Seed for the pseudorandom number generator.
 * \return
 *   Sample
 *  \note If an error ocurrs, it exits by \ref ErrorReporter.
 *******************************************************************************/

double RndGamma2 (double s, long int *seed);

// **** Declarations of public functions **** //

/** 
 * \name Root finding methods 
 *******************************************************************************/
///@{

double brent_root (double (*f)(double value, int n, va_list ap), double a, double b, float epsilon, int verbosity, int n_arg, ...)
{
    
    va_list ap;
    int i=0;
    double prev_b,bprev_b,new_val;
    double fa,fb,fprev_b,fnew_val,initial_a=a, initial_b=b;
    double cond1;
    
    int prev_bis=1;
    
    //printf("Density %lf ",c1); //DBG
    va_start(ap,n_arg);
    fa=(*f)(a,n_arg,ap);
    va_start(ap,n_arg);
    fb=(*f)(b,n_arg,ap);
    
    //printf("Test %lf, %lf\n",fa,fb);

    //Checking a and b points
    if (fa*fb>=0)
    {
        fprintf(stderr,"Root-finding function error, root (0) not bracketed\n");
        ErrorReporter(UNEXPECTED_VALUE);
    }
    
    if (fabs(fa)<fabs(fb))
    {
        prev_b=a;
        a=b;
        b=prev_b;
        fprev_b=fa;
        fa=fb;
        fb=fprev_b;
    }
    
    prev_b=a;
    fprev_b=fa;
    
    //Algorithm
    
    while (((fb!=0) && (fabs(a-b) > epsilon)) || (b==initial_a || b==initial_b)) //Convergence test. Remove "|| (b==initial_a || b==initial_b)" to get the original Brent's method.
    {
        //printf("Iteration %lu\n",i);
        if ((fa!=fprev_b) && (fb!=fprev_b)) //Inverse cuadratic interpolation
        {
            new_val = a*fb*fprev_b/(fa-fb)/(fa-fprev_b)+b*fa*fprev_b/(fb-fa)/(fb-fprev_b)+prev_b*fa*fb/(fprev_b-fa)/(fprev_b-fb);
        }
        else //Secant method
        {
            new_val=b-fb*(b-a)/(fb-fa);
        }
        
        cond1=(3*a+b)/4;
        
        //Bisection method
        if ((!(((new_val > cond1) && (new_val < b)) || ((new_val < cond1) && (new_val > b)))) || (prev_bis && (fabs(new_val - b) >= (fabs(b - prev_b) / 2))) || (!prev_bis && (fabs(new_val - b) >= (fabs(prev_b - bprev_b) / 2))))
        {
            new_val = (a + b) / 2;
            prev_bis = 1;
        }
        else
        {
            if ((prev_bis && (fabs(b - prev_b) < epsilon)) || (!prev_bis && (fabs(prev_b - bprev_b) < epsilon)))
            {
                new_val = (a + b) / 2;
                prev_bis = 1;
            }
            else
                prev_bis = 0;
        }
        
        //New values
        va_start(ap,n_arg);
        fnew_val=(*f)(new_val,n_arg,ap);
                
        bprev_b=prev_b;
        prev_b=b;
        
        if (fa*fnew_val<0)
        {
            b=new_val;
            fb=fnew_val;
        }
        else
        {
            a=new_val;
            fa=fnew_val;
        }
        
        if (fabs(fa)<fabs(fb))
        {
            prev_b=a;
            a=b;
            b=prev_b;
            fprev_b=fa;
            fa=fb;
            fb=fprev_b;
        }
        
        ++i;
        
        if(i>MAX_IT)
            ErrorReporter(LOOP_ERROR);
        
    }
    return(b);
}
//@}

/** 
 * \name Coalescent theory related numerical methods 
 *******************************************************************************/
///@{

double cdf_bounded_coalescent(double w_time, int n_leaves,int pop_size, double bound_time)
{
	int i=0;
    int pn_leaves=n_leaves;
    long double ltime=(long double)w_time, lbound_time=(long double)bound_time; //Increasing precision
    long double cak1=1, cak2=1,sum1=0, sum2=0, f1=0,f2=0,c1=0,leaves_i=0;
    long double error1=0,error2=0, next_sum1=0, next_sum2=0, addend1=0, addend2=0; //Compensated summation
    
    for (i=1;i<pn_leaves;++i)
    {
        c1=expl(-lbound_time*((i*(i-1))/(long double)(2*pop_size))) * (2*i-1);
        f1= c1* cak1;
        f2= c1* cak2;
        leaves_i= i*(i-1)-pn_leaves * (pn_leaves-1);
        
        //Compensated summation
        addend1=  (f1 * (expl(ltime*leaves_i/(2*pop_size))-1)/leaves_i)-error1;
        next_sum1=sum1+addend1;
        error1=(next_sum1-sum1)-addend1;
        addend2= f2-error2;
        next_sum2=sum2+addend2;
        error2=(next_sum2-sum2)-addend2;
        
        sum1=next_sum1;
        sum2=next_sum2;
        cak1*= (i-pn_leaves+1)/(long double)(i+pn_leaves-1);
        cak2*= (i-pn_leaves)/(long double)(i+pn_leaves);
    }
    sum2+=(expl (-lbound_time * pn_leaves * (pn_leaves-1)/(2*pop_size)) * (2 * pn_leaves - 1 ) * cak2)-error2;
    
   	return (double)(pn_leaves/sum2) * ((pn_leaves-1)*sum1);
}

//long double cdf_bounded_coalescent_mpfr(long double time, int n_leaves,int pop_size, long double bound_time, int precision)
//{
//	//Variable initialization
//	
//	//Normal variables
//	long int i=0;
//    long int pn_leaves=n_leaves, leaves_i;
//    long double result_ldouble=0;
//    
//    //Mpfr declaration and initialization
//    mpfr_t cak1,cak2,cak_fact,sum1,sum2,f1,f1_fact,f2,c1,mpfr_time,mpfr_bound_time,result;
//    mpfr_inits2(precision,cak1,cak2,cak_fact,sum1,sum2,f1,f2,c1,f1_fact,mpfr_time,mpfr_bound_time,result,(mpfr_ptr)0);
//    mpfr_set_zero(sum1,1);
//    mpfr_set_zero(sum2,1);
//    mpfr_set_zero(f1,1);
//    mpfr_set_zero(f2,1);
//    mpfr_set_zero(c1,1);
//    mpfr_set_si(cak1,1,MPFR_RNDN);
//    mpfr_set_si(cak2,1,MPFR_RNDN);
//    mpfr_set_ld(mpfr_time,time,MPFR_RNDN);
//    mpfr_set_ld(mpfr_bound_time,bound_time,MPFR_RNDN);
//    //long double error1=0,error2=0, next_sum1=0, next_sum2=0, addend1=0, addend2=0; //Compensated summation
//    
//    for (i=1;i<pn_leaves;++i)
//    {
//    	//C1
//    	mpfr_set_si(c1,i*(i-1),MPFR_RNDN);
//		mpfr_div_ui(c1,c1,(2*pop_size),MPFR_RNDN);
//		mpfr_mul(c1,c1,mpfr_bound_time,MPFR_RNDN);
//		mpfr_mul_si(c1,c1,-1,MPFR_RNDN);
//		mpfr_exp(c1,c1,MPFR_RNDN);
//		mpfr_mul_si(c1,c1,(2*i-1),MPFR_RNDN);
//		
//		//F1 && F2
//        mpfr_mul(f1,c1,cak1,MPFR_RNDN);
//        mpfr_mul(f2,c1,cak2,MPFR_RNDN);
//        
//        leaves_i= i*(i-1)-pn_leaves * (pn_leaves-1);
//        
//        //  //Compensated summation
//        //  addend1=  (f1 * (expl(time*leaves_i/(2*pop_size))-1)/leaves_i)-error1;
//        //  next_sum1=sum1+addend1;
//        //  error1=(next_sum1-sum1)-addend1;
//        //  addend2= f2-error2;
//        //  next_sum2=sum2+addend2;
//        //  error2=(next_sum2-sum2)-addend2;
//        //
//        //  sum1=next_sum1;
//        //  sum2=next_sum2;
//        
// 		//F1 fact
// 		mpfr_set_si(f1_fact,leaves_i,MPFR_RNDN);
// 		mpfr_div_ui(f1_fact,f1_fact,(2*pop_size),MPFR_RNDN);
// 		mpfr_mul(f1_fact,f1_fact,mpfr_time,MPFR_RNDN);
// 		mpfr_exp(f1_fact,f1_fact,MPFR_RNDN);
// 		mpfr_sub_ui(f1_fact,f1_fact,1,MPFR_RNDN);
// 		mpfr_div_si(f1_fact,f1_fact,leaves_i,MPFR_RNDN);
// 		mpfr_mul(f1_fact,f1,f1_fact,MPFR_RNDN);
// 		
// 		//Summatories
// 		mpfr_add(sum1,sum1,f1_fact,MPFR_RNDN);
//		mpfr_add(sum2,sum2,f2,MPFR_RNDN);
//       	//New cak1
//       	mpfr_set_si(cak_fact,i-pn_leaves+1,MPFR_RNDN);
//       	mpfr_div_si(cak_fact,cak_fact,i+pn_leaves-1,MPFR_RNDN);
//       	mpfr_mul(cak1,cak1,cak_fact,MPFR_RNDN);
//       	
//       	//New cak2
//       	mpfr_set_si(cak_fact,i-pn_leaves,MPFR_RNDN);
//       	mpfr_div_si(cak_fact,cak_fact,i+pn_leaves,MPFR_RNDN);
//       	mpfr_mul(cak2,cak2,cak_fact,MPFR_RNDN);
//    }
//    
//    //Last Sum2
//    mpfr_set_si(c1,pn_leaves*(pn_leaves-1),MPFR_RNDN);
//	mpfr_div_ui(c1,c1,(2*pop_size),MPFR_RNDN);
//	mpfr_mul(c1,c1,mpfr_bound_time,MPFR_RNDN);
//	mpfr_mul_si(c1,c1,-1,MPFR_RNDN);
//	mpfr_exp(c1,c1,MPFR_RNDN);
//	mpfr_mul_si(c1,c1,(2*pn_leaves-1),MPFR_RNDN);
//	
//	mpfr_mul(f2,c1,cak2,MPFR_RNDN);
//	
//    mpfr_add(sum2,sum2,f2,MPFR_RNDN);
//    
//    mpfr_set_si(result,pn_leaves,MPFR_RNDN);
//    mpfr_div(result,result,sum2,MPFR_RNDN);
//    mpfr_mul_si(sum1,sum1,pn_leaves-1,MPFR_RNDN);
//    mpfr_mul(result,result,sum1,MPFR_RNDN);
//    
//    result_ldouble=mpfr_get_ld(result,MPFR_RNDN);
//    mpfr_clears(cak1,cak2,cak_fact,sum1,sum2,f1,f1_fact,f2,c1,mpfr_time,mpfr_bound_time,result,(mpfr_ptr)0);
//    mpfr_free_cache();
//    
//    
//   	return result_ldouble;
//}

//double sample_bounded_coalescent(double time, double density, int n_leaves, int pop_size, double bound_time, int precision)
//{
//    long double ltime=(long double)time, ldensity=(long double) density, lbound_time=(long double)bound_time; //Increasing precision
//    
//    if (time <=0)
//    {
//        return time - density; //CDF =0. Way of differenciate all the results of time<=0, the more deviated from 0 the lower resulting value.
//    }
//    else if (time >= bound_time)
//    {
//        return (double)(1 - ldensity + (ltime - lbound_time)); //CDF =1. Way of differenciate all the results of time>=bound_time, the more diviated from bound_time, the higher resulting value.
//    }
//    else
//    {
//        return(double)(cdf_bounded_coalescent_mpfr(ltime, n_leaves, pop_size, lbound_time, precision)-ldensity);
//    }
//}

double ProbCoalFromXtoYLineages(int i_lin,int o_lin,double bound_time,int Ne)
{
    double C=1, sum=0;
    int i=0, pi=0;
    
    for (i=0;i<o_lin;++i)
    {
        C*=(o_lin+i)*(i_lin-i)/(double)(i_lin+i);
    }
    
    sum=exp(-o_lin*(o_lin-1)*bound_time/2/Ne)*C;
    
    for (i=o_lin+1; i<=i_lin;++i)
    {
        pi=i-1;
        C*=(double)(o_lin+pi)*(i_lin-pi)/(i_lin+pi)/(o_lin-i);
        sum+=exp(-i*pi*bound_time/2/Ne)*(2*i-1)/(pi+o_lin)*C;
    }
    
    if (o_lin<GSL_SF_FACT_NMAX)
        return sum/gsl_sf_fact(o_lin);
    else
    {
        fprintf(stderr,"Error using a factorial function. A big numbers factorial library should be used\n");
#ifdef DBG
        fflush(stderr);
#endif
        ErrorReporter(UNEXPECTED_VALUE);
        return -1;
    }
}

double SampleCoalTimeMLCFromXtoYLineages(int i_lin, int o_lin, double bound_time, int pop_size, double brent_epsilon, double density, int verbosity)
{
    double lambda_i=0, lambda_o=0, C0=1, out_sum=0;
    int i=0;
    
    lambda_i= -i_lin*(i_lin-1)/(double)(2*pop_size);
    lambda_o= -o_lin*(o_lin-1)/(double)(2*pop_size);
    for (i=0;i<o_lin;++i)
    {
        C0*=(o_lin+i)*(i_lin-i-1)/(double)(i_lin-1+i);
    }
    
    if (o_lin>GSL_SF_FACT_NMAX)
    {
        fprintf(stderr,"Error using a factorial function. A big numbers factorial library should be used\n");
#ifdef DBG
        fflush(stderr);
#endif
        ErrorReporter(UNEXPECTED_VALUE);
        return -1;
    }
    
    out_sum=(-lambda_i)/(gsl_sf_fact(o_lin)*ProbCoalFromXtoYLineages(i_lin, o_lin, bound_time, pop_size)); /// \todo test if this is the proper order to avoid precision errors.
    
    return brent_root(*(SampleCDFCoalTimeMLCFromXtoYLineages),0.0,bound_time,brent_epsilon,verbosity,9,i_lin,o_lin,bound_time,pop_size,lambda_i,lambda_o,C0,out_sum,density);
}


double SampleCDFCoalTimeMLCFromXtoYLineages(double w_time, int n_arg, va_list ap)
{
    double C=0,sum=0,lambda_k=0,bound_time=0,lambda_i=0,lambda_o=0,out_sum=0,density=0;
    int k=0,k1=0,i=0,i_lin=0,o_lin=0,pop_size=0;
    
    if (n_arg!=9)
        ErrorReporter(UNEXPECTED_VALUE);
    
    for (i=0; i<n_arg; ++i)
    {
        switch (i)
        {
            case 0:
                i_lin=va_arg(ap,int);
                break;
            case 1:
                o_lin=va_arg(ap,int);
                break;
            case 2:
                bound_time=va_arg(ap,double);
                break;
            case 3:
                pop_size=va_arg(ap,int);
                break;
            case 4:
                lambda_i=va_arg(ap,double);
                break;
            case 5:
                lambda_o=va_arg(ap,double);
                break;
            case 6:
                C=va_arg(ap,double);
                break;
            case 7:
                out_sum=va_arg(ap,double);
                break;
            case 8:
                density=va_arg(ap,double);
                break;
            default:
                ErrorReporter(UNEXPECTED_VALUE);
                break;
                
        }
    }
    va_end(ap);

    k1=o_lin-1;
    
    if (w_time<=0)
    {
        return w_time - density;
    }
    else if (w_time>=bound_time)
    {
        return 1-density + (w_time - bound_time);
    }
    else
    {
        sum= exp(lambda_o*bound_time) * (exp((lambda_i-lambda_o)*w_time)-1)/ (lambda_i-lambda_o)*C; //First iteration of the summatory.
        for (k=o_lin+1; k<i_lin; ++k) //K = number of lineages
        {
            ++k1;
            lambda_k=-k*k1/(double)(2*pop_size);
            C=(o_lin+k1)*(i_lin-1-k1)/(double)(i_lin-1+k1)/(o_lin-k)*C;
            sum += exp(lambda_k*bound_time)*(exp((lambda_i-lambda_k)*w_time)-1)/(lambda_i-lambda_k)*(2*k-1)/(k1+o_lin)*C;
        }
        return sum*out_sum-density;
    }
}
//@}

/**
 * \name Gamma related functions
 *******************************************************************************/
///@{

double RndGamma (double s, long int *seed)
{
    
	double  r=0.0;
    int i=0;
    
	
	if (s <= 0.0)
		ErrorReporter(UNEXPECTED_VALUE);
    else if (ceilf(s) == s) //Integer; in this context, natural.
    {
        for (i=0; i<s;++i)//Sampling of s exp(1), appling the alpha addition.
        {
            r+= -log(RandomNumber(seed));
        }
    }
	else if (s < 1.0)
		r = RndGamma1 (s, seed);
	else if (s > 1.0)
		r = RndGamma2 (s, seed);

	return (r);
    
}

double RndGammaE1 (double s, long int *seed)
{
    
	double  r=0.0;
    int i=0;
    
	
	if (s <= 0.0)
		ErrorReporter(UNEXPECTED_VALUE);
    else if (ceilf(s) == s) //Integer; in this context, natural.
    {
        for (i=0; i<s;++i)//Sampling of s exp(1), appling the alpha addition.
        {
            r+= -log(RandomNumber(seed));
        }
    }
	else if (s < 1.0)
		r = RndGamma1 (s, seed);
	else if (s > 1.0)
		r = RndGamma2 (s, seed);
    
	return (r/s);
    
}

double RndGamma1 (double s, long int *seed)
{
    
	double			r, x=0.0, small=1e-37, w;
	static double   a, p, uf, ss=10.0, d;
	
	if (s!=ss)
    {
		a  = 1.0-s;
		p  = a/(a+s*exp(-a));
		uf = p*pow(small/a,s);
		d  = a*log(a);
		ss = s;
    }
	for (;;)
    {
		r = RandomNumber(seed);
		if (r > p)
			x = a-log((1.0-r)/(1.0-p)), w=a*log(x)-d;
		else if (r>uf)
			x = a*pow(r/p,1/s), w=x;
		else
			return (0.0);
		r = RandomNumber(seed);
		if (1.0-r <= w && r > 0.0)
            if (r*(w+1.0) >= 1.0 || -log(r) <= w)
                continue;
		break;
    }
	return (x);
    
}

double RndGamma2 (double s, long int *seed)
{
    
	double		r ,d, f, g, x;
	static double	b, h, ss=0;
	
	if (s!=ss)
    {
		b  = s-1.0;
		h  = sqrt(3.0*s-0.75);
		ss = s;
    }
	for (;;)
    {
		r = RandomNumber(seed);
		g = r-r*r;
		f = (r-0.5)*h/sqrt(g);
		x = b+f;
		if (x <= 0.0)
			continue;
		r = RandomNumber(seed);
		d = 64*r*r*g*g*g;
		if (d*x < x-2.0*f*f || log(d) < 2*(b*log(x/b)-f))
			break;
    }
	return (x);
    
}
///@}

/**
 * \name Miscellanea
 *******************************************************************************/
///@{

double RandomNumber (long int *seed)

{
    
	long int	lo, hi, test;
    
	hi = (*seed) / 127773;
	lo = (*seed) % 127773;
	test = 16807 * lo - 2836 * hi;
	if (test > 0)
		*seed = test;
	else
		*seed = test + 2147483647;
	return (double)(*seed) / (double)2147483647;
    
}

long count_intdigits(long val, int sign)
{
	long d = 1, c;
    
	if (val >= 0) for (c = 10; c <= val; c *= 10) d++;
	else for (c= -10 ; c >= val; c *= 10) d++;
	return (sign && c < 0) ? ++d : d;
}

int Compare_DBL (const void * n1, const void *n2)
{
    double dif=0;
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Absolute difference between n1 and n2 (abs for double hasn't been implemented in libc)
    if(*(double *)n1-*(double *)n2>0)
        dif=*(double *)n1-*(double *)n2;
        else
            dif=*(double *)n2-*(double *)n1;
            
            // *
            /// If this difference is less than DBL_EPSILON, n1=n2 (It returns 0)
            if (dif<DBL_EPSILON)
            {
                return 0;
            }
    // *
    /// Else, it compares the values and returns 1 (n1>n2) or -1(n1<n2)</dl>
            else
            {
                if (*(double *)n1>*(double *)n2)
                    return 1;
                else
                    return -1;
            }
}

int Compare_periods (const void * n1, const void *n2)
{
    double dif=0, v1=0,v2=0;
    
    v1=((period *)n1)->r_bound;
    v2=((period *)n2)->r_bound;
    
    // **
    /// <dl><dt> Function structure </dt><dd>
    
    // *
    /// Absolute difference between n1 and n2 (abs for double hasn't been implemented in libc)
    if(v1-v2>0)
        dif=v1-v2;
    else
        dif=v2-v1;
    
    // *
    /// If this difference is less than DBL_EPSILON, n1=n2 (It returns 0)
    if (dif<DBL_EPSILON)
    {
        return 0;
    }
    // *
    /// Else, it compares the values and returns 1 (n1>n2) or -1(n1<n2)</dl>
    else
    {
        if (v1>v2)
            return 1;
        else
            return -1;
    }
}

inline size_t SampleNDoubles(size_t n, double * array, gsl_rng *seed)
{
    size_t i=0;
    double sum=0, value=0;
    
    for (i=0; i<n;++i)
    {
        sum+=*(array+i);
    }
    
    value=gsl_rng_uniform_pos(seed)*sum;
    
    sum=0;
    for (i=0; i<n;++i)
    {
        sum+=*(array+i);
        if (sum>=value)
            break;
    }
        
    return i;
}


///@}