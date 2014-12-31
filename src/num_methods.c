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

static double RndGamma1 (double s, long int *seed);

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

static double RndGamma2 (double s, long int *seed);

static double chebyshev_eval (double x, const double *a, const int n);
static double biomcmc_log1p(double x);

// **** Declarations of public functions **** //

/** 
 * \name Root finding methods 
 *******************************************************************************/
///@{

long int brent_root (double (*f)(double value, int n, va_list ap), double a, double b, float epsilon, double *result, int verbosity, int n_arg, ...)
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
        //ErrorReporter(UNEXPECTED_VALUE,"Root-finding function error, root (0) not bracketed\n");
        fprintf(stderr, "Root-finding function error, root (0) not bracketed\n");
#ifdef DBG
        fflush(stderr);
#endif
        *result=b;
        return UNEXPECTED_VALUE;
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
        //printf("%lf,%lf,%lf,%lf\n",a,b,fb,fa);
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
        {
            *result=b;
            return UNEXPECTED_VALUE;
        }
        
    }
    *result=b;
    return NO_ERROR;
}
//@}

/** 
 * \name Coalescent theory related numerical methods 
 *******************************************************************************/
///@{

double CdfBoundedCoalescent(double w_time, int n_leaves,int pop_size, double bound_time)
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

#ifdef __MPFR_H
long double cdf_bounded_coalescent_mpfr(long double time, int n_leaves,int pop_size, long double bound_time, int precision)
{
	//Variable initialization
	
	//Normal variables
	long int i=0;
    long int pn_leaves=n_leaves, leaves_i;
    long double result_ldouble=0;
    
    //Mpfr declaration and initialization
    mpfr_t cak1,cak2,cak_fact,sum1,sum2,f1,f1_fact,f2,c1,mpfr_time,mpfr_bound_time,result;
    mpfr_inits2(precision,cak1,cak2,cak_fact,sum1,sum2,f1,f2,c1,f1_fact,mpfr_time,mpfr_bound_time,result,(mpfr_ptr)0);
    mpfr_set_zero(sum1,1);
    mpfr_set_zero(sum2,1);
    mpfr_set_zero(f1,1);
    mpfr_set_zero(f2,1);
    mpfr_set_zero(c1,1);
    mpfr_set_si(cak1,1,MPFR_RNDN);
    mpfr_set_si(cak2,1,MPFR_RNDN);
    mpfr_set_ld(mpfr_time,time,MPFR_RNDN);
    mpfr_set_ld(mpfr_bound_time,bound_time,MPFR_RNDN);
    //long double error1=0,error2=0, next_sum1=0, next_sum2=0, addend1=0, addend2=0; //Compensated summation
    
    for (i=1;i<pn_leaves;++i)
    {
    	//C1
    	mpfr_set_si(c1,i*(i-1),MPFR_RNDN);
		mpfr_div_ui(c1,c1,(2*pop_size),MPFR_RNDN);
		mpfr_mul(c1,c1,mpfr_bound_time,MPFR_RNDN);
		mpfr_mul_si(c1,c1,-1,MPFR_RNDN);
		mpfr_exp(c1,c1,MPFR_RNDN);
		mpfr_mul_si(c1,c1,(2*i-1),MPFR_RNDN);
		
		//F1 && F2
        mpfr_mul(f1,c1,cak1,MPFR_RNDN);
        mpfr_mul(f2,c1,cak2,MPFR_RNDN);
        
        leaves_i= i*(i-1)-pn_leaves * (pn_leaves-1);
        
 		//F1 fact
 		mpfr_set_si(f1_fact,leaves_i,MPFR_RNDN);
 		mpfr_div_ui(f1_fact,f1_fact,(2*pop_size),MPFR_RNDN);
 		mpfr_mul(f1_fact,f1_fact,mpfr_time,MPFR_RNDN);
 		mpfr_exp(f1_fact,f1_fact,MPFR_RNDN);
 		mpfr_sub_ui(f1_fact,f1_fact,1,MPFR_RNDN);
 		mpfr_div_si(f1_fact,f1_fact,leaves_i,MPFR_RNDN);
 		mpfr_mul(f1_fact,f1,f1_fact,MPFR_RNDN);
 		
 		//Summatories
 		mpfr_add(sum1,sum1,f1_fact,MPFR_RNDN);
		mpfr_add(sum2,sum2,f2,MPFR_RNDN);
       	//New cak1
       	mpfr_set_si(cak_fact,i-pn_leaves+1,MPFR_RNDN);
       	mpfr_div_si(cak_fact,cak_fact,i+pn_leaves-1,MPFR_RNDN);
       	mpfr_mul(cak1,cak1,cak_fact,MPFR_RNDN);
       	
       	//New cak2
       	mpfr_set_si(cak_fact,i-pn_leaves,MPFR_RNDN);
       	mpfr_div_si(cak_fact,cak_fact,i+pn_leaves,MPFR_RNDN);
       	mpfr_mul(cak2,cak2,cak_fact,MPFR_RNDN);
    }
    
    //Last Sum2
    mpfr_set_si(c1,pn_leaves*(pn_leaves-1),MPFR_RNDN);
	mpfr_div_ui(c1,c1,(2*pop_size),MPFR_RNDN);
	mpfr_mul(c1,c1,mpfr_bound_time,MPFR_RNDN);
	mpfr_mul_si(c1,c1,-1,MPFR_RNDN);
	mpfr_exp(c1,c1,MPFR_RNDN);
	mpfr_mul_si(c1,c1,(2*pn_leaves-1),MPFR_RNDN);
	
	mpfr_mul(f2,c1,cak2,MPFR_RNDN);
	
    mpfr_add(sum2,sum2,f2,MPFR_RNDN);
    
    mpfr_set_si(result,pn_leaves,MPFR_RNDN);
    mpfr_div(result,result,sum2,MPFR_RNDN);
    mpfr_mul_si(sum1,sum1,pn_leaves-1,MPFR_RNDN);
    mpfr_mul(result,result,sum1,MPFR_RNDN);
    
    result_ldouble=mpfr_get_ld(result,MPFR_RNDN);
    mpfr_clears(cak1,cak2,cak_fact,sum1,sum2,f1,f1_fact,f2,c1,mpfr_time,mpfr_bound_time,result,(mpfr_ptr)0);
    mpfr_free_cache();
    
    
   	return result_ldouble;
}
#endif

double SampleBoundedCoalescent(double w_time, int n_arg, va_list ap)
{
    int i=0,n_leaves=0, pop_size=0;
    double density=0,bound_time=0;
    
    if (n_arg!=4)
        ErrorReporter(UNEXPECTED_VALUE,NULL);
    
    for (i=0; i<n_arg; ++i)
    {
        switch (i)
        {
            case 0:
                density=va_arg(ap,double);
                break;
            case 1:
                n_leaves=va_arg(ap,int);
                break;
            case 2:
                pop_size=va_arg(ap,int);
                break;
            case 3:
                bound_time=va_arg(ap,double);
                break;
            default:
                ErrorReporter(UNEXPECTED_VALUE,NULL);
                break;
                
        }
    }
    va_end(ap);
    
    if (w_time <=0)
    {
        return w_time - density; //CDF =0. Way of differenciate all the results of time<=0, the more deviated from 0 the lower resulting value.
    }
    else if (w_time >= bound_time)
    {
        return (1 - density + (w_time - bound_time)); //CDF =1. Way of differenciate all the results of time>=bound_time, the more diviated from bound_time, the higher resulting value.
    }
    else
    {
        return CdfBoundedCoalescent(w_time, n_leaves, pop_size, bound_time)-density;
        ///\todo CHECK SOMETHING AND USE MPFR IF NECESSARY???
    }
}

double OriginalProbCoalFromXtoYLineages (int i_lin, int o_lin, double bound_time, int Ne)
{
    double sum=0;
    double c=1;
    int i=0,j=0;
    
    for (i=o_lin; i<=i_lin; ++i)
    {
        c=1;
        for (j=0;j<i;++j)
        {
            c*=(o_lin+j)*(i_lin-j)/(double)(i_lin+j);
        }
        if ((i-o_lin)%2==0)
            sum+=exp(-i*(i-1)*bound_time/2.0/Ne) * (2*i-1)/gsl_sf_fact(i-o_lin)/gsl_sf_fact(o_lin)/(i+o_lin-1)*c;
        else
            sum-=exp(-i*(i-1)*bound_time/2.0/Ne) * (2*i-1)/gsl_sf_fact(i-o_lin)/gsl_sf_fact(o_lin)/(i+o_lin-1)*c;
        
    }
    
    return sum;
}

#ifdef __MPFR_H
double MPFROriginalProbcoalFromXtoYLineages (int i_lin, int o_lin, double bound_time, int Ne, int precision)
{
    
    mpfr_t sum,c,c_mult,n_sum,fact;
    int i=0,j=0;
    double sum_ret=0;
    
    mpfr_inits2(precision,sum,c,c_mult,n_sum,fact,(mpfr_ptr)0);
    mpfr_set_si(sum,0,MPFR_RNDN);
    
    for (i=o_lin; i<=i_lin; ++i)
    {
        mpfr_set_si(c, 1, MPFR_RNDN);
        for (j=0;j<i;++j)
        {
            mpfr_set_si(c_mult, (o_lin+j)*(i_lin-j), MPFR_RNDN);
            mpfr_div_si(c_mult, c_mult, i_lin+j, MPFR_RNDN);
            mpfr_mul(c, c, c_mult, MPFR_RNDN);
        }
        mpfr_set_si(n_sum,-i*(i-1),MPFR_RNDN);
        mpfr_mul_d(n_sum, n_sum, bound_time, MPFR_RNDN);
        mpfr_div_si(n_sum,n_sum,2*Ne,MPFR_RNDN);
        mpfr_exp(n_sum, n_sum, MPFR_RNDN);
        mpfr_mul_si(n_sum, n_sum, 2*i-1, MPFR_RNDN);
        mpfr_mul(n_sum,n_sum,c,MPFR_RNDN);
        mpfr_fac_ui(fact, i-o_lin, MPFR_RNDN);
        mpfr_div(n_sum, n_sum, fact, MPFR_RNDN);
        mpfr_fac_ui(fact, o_lin, MPFR_RNDN);
        mpfr_div(n_sum, n_sum, fact, MPFR_RNDN);
        mpfr_div_si(n_sum, n_sum, i+o_lin-1, MPFR_RNDN);
        if ((i-o_lin)%2==0)
        {
            mpfr_add(sum, sum, n_sum, MPFR_RNDN);
        }
        else
        {
            mpfr_sub(sum, sum, n_sum, MPFR_RNDN);
        }
    }
    sum_ret=mpfr_get_ld(sum,MPFR_RNDN);
    mpfr_clears(sum,c,c_mult,n_sum,fact,(mpfr_ptr)0);
    mpfr_free_cache();
    return sum_ret;
}
#endif

//double LogscaleOriginalProbCoalFromXtoYLineages (int i_lin, int o_lin, double bound_time, int Ne) //Intermediate values of sum are negative. The summatory must be splitted like in LogscaleProbCoalFromXtoYLineages
//{
//    double logsum=0;
//    double logc=0;
//    int i=0,j=0;
//    
//    for (i=o_lin; i<=i_lin; ++i)
//    {
//        logc=0;
//        for (j=0;j<i;++j)
//        {
//            logc+=log(o_lin+j)+log(i_lin-j)-log(i_lin+j);
//        }
//        if (o_lin==i)
//            logsum=(-i*(i-1)*bound_time/2.0/Ne)+log(2*i-1)-LogscaleFact(i-o_lin)-LogscaleFact(o_lin)-log(i+o_lin-1)+logc;
//        else if ((i-o_lin)%2==0)
//            logsum=LogscaleAdd(logsum,(-i*(i-1)*bound_time/2.0/Ne)+log(2*i-1)-LogscaleFact(i-o_lin)-LogscaleFact(o_lin)-log(i+o_lin-1)+logc);
//        else
//            logsum=LogscaleSub(logsum,(-i*(i-1)*bound_time/2.0/Ne)+log(2*i-1)-LogscaleFact(i-o_lin)-LogscaleFact(o_lin)-log(i+o_lin-1)+logc);
//    }
//
//    return exp(logsum);
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
    if (sum==0)
    {
        return 0;
    }
    else if (sum==NAN || sum==INFINITY)
    {
        ErrorReporter(UNEXPECTED_VALUE, "");
        return EXIT_FAILURE;
    }
    else if (o_lin<GSL_SF_FACT_NMAX)
        return sum/gsl_sf_fact(o_lin);
    else
    {        
        fprintf(stderr,"Error using a factorial function. A big numbers factorial library should be used\n");
#ifdef DBG
        fflush(stderr);
#endif
        ErrorReporter(UNEXPECTED_VALUE,"");
        return EXIT_FAILURE;
    }
}

double KahanLogscaleProbCoalFromXtoYLineages(int i_lin, int o_lin, double bound_time, int Ne)
{
    double logC=0, logsum_even=0,logsum_odd=0,comp_even=0,comp_odd=0;
    int i=0, pi=0;
    
    
    for (i=0;i<o_lin;++i)
    {
        logC+=log((o_lin+i)*(i_lin-i)/(double)(i_lin+i));
    }
    logsum_even=(-o_lin*(o_lin-1)*bound_time/2/Ne)+logC;
    
    switch (i_lin-o_lin)
    {
        case 0:
            return exp(logsum_even-LogscaleFact(o_lin));
            break;
        default:
            logC+=log((2*o_lin)*(i_lin-o_lin)/(double)(i_lin+o_lin));
            logsum_odd=(-(o_lin+1)*o_lin*bound_time/2/Ne)+biomcmc_log1p(0.5/o_lin)+logC;
            
            for (i=o_lin+2; i<=i_lin;++i)
            {
                pi=i-1;
                logC+=log((o_lin+pi)*(i_lin-pi)/(double)(i_lin+pi)/(i-o_lin));
                if ((i-o_lin)%2==0)
                    logsum_even=KahanLogscaleAdd(logsum_even,(-i*pi*bound_time/2/Ne)+0.693147180559945309417232121458+log((i-0.5)/(pi+o_lin))+logC,&comp_even);
                else
                    logsum_odd=KahanLogscaleAdd(logsum_odd,(-i*pi*bound_time/2/Ne)+0.693147180559945309417232121458+log((i-0.5)/(pi+o_lin))+logC,&comp_odd);
            }
            break;
    }
    if (logsum_odd>=logsum_even)
        return 0;
    else
        return exp(LogscaleSub(logsum_even,logsum_odd)-LogscaleFact(o_lin));
}

double KahanLogscaleLogProbCoalFromXtoYLineages(int i_lin, int o_lin, double bound_time, int Ne)
{
    double logC=0, logsum_even=0,logsum_odd=0,comp_even=0,comp_odd=0;
    int i=0, pi=0;
    
    
    for (i=0;i<o_lin;++i)
    {
        logC+=log((o_lin+i)*(i_lin-i)/(double)(i_lin+i));
    }
    logsum_even=(-o_lin*(o_lin-1)*bound_time/2.0/Ne)+logC;
    
    switch (i_lin-o_lin)
    {
        case 0:
            return exp(logsum_even-LogscaleFact(o_lin));
            break;
        default:
            logC+=log((2*o_lin)*(i_lin-o_lin)/(double)(i_lin+o_lin));
            logsum_odd=(-(o_lin+1)*o_lin*bound_time/2.0/Ne)+biomcmc_log1p(0.5/o_lin)+logC;
            
            for (i=o_lin+2; i<=i_lin;++i)
            {
                pi=i-1;
                logC+=log((o_lin+pi)*(i_lin-pi)/(double)(i_lin+pi)/(i-o_lin));
                if ((i-o_lin)%2==0)
                    logsum_even=KahanLogscaleAdd(logsum_even,(-i*pi*bound_time/2/Ne)+0.693147180559945309417232121458+log((i-0.5)/(pi+o_lin))+logC,&comp_even);
                else
                    logsum_odd=KahanLogscaleAdd(logsum_odd,(-i*pi*bound_time/2/Ne)+0.693147180559945309417232121458+log((i-0.5)/(pi+o_lin))+logC,&comp_odd);
            }
            break;
    }
    if (logsum_odd>=logsum_even)
        return 0;
    else
        return LogscaleSub(logsum_even,logsum_odd)-LogscaleFact(o_lin);
}

double LogscaleProbCoalFromXtoYLineages(int i_lin, int o_lin, double bound_time, int Ne)
{
    double logC=0, logsum_even=0,logsum_odd=0;
    int i=0, pi=0;
    
    
    for (i=0;i<o_lin;++i)
    {
        logC+=log(o_lin+i)+log(i_lin-i)-log(i_lin+i);
    }
    logsum_even=(-o_lin*(o_lin-1)*bound_time/2.0/Ne)+logC;
    
    switch (i_lin-o_lin)
    {
        case 0:
            return exp(logsum_even-LogscaleFact(o_lin));
            break;
        default:
            logC+=log(2*o_lin)+log(i_lin-o_lin)-log(i_lin+o_lin);
            logsum_odd=(-(o_lin+1)*o_lin*bound_time/2.0/Ne)+log(2*o_lin+1)-log(2*o_lin)+logC; //log((2*o_lin)*(i_lin-1-o_lin)/(double)(i_lin-1-o_lin));
            pi=o_lin;
            for (i=o_lin+2; i<=i_lin;++i)
            {
                ++pi;
                logC+=log(o_lin+pi)+log(i_lin-pi)-log(i_lin+pi)-log(i-o_lin);
                if ((i-o_lin)%2==0)
                    logsum_even=LogscaleAdd(logsum_even,(-i*pi*bound_time/2/Ne)+log(2*i-1)-log(pi+o_lin)+logC);
                else
                    logsum_odd=LogscaleAdd(logsum_odd,(-i*pi*bound_time/2/Ne)+log(2*i-1)-log(pi+o_lin)+logC);
            }
            break;
    }
    if (logsum_odd>=logsum_even)
        return 0;
    else
        return exp(LogscaleSub(logsum_even,logsum_odd)-LogscaleFact(o_lin));
}

double SampleCoalTimeMLCFromXtoYLineages(int i_lin, int o_lin, double bound_time, int pop_size, double brent_epsilon, double density, int verbosity)
{
    double lambda_i=0, lambda_o=0, C0=1, out_sum=0, ret_value=0;
    int i=0;
    
    lambda_i= -i_lin*(i_lin-1)/(double)(2*pop_size);
    lambda_o= -o_lin*(o_lin-1)/(double)(2*pop_size);
    for (i=0;i<o_lin;++i)
    {
        C0*=(o_lin+i)*(i_lin-i-1)/(double)(i_lin-1+i);
    }
    
    if (o_lin>GSL_SF_FACT_NMAX)
        return bound_time+1;
    
    out_sum=(-lambda_i)/(gsl_sf_fact(o_lin)*KahanLogscaleProbCoalFromXtoYLineages(i_lin, o_lin, bound_time, pop_size)); /// \todo test if this is the proper order to avoid precision errors.
    
    if (brent_root(*(SampleCDFCoalTimeMLCFromXtoYLineages),0.0,bound_time,brent_epsilon,&ret_value,verbosity,9,i_lin,o_lin,bound_time,pop_size,lambda_i,lambda_o,C0,out_sum,density)!=NO_ERROR)
        return bound_time+1;
    else
        return ret_value;
}

#ifdef __MPFR_H
double MPFRSampleCoalTimeMLCFromXtoYLineages(int i_lin, int o_lin, double bound_time, int pop_size, double brent_epsilon, double density, int verbosity, int precision)
{
    mpfr_t lambda_i,lambda_o,C0,Cx,out_sum,fact;
    int i=0;
    double result=0;
    
    mpfr_inits2(precision,lambda_i,lambda_o,C0,Cx,out_sum,fact,NULL);
    mpfr_set_si(lambda_i, -i_lin*(i_lin-1), MPFR_RNDN);
    mpfr_set_si(lambda_o, -o_lin*(o_lin-1),MPFR_RNDN);
    mpfr_div_si(lambda_i, lambda_i, 2*pop_size, MPFR_RNDN);
    mpfr_div_si(lambda_o, lambda_o, 2*pop_size, MPFR_RNDN);
    mpfr_set_si(C0, 1, MPFR_RNDN);
    
    for (i=0;i<o_lin;++i)
    {
        mpfr_set_si(Cx, (o_lin+i)*(i_lin-i-1), MPFR_RNDN);
        mpfr_div_si(Cx, Cx, i_lin-1+i, MPFR_RNDN);
        mpfr_mul(C0, C0, Cx, MPFR_RNDN);
    }
    
    mpfr_set(out_sum,lambda_i,MPFR_RNDN);
    mpfr_div_d(out_sum, out_sum, MPFROriginalProbcoalFromXtoYLineages(i_lin, o_lin, bound_time,pop_size,precision), MPFR_RNDN);
    mpfr_fac_ui(fact, o_lin, MPFR_RNDN);
    mpfr_div(out_sum,out_sum,fact,MPFR_RNDN);
    mpfr_mul_si(out_sum, out_sum, -1, MPFR_RNDN);
    ErrorReporter(brent_root(*(MPFRSampleCDFCoalTimeMLCFromXtoYLineages),0.0,bound_time,brent_epsilon,&result,verbosity,10,i_lin,o_lin,bound_time,pop_size,(mpfr_ptr)&lambda_i,(mpfr_ptr)&lambda_o,(mpfr_ptr)&C0,(mpfr_ptr)&out_sum,density,precision),"Error rooting the MPFRCDFCoalTimeFromXtoYLineages\n");
    mpfr_clears(Cx,fact,lambda_i,lambda_o,C0,out_sum,NULL);
     return result;
}
#endif

double LogscaleSampleCoalTimeMLCFromXtoYLineages(int i_lin, int o_lin, double bound_time, int pop_size, double brent_epsilon, double density, int verbosity)
{
    double lambda_i=0, lambda_o=0, logC0=0, logoutsum=0, result=0;
    int i=0;
    
    lambda_i= -i_lin*(i_lin-1)/(double)(2*pop_size);
    lambda_o= -o_lin*(o_lin-1)/(double)(2*pop_size);
    for (i=0;i<o_lin;++i)
    {
        logC0+=log((o_lin+i)*(i_lin-i-1)/(double)(i_lin-1+i));
    }
    
    logoutsum=log(-lambda_i)-LogscaleFact(o_lin)-KahanLogscaleLogProbCoalFromXtoYLineages(i_lin, o_lin, bound_time, pop_size);
    
    if(brent_root(*(LogscaleSampleCDFCoalTimeMLCFromXtoYLineages),0.0,bound_time,brent_epsilon,&result,verbosity,9,i_lin,o_lin,bound_time,pop_size,lambda_i,lambda_o,logC0,logoutsum,density)!=NO_ERROR)
        return bound_time+1;
    else
        return result;
}


double SampleCDFCoalTimeMLCFromXtoYLineages(double w_time, int n_arg, va_list ap)
{
    double C=0,sum=0,lambda_k=0,bound_time=0,lambda_i=0,lambda_o=0,out_sum=0,density=0;//,cprob=0;
    int k=0,k1=0,i=0,i_lin=0,o_lin=0,pop_size=0;
    
    if (n_arg!=9)
        ErrorReporter(UNEXPECTED_VALUE,NULL);
    
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
                ErrorReporter(UNEXPECTED_VALUE,NULL);
                break;
                
        }
    }
    va_end(ap);
    
//    printf("Original i_lin %d, o_lin %d, bound_time %lf, pop_size %d, lambda_i %lf, lambda_o %lf, C %lf, out_sum %lf, density %lf\n", i_lin, o_lin, bound_time, pop_size, lambda_i, lambda_o,C,out_sum,density);
//    fflush(stdout);
    
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
//        cprob=sum*out_sum;
//        printf("original cpprob %lf\n",cprob);
//        fflush(stdout);
        return sum*out_sum-density;
    }
}

double LogscaleSampleCDFCoalTimeMLCFromXtoYLineages(double w_time, int n_arg, va_list ap)
{
    double logC=0,logsum_even=0,logsum_odd=0,lambda_k=0,bound_time=0,lambda_i=0,lambda_o=0,logoutsum=0,density=0,cprob=0,comp_even=0,comp_odd=0;
    int k=0,k1=0,i=0,i_lin=0,o_lin=0,pop_size=0;
    
    if (n_arg!=9)
        ErrorReporter(UNEXPECTED_VALUE,NULL);
    
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
                logC=va_arg(ap,double);
                break;
            case 7:
                logoutsum=va_arg(ap,double);
                break;
            case 8:
                density=va_arg(ap,double);
                break;
            default:
                ErrorReporter(UNEXPECTED_VALUE,NULL);
                break;
                
        }
    }
//     printf("logscale i_lin %d, o_lin %d, bound_time %lf, pop_size %d, lambda_i %lf, lambda_o %lf, C %lf, out_sum %lf, density %lf\n", i_lin, o_lin, bound_time, pop_size, lambda_i, lambda_o,exp(logC),exp(logoutsum),density);
//    fflush(stdout);
    va_end(ap);
    
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
        logsum_even= lambda_o*bound_time + LogscaleSub(0,(lambda_i-lambda_o)*w_time) - log(lambda_o-lambda_i)+logC; //First iteration of the summatory.
        if (o_lin==i_lin-1)
        {
            return exp(logsum_even+logoutsum)-density;
        }
        else
        {
            lambda_k=-(o_lin+1)*o_lin/(double)(2*pop_size);
            logC+=log(2*o_lin*(i_lin-1-o_lin)/(double)(i_lin-1+o_lin));
            logsum_odd=lambda_k*bound_time + LogscaleSub(0,(lambda_i-lambda_k)*w_time)-log(lambda_k-lambda_i)+0.693147180559945309417232121458+log((o_lin+0.5)/(2*o_lin))+logC;
            k1=o_lin;
            for (k=o_lin+2; k<i_lin; ++k) //K = number of lineages
            {
                ++k1;
                lambda_k=-k*k1/(double)(2*pop_size);
                logC+=log((o_lin+k1)*(i_lin-1-k1)/(double)(i_lin-1+k1)/(k-o_lin));
                if ((k-o_lin)%2==0)
                    logsum_even=KahanLogscaleAdd(logsum_even, lambda_k*bound_time + LogscaleSub(0,(lambda_i-lambda_k)*w_time)-log(lambda_k-lambda_i)+0.693147180559945309417232121458+log((k-0.5)/(k1+o_lin))+logC,&comp_even);
                else
                    logsum_odd=KahanLogscaleAdd(logsum_odd, lambda_k*bound_time + LogscaleSub(0,(lambda_i-lambda_k)*w_time)-log(lambda_k-lambda_i)+0.693147180559945309417232121458+log((k-0.5)/(k1+o_lin))+logC,&comp_odd);
            }
            if (logsum_odd>=logsum_even)
                return 0-density;
            else
            {
                cprob=exp(LogscaleSub(logsum_even,logsum_odd)+logoutsum);
//                printf("logscale cpprob %lf\n",cprob);
//                fflush(stdout);
                //if (cprob==cprob) //cprob!=NAN
                    return cprob-density;
                //else
                    //return 0-density;
            }
        }
    }
}

#ifdef __MPFR_H
double MPFRSampleCDFCoalTimeMLCFromXtoYLineages(double w_time, int n_arg, va_list ap)
{
    double bound_time=0,density=0,cprob=0;
    mpfr_ptr lambda_i=NULL, lambda_o=NULL, initC=NULL, out_sum=NULL;
    mpfr_t lambda_k,sum,dif_lambdas,nsum,t2,Cx,C;
    int k=0,k1=0,i=0,i_lin=0,o_lin=0,pop_size=0,precision=512;
    
    if (n_arg!=10)
        ErrorReporter(UNEXPECTED_VALUE,NULL);
    
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
                lambda_i=va_arg(ap,mpfr_ptr);
                break;
            case 5:
                lambda_o=va_arg(ap,mpfr_ptr);
                break;
            case 6:
                initC=va_arg(ap,mpfr_ptr);
                break;
            case 7:
                out_sum=va_arg(ap,mpfr_ptr);
                break;
            case 8:
                density=va_arg(ap,double);
                break;
            case 9:
                precision=va_arg(ap,int);
                break;
            default:
                ErrorReporter(UNEXPECTED_VALUE,NULL);
                break;
                
        }
    }
    va_end(ap);
    
//    printf("mpfr i_lin %d, o_lin %d, bound_time %lf, pop_size %d, lambda_i %lf, lambda_o %lf, C %lf, out_sum %lf, density %lf\n", i_lin, o_lin, bound_time, pop_size, mpfr_get_d(lambda_i,MPFR_RNDN), mpfr_get_d(lambda_o,MPFR_RNDN),mpfr_get_d(initC,MPFR_RNDN),mpfr_get_d(out_sum,MPFR_RNDN),density);
//    fflush(stdout);
    
    mpfr_inits2(precision,lambda_k,sum,dif_lambdas,t2,Cx,nsum,C, NULL);
    
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
        mpfr_set(C, initC, MPFR_RNDN);
        mpfr_mul_d(sum, lambda_o, bound_time, MPFR_RNDN);
        mpfr_exp(sum, sum, MPFR_RNDN);
        mpfr_mul(sum,sum,C,MPFR_RNDN);
        mpfr_sub(dif_lambdas,lambda_i,lambda_o,MPFR_RNDN);
        mpfr_div(sum, sum, dif_lambdas, MPFR_RNDN);
        mpfr_mul_d(t2,dif_lambdas,w_time,MPFR_RNDN);
        mpfr_exp(t2,t2,MPFR_RNDN);
        mpfr_sub_si(t2, t2, 1, MPFR_RNDN);
        mpfr_mul(sum,sum,t2,MPFR_RNDN);
        
        for (k=o_lin+1; k<i_lin; ++k) //K = number of lineages
        {
            ++k1;
            mpfr_set_si(lambda_k, -k*k1, MPFR_RNDN);
            mpfr_div_si(lambda_k, lambda_k, 2*pop_size, MPFR_RNDN);
            mpfr_set_si(Cx, (o_lin+k1)*(i_lin-1-k1), MPFR_RNDN);
            mpfr_div_si(Cx, Cx, i_lin-1+k1, MPFR_RNDN);
            mpfr_div_si(Cx, Cx, o_lin-k, MPFR_RNDN);
            mpfr_mul(C, C, Cx, MPFR_RNDN);
            mpfr_mul_d(nsum,lambda_k,bound_time,MPFR_RNDN);
            mpfr_exp(nsum, nsum, MPFR_RNDN);
            mpfr_sub(dif_lambdas,lambda_i,lambda_k,MPFR_RNDN);
            mpfr_mul_d(t2, dif_lambdas, w_time, MPFR_RNDN);
            mpfr_exp(t2,t2,MPFR_RNDN);
            mpfr_sub_si(t2, t2, 1, MPFR_RNDN);
            mpfr_mul(nsum,nsum,t2,MPFR_RNDN);
            mpfr_div(nsum,nsum,dif_lambdas,MPFR_RNDN);
            mpfr_mul_si(nsum, nsum, 2*k-1, MPFR_RNDN);
            mpfr_div_si(nsum, nsum, k1+o_lin, MPFR_RNDN);
            mpfr_mul(nsum,nsum,C,MPFR_RNDN);
            mpfr_add(sum, sum, nsum, MPFR_RNDN);
        }
        mpfr_mul(sum, sum, out_sum, MPFR_RNDN);
        cprob=mpfr_get_d(sum,MPFR_RNDN);
//      printf("mpfr cpprob %lf\n",cprob);
//      fflush(stdout);
        mpfr_clears(lambda_k,sum,dif_lambdas,t2,Cx,nsum,C, NULL);
        
        if (cprob>1 || cprob<0)
        {
            ErrorReporter(UNEXPECTED_VALUE, "MPFR out of precision!!!! That options must be completely out of sense\n");
        }
        return cprob-density;
    }
}
#endif

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
		ErrorReporter(UNEXPECTED_VALUE,NULL);
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
		ErrorReporter(UNEXPECTED_VALUE,NULL);
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
 * \name Logscale
 *******************************************************************************/
///@{
double LogscaleAdd(double logx, double logy)
{
    if (logx > logy)
        return logx + biomcmc_log1p(exp(logy - logx));
    else
        return logy + biomcmc_log1p(exp(logx - logy));
}

double KahanLogscaleAdd(double logx, double logy, double * comp_sum)
{
    double corr=0,nsum=0,sum=0,input=0;

    if (logx > logy)
    {
        sum=logx;
        input=biomcmc_log1p(exp(logy - logx));
        corr=input-*comp_sum;
        nsum=sum+corr;
        *comp_sum=(nsum-sum)-corr;
        return nsum;
    }
    else
    {
        sum=logy;
        input=biomcmc_log1p(exp(logx - logy));
        corr=input-*comp_sum;
        nsum=sum+corr;
        *comp_sum=(nsum-sum)-corr;
        return nsum;
    }
    
}

double LogscaleSub(double logx, double logy)
{
    return logx + biomcmc_log1p(-exp(logy - logx));
}

double LogscaleFact(int n)
{
    double x=0;
    static const double lf[255] =
    {
        0.000000000000000,
        0.000000000000000,
        0.693147180559945,
        1.791759469228055,
        3.178053830347946,
        4.787491742782046,
        6.579251212010101,
        8.525161361065415,
        10.604602902745251,
        12.801827480081469,
        15.104412573075516,
        17.502307845873887,
        19.987214495661885,
        22.552163853123421,
        25.191221182738683,
        27.899271383840894,
        30.671860106080675,
        33.505073450136891,
        36.395445208033053,
        39.339884187199495,
        42.335616460753485,
        45.380138898476908,
        48.471181351835227,
        51.606675567764377,
        54.784729398112319,
        58.003605222980518,
        61.261701761002001,
        64.557538627006323,
        67.889743137181526,
        71.257038967168000,
        74.658236348830158,
        78.092223553315307,
        81.557959456115029,
        85.054467017581516,
        88.580827542197682,
        92.136175603687079,
        95.719694542143202,
        99.330612454787428,
        102.968198614513810,
        106.631760260643450,
        110.320639714757390,
        114.034211781461690,
        117.771881399745060,
        121.533081515438640,
        125.317271149356880,
        129.123933639127240,
        132.952575035616290,
        136.802722637326350,
        140.673923648234250,
        144.565743946344900,
        148.477766951773020,
        152.409592584497350,
        156.360836303078800,
        160.331128216630930,
        164.320112263195170,
        168.327445448427650,
        172.352797139162820,
        176.395848406997370,
        180.456291417543780,
        184.533828861449510,
        188.628173423671600,
        192.739047287844900,
        196.866181672889980,
        201.009316399281570,
        205.168199482641200,
        209.342586752536820,
        213.532241494563270,
        217.736934113954250,
        221.956441819130360,
        226.190548323727570,
        230.439043565776930,
        234.701723442818260,
        238.978389561834350,
        243.268849002982730,
        247.572914096186910,
        251.890402209723190,
        256.221135550009480,
        260.564940971863220,
        264.921649798552780,
        269.291097651019810,
        273.673124285693690,
        278.067573440366120,
        282.474292687630400,
        286.893133295426990,
        291.323950094270290,
        295.766601350760600,
        300.220948647014100,
        304.686856765668720,
        309.164193580146900,
        313.652829949878990,
        318.152639620209300,
        322.663499126726210,
        327.185287703775200,
        331.717887196928470,
        336.261181979198450,
        340.815058870798960,
        345.379407062266860,
        349.954118040770250,
        354.539085519440790,
        359.134205369575340,
        363.739375555563470,
        368.354496072404690,
        372.979468885689020,
        377.614197873918670,
        382.258588773060010,
        386.912549123217560,
        391.575988217329610,
        396.248817051791490,
        400.930948278915760,
        405.622296161144900,
        410.322776526937280,
        415.032306728249580,
        419.750805599544780,
        424.478193418257090,
        429.214391866651570,
        433.959323995014870,
        438.712914186121170,
        443.475088120918940,
        448.245772745384610,
        453.024896238496130,
        457.812387981278110,
        462.608178526874890,
        467.412199571608080,
        472.224383926980520,
        477.044665492585580,
        481.872979229887900,
        486.709261136839360,
        491.553448223298010,
        496.405478487217580,
        501.265290891579240,
        506.132825342034830,
        511.008022665236070,
        515.890824587822520,
        520.781173716044240,
        525.679013515995050,
        530.584288294433580,
        535.496943180169520,
        540.416924105997740,
        545.344177791154950,
        550.278651724285620,
        555.220294146894960,
        560.169054037273100,
        565.124881094874350,
        570.087725725134190,
        575.057539024710200,
        580.034272767130800,
        585.017879388839220,
        590.008311975617860,
        595.005524249382010,
        600.009470555327430,
        605.020105849423770,
        610.037385686238740,
        615.061266207084940,
        620.091704128477430,
        625.128656730891070,
        630.172081847810200,
        635.221937855059760,
        640.278183660408100,
        645.340778693435030,
        650.409682895655240,
        655.484856710889060,
        660.566261075873510,
        665.653857411105950,
        670.747607611912710,
        675.847474039736880,
        680.953419513637530,
        686.065407301994010,
        691.183401114410800,
        696.307365093814040,
        701.437263808737160,
        706.573062245787470,
        711.714725802289990,
        716.862220279103440,
        722.015511873601330,
        727.174567172815840,
        732.339353146739310,
        737.509837141777440,
        742.685986874351220,
        747.867770424643370,
        753.055156230484160,
        758.248113081374300,
        763.446610112640200,
        768.650616799717000,
        773.860102952558460,
        779.075038710167410,
        784.295394535245690,
        789.521141208958970,
        794.752249825813460,
        799.988691788643450,
        805.230438803703120,
        810.477462875863580,
        815.729736303910160,
        820.987231675937890,
        826.249921864842800,
        831.517780023906310,
        836.790779582469900,
        842.068894241700490,
        847.352097970438420,
        852.640365001133090,
        857.933669825857460,
        863.231987192405430,
        868.535292100464630,
        873.843559797865740,
        879.156765776907600,
        884.474885770751830,
        889.797895749890240,
        895.125771918679900,
        900.458490711945270,
        905.796028791646340,
        911.138363043611210,
        916.485470574328820,
        921.837328707804890,
        927.193914982476710,
        932.555207148186240,
        937.921183163208070,
        943.291821191335660,
        948.667099599019820,
        954.046996952560450,
        959.431492015349480,
        964.820563745165940,
        970.214191291518320,
        975.612353993036210,
        981.015031374908400,
        986.422203146368590,
        991.833849198223450,
        997.249949600427840,
        1002.670484599700300,
        1008.095434617181700,
        1013.524780246136200,
        1018.958502249690200,
        1024.396581558613400,
        1029.838999269135500,
        1035.285736640801600,
        1040.736775094367400,
        1046.192096209724900,
        1051.651681723869200,
        1057.115513528895000,
        1062.583573670030100,
        1068.055844343701400,
        1073.532307895632800,
        1079.012946818975000,
        1084.497743752465600,
        1089.986681478622400,
        1095.479742921962700,
        1100.976911147256000,
        1106.478169357800900,
        1111.983500893733000,
        1117.492889230361000,
        1123.006317976526100,
        1128.523770872990800,
        1134.045231790853000,
        1139.570684729984800,
        1145.100113817496100,
        1150.633503306223700,
        1156.170837573242400,
    };
    
    if (n < 0)
    {
        ErrorReporter(UNEXPECTED_VALUE, "Factorial of a negative integer");
        return EXIT_FAILURE;
    }
    else if (n > 254)
    {
        x = n + 1;
        return (x - 0.5)*log(x) - x + 0.5*log(2*M_PI) + 1.0/(12.0*x);
    }
    else
    {
        return lf[n];
    }
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

size_t SampleNDoubles(size_t n, double * array, gsl_rng *seed)
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

double VKahanSum(double * array, int n_elements)
{
    int i=0;
    double comp=0, corr=0, sum=0,nsum=0;
    
    for (i=0; i<n_elements; ++i)
    {
        corr=*(array+i)-comp;
        nsum=sum+corr;
        comp=(nsum-sum)-corr;
        sum=nsum;
    }
    
    return sum;
}

void SKahanSum(double *sum, double input, double *compensation)
{
    double corr=0,nsum=0;
    
    corr=input-*compensation;
    nsum=*sum+corr;
    *compensation=(nsum-*sum)-corr;
    *sum=nsum;
}

double chebyshev_eval (double x, const double *a, const int n)
{/* evaluates the n-term Chebyshev series "a" at "x" */
    double b0, b1, b2, twox;
    int i;
    
    if ((n < 1) || (n > 1000) || (x < -1.1) || (x > 1.1)) return NAN;
    
    twox = x * 2;
    b2 = b1 = 0;
    b0 = 0;
    for (i = 1; i <= n; i++) {
        b2 = b1;
        b1 = b0;
        b0 = twox * b1 - b2 + a[n - i];
    }
    return (b0 - b2) * 0.5;
}

double biomcmc_log1p(double x)
{/* Compute the relative error logarithm [ log(1 + x) ] */
    /* series for biomcmc_log1p on the interval [-.375,.375] with weighted error =  6.35e-32, log weighted error = 31.20,
     * significant figures required = 30.93 and decimal places required = 32.01 */
    static const double alnrcs[43] = {
        +.10378693562743769800686267719098e+1, -.13364301504908918098766041553133e+0,
        +.19408249135520563357926199374750e-1, -.30107551127535777690376537776592e-2,
        +.48694614797154850090456366509137e-3, -.81054881893175356066809943008622e-4,
        +.13778847799559524782938251496059e-4, -.23802210894358970251369992914935e-5,
        +.41640416213865183476391859901989e-6, -.73595828378075994984266837031998e-7,
        +.13117611876241674949152294345011e-7, -.23546709317742425136696092330175e-8,
        +.42522773276034997775638052962567e-9, -.77190894134840796826108107493300e-10,
        +.14075746481359069909215356472191e-10, -.25769072058024680627537078627584e-11,
        +.47342406666294421849154395005938e-12, -.87249012674742641745301263292675e-13,
        +.16124614902740551465739833119115e-13, -.29875652015665773006710792416815e-14,
        +.55480701209082887983041321697279e-15, -.10324619158271569595141333961932e-15,
        +.19250239203049851177878503244868e-16, -.35955073465265150011189707844266e-17,
        +.67264542537876857892194574226773e-18, -.12602624168735219252082425637546e-18,
        +.23644884408606210044916158955519e-19, -.44419377050807936898878389179733e-20,
        +.83546594464034259016241293994666e-21, -.15731559416479562574899253521066e-21,
        +.29653128740247422686154369706666e-22, -.55949583481815947292156013226666e-23,
        +.10566354268835681048187284138666e-23, -.19972483680670204548314999466666e-24,
        +.37782977818839361421049855999999e-25, -.71531586889081740345038165333333e-26,
        +.13552488463674213646502024533333e-26, -.25694673048487567430079829333333e-27,
        +.48747756066216949076459519999999e-28, -.92542112530849715321132373333333e-29,
        +.17578597841760239233269760000000e-29, -.33410026677731010351377066666666e-30,
        +.63533936180236187354180266666666e-31, };
    
    if (x == 0.)  return 0.;
    if (x == -1.) return INFINITY;
    if (x < -1)   return NAN;
    
    if (fabs(x) <= .375) {
        /* Improve on speed (only); again give result accurate to IEEE double precision: */
        if(fabs(x) < .5 * DBL_EPSILON) return x;
        if( ((0 < x) && (x < 1e-8)) || (-1e-9 < x && x < 0)) return x * (1 - .5 * x);
        /* "22" for IEEE double precision where DBL_EPSILON =  2.22044604925031e-16 */
        return x * (1 - x * chebyshev_eval(x / .375, alnrcs, 22));
    }
#ifdef DBG_EXTREME
    if (x < -0.999999985) fprintf (stderr, "biomcmc DEBUG: Answer less than half precision since x=%lf too near -1 in biomcmc_log1p()\n", x);
#endif
    return gsl_log1p(x);
}


///@}


