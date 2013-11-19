/**
 *
 * \file sampling.c
 * Source of the library of bayesian-like sampling to drive simulations
 *******************************************************************************/


#include "sampling.h"

// **** Declarations of public functions **** //

/**
 * \name Public functions
 * Better organized in \ref sampling.h
 *******************************************************************************/
///@{


// *** Sampling variables managing *** //
inline void set_sampling_uint(sampling_unit *variable, unsigned int value)
{
    variable->distribution_code=FIXED;
    variable->value=(ui_d)value;
    variable->params[0]=value;
    variable->vtype=UI;
}

inline void set_sampling_double(sampling_unit *variable, double value)
{
    variable->distribution_code=FIXED;
    variable->value=(ui_d)value;
    variable->params[0]=value;
    variable->vtype=D;
}

long int sample_distr(gsl_rng *r,int n_arg,...)
{
    va_list ap;
    unsigned int i;
    sampling_unit * param;
    
    va_start(ap, n_arg);
    
    for (i=0; i<n_arg; ++i)
    {
        param=va_arg(ap,sampling_unit*);
        switch (param->distribution_code)
        {
            case FIXED:
                set_propsampling(param, *param->params);
                
                break;
            case UNIFORM:
                set_propsampling(param,gsl_ran_flat(r,param->params[0],param->params[1]));
                break;
            case NORMAL:
                set_propsampling(param,(gsl_ran_gaussian_ziggurat(r, param->params[1])+param->params[0]));
                break;
            case EXPONENTIAL:
                set_propsampling(param, gsl_ran_exponential(r, 1/param->params[0]));
                break;
            case GAMMA:
                set_propsampling(param,gsl_ran_gamma(r, param->params[0], param->params[1]));
                break;
            case LOGNORMAL:
                set_propsampling(param, gsl_ran_lognormal(r,param->params[0],param->params[1]));
                break;
            case LOGNORMAL_MULT:
                set_propsampling(param,(gsl_ran_lognormal(r,param->params[0],param->params[1])*param->params[2]));
                break;
                
            default:
                return UNEXPECTED_VALUE;
                break;
        }
    }
    
    va_end(ap);
    
    return NO_ERROR;
}

long int ParseSampling(char * p, sampling_unit * sample)
{
    char * buffer=calloc(strlen(p), sizeof(char));
    unsigned int in_parsing=0;
    unsigned int b_index=0;
    int p_index=-1, is_neg=0;
    
    switch (toupper(*p))
    {
        case 'F':
            sample->distribution_code=FIXED;
            break;
        case 'U':
            sample->distribution_code=UNIFORM;
            break;
        case 'N':
            sample->distribution_code=NORMAL;
            break;
        case 'E':
            sample->distribution_code=EXPONENTIAL;
            break;
        case 'G':
            sample->distribution_code=GAMMA;
            break;
        case 'L':
            sample->distribution_code=LOGNORMAL;
            break;
        case 'S':
            switch (toupper(*(p+1)))
            {
                case 'L':
                sample->distribution_code=LOGNORMAL_MULT;
                break;
                
                default:
                return SETTINGS_ERROR;
                break;
            }
            break;
        default:
            return SETTINGS_ERROR;
            break;
    }
    ++p;
    
    while (*p!='\0')
    {
        if (p_index>4)
        {
            return SETTINGS_ERROR;
        }
        
        if ((isdigit(*p)|| *p=='.'|| *p=='-')&&in_parsing==0)
        {
            ++p_index;
            in_parsing=1;
            b_index=0;
            is_neg=0;
            if (*p=='-')
                is_neg=1;

        }
        else if (!(isdigit(*p)|| *p=='.'))
        {
            in_parsing=0;
            buffer[b_index]='\0';
            if(sscanf(buffer,"%lf",sample->params+p_index)==0)
                return SETTINGS_ERROR;
            if (is_neg==1)
                *(sample->params+p_index)*=-1;
        }
        if (in_parsing==1 && *p!='-')
        {
            buffer[b_index]=*p;
            ++b_index;
        }
        
        
        ++p;
    }
    buffer[b_index]='\0';
    if(sscanf(buffer,"%lf",sample->params+p_index)==0)
        return SETTINGS_ERROR;
    
    if (sample->distribution_code==FIXED)
    {
        set_propsampling(sample, sample->params[0]);
    }
    
    free(buffer);
    
    return NO_ERROR;
    
}

void Print_Sampling(sampling_unit sample, char * buffer)
{

    switch (sample.distribution_code)
    {
        case FIXED:
            sample.vtype==UI?sprintf(buffer,"Fixed %u",(unsigned int)sample.params[0]):sprintf(buffer,"Fixed %lf",sample.params[0]);
            break;
        case UNIFORM:
            sprintf(buffer,"Uniform [%lf,%lf)",sample.params[0],sample.params[1]);
            break;
        case NORMAL:
            sprintf(buffer,"Normal (%lf,%lf)",sample.params[0],sample.params[1]);
            break;
        case EXPONENTIAL:
            sprintf(buffer,"Exp (%lf)",sample.params[0]);
            break;
        case GAMMA:
            sprintf(buffer,"Gamma (%lf,%lf)",sample.params[0],sample.params[1]);
            break;
        case LOGNORMAL:
            sprintf(buffer,"LogN(%lf,%lf)",sample.params[0],sample.params[1]);
            break;
        case LOGNORMAL_MULT:
            sprintf(buffer,"LogN(%lf,%lf)*%lf",sample.params[0],sample.params[1],sample.params[2]);
            break;
        default:
            break;
            
    }
}

unsigned int is_sampling_set(sampling_unit value)
{
    if (value.vtype==UI)
    {
        if (value.distribution_code==FIXED&&value.value.i==0)
        {
            return 0;
        }
        else
            return 1;
    }
    else
    {
        if (value.distribution_code==FIXED&&value.value.d==0)
        {
            return 0;
        }
        else
            return 1;
    }
}
///@}
