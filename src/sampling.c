/**
 *
 * \file sampling.c
 * Source of the library of bayesian-like sampling to drive simulations
 *******************************************************************************/


#include "sampling.h"

// **** Prototipes of private functions **** //
static double GetCastedDoubleUSU(u_idsu var, int kind);
static int GetCastedIntUSU(u_idsu var, int kind);


// **** Declarations of public functions **** //

/**
 * \name Public functions
 * Better organized in \ref sampling.h
 *******************************************************************************/
///@{


// *** Sampling variables managing *** //
inline void set_sampling_uint(sampling_unit *variable, int value)
{
    variable->distribution_code=FIXED;
    variable->dependent_index=0;
    variable->value.i=value;
    variable->params[0].i=value;
    variable->params_type[0]=UI;
    variable->vtype=UI;
}

inline void set_sampling_double(sampling_unit *variable, double value)
{
    variable->distribution_code=FIXED;
    variable->dependent_index=0;
    variable->value.d=value;
    variable->params[0].d=value;
    variable->params_type[0]=D;
    variable->vtype=D;
}

inline void set_sampling_pointeruint(sampling_unit *variable, sampling_unit *dependent)
{
    variable->distribution_code=FIXED;
    variable->dependent_index=-1;
    variable->value.i=0;
    variable->params[0].p=dependent;
    variable->params_type[0]=SU;
    variable->vtype=UI;
}

inline void set_sampling_pointerdouble(sampling_unit *variable, sampling_unit *dependent)
{
    variable->distribution_code=FIXED;
    variable->dependent_index=-1;
    variable->value.d=0;
    variable->params[0].p=dependent;
    variable->params_type[0]=SU;
    variable->vtype=D;
}


int MeasureSampleUnitRecurency(sampling_unit * su, int * it)
{
    int i=0,mdepth=0,depth;
    ++*it;
    if (*it>=MAX_IT)
        ErrorReporter(SETTINGS_ERROR, "Recurrent sampling distribution parameters. Please, check your simulation parameters!\n");
    if (su->dependent_index==-1)
    {
        for (i=0;i<5;++i)
        {
            if (su->params_type[i]==SU)
            {
                depth=MeasureSampleUnitRecurency((su->params[i]).p, it);
                if (depth>mdepth)
                    mdepth=depth;
            }
        }
        su->dependent_index=mdepth+1;
    }
    return su->dependent_index;
}

int Compare_SampleUnits(const void *su1,const void *su2)
{
    int rec_test=0;
    sampling_unit *wsu1=NULL,*wsu2=NULL;
    
    wsu1=*(sampling_unit**)su1;
    wsu2=*(sampling_unit**)su2;
    if (wsu1->dependent_index==-1) //measure dependency index
    {
        wsu1->dependent_index=MeasureSampleUnitRecurency(wsu1,&rec_test);
    }
    if (wsu2->dependent_index==-1) //measure dependency index
    {
        wsu2->dependent_index=MeasureSampleUnitRecurency(wsu2,&rec_test);
    }
    return wsu1->dependent_index-wsu2->dependent_index;
}

long int sample_distr(gsl_rng *r,int n_arg,...)
{
    va_list ap;
    int i;
    sampling_unit *variable,**variables;
    variables=calloc(n_arg, sizeof(sampling_unit*));
    va_start(ap, n_arg);
    
    for (i=0;i<n_arg;++i)
    {
        *(variables+i)=va_arg(ap,sampling_unit*);
    }
    va_end(ap);
    
    qsort(variables, n_arg, sizeof(sampling_unit*), Compare_SampleUnits);
    
    for (i=0; i<n_arg; ++i)
    {
        variable=variables[i];
        switch (variable->distribution_code)
        {
            case FIXED:
                switch (variable->vtype)
                {
                    case D:
                        variable->value.d=GetCastedDoubleUSU(variable->params[0], variable->params_type[0]);
                        break;
                    default:
                        variable->value.i=GetCastedIntUSU(variable->params[0], variable->params_type[0]);
                        break;
                }
                break;
            case UNIFORM:
                set_propsampling(variable,gsl_ran_flat(r,GetCastedDoubleUSU(variable->params[0], variable->params_type[0]),variable->vtype==UI?GetCastedDoubleUSU(variable->params[1], variable->params_type[1])+1:GetCastedDoubleUSU(variable->params[1], variable->params_type[1])));
                break;
            case NORMAL:
                set_propsampling(variable,(gsl_ran_gaussian_ziggurat(r, GetCastedDoubleUSU(variable->params[1], variable->params_type[1]))+GetCastedDoubleUSU(variable->params[0], variable->params_type[0])));
                break;
            case EXPONENTIAL:
                set_propsampling(variable, gsl_ran_exponential(r, 1/GetCastedDoubleUSU(variable->params[0], variable->params_type[0])));
                break;
            case GAMMA:
                set_propsampling(variable,gsl_ran_gamma(r, GetCastedDoubleUSU(variable->params[0], variable->params_type[0]), GetCastedDoubleUSU(variable->params[1], variable->params_type[1])));
                break;
            case LOGNORMAL:
                set_propsampling(variable, gsl_ran_lognormal(r,GetCastedDoubleUSU(variable->params[0], variable->params_type[0]),GetCastedDoubleUSU(variable->params[1], variable->params_type[1])));
                break;
            case LOGNORMAL_MULT:
                set_propsampling(variable,(gsl_ran_lognormal(r,GetCastedDoubleUSU(variable->params[0], variable->params_type[0]),GetCastedDoubleUSU(variable->params[1], variable->params_type[1]))*GetCastedDoubleUSU(variable->params[2], variable->params_type[2])));
                break;
                
            default:
                return UNEXPECTED_VALUE;
                break;
        }
    }
    
    va_end(ap);
    free(variables);
    
    return NO_ERROR;
}

long int ParseSampling(char * p, sampling_unit * sample, const sampling_table sampling_vars)
{
    char * buffer=calloc(strlen(p)+1, sizeof(char)),*last_parsed=NULL;
    int in_parsing=0,is_pointer=0,b_index=0,p_index=0,i=0,j=0,p_found=0,n_p=0;
    long int temp_int=0;
    
    switch (toupper(*p))
    {
        case 'F':
            sample->distribution_code=FIXED;
            n_p=1;
            break;
        case 'U':
            sample->distribution_code=UNIFORM;
            n_p=2;
            break;
        case 'N':
            sample->distribution_code=NORMAL;
            n_p=2;
            break;
        case 'E':
            sample->distribution_code=EXPONENTIAL;
            n_p=1;
            break;
        case 'G':
            sample->distribution_code=GAMMA;
            n_p=2;
            break;
        case 'L':
            sample->distribution_code=LOGNORMAL;
            n_p=2;
            break;
        case 'S':
            switch (toupper(*(p+1)))
            {
                case 'L':
                    sample->distribution_code=LOGNORMAL_MULT;
                    n_p=3;
                    ++p;
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
    
    switch (*p)
    {
        case ':':
            break;
        default:
            return SETTINGS_ERROR;
            break;
    }

    do
    {
        ++i;
        ++p;
        if (p_index>4)
        {
            return SETTINGS_ERROR;
        }
        
        switch (in_parsing)
        {
            case 0:
                in_parsing=1;
                b_index=1;
                buffer[0]=*p;
                if (!(isdigit(*p)|| *p=='.'|| *p=='-'|| *p=='+')) //Pointer
                    is_pointer=1;
                break;
                
            default:
                switch (*p)
                {
                    case ',': //stop buffering and parse
                    case '\0':
                    case '/':
                    case '\n':
                        in_parsing=0;
                        buffer[b_index]='\0';
                        switch (is_pointer)
                        {
                            default: //value
                                temp_int=strtol(buffer, &last_parsed,10);
                                if ((temp_int>=INT_MIN && temp_int<=INT_MAX) && (*last_parsed=='\0' || *last_parsed==',' || *last_parsed=='/'|| *last_parsed=='\n'))
                                {
                                    sample->params[p_index].i=(int)temp_int;
                                    sample->params_type[p_index]=UI;
                                }
                                else
                                {
                                    sample->params[p_index].d=strtod(buffer,&last_parsed);
                                    if(!(*last_parsed=='\0' || *last_parsed==',' || *last_parsed=='/'|| *last_parsed=='\n'))
                                        return SETTINGS_ERROR;
                                    else
                                        sample->params_type[p_index]=D;
                                }
                                ++p_index;
                                break;
                            case 1: //pointer
                                if(b_index>2)
                                    return SETTINGS_ERROR;
                                p_found=0;
                                for (j=0; j<sampling_vars.n_duples; ++j)
                                {
                                    if (strcasecmp(buffer, (sampling_vars.table+j)->name)==0)
                                    {
                                        sample->params[p_index].p=(sampling_vars.table+j)->p;
                                        sample->params_type[p_index]=SU;
                                        sample->dependent_index=-1;
                                        p_found=1;
                                        break;
                                    }
                                }
                                if(p_found==0)
                                    return SETTINGS_ERROR;
                                ++p_index;
                                break;
                        }
                        is_pointer=0;
                        break;
                    case ' ':
                        break;
                    default:
                        buffer[b_index]=*p;
                        ++b_index;
                        break;
                }
                break;
        }
    }
    while (*p!='\0' && *p!='/' && i<= MAX_IT);
    if (i>=MAX_IT)
        return MAX_IT;
    else if (p_index != n_p || (*p=='/' && *(p+1)!='/'))
        return SETTINGS_ERROR;
    if (sample->distribution_code==FIXED)
    {
        if(sample->params_type[0]==UI)
            set_propsampling(sample, sample->params[0].i);
        else if (sample->params_type[0]==D)
            set_propsampling(sample, sample->params[0].d);
    }
    
    free(buffer);
    
    return NO_ERROR;
    
}

void Print_Sampling(sampling_unit *sample, char * buffer, const sampling_table stable)
{
    int n_param=0,i=0,j=0;
    
    switch (sample->distribution_code)
    {
        case FIXED:
            sprintf(buffer,"Fixed ");
            n_param=1;
            break;
        case UNIFORM:
            sprintf(buffer,"Uniform[");
            n_param=2;
            break;
        case NORMAL:
            sprintf(buffer,"Normal(");
            n_param=2;
            break;
        case EXPONENTIAL:
            sprintf(buffer,"Exp(");
            n_param=1;
            break;
        case GAMMA:
            sprintf(buffer,"Gamma(");
            n_param=2;
            break;
        case LOGNORMAL:
            sprintf(buffer,"LogN(");
            n_param=2;
            break;
        case LOGNORMAL_MULT:
            sprintf(buffer,"LogN(");
            n_param=3;
            break;
        default:
            break;
            
    }

    for(i=0;i<n_param;++i)
    {
        switch (*(sample->params_type+i))
        {
            case UI:
                sprintf(buffer+strlen(buffer),"%u,",(sample->params+i)->i);
                break;
            case D:
                sprintf(buffer+strlen(buffer),"%e,",(sample->params+i)->d);
                break;
            case SU:
                for (j=0;j<=stable.n_duples;++j)
                {
                    if((stable.table+j)->p==(sample->params+i)->p)
                        break;
                }
                sprintf(buffer+strlen(buffer),"%s,",(stable.table+j)->name);
                break;
            default:
                ErrorReporter(UNEXPECTED_VALUE, "");
                break;
        }
        
    }

    switch (sample->distribution_code)
    {
        default:
            sprintf(buffer+strlen(buffer)-1, ")");
            break;
        case UNIFORM:
            switch (sample->vtype)
            {
                case UI:
                    sprintf(buffer+strlen(buffer)-1, "]");
                    break;
                default:
                    sprintf(buffer+strlen(buffer)-1, ")");
                    break;
            }
        case FIXED:
            break;
    }
}

int is_variable(sampling_unit value)
{
    switch (value.distribution_code)
    {
        case FIXED:
            switch (value.params_type[0])
            {
                case SU:
                    return is_variable(*value.params[0].p);
                    break;
                case UI:
                    switch (value.params[0].i)
                    {
                        case 0:
                            return 0;
                            break;
                            
                        default:
                            return 1;
                            break;
                    }
                    break;
                case D:
                    if (value.params[0].d==0)
                        return 0;
                    else
                        return 1;
                    break;
            }
            break;
            
        default:
            return 1;
            break;
    }
    return 0;
}

int is_sampling_set(sampling_unit value)
{
    if (value.vtype==UI)
    {
        if ((value.distribution_code==FIXED)&&(value.value.i==0)&&(value.dependent_index==0))
        {
            return 0;
        }
        else
            return 1;
    }
    else
    {
        if ((value.distribution_code==FIXED)&&(value.value.d==0)&&(value.dependent_index==0))
        {
            return 0;
        }
        else
            return 1;
    }
}
///@}

/**
 * \name Private functions
 *******************************************************************************/
///@{
inline double GetCastedDoubleUSU(u_idsu var, int kind)
{
    if (kind==SU)
        return (double) get_sampling((*var.p));
    else if (kind==D)
        return var.d;
    else
        return (double) var.i;
}

inline int GetCastedIntUSU(u_idsu var, int kind)
{
    if (kind==SU)
        return (int) get_sampling((*var.p));
    else if (kind==D)
        return (int) var.d;
    else
        return var.i;
}
///@}

