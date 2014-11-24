

/**
 *
 * \mainpage SimPhy
 * Gene tree simulator from a given species tree,
 * taking into account ILS,GDL and HGT.
 * \author Mallo D.
 * \date October-2013/...
 * \section usage User Guide
 *
 * <B>Usage: ./SimPhy -[Parameter code] value ...</B>
 *
 * Essential parameters\n
 * ------------------
 * <ul><li>-Guide tree (\"or\" list):</li>
 *       <ul><li>-Fixed species tree:\n</li>
 *           <ul><li>-S newick_tree string:(\"(...);\")\n</li></ul>
 *       <li>-Species tree simulation parameters:\n</li>
 *           <ul><li>-Sb birth rate (relative to time if -G != 1) rational:(0,...);\n</li>
 *           <li>-Sd death rate (relative to time if -G != 1) rational:[birth_rate,...);\n</li>
 *           <li>-Sl number of leaves natural:[2,...);\n</li>
 *           <li>-St tree max time rational:(0,...);\n</li>
 *           <li>-So internal branch length deviation from the half of the height of the ingroup:[0,...);\n</li></ul>
 *       <li>-Fixed locus tree:\n</li>
 *           <ul><li>-L newick_tree string:(\"(...);\")\n</li></ul></ul>
 * <li>-Number of trees:</li>
 *       <ul><li>-R general replicates natural:[1,...) (number of locus trees, with one gene tree from each, instead if there is no locus tree, which implies using it as number of gene trees);\n</li>
 *       <li>-Rl locus trees natural:[1,...);\n</li>
 *       <li>-Rg number of gene trees from each locus tree natural:[1,...);\n</li></ul>
 * <li>-Effective haploid population size:\n</li>
 *       <ul><li>-P natural:[1,...);\n</li></ul></ul>
 *
 * Optional parameters\n
 * -----------------
 * <ul><li>-Locus tree simulation parameters:\n</li>
 *       <ul><li>-Lb birth rate (relative to generations) rational:(0,...);\n</li>
 *       <li>-Ld death rate (relative to generations) rational:(0,...);\n</li>
 *       <li>-Lt transfer rate (relative to generations) rational:(0,...);\n</li>
 *       <li>-Ll minimum number of tips natural:[2,species_tree_n_leaves];</li>
 *       <li>-Ls minimum number of tips from different species natural:[2,species_tree_n_leaves];</li></ul>
 * <li>-Number of individuals per each species (common):\n
 *      <ul><li>-I number of individuals per species natural:[1,...);\n</li></ul>
 * <li>-Generation time:\n
 *       <ul><li>-G generation time rational:(0,...) (if it is not set, the species tree is considered with time in number of generations);\n</li></ul>
 * <li>-Substitution rate:\n
 *      <ul><li>-U substitution time (substitutions per time) rational:(0,...);\n</li></ul>
 * <li>-Substitution rate heterogeneity:\n
 *      <ul><li>-Hs Alpha parameter of the gamma distribution with mean=1 to sample rate heterogeneity mulipliers. Lineage specific rate heterogeneity. rational:[0,...) (0=>no heterogeneity);\n</li>
 *      <li>-Hl Alpha parameter of the gamma distribution with mean=1 to sample rate heterogeneity mulipliers. Gene family (gene tree) specific rate heterogeneity. rational:[0,...) (0=>no heterogeneity);\n</li>
 *      <li>-Hg Alpha parameter of the gamma distribution with mean=1 to sample rate heterogeneity mulipliers. Gene tree branch specific rate heterogeneity. rational:[0,...) (0=>no heterogeneity);\n</li></ul>
 * <li>-Bounded Multispecies Coalescent sampling related parameters:\n
 *  <ul><li>-E Brent root epsilon (maximum difference between iterations to consider them convergent): rational: (0,...);\n</li>
 *      <li>-Ep Minimum number of coalescent units to the bound to use a politomy instead of sampling the bounded coalescent. rational: (0,...) (0 = no politomies);\n</li></ul>
 * <li>-Verbosity:\n</li>
 *       <ul><li>-V natural:[0-6];\n</li></ul>
 * <li>-General options:\n</li>
 *       <ul><li>-S seed rational:(0,...);\n</li></ul>
 * <li>-Output related options:\n</li>
 *       <ul><li>-O prefix for all the output files string:*;\n</li>
 *       <li>-Os logical flag to activate the csv *.stats file with some summary statistics binary;\n</li>
 *       <li>-Or reconciliation output version. natural: (0, 3)\n
 *              <ul><li>0: No reconciliation.\n</li>
 *              <li>1: *.reconsl and *.reconlg files with the reconciliations between species tree and locus tree and between the locus tree and the gene tree respectively.\n</li>
 *              <li>2: *.recon files with the reconciliations between the species tree the gene trees, using losses (locus tree level) as proper gene tree leaves. NOT FOR GENERAL USE.\n</li>
 *              <li>3: Both reconciliation formats</li></ul>
 *      <li>-Od logical flag to activate the output of an SQLite database with 3 linked tables (Species_Trees, Locus_Trees and Gene_Trees) with some characteristics of each simulated tree.</li>
 *      <li>-Op logical flag to activate an output file (prefix.params) with the general parameters of the simulation</li></ul></ul>
 * \n<B>Example\n</B>
 *       ./SimPhy -Rl 10 -Rg 10 -Lb 0.001 -Ld 0.0011 -P 200 -V 2 -O test1 -S \"((A:250#2000/2,B:250):250,(C:300,D:300):200);\"\n
 *
 * <B> Description:</B>
 *
 * This software generates a locus tree using a birth-death model which grows inside
 * a given species tree. This locus tree is used in a second step to simulate the
 * gene tree using a multispecies coalescent process (and an special bounded multispecies
 * coalescent process for one of the linage resulting of a duplication).
 *
 * This simulator also can simulate the species tree using a birth-death/conditioned birth-death approach* (or special cases like the Yule model), or only simulate gene trees from a fixed locus tree.
 * The trees should be specified in a modified newick format, which allows some extra tags to add useful info:\n
 * # : Private effective haploid population size\n
 * / : Number of individuals (applied to the gene tree, tagged in the species tree or locus tree).\n
 * % : Duplication state (0:speciation node, 1:duplication node), only appliable to a fixed locus tree.\n
 * * : Lineage (branch) specific substitution rate multiplier.\n
 *\note The species tree conditioned birth-death (fixed number of leaves) requires a birth rate equal or bigger than the death rate.
 *\todo Include in the manual a detailed explanation on the epsilon values.
 *******************************************************************************/

/**
 *
 * \file main.c
 * Main file of the project.
 *
 *******************************************************************************/
#ifndef trees_h
#include "trees.h"
#endif

#ifndef sampling_h
#include "sampling.h"
#endif

#ifndef num_methods_h
#include "num_methods.h"
#endif

#ifndef sql_managing_h
#include "sql_managing.h"
#endif

#include <time.h>
#include <limits.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include <stdlib.h>
extern int errno;

int MAX_IT=10000;
int MAX_NAME=100;
int MAX_LEAVES=1000;
int NUM_BUFFER=20;
int IO_BUFFER=1000;
char TEST_CHAR=7;

//#define DBG

#undef NO_VAR ///< If NO_VAR is defined, the random number generator allways leads the same numbers.
#undef SORTHOLOGS

#ifdef DBG
#define NO_OUT
#endif

//\todo Implement branch heterogeneity to the locus tree (different speed for each paralog)//


// **** Function prototypes **** //

// NOTE: Function documentation is in the bottom of the text, in the function declaration section.

// ** I/O ** //
static inline long int GetSettings(int argc, char **argv, int *ns_trees, sampling_unit *nl_trees, int *ng_trees, char ** newick_stree, FILE **stree_ifile, int *n_istrees, char **stree_iname, sampling_unit * gen_time, char ** newick_ltree, FILE **ltree_ifile, int *n_iltrees, char **ltree_iname,sampling_unit * b_rate, sampling_unit * d_rate, sampling_unit *t_rate, sampling_unit *gc_rate,sampling_unit * lb_rate, sampling_unit * ld_rate, sampling_unit *lt_rate, sampling_unit *lgc_rate, int *t_kind, int *min_lleaves, int *min_lsleaves, sampling_unit *ind_per_sp,sampling_unit * sb_rate, sampling_unit *sd_rate, sampling_unit * bds_leaves, sampling_unit * bds_length, sampling_unit *outgroup, sampling_unit *Ne,sampling_unit *mu, sampling_unit *alpha_s, sampling_unit *alpha_l, sampling_unit *alpha_g,sampling_unit *lalpha_g,sampling_unit *salpha_g, float *epsilon,int *verbosity,char **out_name, int *stats, int *recon, int *db, int *params, int *commands, float *u_seed, char **confile_name,const sampling_table sampling_vars);
static inline long int GetSettingsFromFile(FILE *input_file,int *ns_trees, sampling_unit *nl_trees, int *ng_trees, char ** newick_stree, FILE ** stree_ifile, int *n_istrees,char **stree_iname,sampling_unit * gen_time, char ** newick_ltree, FILE **ltree_ifile, int *n_iltrees,char **ltree_iname, sampling_unit * b_rate, sampling_unit * d_rate, sampling_unit * t_rate, sampling_unit *gc_rate,sampling_unit * lb_rate, sampling_unit * ld_rate, sampling_unit * lt_rate, sampling_unit *lgc_rate, int *t_kind, int *min_lleaves, int *min_lsleaves, sampling_unit *ind_per_sp,sampling_unit * sb_rate, sampling_unit * sd_rate, sampling_unit * bds_leaves, sampling_unit * bds_length, sampling_unit * outgroup, sampling_unit *Ne,sampling_unit *mu,sampling_unit *alpha_s, sampling_unit *alpha_l, sampling_unit *salpha_g, sampling_unit *lalpha_g,sampling_unit *alpha_g,float *epsilon, int *verbosity,char **out_name, int *stats, int *recon, int *db, int *params, int *commands, float *u_seed, const sampling_table sampling_vars);
static inline void PrintXStringError(char **string, int x, char *errormsg);
void PrintUsage(void);
static inline void PrintGlobalSettings(FILE *output_file,char * s_tree_newick, char *stree_ifile, int n_istrees,int n_strees,sampling_unit *gen_time, char * l_tree_newick, char *ltree_ifile, int n_iltrees, sampling_unit *sb_rate, sampling_unit *sd_rate, sampling_unit *s_leaves, sampling_unit *s_time, sampling_unit *outgroup, sampling_unit *b_rate, sampling_unit *d_rate, sampling_unit *t_rate, sampling_unit *gc_rate,sampling_unit *lb_rate, sampling_unit *ld_rate, sampling_unit *lt_rate, sampling_unit *lgc_rate, int t_kind, int min_lleaves, int min_lsleaves, sampling_unit *ind_per_sp,sampling_unit *nl_trees,int ng_trees,sampling_unit *Ne,sampling_unit *mu,sampling_unit *alpha_s, sampling_unit *alpha_l, sampling_unit *alpha_g,sampling_unit *lalpha_g,sampling_unit *salpha_g,float epsilon, int verbosity, char * out_file,float u_seed,int stats,int recon,int db, int params, int command, const sampling_table stable);
static inline void PrintSettingsSloop(char * s_tree_newick, sampling_unit gen_time, char *l_tree_newick, char *ltree_ifile, int n_iltrees,sampling_unit sb_rate, sampling_unit sd_rate, sampling_unit s_leaves, sampling_unit s_time,sampling_unit lb_rate, sampling_unit ld_rate, sampling_unit lt_rate, sampling_unit lgc_rate,int t_kind,sampling_unit outgroup, int min_lleaves, int min_lsleaves, sampling_unit ind_per_sp,sampling_unit nl_trees,int ng_trees,sampling_unit Ne,sampling_unit mu,sampling_unit alpha_s, sampling_unit alpha_l, sampling_unit alpha_g,float epsilon, int verbosity, char * out_file, int n_replicate);
static inline void PrintSettingsLloop(int curr_ltree,sampling_unit lb_rate,sampling_unit ld_rate,sampling_unit lt_rate,sampling_unit lgc_rate,double gamma_l,sampling_unit lalpha_g);
static inline void PrintSettingsGloop(int curr_gtree,sampling_unit alpha_g);
static inline long int CheckSampledSettingsSloop(sampling_unit bds_leaves, sampling_unit bds_length, sampling_unit sb_rate, sampling_unit sd_rate, sampling_unit outgroup, sampling_unit ind_per_sp, sampling_unit nl_trees, sampling_unit b_rate, sampling_unit d_rate, sampling_unit t_rate, sampling_unit gc_rate, sampling_unit Ne, sampling_unit alpha_s, sampling_unit alpha_l, sampling_unit alpha_g, sampling_unit mu, sampling_unit gen_time, int min_lleaves);
static inline long int CheckSampledSettingsLloop(sampling_unit lb_rate, sampling_unit ld_rate, sampling_unit lt_rate, sampling_unit lgc_rate, sampling_unit lalpha_g);
static inline long int CheckSampledSettingsGloop(sampling_unit alpha_g);

/**
 * Main function
 *
 * Mainly it manages all the simulation using the libraries trees.h and num_methods.h,
 * as a few more I/O-related functions. This version has been written thinking on
 * performance, so it is a bit long with less function calls than expected.
 *******************************************************************************/

int main (int argc, char **argv)
{
    
    // ********
    /// <dl><dt> Declaration of variables </dt><dd>
    const gsl_rng_type * T;
    gsl_rng * r;
    
    /* create a generator chosen by the
     environment variable GSL_RNG_TYPE */
    
    gsl_rng_env_setup();
    
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    
    // *******
    /// <dl><dt>General variables</dt><dd>
    
    // ******
    /// Structures
    s_tree *sp_tree=NULL;
    l_tree *locus_tree=NULL;
    g_tree *gene_tree=NULL;
    name_c * names=NULL;
    l_node **node_ptrs=NULL;
    
    // ******
    /// Configuration variables
    int ns_trees=1,ng_trees=1,verbosity=1,min_lleaves=2,min_lsleaves=3, stats=0, recon=0, db=0, params=0, command=1, weirdness=0, t_kind=1;
    float epsilon_brent=0.000001, u_seed=(time (NULL) * clock());
    char *species_tree_str=NULL,*locus_tree_str=NULL,*out_name=NULL,*confile_name=NULL;
    char *buffer=NULL;
    
    buffer=getenv("SIMPHY_MAXIT");
    if (buffer!=NULL)
    {
        sscanf(buffer,"%u",&MAX_IT);
        buffer=NULL;
    }
    buffer=getenv("SIMPHY_MAXNAME");
    if (buffer!=NULL)
    {
        sscanf(buffer,"%u",&MAX_NAME);
        buffer=NULL;
    }
    buffer=getenv("SIMPHY_MAXLEAVES");
    if (buffer!=NULL)
    {
        sscanf(buffer,"%u",&MAX_LEAVES);
        buffer=NULL;
    }
    buffer=getenv("SIMPHY_NUMBUFFER");
    if (buffer!=NULL)
    {
        sscanf(buffer,"%u",&NUM_BUFFER);
        buffer=NULL;
    }
    buffer=getenv("SIMPHY_TESTCHAR");
    if (buffer!=NULL)
    {
        sscanf(buffer,"%c",&TEST_CHAR);
        buffer=NULL;
    }
    buffer=getenv("SIMPHY_IOBUFFER");
    if (buffer!=NULL)
    {
        sscanf(buffer,"%u",&IO_BUFFER);
        buffer=NULL;
    }
    
    
    // ******
    /// Sampling variables
    sampling_unit Ne, bds_leaves, ind_per_sp, nl_trees, lb_rate, ld_rate, lt_rate, lgc_rate, b_rate, d_rate, t_rate, gc_rate, sb_rate, sd_rate, bds_length, outgroup, gen_time, mu,alpha_s, alpha_l, salpha_g,lalpha_g,alpha_g;
    
    //**Sampled once
    set_sampling_double(&alpha_s,0);
    
    //**Sampled in the species loop
        //*Global
        set_sampling_uint(&Ne,0);
        set_sampling_double(&mu,1);
        set_sampling_double(&gen_time,1);
        set_sampling_uint(&nl_trees,1);
        //*ST_related
        set_sampling_double(&sb_rate,0);
        set_sampling_double(&sd_rate,0);
        set_sampling_uint(&ind_per_sp,1);
        set_sampling_double(&outgroup,0);
        set_sampling_double(&bds_length,0);
        set_sampling_uint(&bds_leaves,0);
        //*LT_related
        set_sampling_double(&b_rate,0);
        set_sampling_double(&d_rate,0);
        set_sampling_double(&t_rate,0);
        set_sampling_double(&gc_rate,0);
        set_sampling_double(&alpha_l,0);
        //*GT_related
        set_sampling_double(&salpha_g, 0);
    
    //**Sampled in the locus loop
        //*LT_related
        set_sampling_pointerdouble(&lb_rate,&b_rate);
        set_sampling_pointerdouble(&ld_rate,&d_rate);
        set_sampling_pointerdouble(&lt_rate,&t_rate);
        set_sampling_pointerdouble(&lgc_rate,&gc_rate);
        //*GT_related
        set_sampling_pointerdouble(&lalpha_g, &salpha_g);
    
    //**Sampled in the gene loop
        set_sampling_pointerdouble(&alpha_g,&lalpha_g);
    
    //HARDCODED sampling_douple table
    const sampling_table sampling_vars=
    {
    (sampling_duple[]){
        {"P",&Ne},
        {"SL",&bds_leaves},
        {"SI",&ind_per_sp},
        {"RL",&nl_trees},
        {"PB",&b_rate},
        {"PD",&d_rate},
        {"PT",&t_rate},
        {"PG",&gc_rate},
        {"PH",&lalpha_g},
        {"PP",&salpha_g},
        {"LB",&lb_rate},
        {"LD",&ld_rate},
        {"LT",&lt_rate},
        {"LG",&lgc_rate},
        {"SB",&sb_rate},
        {"SD",&sd_rate},
        {"ST",&bds_length},
        {"SO",&outgroup},
        {"G",&gen_time},
        {"U",&mu},
        {"HS",&alpha_s},
        {"HL",&alpha_l},
        {"HG",&alpha_g}
    },23};
    
    // ******
    /// Loop related variables</dd></dl>
    int curr_stree=1,curr_ltree=1,curr_gtree=1, i=0,prev_ltree_eq_stree=0;
    
    // *******
    /// <dl><dt>I/O variables</dt><dd></dd></dl>
    char g_prefix[8]="g_trees", command_sufix[8]=".comand", recon_sufix2[8]="g.recon", reconsl_sufix2[9]=".reconsl", reconlg_sufix2[10]="g.reconlg", db_sufix[4]=".db", params_sufix[8]=".params", tree_sufix[7]=".trees", weirdg_sufix[8]=".ralpha", s_outname[13]="s_tree.trees",l_outname[14]="l_trees.trees", *g_outname=NULL, stat_outname[10]="stats.txt",weirds_outname[14]="s_tree.ralpha",confile_sufix[6]=".conf",c=0;
    char  *recon_outname=NULL, *reconsl_outname=NULL, *reconlg_outname=NULL, *db_outname=NULL, *params_outname=NULL, *command_outname=NULL, *curr_outdir=NULL,  *weirdg_outname=NULL, *stree_iname=NULL, *ltree_iname=NULL, *confile_outname=NULL;
    FILE *stree_ifile=NULL, *ltree_ifile=NULL, *s_outfile=NULL,*l_outfile=NULL,*g_outfile=NULL,*stat_outfile=NULL, *params_outfile=NULL, *command_outfile=NULL, *weirds_outfile=NULL, *weirdg_outfile=NULL, *confile_infile=NULL, *confile_outfile=NULL;
    sqlite3 *database;
    int n_sdigits=0, n_ldigits=0, n_gdigits=0, n_istrees=0, n_iltrees=0, error=0;

#ifdef SORTHOLOGS
    //Sorthologs
    int sorthologs_test=1;
    char *sorthologs_outname=NULL, sorthologs_prefix[7]="g_tree", sorthologs_midfix[7]="l_tree", sorthologs_sufix[12]=".sorthologs";;
    long double t_ddup=0,t_dgc=0, t_dtrfr=0;
    unsigned long long t_ndup=0, t_ngc=0, t_ntrfr=0;
    FILE *sorthologs_outfile=NULL;
    double *dupdistances=NULL,*gcdistances=NULL,*trfrdistances=NULL;
    //
#endif
    // *******
    /// <dl><dt>Statistical variables</dt><dd></dd></dl></dd></dl>
    long double t_height_cu=0,t_height_bl=0,height_cu=0,height_bl=0,st_height=0,st_length=0,length_bl=0;
    double gamma_l=0;
    int st_lcoals=0,st_dups=0,st_losses=0,st_transf=0,st_gc=0,st_leaves=0,st_gleaves=0;
    unsigned long long t_lcoals=0,t_n_ltree=0;
    
    // ********
    /// <dl><dt>Main program</dt><dd>
    
    printf("\nSimPhy1.0\n--------------------\n"); /// Welcome
    
    // *******
    /// <dl><dt>Setting-dependent I/O</dt><dd>
    
    
    
    if (NUM_BUFFER<20 || MAX_LEAVES<100 || MAX_NAME<10 || MAX_IT<10 || IO_BUFFER<100)
    {
        fprintf(stderr,"Environment variables have been modified to unproper values, please, do not change them unless it is completely necessary\n");
        return (EXIT_FAILURE);
    }
    
    // ******
    /// Settings obtaining from command line

    printf("\nGetting settings from command line...");
#ifdef DBG
    fflush(stdout);
#endif
    
    ErrorReporter(GetSettings(argc, argv, &ns_trees, &nl_trees, &ng_trees, &species_tree_str, &stree_ifile, &n_istrees, &stree_iname, &gen_time, &locus_tree_str, &ltree_ifile, &n_iltrees, &ltree_iname, &b_rate, &d_rate, &t_rate, &gc_rate, &lb_rate, &ld_rate, &lt_rate, &lgc_rate, &t_kind, &min_lleaves, &min_lsleaves, &ind_per_sp, &sb_rate, &sd_rate, &bds_leaves, &bds_length, &outgroup, &Ne, &mu, &alpha_s, &alpha_l, &alpha_g, &lalpha_g, &salpha_g, &epsilon_brent, &verbosity, &out_name, &stats, &recon, &db, &params, &command, &u_seed, &confile_name, sampling_vars), ": getting settings");

#ifdef NO_VAR
    gsl_rng_set(r,5);
#endif
#ifndef NO_VAR
    gsl_rng_set(r,u_seed);
#endif
    
    if(out_name==NULL)
    {
        out_name=malloc(sizeof("SimPhy_outfiles"));
        strcpy(out_name,"SimPhy_outfiles"); //Default value
    }
    
    printf("Done \n");
#ifdef DBG
    fflush(stdout);
#endif
    
    
    // ******
    /// Printing of recognized settings
    if (verbosity)
    {
        PrintGlobalSettings(stdout,species_tree_str,stree_iname,n_istrees,ns_trees,&gen_time,locus_tree_str,ltree_iname,n_iltrees,&sb_rate,&sd_rate,&bds_leaves,&bds_length,&outgroup,&b_rate,&d_rate,&t_rate,&gc_rate,&lb_rate, &ld_rate, &lt_rate, &lgc_rate,t_kind,min_lleaves,min_lsleaves,&ind_per_sp,&nl_trees,ng_trees,&Ne,&mu,&alpha_s,&alpha_l,&alpha_g,&lalpha_g,&salpha_g,epsilon_brent,verbosity,out_name, u_seed,stats,recon,db,params,command,sampling_vars);
    }
    
    
    // ******
    /// Initialization of file-I/O (output dir, and opening of some files)</dd></dl>
#ifndef NO_OUT
    n_sdigits=(int)count_intdigits((long)ns_trees,0);
    curr_outdir=malloc((n_sdigits+1)*sizeof(char));
    
    if (verbosity>4)
    {
        printf("\n\nI/O directories managing\n-------------------------\nGenerating the output directory %s... ",out_name);
#ifdef DBG
        fflush(stdout);
#endif
    }
    
    umask(0);
    
    if(mkdir(out_name, S_IRWXU | S_IRWXG | S_IRWXO)!=0)
    {
        if (errno==EEXIST)
        {
            if (verbosity>4)
            {
                printf("\n\t WARNING, The output directory already exists\n");
#ifdef DBG
            fflush(stdout);
#endif
            }
        }
        else
        {
            fprintf(stderr,"\n\t ERROR, The output folder, is not accesible. Errno %d\n",errno);
#ifdef DBG
            fflush(stderr);
#endif
            return (IO_ERROR);
        }
    }
    
    if (confile_name != NULL && (confile_infile=fopen(confile_name,"r"))==NULL)
    {
        perror("Error opening input .conf file:");
        ErrorReporter(IO_ERROR,NULL);
    }
    
    error=chdir(out_name);
    
    if(error!=0)
    {
        fprintf(stderr,"\n\t ERROR, The output folder, is not accesible\n");
#ifdef DBG
        fflush(stderr);
#endif
        return IO_ERROR;
    }
    if (verbosity>4)
    {
        printf("Done\n");
    }
    
    if (db>0)
    {
        db_outname=malloc((strlen(out_name)+strlen(db_sufix)+1)*sizeof(char));
        strcpy(db_outname,out_name);
        strcat(db_outname,db_sufix);
        if (verbosity>4)
        {
            printf("Initiating the SQLite database %s... ",db_outname);
#ifdef DBG
            fflush(stdout);
#endif
        }
        ErrorReporter(InitDB(&database, db_outname), NULL);
        if (verbosity>4)
        {
            printf("Done\n");
#ifdef DBG
            fflush(stdout);
#endif
        }
        
    }
    if (params>0)
    {
        params_outname=malloc((strlen(out_name)+strlen(params_sufix)+1)*sizeof(char));
        strcpy(params_outname,out_name);
        strcat(params_outname,params_sufix);
        
        if (verbosity>4)
        {
            printf("Saving the global settings in %s... ",params_outname);
#ifdef DBG
            fflush(stdout);
#endif
        }
        if ((params_outfile=fopen(params_outname,"w"))==NULL)
        {
            perror("Error opening .params file:");
            ErrorReporter(IO_ERROR,NULL);
        }
        PrintGlobalSettings(params_outfile,species_tree_str,stree_iname,n_istrees,ns_trees,&gen_time,locus_tree_str,ltree_iname,n_iltrees,&sb_rate,&sd_rate,&bds_leaves,&bds_length,&outgroup,&b_rate,&d_rate,&t_rate,&gc_rate,&lb_rate, &ld_rate, &lt_rate, &lgc_rate,t_kind,min_lleaves,min_lsleaves,&ind_per_sp,&nl_trees,ng_trees,&Ne,&mu,&alpha_s,&alpha_l,&alpha_g,&lalpha_g,&salpha_g,epsilon_brent,verbosity,out_name, u_seed,stats,recon,db,params,command,sampling_vars);
        fclose(params_outfile);
        free(params_outname);
        
        if (verbosity>4)
        {
            printf("Done\n");
#ifdef DBG
            fflush(stdout);
#endif
        }
        
    }
    if (command>0)
    {
        command_outname=malloc((strlen(out_name)+strlen(command_sufix)+1)*sizeof(char));
        strcpy(command_outname,out_name);
        strcat(command_outname,command_sufix);
        
        if (confile_name!=NULL)
        {
            confile_outname=malloc((strlen(out_name)+strlen(confile_sufix)+1)*sizeof(char));
            strcpy(confile_outname,out_name);
            strcat(confile_outname,confile_sufix);
            if (verbosity>4)
            {
                printf("Saving the command line settings in %s and the configuration file in %s... ",command_outname,confile_outname);
#ifdef DBG
                fflush(stdout);
#endif
            }
        }
        else if (verbosity>4)
        {
            printf("Saving the original command line settings in %s... ",command_outname);
#ifdef DBG
            fflush(stdout);
#endif
        }
        if ((command_outfile=fopen(command_outname,"w"))==NULL)
        {
            perror("Error opening .command file:");
            ErrorReporter(IO_ERROR, NULL);
        }
        for (i=0; i<argc;++i)
            fprintf(command_outfile,"%s ",argv[i]);
        fclose(command_outfile);
        free(command_outname);

        if(confile_name!=NULL)
        {
            if ((confile_outfile=fopen(confile_outname,"w"))==NULL)
            {
                perror("Error opening output .conf file:");
                ErrorReporter(IO_ERROR,NULL);
            }
            while ((c=fgetc(confile_infile)) && c!=EOF)
            {
                fputc(c,confile_outfile);
                if (ferror(confile_outfile))
                {
                    ErrorReporter(IO_ERROR,NULL);
                }
            }
            fclose(confile_infile);
            fclose(confile_outfile);
            free(confile_outname);
        }
        
        if (verbosity>4)
        {
            printf("Done\n");
#ifdef DBG
            fflush(stdout);
#endif
        }
    }
    
    if (n_istrees>0)
    {
        ErrorReporter(InitNexusParser(stree_ifile), ": parsing Nexus file");
    }
    
#endif
    
    // *******
    /// <dl><dt>Replication loop</dt><dd>

    if (verbosity>0 && verbosity<3)
        printf("\n\nSimulation:\n-----------\n");
    
#ifdef DBG
    fflush(stdout);
#endif
    
    for (curr_stree=1;curr_stree<=ns_trees;++curr_stree)
    {
        
        // ******
        /// <dl><dt>Species-level parameter sampling</dt></dl>

        ErrorReporter(sample_distr(r,17,&bds_leaves,&bds_length,&sb_rate,&sd_rate,&outgroup,&ind_per_sp,&nl_trees,&b_rate,&d_rate,&t_rate,&gc_rate,&Ne,&alpha_s,&alpha_l,&salpha_g,&mu,&gen_time),": sampling distributions");
        ErrorReporter(CheckSampledSettingsSloop(bds_leaves,bds_length,sb_rate,sd_rate,outgroup,ind_per_sp,nl_trees,b_rate,d_rate,t_rate,gc_rate,Ne,alpha_s,alpha_l,salpha_g,mu,gen_time,min_lleaves), ": improper sampled values");
        
        if (n_istrees>0)
        {
            ErrorReporter(NextNexusTree(stree_ifile,&species_tree_str)," : getting next Nexus tree");
        }
        else if (n_iltrees>0)
        {
            ErrorReporter(InitNexusParser(ltree_ifile),": Parsing Nexus tree");
        }
        
        if (verbosity>2)
            PrintSettingsSloop(species_tree_str,gen_time,locus_tree_str,ltree_iname,n_iltrees,sb_rate,sd_rate,bds_leaves,bds_length,b_rate,d_rate,t_rate,gc_rate,t_kind,outgroup,min_lleaves,min_lsleaves,ind_per_sp,nl_trees,ng_trees,Ne,mu,alpha_s,alpha_l,salpha_g,epsilon_brent,verbosity,out_name,curr_stree);
        
        // ******
        /// Initialization of setting-dependent variables
        n_ldigits=(int)count_intdigits((long)get_sampling(nl_trees),0);
        n_gdigits=(int)count_intdigits((long)ng_trees,0);
        
#ifndef NO_OUT
#ifdef SORTHOLOGS
        //Sorthologs
        sorthologs_outname=realloc(sorthologs_outname,(strlen(sorthologs_prefix)+n_gdigits+strlen(sorthologs_midfix)+n_ldigits+strlen(sorthologs_sufix)+1)*sizeof(char));
        t_ndup=0;
        t_ngc=0;
        t_ntrfr=0;
        t_ddup=0;
        t_dgc=0;
        t_dtrfr=0;
        //
#endif
#endif
        // ******
        /// Initialization of file-I/O variables and s_tree output opening</dd></dl>
#ifndef NO_OUT
        if (verbosity>4)
        {
            printf("Generating the output directory %.*d... ",n_sdigits,curr_stree);
#ifdef DBG
            fflush(stdout);
#endif
        }
        sprintf(curr_outdir,"%.*d",n_sdigits,curr_stree);
        
        if(mkdir(curr_outdir, S_IRWXU | S_IRWXG | S_IRWXO)!=0)
        {
            if (errno==EEXIST)
            {
                if (verbosity>4)
                {
                    printf("\n\t WARNING, The output directory for the replicate %d already exists\n",curr_stree);
#ifdef DBG
                    fflush(stdout);
#endif
                }
            }
            else
            {
                fprintf(stderr,"\n\t ERROR, with the output folder for the replicate %d: %d\n",curr_stree,errno);
#ifdef DBG
                fflush(stderr);
#endif
                return (IO_ERROR);
            }
        }
        error=chdir(curr_outdir);
        if(error!=0)
        {
            fprintf(stderr,"\n\t ERROR, The output folder %s, is not accesible\n",curr_outdir);
#ifdef DBG
            fflush(stderr);
#endif
            return IO_ERROR;
        }
        
        g_outname=malloc((strlen(g_prefix)+n_ldigits+strlen(tree_sufix)+1)*sizeof(char));
        
        if (recon>0)
        {
            if (recon!=2)
            {
                reconsl_outname=malloc((n_ldigits+strlen(reconsl_sufix2)+1)*sizeof(char));
                reconlg_outname=malloc((n_ldigits+1+n_gdigits+1+strlen(reconlg_sufix2)+1)*sizeof(char));
            }
            if (recon!=1)
            {
                recon_outname=malloc((n_ldigits+1+n_gdigits+strlen(recon_sufix2)+1)*sizeof(char));
            }
        }
        
        if (locus_tree_str==NULL && n_iltrees==0)
        {
            if (verbosity>4)
            {
                printf("Opening the s_trees.tree file... ");
#ifdef DBG
                fflush(stdout);
#endif
            }
            if ((s_outfile=fopen(s_outname, "w"))==NULL)
            {
                perror("Error opening s_trees.tree file..");
                ErrorReporter(IO_ERROR,NULL);
            }
            if (verbosity>4)
            {
                printf("Done\n");
#ifdef DBG
                fflush(stdout);
#endif
            }
        }
        
        if (stats==1)
        {
            if (verbosity>4)
            {
                printf("Opening and initializing stats.txt... ");
#ifdef DBG
                fflush(stdout);
#endif
            }
            if ((stat_outfile=fopen(stat_outname, "w"))==NULL)
            {
                perror("Error opening stats.txt:");
                ErrorReporter(IO_ERROR,NULL);
            }
            fprintf(stat_outfile,"L_tree;N_losses;N_duplications;N_transfers;N_lt_totaltips;N_gt_presenttips;Mean_gt_height(cu);Mean_gt_height(ec);Mean_extra_lineages\n");
            if (verbosity>4)
            {
                printf("Done\n");
#ifdef DBG
                fflush(stdout);
#endif
            }
        }
        
        if (verbosity>4)
        {
            printf("Opening the l_trees.tree file... ");
#ifdef DBG
            fflush(stdout);
#endif
        }
        
        if ((l_outfile=fopen(l_outname,"w"))==NULL)
        {
            perror("Error opening l_trees.tree file: ");
            ErrorReporter(IO_ERROR,NULL);
        }
        if (verbosity>4)
        {
            printf("Done\n");
#ifdef DBG
            fflush(stdout);
#endif
        }
        
        if (weirdness!=0)
        {
            if ((weirds_outfile=fopen(weirds_outname, "w"))==NULL)
            {
                perror("Error opening weirds file:");
                ErrorReporter(IO_ERROR,NULL);
            }
            weirdg_outname=malloc((strlen(g_prefix)+n_ldigits+strlen(weirdg_sufix)+1)*sizeof(char));
        }
#endif
        
        // ******
        /// <dl><dt>Species tree and related structures preparation</dt><dd>
        
        // ******
        /// <dl><dt>If the locus tree is not fixed</dt><dd>
        if (locus_tree_str==NULL && n_iltrees==0)
        {
            if(species_tree_str!=NULL || n_istrees>0)
            {
                // *****
                /// Species tree allocation (\ref ParseNexusSTree) and collapse-reallocation (\ref CollapseSTree) if it is fixed. Pre-order if the locus tree will be simulated, post-order if it will not.
                sp_tree=ParseNexusSTree(species_tree_str,&names,verbosity,get_sampling(gen_time)!=1?get_sampling(gen_time):1,get_sampling(Ne),get_sampling(mu),get_sampling(ind_per_sp));
                CollapseSTree(sp_tree,(get_sampling(b_rate) == get_sampling(d_rate)) && get_sampling(b_rate)==0?1:0); //Memory reallocation in post-order
                
                if (sp_tree->n_leaves<min_lleaves)
                    ErrorReporter(SETTINGS_ERROR,"\n\t\tERROR:The minimum number of locus tree leaves is bigger than the number of leaves of the fixed species tree. Please, check the -Ll parameter and the input species tree.\n");
                
            }
            else
            {
                if (verbosity>2)
                {
                    printf("\nSimulating the species tree...");
#ifdef DBG
                    fflush(stdout);
#endif
                }
                // ******
                /// Species tree simulation using a birth-death process (\ref NewBDSTree) and collapse-reallocation (\ref CollapseSTree) if it is not fixed. Pre-order if the locus tree will be simulated, post-order if it will not.
                ErrorReporter(NewBDSTree(&sp_tree,get_sampling(bds_leaves), get_sampling(bds_length), get_sampling(sb_rate), get_sampling(sd_rate),get_sampling(gen_time)!=1?get_sampling(gen_time):1,get_sampling(Ne),get_sampling(mu),get_sampling(ind_per_sp),get_sampling(outgroup),0,1,r, verbosity), ": simulating the species tree");//complete 0 and SSA simulation at least.
                CollapseSTree(sp_tree,(get_sampling(b_rate) == get_sampling(d_rate)) && get_sampling(b_rate)==0?1:0);
                
                if (sp_tree->n_leaves<min_lleaves)
                {
                    if (verbosity>2)
                    {
                        printf("\n\t\tThe number of minimum locus tree leaves will be changed to the number of species tree leaves to avoid infinite loops.\n\t\t");
                    }
#ifdef DBG
                    fflush(stdout);
#endif
                    min_lleaves=sp_tree->n_leaves;
                }
                
            }
            
            // *****
            /// Species tree modification (substitution rate and generation time heterogeneities)
            if (get_sampling(alpha_s)!=0)
            {
                if (verbosity>2)
                {
                    printf("\tGenerating lineage specific substitution rate heterogeneity...");
#ifdef DBG
                    fflush(stdout);
#endif
                }
                Rateheter_lineagespec(sp_tree, get_sampling(alpha_s), r, weirds_outfile);
                if (verbosity>2)
                {
                    printf(" Done\n");
#ifdef DBG
                    fflush(stdout);
#endif
                }
                
            }
                //\todo Implement lineage specific generation time heterogeneity
            
//            if (get_sampling(alpha_X)!=0)
//            {
//                if (verbosity>2)
//                {
//                    printf("\tGenerating lineage specific generation time heterogeneity...");
//#ifdef DBG
//                    fflush(stdout);
//#endif
//                }
//                GenTimeheter_lineagespec(sp_tree, get_sampling(alpha_X?), r, weirdx?);
//                if (verbosity>2)
//                {
//                    printf(" Done\n\t");
//#ifdef DBG
//                    fflush(stdout);
//#endif
//                }
//
//            }
            
            if (verbosity>2)
            {
                printf("Done\n");
#ifdef DBG
                fflush(stdout);
#endif
            }
            
            // *****
            /// <dl><dt>Pointer allocation for locus tree simulation</dt><dd>
            if (is_variable(lb_rate) || is_variable(ld_rate) || is_variable(lt_rate) || is_variable(lgc_rate))
            {
                // ****
                /// Intermediate node pointers allocation </dd></dl></dd></dl>
                node_ptrs=calloc(MAX_LEAVES,sizeof(l_node *));
                ErrorReporter((long int)node_ptrs, NULL);
            }
        }
        
        // *******
        /// <dl><dt>SIMULATION</dt><dd>
        
        switch (verbosity)
        {
            case 0:
                break;
            case 1:
                printf("Replicate %u of %u: Simulating %u gene trees from %d locus trees... ",curr_stree,ns_trees,nl_trees.value.i*ng_trees,nl_trees.value.i);
                break;
            case 2:
                printf("\nReplicate %u of %u: Simulating %u gene trees from %d locus trees...\n",curr_stree,ns_trees,nl_trees.value.i*ng_trees,nl_trees.value.i);
                break;
            default:
                break;
        }
#ifdef DBG
        fflush(stdout);
#endif
        
#ifndef NO_OUT
        if (s_outfile!=NULL)
        {
            WriteSTreeFile(s_outfile,sp_tree,names,0);
        }
        if (db>0) //There is an species tree
        {
            if (locus_tree_str!=NULL || n_iltrees!=0)
            {
                ErrorReporter(WriteSTreeDB(&database,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), ": writting the species tree table");
            }
            else
            {
                Measure_ST_height(sp_tree, &st_height, CU);
                Measure_ST_length(sp_tree, &st_length, CU);
                ErrorReporter(WriteSTreeDB(&database, sp_tree->n_leaves, get_sampling(sb_rate),get_sampling(sd_rate),get_sampling(bds_leaves),get_sampling(bds_length),st_height,st_length,(*sp_tree->root->children)->gen_length, get_sampling(ind_per_sp), get_sampling(nl_trees),get_sampling(b_rate),get_sampling(d_rate),get_sampling(t_rate),get_sampling(gc_rate),get_sampling(alpha_s), get_sampling(alpha_l), get_sampling(salpha_g), get_sampling(Ne), get_sampling(mu), get_sampling(gen_time)), ": writting the species tree table");
            }

        }
#endif
        
        // ********
        /// <dl><dt> Main loop for each locus tree </dt><dd>
        prev_ltree_eq_stree=0;
        
        for (curr_ltree=1; curr_ltree<=get_sampling(nl_trees);++curr_ltree)
        {
            ErrorReporter(sample_distr(r,5,&lb_rate,&ld_rate,&lt_rate,&lgc_rate,&lalpha_g),": sampling distributions");
            ErrorReporter(CheckSampledSettingsLloop(lb_rate,ld_rate,lt_rate,lgc_rate,lalpha_g), ": improper sampled values");
            
            if(get_sampling(alpha_l)!=0)
                gamma_l=gsl_ran_gamma(r, get_sampling(alpha_l), 1/(get_sampling(alpha_l)));
            if (verbosity>2)
                PrintSettingsLloop(curr_ltree,lb_rate,ld_rate,lt_rate,lgc_rate,gamma_l,lalpha_g);
            
            // *******
            /// <dl><dt>Locus tree simulation (if it is necessary)</dt><dd>
            if (get_sampling(lb_rate) !=0 || get_sampling(ld_rate)!= 0 || get_sampling(lt_rate)!=0 || get_sampling(lgc_rate)!=0) // l_tree != s_tree
            {
                switch (verbosity)
                {
                    case 0:
                    case 1:
                    case 2:
                        break;
                    default:
                        printf("\tSimulating locus tree %d ... ",curr_ltree);
                        break;
                }
                    
#ifdef DBG
                fflush(stdout);
#endif

                // ******
                /// Locus tree simulation
                ErrorReporter(SimBDLHTree(sp_tree, &locus_tree, node_ptrs, get_sampling(lb_rate), get_sampling(ld_rate), get_sampling(lt_rate), get_sampling(lgc_rate),t_kind,r, min_lleaves, min_lsleaves, verbosity, &st_losses, &st_dups, &st_transf, &st_gc, &st_leaves, &st_gleaves), ": simulating a locus tree");
                
                // ******
                /// Species tree reindexation in post-order
                i=0;
                PostReorderSNodes(sp_tree->root,&i);

                // ******
                /// Reallocation of locus tree in post-order and reconciliation output (if necessary)
                if (recon>1) //If I used includelosses (more than one extra lineage per loss), this should be st_losses*includelosses
                    locus_tree->n_gleaves+=st_losses;
                
                ErrorReporter(CollapseLTree(locus_tree,1,0,(st_dups==0 && st_transf==0 && st_gc==0)?0:1), NULL);
#ifndef NO_OUT
                if (recon>0 && recon!=2)
                {
                    sprintf(reconsl_outname,"%.*d%s",n_ldigits,curr_ltree,reconsl_sufix2);
                    WriteReconSL(sp_tree, locus_tree, names, reconsl_outname);
                }
#endif
                
                // ******
                /// Gene tree allocation to perform its posterior simulation</dd></dl>
                if (gene_tree!=NULL)
                    FreeGTree(&gene_tree, 1);
                gene_tree=NewGTree((locus_tree->n_gleaves*2)-1,locus_tree->max_children, locus_tree->gen_time); //Maximum number of nodes.
                
            }
            // ******
            /// <dl><dt>Else obtains the fixed locus tree, either a fixed one or by copying it from the species tree</dt><dd>
            else if (locus_tree_str!=NULL || n_iltrees>0)
            {
                // ***
                /// Fixed case
                if (verbosity>2)
                {
                    printf("\tObtaining the preset locus tree (there is no birth-death process)... ");
#ifdef DBG
                    fflush(stdout);
#endif
                }
                if (n_iltrees>0)
                {
                    ErrorReporter(NextNexusTree(ltree_ifile,&locus_tree_str),": parsing the next Nexus tree");
                }
                
                // ***
                /// Locus tree allocation and collapse-reallocation (post-order)
                locus_tree=ParseNexusLTree(locus_tree_str, &names, verbosity,get_sampling(gen_time)!=1?get_sampling(gen_time):1,get_sampling(Ne),get_sampling(mu),get_sampling(ind_per_sp));
                CollapseLTree(locus_tree,1,0,1);//Memory reallocation in post-order
                
                // ***
                /// Statistical measurements
                if (stats==1)
                {
                    st_leaves=locus_tree->n_leaves;
                    st_gleaves=locus_tree->n_gleaves;
                    ErrorReporter(Count_duplications(locus_tree,&st_dups),NULL);
                }
                
                // ****
                /// Gene tree allocation
                gene_tree=NewGTree((locus_tree->n_gleaves*2)-1,locus_tree->max_children,locus_tree->gen_time);
            }
            else if (prev_ltree_eq_stree==0)
            {
                // ****
                /// From the species tree
                
                prev_ltree_eq_stree=1;
                
                // ****
                /// Locus tree and working locus tree allocation
                locus_tree=NewLTree(sp_tree->n_nodes, sp_tree->n_leaves, sp_tree->n_gleaves, sp_tree->max_children,sp_tree->gen_time, sp_tree->Ne, sp_tree->mu,0);
                
                // ****
                /// Gene tree allocation
                gene_tree=NewGTree((sp_tree->n_gleaves*2)-1,locus_tree->max_children, locus_tree->gen_time); //Gene tree does not allow polytomies.
                
                if (verbosity>2)
                {
                    printf("\tInvariable locus tree: Copying the species tree as locus tree... ");
#ifdef DBG
                    fflush(stdout);
#endif
                }
                
                // ****
                /// Direct conversion of s_tree into l_tree.
                CopyStoLTree(sp_tree, locus_tree);
                
                // ***
                /// Statistical measurements </dd></dl>
                if (stats==1)
                {
                    st_leaves=locus_tree->n_leaves;
                    st_gleaves=locus_tree->n_gleaves;
                }
            }
            else
            {
                if (verbosity>2)
                {
                    printf("\tInvariable locus tree: Reusing the previous locus tree... ");
#ifdef DBG
                    fflush(stdout);
#endif
                }
            }
            
            switch (verbosity)
            {
                case 0:
                case 1:
                case 2:
                    break;
                case 3:
                    printf(" Done");
                    break;
                default:
                    printf("\n\tDone: ");
                    WriteLTree(locus_tree,names,0);
                    break;
            }
#ifdef DBG
            fflush(stdout);
#endif
            
            // *****
            /// Locus tree modification (substitution rate heterogeneity)
            
            if (get_sampling(alpha_l)!=0)
            {
                if (verbosity>2)
                {
                    printf("\n\tGenerating gene family (locus tree) specific substitution rate heterogeneity...");
#ifdef DBG
                    fflush(stdout);
#endif
                }
                Rateheter_genespec(locus_tree, sp_tree->mu*gamma_l); //gsl_rang_gamma(r,shape,scale) -> Mean= shape*scale
                if (verbosity>2)
                {
                    printf(" Done\n");
#ifdef DBG
                    fflush(stdout);
#endif
                }
            }
            
            // *******
            // Locus tree output
#ifndef NO_OUT
            WriteLTreeFile(l_outfile,locus_tree,names,0);

            if (db>0)
            {
                ErrorReporter(WriteLTreeDB(&database, curr_ltree, curr_stree, get_sampling(lb_rate), get_sampling(ld_rate), get_sampling(lt_rate), get_sampling(lgc_rate),locus_tree->n_leaves, st_dups, st_losses,st_transf,st_gc,gamma_l,get_sampling(lalpha_g)),": writting the locus tree rable");
                t_n_ltree++;
            }
            
            if (stats==1)
            {
                t_height_bl=0;
                t_height_cu=0;
                t_lcoals=0;
            }
            
            // *******
            /// <dl><dt>Gene tree simulation</dt><dd>
            
            // ******
            // Gene tree output opening
            
            //Gene tree I/O
            sprintf(g_outname,"%s%.*d%s",g_prefix,n_ldigits,curr_ltree,tree_sufix);
            if (weirdness!=0)
                sprintf(weirdg_outname,"%s%.*d%s",g_prefix,n_ldigits,curr_ltree,weirdg_sufix);
            
            //I/O s_trees and l_trees opening
            if (verbosity>4)
            {
                printf("\n\t\tOpening and initializing gene tree output file for this locus and replicate... ");
#ifdef DBG
                fflush(stdout);
#endif
            }
            if ((g_outfile=fopen(g_outname, "w"))==NULL)
            {
                perror("Error opening g_tree file: ");
                ErrorReporter(IO_ERROR,NULL);
            }
            if (weirdness!=0 && (weirdg_outfile=fopen(weirdg_outname, "w"))==NULL)
            {
                perror("Error opening weirdg_tree file: ");
                ErrorReporter(IO_ERROR,NULL);
            }
            if (verbosity>4)
            {
                printf("Done");
#ifdef DBG
                fflush(stdout);
#endif
            }
#endif
            // ******
            /// <dl><dt> Main loop for each gene tree </dt><dd>
            for (curr_gtree=1;curr_gtree<=ng_trees; ++curr_gtree)
            {
                ErrorReporter(sample_distr(r,1,&alpha_g),": sampling distributions");
                ErrorReporter(CheckSampledSettingsGloop(alpha_g), ": improper sampled values");
                if (verbosity>2)
                    PrintSettingsGloop(curr_gtree, alpha_g);
                
                // ****
                /// Reconciliation output file name generation
#ifndef NO_OUT
                if (recon>0)
                {
                    if (recon!=2)
                        sprintf(reconlg_outname,"%.*d%s%.*d%s",n_ldigits,curr_ltree,"l",n_gdigits,curr_gtree,reconlg_sufix2);
                    if (recon!=1)
                        sprintf(recon_outname,"%.*d%s%.*d%s",n_ldigits,curr_ltree,"l",n_gdigits,curr_gtree,recon_sufix2);
                }
                
#endif
                
                // IO
                
                switch (verbosity)
                {
                    case 0:
                    case 1:
                        break;
                    case 2:
                        printf("\tSimulating gene tree %d... ",(curr_ltree-1)*ng_trees+curr_gtree);
                        break;
                    default:
                        printf("\n\t\tSimulating gene tree %d using a %s process...",curr_gtree,(st_dups==0 && st_transf==0 && st_gc==0)?"multispecies coalescent":"multilocus coalescent");
                        break;
                }
#ifdef DBG
                fflush(stdout);
#endif

                // ****
                /// Gene tree simulation
                if (st_dups==0 && st_transf==0 && st_gc==0)
                    ErrorReporter(SimMSCGTree(locus_tree,&gene_tree,names,epsilon_brent,r,&st_lcoals,recon>1?1:0,verbosity,get_sampling(gen_time)),": simulating a gene tree");
                else
                    ErrorReporter(SimMLCGTree(locus_tree,&gene_tree,names,epsilon_brent,r,&st_lcoals,recon>1?1:0,verbosity,get_sampling(gen_time)),": simulating a gene tree");
                
                // ****
                /// <dl><dt>Gene tree bl modifications</dt><dd>
                
                // *****
                /// Branch length heterogeneity generation</dd></dl>
                
                if (get_sampling(alpha_g)!=0)
                {
                    if (verbosity>2)
                    {
                        printf("\n\t\tGenerating gene tree branch specific substitution rate heterogeneity...");
#ifdef DBG
                        fflush(stdout);
#endif
                    }
                    Rateheter_GTbranchspec(gene_tree,get_sampling(alpha_g),r, weirdg_outfile);
                    if (verbosity>3)
                    {
                        printf(" Done");
#ifdef DBG
                        fflush(stdout);
#endif
                    }
                    else if (verbosity>2)
                    {
                        printf(" Done\n\t");
                    }
                }
                
                // ****
                /// Writes the current gene tree in Newick format
                
#ifndef NO_OUT
                WriteGTreeFile(g_outfile,gene_tree, names);
                
#endif
                switch (verbosity)
                {
                    case 0:
                    case 1:
                        break;
                    case 3:
                    case 2:
                        printf("\tDone\n");
                        break;
                    default:
                        printf("\n\t\tDone: ");
                        WriteGTree(gene_tree,names);
                        printf("\n");
                        break;
                }
#ifdef DBG
                fflush(stdout);
#endif

#ifndef NO_OUT
                // ****
                /// Reconciliation output and closing
                if (recon>0)
                {
                    if (recon!=2)
                        ErrorReporter(WriteReconLG(gene_tree, names, reconlg_outname),": obtaining LG reconciliation");
                }
#endif
                // ****
                /// Statistical measurements</dd></dl></dd></dl>
                
#ifdef SORTHOLOGS
                //This function is not intended for general usage.
                if (sorthologs_test==1)
                {
#ifndef NO_OUT
                    sprintf(sorthologs_outname,"%s%.*d%s%.*d%s",sorthologs_prefix,n_gdigits,curr_gtree,sorthologs_midfix,n_ldigits,curr_ltree,sorthologs_sufix);

                    if ((sorthologs_outfile=fopen(sorthologs_outname, "w"))==NULL)
                    {
                        perror("Error opening sorthologs file:");
                        ErrorReporter(IO_ERROR);
                    }
#endif
                    if (st_dups>0)
                    {
                        ErrorReporter(MeasureMRCAEVdistance(gene_tree,DUP,&dupdistances, st_dups, GL));

                        for (i=0; i<st_dups; ++i)
                        {
#ifndef NO_OUT
                            fprintf(sorthologs_outfile,"DUP\t%lf\n",*(dupdistances+i));
#endif
                            if (*(dupdistances+i)!=-1)
                            {
                                t_ddup+=*(dupdistances+i);
                                ++t_ndup;
                            }
                        }
                    }
                    if (st_gc>0)
                    {
                        ErrorReporter(MeasureMRCAEVdistance(gene_tree,GC,&gcdistances, st_gc, GL));

                        for (i=0; i<st_gc; ++i)
                        {
#ifndef NO_OUT
                            fprintf(sorthologs_outfile,"GC\t%lf\n",*(gcdistances+i));
#endif
                            if (*(gcdistances+i)!=-1)
                            {
                                t_dgc+=*(gcdistances+i);
                                ++t_ngc;
                            }
                        }

                    }
                    if (st_transf>0)
                    {
                        ErrorReporter(MeasureMRCAEVdistance(gene_tree,TRFR,&trfrdistances, st_transf, GL));

                        for (i=0; i<st_transf; ++i)
                        {
#ifndef NO_OUT
                            fprintf(sorthologs_outfile,"TRFR\t%lf\n",*(trfrdistances+i));
#endif
                            if (*(trfrdistances+i)!=-1)
                            {
                                t_dtrfr+=*(trfrdistances+i);
                                ++t_ntrfr;
                            }
                        }

                    }
#ifndef NO_OUT
                    fclose(sorthologs_outfile);
#endif

                }
#endif
                if (db>0 || stats>0)
                {
                    Measure_GT_height(gene_tree,&height_cu,CU);
                    Measure_GT_length(gene_tree,&length_bl,BL);
                    
                }
                if (stats>0)
                {
                    Measure_GT_height(gene_tree,&height_bl,1);
                    t_height_cu+=height_cu;
                    t_height_bl+=height_bl;
                    t_lcoals+=st_lcoals;
                }
#ifndef NO_OUT
                if (db>0)
                {
                    ErrorReporter(WriteGTreeDB(&database, curr_gtree, t_n_ltree,curr_ltree, curr_stree, get_sampling(alpha_g),locus_tree->n_gleaves, st_lcoals, height_cu, length_bl),": writting gene tree table");
                }
                
#endif
                
            }
            
#ifndef NO_OUT
            // *******
            /// Statistical output </dd></dl>
            if (stats==1)
            {
                fprintf(stat_outfile,"%.*d;%d;%d;%d;%d;%d;%Le;%Le;%Le\n",n_ldigits,curr_ltree,st_losses,st_dups,st_transf,st_leaves,st_gleaves,t_height_cu/ng_trees,t_height_bl/ng_trees,((long double)t_lcoals)/ng_trees);
            }
            // *******
            /// Gene tree I/O close</dd></dl>
            fclose(g_outfile);
            if (weirdness!=0)
                fclose(weirdg_outfile);
            
#endif
        }
        
#ifdef SORTHOLOGS
        //Sorthologs
        printf("\n\nSorthologs' paper related statistics: Duplication %Lf, Tranfers %Lf, Gene conversions %Lf\n", t_ddup/t_ndup, t_dtrfr/t_ntrfr, t_dgc/t_ngc);
        //
#endif
        error=chdir("..");
        if(error!=0)
        {
            fprintf(stderr,"\n\t ERROR, impossible to come back to the parent directory\n");
            return IO_ERROR;
        }
        
        // Freeing allocated memory
        
        if (names!=NULL)
            FreeNames(&names);
        if (sp_tree!=NULL)
            FreeSTree(&sp_tree);
        if (locus_tree!=NULL)
            FreeLTree(&locus_tree);
        if (gene_tree!=NULL)
            FreeGTree(&gene_tree, 1);
        
#ifndef NO_OUT
        if (g_outname!=NULL)
            free(g_outname);
        if (weirdg_outname!=NULL)
            free(weirdg_outname);
        if (recon_outname!=NULL)
            free(recon_outname);
        if (reconlg_outname!=NULL)
            free(reconlg_outname);
        if (reconsl_outname!=NULL)
            free(reconsl_outname);
        if (l_outfile!=NULL)
            fclose(l_outfile);
        if (stat_outfile!=NULL)
            fclose(stat_outfile);
        if (s_outfile!=NULL)
            fclose(s_outfile);
        if (weirds_outfile!=NULL)
            fclose(weirds_outfile);
#endif
        switch (verbosity)
        {
            case 1:
            case 2:
            case 3:
            case 4:
                printf("Done\n");
                break;
            default:
                break;
        }
        
    }
    
    switch (verbosity)
    {
        case 0:
            printf("Done\n");
            break;
        case 1:
        case 2:
            break;
        default:
            printf("\nDone \n");
            break;
    }

#ifdef DBG
    fflush(stdout);
#endif
    
    // ********
    /// Freeing of allocated memory and I/O closing </dd></dl>
    if (species_tree_str!=NULL)
        free(species_tree_str);
    if (locus_tree_str!=NULL)
        free(locus_tree_str);
    if (out_name!=NULL)
        free(out_name);
    if (db_outname!=NULL)
        free(db_outname);
    if (stree_iname!=NULL)
        free(stree_iname);
    if (ltree_iname!=NULL)
        free(ltree_iname);
#ifndef NO_OUT
    if (db>0)
        ErrorReporter(CloseDB(&database),NULL);
#endif
#ifdef SORTHOLOGS
    if (dupdistances!=NULL)
        free(dupdistances);
    if (gcdistances!=NULL)
        free(gcdistances);
    if (trfrdistances!=NULL)
        free(trfrdistances);
    if (sorthologs_outname!=NULL)
        free(sorthologs_outname);
#endif
    
    gsl_rng_free (r);
    
    fflush(stdout);
    return (EXIT_SUCCESS);
    /// </dd></dl>
}


// *********************** Declaration of functions ***************************** //

/**
 * \name I/O
 * Input and output functions
 *******************************************************************************/

///@{

/**
 * Writes a short description of program usage.
 * \param void
 *******************************************************************************/

void PrintUsage(void)
{
    printf("\nUsage: ./SimPhy -[Parameter code] value(i|d|c|b|s) ...\n\nValue kinds\n\ti=integer\n\td=decimal\n\tc=character string\n\tb=boolean(1 or 0)\n\ts=sampling notation\n\nPARAMETERS\n_________________\n\n-RX→ Replicates\n_________________\n\n\t-RS i: Number of species tree replicates (study replicates).\n\t-RL s: Distribution of number of locus trees per species tree.\n\t-RG i: Number of gene trees per locus tree.\n\n-HX→ Heterogeneity generators\n_________________\n\n\t-HS s: Distribution of species-specific branch rate heterogeneity modifiers.\n\t-HL s: Distribution of locus-specific rate heterogeneity modifiers.\n\t-HG s: Distribution of gene-tree-branch-specific rate heterogeneity modifiers.\n\t-P s: Distribution of global population sizes.\n\t-U s: Distribution of global substitution rates\n\t-E d: Precision of the Brent’s method for root-finding when sampling the multilocus coalescent. (Devs advice: do not touch if you do not really need it) \n\t-G s: Global generation time.\n\t-V [0,6]: Verbosity. The bigger the verbosity the slower SimPhy becomes. Levels over 3 may only make sense for debugging or study the way SimPhy works.\n\n-OX→ Output\n_________________\n\n\t-O c: Common output prefix.\n\t-OS b: Activates the statistical output.\n\t-OR b: Activates the reconciliation output.\n\t-OD b: Activates the SQLite database output.\n\t-OP b: Activates the output of a file with a summary of the sampled options .\n\t-OC b: Activates the output of original command line parameters and input configuration file if used.\n\n-SX → Species tree\n_________________\n\n\t-S c: Fixed species tree string in modified Newick format . \n\t-SR c: Nexus file containing species trees.\n\t-SB s: Species tree simulation parameter: Distribution of birth rates.\n\t-SD s: Species tree simulation parameter: Distribution of death rates.\n\t-ST s: Species tree simulation parameter: Distribution of maximum tree lengths.\n\t-SL s: Species tree simulation parameter: Distribution of number of species (species tree leaves).\n\t-SO s: Species tree simulation parameter: Distribution of multipliers of half of the tree length to be applied as the branch length of the internal node generated when introducing an outgroup. If this parameter is not set the outgroup is not simulated. \n\t-SI s: Species tree simulation parameter: Distribution of number of individuals per species.\n\n-LX → Locus tree\n_________________\n\n\t-L c: Fixed locus tree string in modified Newick format .\n\t-LR c: Nexus file containing locus trees .\n\t-LB s: Locus tree simulation parameter: Distribution of duplication rates.\n\t-LD s: Locus tree simulation parameter: Distribution of loss rates.\n\t-LT s: Locus tree simulation parameter: Distribution of transfer rates. \n\t-LG s: Locus tree simulation parameter: Distribution of gene conversion rates.\n\t-LK b: Locus tree simulation parameter: Determines whether the sampling of acceptors for both transfers and gene conversions depends (1) or not (0) on the evolutionary distance between candidate acceptors and donors.\n\t-LL s: Locus tree simulation parameter: Distribution of minimum number of leaves to keep the simulated locus tree.\n\t-LS s: Locus tree simulation parameter: Distribution of minimum number of leaves from different species to keep the simulated locus tree.\n\n-CX → Global options\n_________________\n\n\t-CS d: Random number generator seed.\n\t-I c: Input config file (only as a command line parameter)\n\n\nSampling notation\n___________________\n\t[.....]\n\nExample\n___________\n\nSimPhy -rs 100 -rl U10-100 -sb f0.000001 -sl U20-50 -ld f0.0000005 -lb f0.0000005 -lt f0.0000005 -p f10000 -g f1 -o test -os 1 -oc 1 -or 1 -od 1 -v 2\n\n");
    
    fflush(stdout);
    
}

/**
 * Gets settings from command line.
 *
 * \param argc
 *   Number of arguments passed by command line.
 * \param argv
 *   Array of arguments
 * \param ns_trees
 *   Number of species trees to simulate.
 * \param nl_trees
 *   Number of locus trees to simulate.
 * \param ng_trees
 *   Number of gene trees from each locus tree to simulate.
 * \param newick_stree
 *   Species tree string.
 * \param stree_ifile
 *   Species tree nexus filename.
 * \param n_istrees
 *   Number of input (fixed) species trees present in stree_ifile.
 * \param stree_ifile
 *   Species tree nexus file.
 * \param gen_time
 *   Generation time.
 * \param newick_ltree
 *   Locus tree string.
 * \param ltree_ifile
 *   Locus tree nexus file.
 * \param n_iltrees
 *   Number of input (fixed) locus trees present in ltree_ifile.
 * \param ltree_ifile
 *   Locus tree nexus filename.
 * \param b_rate
 *   Birth rate for the locus_tree simulation.
 * \param d_rate
 *   Death rate for the locus_tree simulation.
 * \param t_rate
 *   Transfer rate for the locus_tree simulation.
 * \param gc_rate
 *   Gene conversion rate for the locus_tree simulation.
 * \param lb_rate
 *   Birth rate for the locus_tree simulation.
 * \param ld_rate
 *   Death rate for the locus_tree simulation.
 * \param lt_rate
 *   Transfer rate for the locus_tree simulation.
 * \param lgc_rate
 *   Gene conversion rate for the locus_tree simulation.
 * \param t_kind
 *   Logical flag. 0=> RTRFR randomly sampled. 1=> RTRFR sampled with probability inversely related to distance (generations).
 * \param min_lleaves
 *   Minimum number of leaves for each locus tree.
 * \param min_lsleaves
 *   Minimum number of leaves from different species for each locus tree.
 * \param ind_per_sp
 *   Number of individuals per species.
 * \param sb_rate
 *   Birth rate for the species_tree simulation.
 * \param sd_rate
 *   Death rate for the species_tree simulation.
 * \param bds_leaves
 *   Max number of nodes in the b-d simulation of the species tree.
 * \param bds_length
 *   Max lenght in the b-d simulation of the species tree.
 * \param Outgroup
 *   Internal branch length deviation from a half of the ingroup height,for outgroup addition.
 * \param Ne
 *   Efective population size.
 * \param mu
 *   Substitution rate.
 * \param alpha_s
 *   Alpha parameter of the gamma distribution with mean=1 to sample rate heterogeneity mulipliers. Lineage specific rate heterogeneity.
 * \param alpha_l
 *   Alpha parameter of the gamma distribution with mean=1 to sample rate heterogeneity mulipliers. Gene family specific rate heterogeneity.
 * \param alpha_g
 *   Alpha parameter of the gamma distribution with mean=1 to sample rate heterogeneity mulipliers. Gene tree branch specific rate heterogeneity.
 * \param lalpha_g
 *   Alpha_g prior. Gene tree branch specific rate heterogeneity prior.
 * \param salpha_g
 *   Alpha_l prior (alpha g hyperprior). Gene tree branch specific rate heterogeneity hyperprior.
 * \param epsilon
 *  Epsilon of the convergence of the brent method for sampling Bounded multispecies coalescent
 * \param verbosity
 *   Config about verbosity.
 * \param out_name
 *   Name of the file with the gene trees in Newick format.
 * \param stats
 *   Logical flag, if ==1 it generates a csv with some stats.
 * \param recon
 *   Logical flag, if ==1 it generates a reconciliation between species and gene trees.
 * \return Error-control code.
 * \param db
 *   Logical flag, if ==1 it generates a SQLite database with 3 linked tables (Species_Trees, Locus_Trees and Gene_Trees) with some characteristics of each tree.
 * \param params
 *  Logical flag to activate an output file (prefix.params) with the general parameters of the simulation.
 * \param commands
 *  Logical flag to activate the output of both the original command-line arguments (prefix.command) and the input configuration file (prefix.conf)
 * \param u_seed
 *  Seed for the random number generator.
 * \param confile_name
 *  Configuration file name.
 * \return Error-control code.
 *
 * \note n_trees and newick_stree/newick_ltree are essential.
 * \todo Think about the idea of using this library https://github.com/Cofyc/argparse to get the options, because it is much more flexible.
 *******************************************************************************/

long int GetSettings(int argc, char **argv, int *ns_trees, sampling_unit *nl_trees, int *ng_trees, char ** newick_stree, FILE **stree_ifile, int *n_istrees, char **stree_iname, sampling_unit * gen_time, char ** newick_ltree, FILE **ltree_ifile, int *n_iltrees, char **ltree_iname,sampling_unit * b_rate, sampling_unit * d_rate, sampling_unit *t_rate, sampling_unit *gc_rate,sampling_unit * lb_rate, sampling_unit * ld_rate, sampling_unit *lt_rate, sampling_unit *lgc_rate, int *t_kind, int *min_lleaves, int *min_lsleaves, sampling_unit *ind_per_sp,sampling_unit * sb_rate, sampling_unit *sd_rate, sampling_unit * bds_leaves, sampling_unit * bds_length, sampling_unit *outgroup, sampling_unit *Ne,sampling_unit *mu, sampling_unit *alpha_s, sampling_unit *alpha_l, sampling_unit *alpha_g,sampling_unit *lalpha_g,sampling_unit *salpha_g, float *epsilon,int *verbosity,char **out_name, int *stats, int *recon, int *db, int *params, int *commands, float *u_seed, char **confile_name,const sampling_table sampling_vars)
{
    int i=0,read_file=0;
    char code=0,sub_code=0,buff_char=0;
    FILE *input_file=NULL;
    
    // ******
    /// <dl><dt> Function structure </dt><dd>
    
    // *****
    /// Control of number of arguments
//    if (argc<3 || argc%2==0)
//    {
//        fprintf(stderr,"Improper number of arguments, they have to be even and more than 3\n");
//        return (SETTINGS_ERROR);
//    }
    
    // *****
    /// <dl><dt> Arguments loop </dt><dd>
    for (i=1; i<argc;++i)
    {
        // ****
        /// Option/value detection. The first character has to be a '-' (values will be jumped during its parsing)
        buff_char=*argv[i]; //Reads a character
        if (buff_char!='-' && code==0) //Is not a - and is not placed after a recogniced code.
        {
            PrintXStringError(argv,i, "|<– ERROR in this parameter\n");
            return (SETTINGS_ERROR);
        }
        // ****
        /// <dl><dt>Option recognition and value parsing, taking into account 2 char options and skipping the next argument in the loop (value)</dt><dd>
        else if (buff_char == '-') //If the argument starts with and - reads the character to recognice the option.
        {
            code=*(argv[i]+1);
            sub_code=*(argv[i]+2);
            ++i;
            switch (toupper(code)) //If it is a lowercase converts it in the uppercase.
            {
                    // ***
                    /// <dl><dt>-Rx. Replicates</dt><dd>
                case 'R':
                    switch (toupper(sub_code))
                    {
                        // **
                        /// -RS. Number of study replicates (species trees replicates).
                    case 'S':
                        if(sscanf(argv[i],"%u",ns_trees)==0 || *ns_trees<1)
                        {
                            PrintXStringError(argv,i, "|<– ERROR in this parameter\n");
                            return(SETTINGS_ERROR);
                        }
                        break;
                        
                        // **
                        /// -RG. Number of gene trees for each locus tree.
                    case 'G':
                        if(sscanf(argv[i],"%u",ng_trees)==0 || *ng_trees<1)
                        {
                            PrintXStringError(argv,i, "|<– ERROR in this parameter\n");
                            return(SETTINGS_ERROR);
                        }
                        break;
                        
                        // **
                        /// -RL. Number of locus trees</dd></dl>
                    case 'L':
                        if(ParseSampling(argv[i],nl_trees,sampling_vars)!=NO_ERROR)
                        {
                            PrintXStringError(argv,i, "|<– ERROR in this parameter\n");
                            return SETTINGS_ERROR;
                        }
                        break;
                        
                    default:
                        PrintXStringError(argv,i, "|<– Unrecognized parameter\n");
                        return(SETTINGS_ERROR);
                        break;
                }
                    break;
                    // ***
                    /// <dl><dt>-Hx. Rate heterogeneity parameters</dt><dd>
                case 'H':
                    switch (toupper(sub_code))
                    {
                        // **
                        /// -HS. Alpha parameter for the heterogeneity gamma sampling of lineage specific substitution rates.
                    case 'S':
                        if(ParseSampling(argv[i],alpha_s,sampling_vars)!=NO_ERROR)
                        {
                            PrintXStringError(argv,i, "|<– ERROR in this parameter\n");
                            return SETTINGS_ERROR;
                        }
                        break;
                        // **
                        /// -HB. Reserved, it is not implemented yet
                        // **
                        /// -HL. Alpha parameter for the heterogeneity gamma sampling of gene family (gene tree) specific substitution rates
                    case 'L':
                        if(ParseSampling(argv[i],alpha_l,sampling_vars)!=NO_ERROR)
                        {
                            PrintXStringError(argv,i, "|<– ERROR in this parameter\n");
                            return SETTINGS_ERROR;
                        }
                        break;
                        // **
                        /// -HG. Alpha parameter for the heterogeneity gamma sampling of gene family branch specific substitution rates</dd></dl>
                    case 'G':
                        if(ParseSampling(argv[i],alpha_g,sampling_vars)!=NO_ERROR)
                        {
                            PrintXStringError(argv,i, "|<– ERROR in this parameter\n");
                            return SETTINGS_ERROR;
                        }
                        break;
                    default:
                        PrintXStringError(argv,i, "|<– Unrecognized parameter\n");
                        return(SETTINGS_ERROR);
                        break;
                }
                    break;
                    // ***
                    /// -P. Haploid population size
                case 'P':
                    if(ParseSampling(argv[i],Ne,sampling_vars)!=NO_ERROR)
                    {
                        PrintXStringError(argv,i, "|<– ERROR in this parameter\n");
                        return SETTINGS_ERROR;
                    }
                    break;
                    // ***
                    /// -U. Global substitution rate
                case 'U':
                    if(ParseSampling(argv[i],mu,sampling_vars)!=NO_ERROR)
                    {
                        PrintXStringError(argv,i, "|<– ERROR in this parameter\n");
                        return SETTINGS_ERROR;
                    }
                    break;
                    // ***
                    /// <dl><dt>-Ex. Epsilons for the bounded multispecies coalescent sampling </dt><dd>
                case 'E':
                    switch (toupper(sub_code))
                    {
                        // **
                        /// -E. Brent method </dd></dt>
                    case '\0':
                        if(sscanf(argv[i],"%f",epsilon)==0 || *epsilon<0)
                        {
                            PrintXStringError(argv,i, "|<– ERROR in this parameter\n");
                            return (SETTINGS_ERROR);
                        }
                        break;

                    default:
                        PrintXStringError(argv,i, "|<– Unrecognized parameter\n");
                        return (SETTINGS_ERROR);
                        break;
                    }
                    
                    break;
                    // ***
                    /// -G. Generation time (branch lengths as time instead of generations, input and output trees, inside the software, the trees are in generations)
                case 'G':
                    if(ParseSampling(argv[i],gen_time,sampling_vars)!=NO_ERROR)
                    {
                        PrintXStringError(argv,i, "|<– ERROR in this parameter\n");
                        return SETTINGS_ERROR;
                    }
                    break;
                    // ***
                    /// -V. Verbosity
                case 'V':
                    if(sscanf(argv[i],"%u",verbosity)==0)
                    {
                        PrintXStringError(argv,i, "|<– ERROR in this parameter\n");
                        return (SETTINGS_ERROR);
                    }
                    break;
                    // ***
                    /// <dl><dt>-Ox. Output related options</dt><dd>
                case 'O':
                    switch (toupper(sub_code))
                {
                        // **
                        /// -O. Output common prefix
                    case '\0':
                        *out_name=calloc(strlen(argv[i])+1,sizeof(char));
                        sscanf(argv[i],"%s",*out_name);
                        break;
                        // **
                        /// -OS. Statistics output
                    case 'S':
                        if(sscanf(argv[i],"%u",stats)==0)
                        {
                            PrintXStringError(argv,i, "|<– ERROR in this parameter\n");
                            return (SETTINGS_ERROR);
                        }
                        break;
                        // **
                        /// -OR. Reconciliation output
                    case 'R':
                        if(sscanf(argv[i],"%u",recon)==0)
                        {
                            PrintXStringError(argv,i, "|<– ERROR in this parameter\n");
                            return (SETTINGS_ERROR);
                        }
                        break;
                        // **
                        /// -OD. SQLite relational database output
                    case 'D':
                        if(sscanf(argv[i],"%u",db)==0)
                        {
                            PrintXStringError(argv,i, "|<– ERROR in this parameter\n");
                            return (SETTINGS_ERROR);
                        }
                        break;
                        // **
                        /// -OP. Params output
                    case 'P':
                        if(sscanf(argv[i],"%u",params)==0)
                        {
                            PrintXStringError(argv,i, "|<– ERROR in this parameter\n");
                            return (SETTINGS_ERROR);
                        }
                        break;
                        // **
                        /// -OC. Commands output </dd></dl>
                    case 'C':
                        if(sscanf(argv[i],"%u",commands)==0)
                        {
                            PrintXStringError(argv,i, "|<– ERROR in this parameter\n");
                            return (SETTINGS_ERROR);
                        }
                        break;
                    default:
                        PrintXStringError(argv,i, "|<– Unrecognized parameter\n");
                        return (SETTINGS_ERROR);
                        break;
                }
                    break;
                    // ***
                    /// <dl><dt>-Sx. Species tree</dt><dd>
                case 'S':
                    switch (toupper(sub_code))
                {
                        // **
                        /// -S. One fixed species tree. Parses the string as a species tree
                    case '\0':
                        *newick_stree=calloc(strlen(argv[i])+1,sizeof(char));
                        sscanf(argv[i],"%s",*newick_stree);
                        break;
                        // **
                        /// -SR. Fixed species trees. Parses the string as the Nexus file containing the species tree/s
                    case 'R':
                        *stree_iname=calloc(strlen(argv[i])+1,sizeof(char));
                        sscanf(argv[i],"%s",*stree_iname);
                        if ((*stree_ifile=fopen(*stree_iname, "r"))==NULL)
                        {
                            PrintXStringError(argv,i, "|<- ERROR in this parameter\n");
                            perror("Error opening input_file: ");
                            ErrorReporter(IO_ERROR,NULL);
                        }
                        else
                        {
                            ErrorReporter(NNexusTrees(*stree_ifile,n_istrees),NULL);
                        }
                        break;
                        // **
                        /// <dl><dt>-Sx. Birth-death species tree</dt><dd>
                        // *
                        /// -SB. Birth rate
                    case 'B':
                        if(ParseSampling(argv[i],sb_rate,sampling_vars)!=NO_ERROR)
                        {
                            PrintXStringError(argv,i, "|<– ERROR in this parameter\n");
                            return SETTINGS_ERROR;
                        }
                        break;
                        // *
                        /// -SD. Death rate
                    case 'D':
                        if(ParseSampling(argv[i],sd_rate,sampling_vars)!=NO_ERROR)
                        {
                            PrintXStringError(argv,i, "|<– ERROR in this parameter\n");
                            return SETTINGS_ERROR;
                        }
                        break;
                        // *
                        /// -ST. Maximum time to stop the birth-death process.
                    case 'T':
                        if(ParseSampling(argv[i],bds_length,sampling_vars)!=NO_ERROR)
                        {
                            PrintXStringError(argv,i, "|<– ERROR in this parameter\n");
                            return SETTINGS_ERROR;
                        }
                        break;
                        // *
                        /// -SL. Desired number of leaves to stop the birth-death process.
                    case 'L':
                        if(ParseSampling(argv[i],bds_leaves,sampling_vars)!=NO_ERROR)
                        {
                            PrintXStringError(argv,i, "|<– ERROR in this parameter\n");
                            return SETTINGS_ERROR;
                        }
                        break;
                        // *
                        /// -SO. Deviation of a half of the ingroup heigth for the new internal branch due to the outgroup addition
                    case 'O':
                        if(ParseSampling(argv[i],outgroup,sampling_vars)!=NO_ERROR)
                        {
                            PrintXStringError(argv,i, "|<– ERROR in this parameter\n");
                            return SETTINGS_ERROR;
                        }
                        break;
                        // ***
                        /// -SI. Number of individuals per species</dd></dl></dd></dl>
                    case 'I':
                        if(ParseSampling(argv[i],ind_per_sp,sampling_vars)!=NO_ERROR)
                        {
                            PrintXStringError(argv,i, "|<– ERROR in this parameter\n");
                            return SETTINGS_ERROR;
                        }
                        break;
                        
                    default:
                        PrintXStringError(argv,i, "|<– Unrecognized parameter\n");
                        return (SETTINGS_ERROR);
                        break;
                }
                    break;
                    // ***
                    /// <dl><dt>-Lx. Locus tree.</dt><dd>
                case 'L':
                    switch (toupper(sub_code))
                {
                        // **
                        /// -L. Fixed locus tree. Parses the string as a locus tree
                    case'\0':
                        *newick_ltree=calloc(strlen(argv[i])+1,sizeof(argv[i]));
                        sscanf(argv[i],"%s",*newick_ltree);
                        break;
                        // **
                        /// -LR. Fixed locus trees. Parses the string as the Nexus file containing the locus tree/s
                    case 'R':
                        *ltree_iname=calloc(strlen(argv[i])+1,sizeof(char));
                        sscanf(argv[i],"%s",*ltree_iname);
                        if ((*ltree_ifile=fopen(*ltree_iname, "r"))==NULL)
                        {
                            PrintXStringError(argv,i, "|<- ERROR in this parameter\n");
                            perror("Error opening input_file: ");
                            ErrorReporter(IO_ERROR,NULL);
                        }
                        else
                        {
                            ErrorReporter(NNexusTrees(*ltree_ifile,n_iltrees),NULL);
                        }
                        break;
                        // **
                        /// -LB. Simulated locus tree. Birth rate.
                    case 'B':
                        if(ParseSampling(argv[i],b_rate,sampling_vars)!=NO_ERROR)
                        {
                            PrintXStringError(argv,i, "|<– ERROR in this parameter\n");
                            return SETTINGS_ERROR;
                        }
                        break;
                        // **
                        /// -LD. Simulated locus tree. Death rate.
                    case 'D':
                        if(ParseSampling(argv[i],d_rate,sampling_vars)!=NO_ERROR)
                        {
                            PrintXStringError(argv,i, "|<– ERROR in this parameter\n");
                            return SETTINGS_ERROR;
                        }
                        break;
                        // **
                        /// -LT. Simulated locus tree. Transfer rate.
                    case 'T':
                        if(ParseSampling(argv[i],t_rate,sampling_vars)!=NO_ERROR)
                        {
                            PrintXStringError(argv,i, "|<– ERROR in this parameter\n");
                            return SETTINGS_ERROR;
                        }
                        break;
                        // **
                        /// -LG. Simulated locus tree. Gene conversion rate.
                    case 'G':
                        if(ParseSampling(argv[i],gc_rate,sampling_vars)!=NO_ERROR)
                        {
                            PrintXStringError(argv,i, "|<– ERROR in this parameter\n");
                            return SETTINGS_ERROR;
                        }
                        break;
                        // **
                        /// -LK. Simulated locus tree. Aceptor sampling dependent (1) or independent (0) to the evolutionary distance.
                    case 'K':
                        if (sscanf(argv[i],"%u",t_kind)==0 || *t_kind<0 || *t_kind>1)
                        {
                            PrintXStringError(argv,i, "|<- ERROR in this parameter\n");
                            return (SETTINGS_ERROR);
                        }
                        break;
                        // **
                        /// -LL. Simulated locus tree. Minimum number of leaves of simulated locus trees.
                    case 'L':
                        if(sscanf(argv[i],"%u",min_lleaves)==0 || *min_lleaves<2)
                        {
                            PrintXStringError(argv,i, "|<– ERROR in this parameter\n");
                            return (SETTINGS_ERROR);
                        }
                        break;
                        // **
                        /// -LS. Simulated locus tree. Minimum number of leaves from different species of simulated locus trees.</dd></dl>
                    case 'S':
                        if(sscanf(argv[i],"%u",min_lsleaves)==0 || *min_lsleaves<2)
                        {
                            PrintXStringError(argv,i, "|<– ERROR in this parameter\n");
                            return (SETTINGS_ERROR);
                        }
                        break;
                    default:
                        PrintXStringError(argv,i, "|<– Unrecognized parameter\n");
                        return (SETTINGS_ERROR);
                        break;
                }
                    break;
                    // ***
                    /// <dl><dt>-Cx. Global options.</dt><dd>
                case 'C':
                    switch (toupper(sub_code))
                {       // **
                        /// -Cs. Seed for the random number generator.</dd></dl>
                    case 'S':
                    {
                        if(sscanf(argv[i],"%f",u_seed)==0)
                        {
                            PrintXStringError(argv,i, "|<– ERROR in this parameter\n");
                            return (SETTINGS_ERROR);
                        }
                        
                    }
                        break;
                    default:
                        PrintXStringError(argv,i, "|<– Unrecognized parameter\n");
                        return (SETTINGS_ERROR);
                        break;
                        
                        
                }
                    break;
                    // ***
                    /// -I. Input file</dd></dl>
                case 'I':
                    if ((input_file=fopen(argv[i], "r"))==NULL)
                    {
                        PrintXStringError(argv,i, "|<- ERROR in this parameter\n");
                        perror("Error opening input_file: ");
                        ErrorReporter(IO_ERROR,NULL);
                    }
                    else
                    {
                        *confile_name=argv[i];
                        if (argc>3)
                        {
                            printf("\n\tWARNING! You are using both the input file option -I and command line parameters. If you define a parameter using both input options, input file options will have preference\n");
                        }
                        read_file=1;
                    }
                    break;
                    
                default:
                    PrintXStringError(argv,i, "|<– Unrecognized parameter\n");
                    return (SETTINGS_ERROR);
                    break;
            }
            code=0;
            
        }
    }
    
    // *****
    /// <dl><dt>Getting input parameters from input file if it is present</dt><dd></dd></dl>
    switch (read_file)
    {
        case 1:
            
            ErrorReporter(GetSettingsFromFile(input_file, ns_trees, nl_trees, ng_trees, newick_stree, stree_ifile, n_istrees, stree_iname, gen_time, newick_ltree, ltree_ifile, n_iltrees, ltree_iname, b_rate, d_rate, t_rate, gc_rate, lb_rate, ld_rate, lt_rate, lgc_rate, t_kind, min_lleaves, min_lsleaves, ind_per_sp, sb_rate, sd_rate, bds_leaves, bds_length, outgroup, Ne, mu, alpha_s, alpha_l, salpha_g, lalpha_g, alpha_g, epsilon, verbosity, out_name, stats, recon, db, params, commands, u_seed, sampling_vars),"");
            
            fclose(input_file);
            break;
    }
    
    // *****
    /// <dl><dt>Error detection</dt><dd>
    
    // *****
    /// <dl><dt>Input trees/b-d parameters</dt><dd>
    
    // ****
    /// Locus tree and species tree fixed (not allowed)
    if (*newick_ltree!=NULL && *newick_stree!=NULL)
    {
        fprintf(stderr,"\n\tERROR!!! Using -S and -L (fixed species tree and locus tree) at the same time is not allowed. Please, use only a fixed locus tree (-L) or species tree (-S). In the second option, the species tree will be converted directly into a locus tree without duplication and losses. The first option allow you to add duplications (losses does not affect to the gene tree simulation.\n");
        return (SETTINGS_ERROR);
    }
    
    if (*newick_stree!=NULL && *n_istrees>0)
    {
        fprintf(stderr, "\n\tWARNING!!! Using a directly given fixed tree and a species tree input file is not allowed. Only the file will be used \n");
        free(*newick_stree);
        newick_stree=NULL;
    }
    
    if (*newick_ltree!=NULL && *n_iltrees>0)
    {
        fprintf(stderr, "\n\tWARNING!!! Using a directly given fixed tree and a locus tree input file is not allowed. Only the file will be used \n");
        free(*newick_ltree);
        newick_ltree=NULL;
    }
    
    // ****
    /// <dl><dt>No fixed trees</dt><dd>
    if ((*newick_stree==NULL && n_istrees==0) && (*newick_ltree==NULL && *n_iltrees==0))
    {
        // ***
        /// Species tree simulation parameters test</dd></dl>
        if ((is_sampling_set(*bds_length)|| is_sampling_set(*bds_leaves))==0)
        {
            fprintf(stderr,"\n\tERROR!!! You are not using neither a fixed locus tree (-L) nor a correctly specified species tree (fixed or simulated). To correctly specify a simulated species tree, you need at least one stopping rule, either bds_length or bds_leaves (-St, -Sl)\n");
            return (SETTINGS_ERROR);
        }
    }
    // ****
    /// <dl><dt>Fixed tree (locus or species)</dt><dd>
    if (((*newick_stree!=NULL || *n_istrees>0) && (*newick_ltree==NULL && *n_iltrees==0))||((*newick_ltree!=NULL || *n_iltrees>0) && (*newick_stree==NULL && *n_istrees==0)))
    {
        // ***
        /// Superfluous species tree simulation parameters check.
        if (is_sampling_set(*sb_rate) || is_sampling_set(*sd_rate) || is_sampling_set(*bds_length) || is_sampling_set(*bds_leaves))
        {
            fprintf(stderr, "\n\tWARNING!!! Using a fixed tree and species tree simulation parameters has no sense, and these simulation parameters will be ignored\n");
            ParseSampling("f0", sb_rate,sampling_vars);
            ParseSampling("f0", sd_rate,sampling_vars);
            ParseSampling("f0", bds_length,sampling_vars);
            ParseSampling("f0", bds_leaves,sampling_vars);
        }
        // ***
        /// Superfluous locus tree simulation parameters check if the locus tree is fixed.</dd></dl>
        if(*newick_ltree!=NULL || *n_iltrees>0)
        {
            if ((is_sampling_set(*b_rate) || is_sampling_set(*d_rate) || is_sampling_set(*t_rate) || is_sampling_set(*gc_rate) || *min_lleaves>1|| *min_lsleaves>1))
            {
                fprintf(stderr,"\n\tWARNING!!! Using locus tree birth-death parameters with a fixed locus tree has no sense, and these parameters (birth, death and minimum number of leaves) will be ignored\n");
                ParseSampling("f0", b_rate,sampling_vars);
                ParseSampling("f0",d_rate,sampling_vars);
                ParseSampling("f0", t_rate,sampling_vars);
                ParseSampling("f0",gc_rate,sampling_vars);
                *min_lleaves=0;
                *min_lsleaves=0;
            }
            if (is_sampling_set(*alpha_s))
            {
                fprintf(stderr,"\n\tWARNING!!! Using the lineage specific substitution rate heterogeneity alpha with a fixed locus tree has no sense, and it will be ignored\n");
                ParseSampling("f0", alpha_s,sampling_vars);
            }
            if (recon!=0)
            {
                fprintf(stderr,"\n\tWARNING!!! Reconciliation outputs will not be generated as there is no species tree when a fixed locus tree is provided\n");
                recon=0;
            }
            fprintf(stderr,"\n\tWARNING!!! Using a fixed user–specified locus tree is an advanced option, and it should not be used by general users. Some locus–tree–specific parameters may not be checked and could induce biologically senseless simulation scenarios. Examples: Different generation time for paralogs in the same species tree branch.\n");
        }
        
    }
    // ****
    /// Check of species tree simulation parameters</dd></dl></dd></dl>
    else if (*newick_stree==NULL && *n_istrees==0)
    {
        
        if(is_sampling_set(*bds_length)&& *min_lleaves>2)
        {
            fprintf(stderr,"\n\tWARNING!!! Using the minimum number of locus trees (-Ll) and a time-guided simulation of the species tree (-St) could produce incompatible values (minimum locus tree leaves > species tree leaves). If this happens, the number of minimum locus tree leaves will be changed to the number of species tree leaves.\n");
        }
        
    }
    
    if (is_sampling_set(*gc_rate) && is_sampling_set(*b_rate)==0)
    {
        fprintf(stderr,"\n\tWARNING!!! Using locus tree gene conversion simulation with no paralogs (without births) has no sense, and this parameter will be ignored\n");
        ParseSampling("f0", gc_rate,sampling_vars);
    }
    
    //Provisional

    if (*n_istrees>0 && *ns_trees>1)
    {
        fprintf(stderr,"\n\tWARNING!!! Setting both a number of species trees and a species tree input file is not allowed, and therefore the number of replicates will be set equal to the number of species trees present in the file\n");
        *ns_trees=*n_istrees;
    }
    if (*n_iltrees>0 && nl_trees->value.i>0)
    {
        fprintf(stderr,"\n\tWARNING!!! Setting both a number of locus trees and a locus tree input file is not allowed, and therefore the number of replicates will be set equal to the number of locus trees present in the file\n");
        nl_trees->distribution_code=0;
        nl_trees->params[0].i=*n_iltrees;
        nl_trees->params_type[0]=UI;
        nl_trees->value.i=*n_iltrees;
        nl_trees->vtype=0;
    }
    //
    
    fflush(stdout);
    fflush(stderr);
    return(NO_ERROR);
}

/**
 * Gets settings from command line.
 *
 * \param file
 *   Opened file pointer.
 * \param ns_trees
 *   Number of species trees to simulate.
 * \param nl_trees
 *   Number of locus trees to simulate.
 * \param ng_trees
 *   Number of gene trees from each locus tree to simulate.
 * \param newick_stree
 *   Species tree string.
 * \param stree_ifile
 *   Species tree nexus file.
 * \param n_istrees
 *   Number of input (fixed) species trees present in stree_ifile.
 * \param stree_iname
 *   Species tree nexus file name.
 * \param gen_time
 *   Generation time.
 * \param newick_ltree
 *   Locus tree string.
 * \param ltree_ifile
 *   Locus tree nexus file.
 * \param n_iltrees
 *   Number of input (fixed) locus trees present in ltree_ifile.
 * \param stree_iname
 *   Locus tree nexus file name.
 * \param b_rate
 *   Birth prior rate for the locus_tree simulation.
 * \param d_rate
 *   Death prior rate for the locus_tree simulation.
 * \param t_rate
 *   Transfer prior rate for the locus_tree simulation.
 * \param gc_rate
 *   Gene conversion prior rate for the locus_tree simulation.
 * \param lb_rate
 *   Birth rate for the locus_tree simulation.
 * \param ld_rate
 *   Death rate for the locus_tree simulation.
 * \param lt_rate
 *   Transfer rate for the locus_tree simulation.
 * \param lgc_rate
 *   Gene conversion rate for the locus_tree simulation.
 * \param t_kind
 *   Logical flag. 0=> RTRFR randomly sampled. 1=> RTRFR sampled with probability inversely related to distance (generations).
 * \param min_lleaves
 *   Minimum number of leaves for each locus tree.
 * \param min_lsleaves
 *   Minimum number of leaves from different species for each locus tree.
 * \param ind_per_sp
 *   Number of individuals per species.
 * \param sb_rate
 *   Birth rate for the species_tree simulation.
 * \param sd_rate
 *   Death rate for the species_tree simulation.
 * \param bds_leaves
 *   Max number of nodes in the b-d simulation of the species tree.
 * \param bds_length
 *   Max lenght in the b-d simulation of the species tree.
 * \param Outgroup
 *   Internal branch length deviation from a half of the ingroup height,for outgroup addition.
 * \param Ne
 *   Efective population size.
 * \param mu
 *   Substitution rate.
 * \param alpha_s
 *   Alpha parameter of the gamma distribution with mean=1 to sample rate heterogeneity mulipliers. Lineage specific rate heterogeneity.
 * \param alpha_l
 *   Alpha parameter of the gamma distribution with mean=1 to sample rate heterogeneity mulipliers. Gene family specific rate heterogeneity.
 * \param alpha_g
 *   Alpha parameter of the gamma distribution with mean=1 to sample rate heterogeneity mulipliers. Gene tree branch specific rate heterogeneity.
 * \param lalpha_g
 *   Alpha_g prior. Gene tree branch specific rate heterogeneity prior.
 * \param salpha_g
 *   Alpha_l prior (alpha g hyperprior). Gene tree branch specific rate heterogeneity hyperprior.
 * \param epsilon
 *  Epsilon of the convergence of the brent method for sampling Bounded multispecies coalescent
 * \param verbosity
 *   Config about verbosity.
 * \param out_name
 *   Name of the file with the gene trees in Newick format.
 * \param stats
 *   Logical flag, if ==1 it generates a csv with some stats.
 * \param recon
 *   Logical flag, if ==1 it generates a reconciliation between species and gene trees.
 * \return Error-control code.
 * \param db
 *   Logical flag, if ==1 it generates a SQLite database with 3 linked tables (Species_Trees, Locus_Trees and Gene_Trees) with some characteristics of each tree.
 * \param params
 *  Logical flag to activate an output file (prefix.params) with the general parameters of the simulation.
 * \param commands
 *  Logical flag to activate the output of both the original command-line arguments (prefix.command) and the input configuration file (prefix.conf)
 * \param u_seed
 *  Seed for the random number generator.
 
 * \return Error-control code.
 *******************************************************************************/

long int GetSettingsFromFile(FILE *input_file,int *ns_trees, sampling_unit *nl_trees, int *ng_trees, char ** newick_stree, FILE ** stree_ifile, int *n_istrees,char **stree_iname,sampling_unit * gen_time, char ** newick_ltree, FILE **ltree_ifile, int *n_iltrees,char **ltree_iname, sampling_unit * b_rate, sampling_unit * d_rate, sampling_unit * t_rate, sampling_unit *gc_rate,sampling_unit * lb_rate, sampling_unit * ld_rate, sampling_unit * lt_rate, sampling_unit *lgc_rate, int *t_kind, int *min_lleaves, int *min_lsleaves, sampling_unit *ind_per_sp,sampling_unit * sb_rate, sampling_unit * sd_rate, sampling_unit * bds_leaves, sampling_unit * bds_length, sampling_unit * outgroup, sampling_unit *Ne,sampling_unit *mu,sampling_unit *alpha_s, sampling_unit *alpha_l, sampling_unit *salpha_g, sampling_unit *lalpha_g,sampling_unit *alpha_g,float *epsilon, int *verbosity,char **out_name, int *stats, int *recon, int *db, int *params, int *commands, float *u_seed, const sampling_table sampling_vars)
{
    int i=0,offset=0;
    unsigned int LENGTH=10000;
    char code=0,sub_code=0, *buffer, test_char=7; //BEL
    char arg[4];
    
    // ******
    /// <dl><dt> Function structure </dt><dd>
    
    buffer=calloc(LENGTH,sizeof(char));
    *(buffer+LENGTH-2)=test_char;
    
    // *****
    /// <dl><dt> Line's loop </dt><dd>
    while (fgets(buffer,LENGTH,input_file)!=NULL && i<=MAX_IT)
    {
        if (*buffer=='#' || (*buffer=='/' && *buffer=='/'))
            continue;
        else if (*buffer!='-' || ferror(input_file)!=0)
        {
            fprintf(stderr,"Error in the parameter: %s\n",buffer);
            return (SETTINGS_ERROR);
        }
        else if (*(buffer+LENGTH-2)!=test_char) // It reached the buffer's limit
        {
            
#ifdef DBG
            fprintf(stderr,"\n\tWARNING, the buffer used reading the parameters file is not enough. Using a bigger size\n");
            fflush(stderr);
#endif
            LENGTH*=10;
            buffer=realloc(buffer, LENGTH*sizeof(char));
            *(buffer+LENGTH-2)=test_char;
            rewind(input_file);
            continue;
        }
        // ****
        /// <dl><dt>Getting the value depending on the 2letter code (-XX)</dt><dd>
        else
        {
            if (sscanf(buffer,"-%3s%n",arg,&offset)!=1)
            {
                fprintf(stderr,"Error in the parameter: %s\n",buffer);
                return (SETTINGS_ERROR);
            }
            else
            {
                code=toupper(arg[0]);
                sub_code=toupper(arg[1]);
                
                switch (code) //If it is a lowercase converts it in the uppercase.
                {
                        // ***
                        /// <dl><dt>-Rx. Replicates</dt><dd>
                    case 'R':
                        switch (sub_code)
                    {
                            // **
                            /// -RS. Number of study replicates (species trees replicates).
                        case 'S':
                            if(sscanf(buffer+offset,"%u",ns_trees)==0 || *ns_trees<1)
                            {
                                
                                fprintf(stderr,"Error in the parameter: %s\n",buffer);
                                return(SETTINGS_ERROR);
                            }
                            break;
                            
                            // **
                            /// -RG. Number of gene trees for each locus tree.
                        case 'G':
                            if(sscanf(buffer+offset,"%u",ng_trees)==0 || *ng_trees<1)
                            {
                                
                                fprintf(stderr,"Error in the parameter: %s\n",buffer);
                                return(SETTINGS_ERROR);
                            }
                            break;
                            
                            // **
                            /// -RL. Number of locus trees</dd></dl>
                        case 'L':
                            if(ParseSampling(buffer+offset+1,nl_trees,sampling_vars)!=NO_ERROR)
                            {
                                
                                fprintf(stderr,"Error in the parameter: %s\n",buffer);
                                return SETTINGS_ERROR;
                            }
                            break;
                            
                        default:
                            
                            fprintf(stderr,"Unrecognized parameter: %s\n",buffer);
                            return(SETTINGS_ERROR);
                            break;
                    }
                        break;
                        // ***
                        /// <dl><dt>-Hx. Rate heterogeneity parameters</dt><dd>
                    case 'H':
                        switch (sub_code)
                    {
                            // **
                            /// -HS. Alpha parameter for the heterogeneity gamma sampling of lineage specific substitution rates.
                        case 'S':
                            if(ParseSampling(buffer+offset+1,alpha_s,sampling_vars)!=NO_ERROR)
                            {
                                
                                fprintf(stderr,"Error in the parameter: %s\n",buffer);
                                return SETTINGS_ERROR;
                            }
                            break;
                            // **
                            /// -HB. Reserved, it is not implemented yet
                            // **
                            /// -HL. Alpha parameter for the heterogeneity gamma sampling of gene family (gene tree) specific substitution rates
                        case 'L':
                            if(ParseSampling(buffer+offset+1,alpha_l,sampling_vars)!=NO_ERROR)
                            {
                                
                                fprintf(stderr,"Error in the parameter: %s\n",buffer);
                                return SETTINGS_ERROR;
                            }
                            break;
                            // **
                            /// -HG. Alpha parameter for the heterogeneity gamma sampling of gene family branch specific substitution rates</dd></dl>
                        case 'G':
                            if(ParseSampling(buffer+offset+1,alpha_g,sampling_vars)!=NO_ERROR)
                            {
                                
                                fprintf(stderr,"Error in the parameter: %s\n",buffer);
                                return SETTINGS_ERROR;
                            }
                            break;
                        default:
                            
                            fprintf(stderr,"Unrecognized parameter: %s\n",buffer);
                            return(SETTINGS_ERROR);
                            break;
                    }
                        break;
                        // ***
                        /// -P. Haploid population size
                    case 'P':
                        if(ParseSampling(buffer+offset+1,Ne,sampling_vars)!=NO_ERROR)
                        {
                            
                            fprintf(stderr,"Error in the parameter: %s\n",buffer);
                            return SETTINGS_ERROR;
                        }
                        break;
                        // ***
                        /// -U. Global substitution rate
                    case 'U':
                        if(ParseSampling(buffer+offset+1,mu,sampling_vars)!=NO_ERROR)
                        {
                            
                            fprintf(stderr,"Error in the parameter: %s\n",buffer);
                            return SETTINGS_ERROR;
                        }
                        break;
                        // ***
                        /// <dl><dt>-Ex. Epsilons for the bounded multispecies coalescent sampling </dt><dd>
                    case 'E':
                        switch (sub_code)
                    {
                            // **
                            /// -E. Brent method </dd></dt>
                        case '\0':
                            if(sscanf(buffer+offset,"%f",epsilon)==0 || *epsilon<0)
                            {
                                
                                fprintf(stderr,"Error in the parameter: %s\n",buffer);
                                return (SETTINGS_ERROR);
                            }
                            break;
                        default:
                            
                            fprintf(stderr,"Unrecognized parameter: %s\n",buffer);
                            return (SETTINGS_ERROR);
                            break;
                    }
                        break;
                        // ***
                        /// -G. Generation time (branch lengths as time instead of generations, input and output trees, inside the software, the trees are in generations)
                    case 'G':
                        if(ParseSampling(buffer+offset+1,gen_time,sampling_vars)!=NO_ERROR)
                        {
                            
                            fprintf(stderr,"Error in the parameter: %s\n",buffer);
                            return SETTINGS_ERROR;
                        }
                        break;
                        // ***
                        /// -V. Verbosity
                    case 'V':
                        if(sscanf(buffer+offset,"%u",verbosity)==0)
                        {
                            
                            fprintf(stderr,"Error in the parameter: %s\n",buffer);
                            return (SETTINGS_ERROR);
                        }
                        break;
                        // ***
                        /// <dl><dt>-Ox. Output related options</dt><dd>
                    case 'O':
                        switch (sub_code)
                    {
                            // **
                            /// -O. Output common prefix
                        case '\0':
                            *out_name=calloc(strlen(buffer+offset)+1,sizeof(char));
                            sscanf(buffer+offset,"%s",*out_name);
                            break;
                            // **
                            /// -OS. Statistics output
                        case 'S':
                            if(sscanf(buffer+offset,"%u",stats)==0)
                            {
                                
                                fprintf(stderr,"Error in the parameter: %s\n",buffer);
                                return (SETTINGS_ERROR);
                            }
                            break;
                            // **
                            /// -OR. Reconciliation output
                        case 'R':
                            if(sscanf(buffer+offset,"%u",recon)==0)
                            {
                                
                                fprintf(stderr,"Error in the parameter: %s\n",buffer);
                                return (SETTINGS_ERROR);
                            }
                            break;
                            // **
                            /// -OD. SQLite relational database output
                        case 'D':
                            if(sscanf(buffer+offset,"%u",db)==0)
                            {
                                
                                fprintf(stderr,"Error in the parameter: %s\n",buffer);
                                return (SETTINGS_ERROR);
                            }
                            break;
                            // **
                            /// -OP. Params output
                        case 'P':
                            if(sscanf(buffer+offset,"%u",params)==0)
                            {
                                
                                fprintf(stderr,"Error in the parameter: %s\n",buffer);
                                return (SETTINGS_ERROR);
                            }
                            break;
                            // **
                            /// -OC. Commands output </dd></dl>
                        case 'C':
                            if(sscanf(buffer+offset,"%u",commands)==0)
                            {
                                
                                fprintf(stderr,"Error in the parameter: %s\n",buffer);
                                return (SETTINGS_ERROR);
                            }
                            break;
                        default:
                            
                            fprintf(stderr,"Unrecognized parameter: %s\n",buffer);
                            return (SETTINGS_ERROR);
                            break;
                    }
                        break;
                        // ***
                        /// <dl><dt>-Sx. Species tree</dt><dd>
                    case 'S':
                        switch (sub_code)
                    {
                            // **
                            /// -S. Fixed species tree. Parses the string as a species tree
                        case '\0':
                            *newick_stree=calloc(strlen(buffer+offset)+1,sizeof(char));
                            sscanf(buffer+offset,"%s",*newick_stree);
                            break;
                            // **
                            /// -SR. Fixed species trees. Parses the string as the Nexus file containing the species tree/s
                        case 'R':
                            *stree_iname=calloc(strlen(buffer+offset)+1,sizeof(char));
                            sscanf(buffer+offset,"%s",*stree_iname);
                            if ((*stree_ifile=fopen(*stree_iname, "r"))==NULL)
                            {
                                fprintf(stderr,"Error in the parameter: %s\n",buffer);
                                perror("Error opening input_file: ");
                                ErrorReporter(IO_ERROR,NULL);
                            }
                            else
                            {
                                ErrorReporter(NNexusTrees(*stree_ifile,n_istrees),NULL);
                            }
                            break;
                            // **
                            /// <dl><dt>-Sx. Birth-death species tree</dt><dd>
                            // *
                            /// -SB. Birth rate
                        case 'B':
                            if(ParseSampling(buffer+offset+1,sb_rate,sampling_vars)!=NO_ERROR)
                            {
                                
                                fprintf(stderr,"Error in the parameter: %s\n",buffer);
                                return SETTINGS_ERROR;
                            }
                            break;
                            // *
                            /// -SD. Death rate
                        case 'D':
                            if(ParseSampling(buffer+offset+1,sd_rate,sampling_vars)!=NO_ERROR)
                            {
                                
                                fprintf(stderr,"Error in the parameter: %s\n",buffer);
                                return SETTINGS_ERROR;
                            }
                            break;
                            // *
                            /// -ST. Maximum time to stop the birth-death process.
                        case 'T':
                            if(ParseSampling(buffer+offset+1,bds_length,sampling_vars)!=NO_ERROR)
                            {
                                
                                fprintf(stderr,"Error in the parameter: %s\n",buffer);
                                return SETTINGS_ERROR;
                            }
                            break;
                            // *
                            /// -SL. Desired number of leaves to stop the birth-death process.
                        case 'L':
                            if(ParseSampling(buffer+offset+1,bds_leaves,sampling_vars)!=NO_ERROR)
                            {
                                
                                fprintf(stderr,"Error in the parameter: %s\n",buffer);
                                return SETTINGS_ERROR;
                            }
                            break;
                            // *
                            /// -SO. Deviation of a half of the ingroup heigth for the new internal branch due to the outgroup addition
                        case 'O':
                            if(ParseSampling(buffer+offset+1,outgroup,sampling_vars)!=NO_ERROR)
                            {
                                
                                fprintf(stderr,"Error in the parameter: %s\n",buffer);
                                return SETTINGS_ERROR;
                            }
                            break;
                            // ***
                            /// -SI. Number of individuals per species</dd></dl></dd></dl>
                        case 'I':
                            if(ParseSampling(buffer+offset+1,ind_per_sp,sampling_vars)!=NO_ERROR)
                            {
                                
                                fprintf(stderr,"Error in the parameter: %s\n",buffer);
                                return SETTINGS_ERROR;
                            }
                            break;
                            
                        default:
                            
                            fprintf(stderr,"Unrecognized parameter: %s\n",buffer);
                            return (SETTINGS_ERROR);
                            break;
                    }
                        break;
                        // ***
                        /// <dl><dt>-Lx. Locus tree.</dt><dd>
                    case 'L':
                        switch (sub_code)
                    {
                            // **
                            /// -L. Fixed locus tree. Parses the string as a locus tree
                        case'\0':
                            *newick_ltree=calloc(strlen(buffer+offset)+1,sizeof(char));
                            sscanf(buffer+offset,"%s",*newick_ltree);
                            break;
                            // **
                            /// -LR. Fixed locus trees. Parses the string as the Nexus file containing the locus tree/s
                        case 'R':
                            *ltree_iname=calloc(strlen(buffer+offset)+1,sizeof(char));
                            sscanf(buffer+offset,"%s",*ltree_iname);
                            if ((*ltree_ifile=fopen(*ltree_iname, "r"))==NULL)
                            {
                                fprintf(stderr,"Error in the parameter: %s\n",buffer);
                                perror("Error opening input_file: ");
                                ErrorReporter(IO_ERROR,NULL);
                            }
                            else
                            {
                                ErrorReporter(NNexusTrees(*ltree_ifile,n_iltrees),NULL);
                            }
                            break;
                            // **
                            /// -LB. Simulated locus tree. Birth rate.
                        case 'B':
                            if(ParseSampling(buffer+offset+1,b_rate,sampling_vars)!=NO_ERROR)
                            {
                                
                                fprintf(stderr,"Error in the parameter: %s\n",buffer);
                                return SETTINGS_ERROR;
                            }
                            break;
                            // **
                            /// -LD. Simulated locus tree. Death rate.
                        case 'D':
                            if(ParseSampling(buffer+offset+1,d_rate,sampling_vars)!=NO_ERROR)
                            {
                                
                                fprintf(stderr,"Error in the parameter: %s\n",buffer);
                                return SETTINGS_ERROR;
                            }
                            break;
                            // **
                            /// -LT. Simulated locus tree. Transfer rate.
                        case 'T':
                            if(ParseSampling(buffer+offset+1,t_rate,sampling_vars)!=NO_ERROR)
                            {
                                
                                fprintf(stderr,"Error in the parameter: %s\n",buffer);
                                return SETTINGS_ERROR;
                            }
                            break;
                            // **
                            /// -LG. Simulated locus tree. Gene conversion rate.
                        case 'G':
                            if(ParseSampling(buffer+offset+1,gc_rate,sampling_vars)!=NO_ERROR)
                            {
                                
                                fprintf(stderr,"Error in the parameter: %s\n",buffer);
                                return SETTINGS_ERROR;
                            }
                            break;
                            // **
                            /// -LK. Simulated locus tree. Aceptor sampling dependent (1) or independent (0) to the evolutionary distance.
                        case 'K':
                            if (sscanf(buffer+offset,"%u",t_kind)==0 || *t_kind<0 || *t_kind>1)
                            {
                                
                                fprintf(stderr,"Error in the parameter: %s\n",buffer);
                                return (SETTINGS_ERROR);
                            }
                            break;
                            // **
                            /// -LL. Simulated locus tree. Minimum number of leaves of simulated locus trees.
                        case 'L':
                            if(sscanf(buffer+offset,"%u",min_lleaves)==0 || *min_lleaves<2)
                            {
                                
                                fprintf(stderr,"Error in the parameter: %s\n",buffer);
                                return (SETTINGS_ERROR);
                            }
                            break;
                            // **
                            /// -LS. Simulated locus tree. Minimum number of leaves from different species of simulated locus trees.</dd></dl>
                        case 'S':
                            if(sscanf(buffer+offset,"%u",min_lsleaves)==0 || *min_lsleaves<2)
                            {
                                
                                fprintf(stderr,"Error in the parameter: %s\n",buffer);
                                return (SETTINGS_ERROR);
                            }
                            break;
                        default:
                            
                            fprintf(stderr,"Unrecognized parameter: %s\n",buffer);
                            return (SETTINGS_ERROR);
                            break;
                    }
                        break;
                        // ***
                        /// <dl><dt>-Cx. Global options.</dt><dd>
                    case 'C':
                        switch (sub_code)
                    {       // **
                            /// -Cs. Seed for the random number generator.</dd></dl></dd></dl>
                        case 'S':
                        {
                            if(sscanf(buffer+offset,"%f",u_seed)==0)
                            {
                                
                                fprintf(stderr,"Error in the parameter: %s\n",buffer);
                                return (SETTINGS_ERROR);
                            }
                            
                        }
                            break;
                        default:
                            
                            fprintf(stderr,"Unrecognized parameter: %s\n",buffer);
                            return (SETTINGS_ERROR);
                            break;
                            
                            
                    }
                        break;
                        
                    default:
                        
                        fprintf(stderr,"Unrecognized parameter: %s\n",buffer);
                        return (SETTINGS_ERROR);
                        break;
                }
            }
        }
        
        *(buffer+LENGTH-2)=test_char;
        ++i;
    }
    free(buffer);
    fflush(stdout);
    fflush(stderr);
    return(NO_ERROR);
}


/**
 * Writes the config of the program in a file (or standard outputs).
 *
 * \param file
 *  Output file
 * \param s_tree_newick
 *  Newick string of a fixed species tree.
 * \param stree_ifile
 *  Input stree file name/route
 * \param n_istrees
 *  Number of input trees present in the stree_ifile.
 * \param n_strees
 *  Number of species tree simulations.
 * \param gen_time
 *  Generation time
 * \param l_tree_newick
 *  Newick string of a fixed locus tree.
 * \param ltree_ifile
 *  Input ltree file name/route
 * \param n_iltrees
 *  Number of input trees present in the ltree_ifile.
 * \param sb_rate
 *  Birth rate of the species tree simulation.
 * \param sd_rate
 *  Death rate of the species tree simulation.
 * \param s_leaves
 *  Number of leaves of the desired species tree (stopping rule).
 * \param s_time
 *  Maximum length of the simulated species tree.
 * \param Outgroup
 *  Internal branch length deviation from a half of the ingroup height,for outgroup addition.
 * \param b_rate
 *  Birth rate hyperparameter (species tree related) of the locus tree simulation.
 * \param d_rate
 *  Death rate hyperparameter (species tree related) of the locus tree simulation.
 * \param t_rate
 *  Transfer rate hyperparameter (species tree related) of the locus tree simulation.
 * \param gc_rate
 *  Gene conversion rate hyperparameter (species tree related) of the locus tree simulation.
 * \param lb_rate
 *  Birth rate of the locus tree simulation.
 * \param ld_rate
 *  Death rate of the locus tree simulation.
 * \param lt_rate
 *  Transfer rate of the locus tree simulation.
 * \param lgc_rate
 *  Gene conversion rate of the locus tree simulation.
 * \param t_kind
 *  Gene transferences (conversion and transfer) depending (1) or not (0) on the distance.
 * \param min_lleaves
 *  Minimum number of leaves for each locus tree.
 *  \param min_lsleaves
 *  Minimum number of tips from different species for each locus tree.
 * \param ind_per_sp
 *  Number of individuals per species.
 * \param nl_trees
 *  Number of locus trees to simulate.
 * \param ng_trees
 *  Number of gene trees to simulate from each locus tree.
 * \param Ne
 *  Efective population size.
 * \param mu
 *  Substitution rate.
 * \param alpha_s
 *  Alpha parameter of the gamma distribution with mean=1 to sample rate heterogeneity mulipliers. Lineage specific rate heterogeneity.
 * \param alpha_l
 *  Alpha parameter of the gamma distribution with mean=1 to sample rate heterogeneity mulipliers. Gene family specific rate heterogeneity.
 * \param alpha_g
 *  Alpha parameter of the gamma distribution with mean=1 to sample rate heterogeneity mulipliers. Gene tree branch specific rate heterogeneity.
 * \param lalpha_g
 *  Hyperparemeter (locus-tree related) for the Alpha parameter of the gamma distribution with mean=1 to sample rate heterogeneity mulipliers. Gene tree branch specific rate heterogeneity.
 * \param salpha_g
 *  Hyperhyperparameter (species-tree related) for the Alpha parameter of the gamma distribution with mean=1 to sample rate heterogeneity mulipliers. Gene tree branch specific rate heterogeneity.
 * \param epsilon
 *  Epsilon of the convergence of the brent method for sampling Bounded multispecies coalescent
 * \param verbosity
 *  Config about verbosity.
 * \param out_file
 *  Prefix of the output files.
 * \param u_seed
 *  Seed for the random number generator.
 *******************************************************************************/

static void PrintGlobalSettings(FILE *file,char * s_tree_newick, char *stree_ifile, int n_istrees, int n_strees,sampling_unit *gen_time, char * l_tree_newick, char *ltree_ifile, int n_iltrees, sampling_unit *sb_rate, sampling_unit *sd_rate, sampling_unit *s_leaves, sampling_unit *s_time, sampling_unit *outgroup, sampling_unit *b_rate, sampling_unit *d_rate, sampling_unit *t_rate, sampling_unit *gc_rate,sampling_unit *lb_rate, sampling_unit *ld_rate, sampling_unit *lt_rate, sampling_unit *lgc_rate, int t_kind, int min_lleaves, int min_lsleaves, sampling_unit *ind_per_sp,sampling_unit *nl_trees,int ng_trees,sampling_unit *Ne,sampling_unit *mu,sampling_unit *alpha_s, sampling_unit *alpha_l, sampling_unit *alpha_g, sampling_unit *lalpha_g, sampling_unit *salpha_g,float epsilon, int verbosity, char * out_file,float u_seed,int stats,int recon,int db, int params, int command, const sampling_table stable)
{
    char * buffer=calloc(100,sizeof(char));
    
    fprintf(file,"\nGlobal settings:\n--------------------\n\nTrees:");
#ifdef DBG
    fflush(stdout);
#endif
    
    if (l_tree_newick!=NULL || n_iltrees>0)
    {
        fprintf(file,"\n\t-Species trees: Not taken into consideration\n");
        if (ltree_ifile!=NULL)
            fprintf(file,"\t-Locus trees: %d fixed trees found in %s",n_iltrees,ltree_ifile);
        else
            fprintf(file,"\t-Locus tree: Fixed\n\t\t--> %s",l_tree_newick);
#ifdef DBG
        fflush(stdout);
#endif

    }
    else
    {
        if (s_tree_newick!=NULL || n_istrees>0)
        {
            if (stree_ifile!=NULL)
                fprintf(file,"\n\t-Species trees: %d fixed trees found in %s",n_istrees,stree_ifile);
            else
                fprintf(file,"\n\t-Species tree: Fixed\n\t\t--> %s",s_tree_newick);
#ifdef DBG
            fflush(stdout);
#endif
        }
        else
        {
            Print_Sampling(sb_rate,buffer,stable);
            fprintf(file,"\n\t-Species trees: %d Birth-death simulations\n\t\t-Birth rate: %s",n_strees,buffer);
            Print_Sampling(sd_rate,buffer,stable);
            fprintf(file,"\n\t\t-Death rate: %s",buffer);
            
            if (get_sampling((*outgroup))>0||outgroup->distribution_code!=0)
            {
                Print_Sampling(outgroup,buffer,stable);
                fprintf(file,"\n\t\t-Outgroup addition:\n\t\t\t-Internal branch length deviation from the half of the height of the ingroup: %s",buffer);
            }
            else
                fprintf(file,"\n\t\t-Outgroup addition: No addition");
            
            if (get_sampling((*s_leaves))>0||s_leaves->distribution_code!=0)
                Print_Sampling(s_leaves,buffer,stable);
            else
                sprintf(buffer,"%s","Not used");
            fprintf(file,"\n\t\t-Stopping rules:\n\t\t\t-Number of leaves: %s",buffer);
            
            if (get_sampling((*s_time))>0||s_time->distribution_code!=0)
                Print_Sampling(s_time,buffer,stable);
            else
                sprintf(buffer,"%s","Not used");
            fprintf(file,"\n\t\t\t-Generations: %s",buffer);
        }
        if ((is_variable(*lb_rate) && is_variable(*ld_rate) && is_variable(*lt_rate) && is_variable(*lgc_rate)) == 0)
        {
            Print_Sampling(nl_trees,buffer,stable);
            fprintf(file,"\n\t-Locus trees: %s directly obtained from each species tree (no birth-death process)\n",buffer);
        }
        else
        {
            Print_Sampling(nl_trees,buffer,stable);
            fprintf(file,"\n\t-Locus trees: %s locus simulations per species tree \n\t\t",buffer);
            Print_Sampling(b_rate,buffer,stable);
            fprintf(file,"-Birth rate: Prior (species-tree related): %s",buffer);
            Print_Sampling(lb_rate,buffer,stable);
            fprintf(file, ", Rate: %s",buffer);
            Print_Sampling(d_rate,buffer,stable);
            fprintf(file,"\n\t\t-Death rate: Prior (species-tree related): %s",buffer);
            Print_Sampling(ld_rate,buffer,stable);
            fprintf(file, ", Rate: %s",buffer);
            Print_Sampling(t_rate,buffer,stable);
            fprintf(file,"\n\t\t-Transference rate: Prior (species-tree related): %s",buffer);
            Print_Sampling(lt_rate,buffer,stable);
            fprintf(file, ", Rate: %s",buffer);
            Print_Sampling(gc_rate,buffer,stable);
            fprintf(file,"\n\t\t-Gene conversion rate: Prior (species-tree related): %s",buffer);
            Print_Sampling(lgc_rate,buffer,stable);
            fprintf(file, ", Rate: %s",buffer);
            fprintf(file,"\n\t\t-Transference/Conversion sampling strategy: %s",t_kind==0?"Random":"Inversely proportional to the distance");
            fprintf(file,"\n\t\t-Minimum number of leaves: %u\n\t\t-Minimum number of leaves from different species: %u\n",min_lleaves,min_lsleaves);
        }
        
    }
    fprintf(file,"\n\t-Gene trees: %d multilocus coalescent simulations\n",ng_trees);
    Print_Sampling(Ne,buffer,stable);
    fprintf(file,"\nParameters:\n\t-Haploid efective population size: %s",buffer);
    Print_Sampling(gen_time,buffer,stable);
    fprintf(file,"\n\t-Generation time: %s",buffer);
    Print_Sampling(mu,buffer,stable);
    fprintf(file,"\n\t-Global substitution rate: %s",buffer);
    Print_Sampling(alpha_s,buffer,stable);
    fprintf(file,"\n\t-Substitution rate heterogeneities\n\t\t-Lineage (species) specific rate heterogeneity gamma shape: %s",(get_sampling((*alpha_s))==0&&is_variable(*alpha_s)==0)?"No heterogeneity":buffer);
    Print_Sampling(alpha_l,buffer,stable);
    fprintf(file,"\n\t\t-Gene family (locus tree) specific rate heterogeneity gamma shape: %s",(get_sampling((*alpha_l))==0&&is_variable(*alpha_l)==0)?"No heterogeneity":buffer);
    Print_Sampling(salpha_g,buffer,stable); ///To check wether it is pure fixed or fixingly depending on another parameter and print no heter just in the first case
    fprintf(file,"\n\t\t-Gene tree branch specific rate heterogeneity gamma shape: HyperPrior(species-tree related) %s",is_sampling_set(*salpha_g)?buffer:"No heterogeneity");
    Print_Sampling(lalpha_g,buffer,stable);
    fprintf(file,", Prior(locus-tree related) %s",(get_sampling((*lalpha_g))!=0||is_variable(*lalpha_g)!=0)?buffer:"No heterogeneity");
    Print_Sampling(alpha_g,buffer,stable);
    fprintf(file,", Gamma parameter(gene-tree related) %s",(get_sampling((*alpha_g))!=0 ||is_variable(*alpha_g)!=0)?buffer:"No heterogeneity");
    Print_Sampling(ind_per_sp,buffer,stable);
    fprintf(file,"\n\t-Individuals per species: %s",buffer);
    Print_Sampling(nl_trees,buffer,stable);
    fprintf(file,"\n\nMisc parameters:\n\t-Rooting method epsilon: %f\n\t-Seed: %f",epsilon,u_seed);
    fprintf(file,"\n\nI/O options:\n\t-Output files prefix: %s\n\t-Verbosity: %d\n\t-Stats file: %s\n\t-Reconciliation: %s\n\t-Database: %s\n\t-Parameterization: %s \n\t-Command-line arguments: %s\n",out_file,verbosity,stats==0?"OFF":"ON",recon==0?"OFF":"ON",db==0?"OFF":"ON",params==0?"OFF":"ON",command==0?"OFF":"ON");
    
    fflush(stdout);
    
    free(buffer);
    
}

/**
 * Writes the config of the program in stdout for each sample (species tree).
 *
 * \param s_tree_newick
 *  Newick string of a fixed species tree.
 * \param gen_time
 *  Generation time
 * \param l_tree_newick
 *  Newick string of a fixed locus tree.
 * \param ltree_ifile
 *  Input ltree file name/route
 * \param n_iltrees
 *  Number of input trees present in the ltree_ifile.
 * \param sb_rate
 *  Birth rate of the species tree simulation.
 * \param sd_rate
 *  Death rate of the species tree simulation.
 * \param s_leaves
 *  Number of leaves of the desired species tree (stopping rule).
 * \param s_time
 *  Maximum length of the simulated species tree.
 * \param lb_rate
 *  Birth rate of the locus tree simulation.
 * \param ld_rate
 *  Death rate of the locus tree simulation.
 * \param lt_rate
 *  Transfer rate of the locus tree simulation.
 * \param lgc_rate
 *  Gene conversion rate of the locus tree simulation.
 * \param t_kind
 *  Gene transferences (conversion and transfer) depending (1) or not (0) on the distance.
 * \param Outgroup
 *  Parameter to add an outgroup with a given/sampled certain length of the new internal branch.
 * \param min_lleaves
 *  Minimum number of leaves for each locus tree.
 *  \param min_lsleaves
 *  Minimum number of tips from different species for each locus tree.
 * \param ind_per_sp
 *  Number of individuals per species.
 * \param nl_trees
 *  Number of locus trees to simulate.
 * \param ng_trees
 *  Number of gene trees to simulate from each locus tree.
 * \param Ne
 *  Efective population size.
 * \param mu
 *  Substitution rate.
 * \param alpha_s
 *  Alpha parameter of the gamma distribution with mean=1 to sample rate heterogeneity mulipliers. Lineage specific rate heterogeneity.
 * \param alpha_l
 *  Alpha parameter of the gamma distribution with mean=1 to sample rate heterogeneity mulipliers. Gene family specific rate heterogeneity.
 * \param alpha_g
 *  Alpha parameter of the gamma distribution with mean=1 to sample rate heterogeneity mulipliers. Gene tree branch specific rate heterogeneity.
 * \param lalpha_g
 *  Hyperparemeter (locus-tree related) for the Alpha parameter of the gamma distribution with mean=1 to sample rate heterogeneity mulipliers. Gene tree branch specific rate heterogeneity.
 * \param salpha_g
 *  Hyperhyperparameter (species-tree related) for the Alpha parameter of the gamma distribution with mean=1 to sample rate heterogeneity mulipliers. Gene tree branch specific rate heterogeneity.
 * \param epsilon
 *  Epsilon of the convergence of the brent method for sampling Bounded multispecies coalescent
 * \param verbosity
 *  Config about verbosity.
 * \param out_file
 *  Prefix of the output files.
 * \param n_replicate
 *  Number of this replicate.
 *******************************************************************************/

void PrintSettingsSloop(char * s_tree_newick, sampling_unit gen_time,char *l_tree_newick, char *ltree_ifile, int n_iltrees,sampling_unit sb_rate, sampling_unit sd_rate, sampling_unit s_leaves, sampling_unit s_time,sampling_unit lb_rate, sampling_unit ld_rate, sampling_unit lt_rate, sampling_unit lgc_rate,int t_kind,sampling_unit outgroup, int min_lleaves, int min_lsleaves, sampling_unit ind_per_sp,sampling_unit nl_trees,int ng_trees,sampling_unit Ne,sampling_unit mu,sampling_unit alpha_s, sampling_unit alpha_l, sampling_unit salpha_g,float epsilon, int verbosity, char * out_file, int n_replicate)
{
    char * buffer=calloc(100,sizeof(char));
    
    printf("\n\nReplicate %u parameters:\n------------------------\n\nTrees:",n_replicate);
#ifdef DBG
    fflush(stdout);
#endif
    
    if (l_tree_newick!=NULL || n_iltrees>0)
    {
        printf("\n\t-Species tree: Not taken into consideration\n");
        if (ltree_ifile!=NULL)
            printf("\t-Locus trees: %d fixed trees found in %s",n_iltrees,ltree_ifile);
        else
            printf("\t-Locus tree: Fixed\n\t\t--> %s",l_tree_newick);
#ifdef DBG
        fflush(stdout);
#endif
        if (get_sampling(gen_time)!=1)
        {
            printf("\n\t\t-Time units: time\n\t\t-Generation time %e\n",get_sampling(gen_time));
#ifdef DBG
            fflush(stdout);
#endif
        }
        else
        {
            printf("\n\t\t-Time units: generations\n");
#ifdef DBG
            fflush(stdout);
#endif
        }
    }
    else
    {
        if (s_tree_newick!=NULL)
        {
            printf("\n\t-Species tree: Fixed\n\t\t--> %s",s_tree_newick);
#ifdef DBG
            fflush(stdout);
#endif
            if (get_sampling(gen_time)!=1)
            {
                printf("\n\t\t-Time units: time\n\t\t-Generation time %e\n",get_sampling(gen_time));
#ifdef DBG
                fflush(stdout);
#endif
            }
            else
            {
                printf("\n\t\t-Time units: generations\n");
#ifdef DBG
                fflush(stdout);
#endif
            }
        }
        else
        {
            if (get_sampling(gen_time)!=1)
            {
                printf("\n\t-Species tree: Birth-death simulation\n\t\t-Birth rate (per time unit): %e",get_sampling(sb_rate));
                printf("\n\t\t-Death rate (per time unit): %e",get_sampling(sd_rate));
                if (get_sampling(outgroup)>0||outgroup.distribution_code!=0)
                    printf("\n\t\t-Outgroup addition:\n\t\t\t-Internal branch length deviation from the half of the height of the ingroup: %e",get_sampling(outgroup));
                else
                    printf("\n\t\t-Outgroup addition: No addition");
                
                if (get_sampling(s_leaves)>0||s_leaves.distribution_code!=0)
                    printf("\n\t\t-Stopping rules:\n\t\t\t-Number of leaves: %lf",get_sampling(s_leaves));
                else
                {
                    sprintf(buffer,"%s","Not used");
                    printf("\n\t\t-Stopping rules:\n\t\t\t-Number of leaves: %s",buffer);
                }
                
                
                if (get_sampling(s_time)>0||s_time.distribution_code!=0)
                    printf("\n\t\t\t-Time: %e",get_sampling(s_time));
                else
                {
                    sprintf(buffer,"%s","Not used");
                    printf("\n\t\t\t-Time: %s",buffer);
                }
                
                printf("\n\t\t-Time units: time\n\t\t-Generation time %e\n",get_sampling(gen_time));
                printf("\t\tNOTE: The intensity of the b-d process corresponds to the time units instead of the number of generations\n");
            }
            else
            {
                printf("\n\t-Species tree: Birth-death simulation\n\t\t-Birth rate: %e",get_sampling(sb_rate));
                printf("\n\t\t-Death rate: %e",get_sampling(sd_rate));
                if (get_sampling(outgroup)>0||outgroup.distribution_code!=0)
                    printf("\n\t\t-Outgroup addition:\n\t\t\t-Internal branch length deviation from the half of the height of the ingroup: %e",get_sampling(outgroup));
                else
                    printf("\n\t\t-Outgroup addition: No addition");
                
                if (get_sampling(s_leaves)>0||s_leaves.distribution_code!=0)
                    printf("\n\t\t-Stopping rules:\n\t\t\t-Number of leaves: %lf",get_sampling(s_leaves));
                else
                {
                    sprintf(buffer,"%s","Not used");
                    printf("\n\t\t-Stopping rules:\n\t\t\t-Number of leaves: %s",buffer);
                    
                }
                
                if (get_sampling(s_time)>0||s_time.distribution_code!=0)
                    printf("\n\t\t\t-Generations: %e",get_sampling(s_time));
                else
                {
                    sprintf(buffer,"%s","Not used");
                    printf("\n\t\t\t-Generations: %s",buffer);
                }
                printf("\n\t\t-Time units: generations\n");
            }
        }
        if (get_sampling(lb_rate)==0 && get_sampling(ld_rate)==0 && get_sampling(lt_rate)==0 && get_sampling(lgc_rate)==0)
        {
            printf("\n\t-%f Locus trees: directly obtained from the species tree (no birth-death process)\n",get_sampling(nl_trees));
        }
        else
        {
            printf("\n\t-%f Locus trees: birth-death simulation\n\t\t-Birth rate prior (per generation): %e",get_sampling(nl_trees),get_sampling(lb_rate));
            printf("\n\t\t-Death rate prior (per generation): %e",get_sampling(ld_rate));
            printf("\n\t\t-Transference rate prior (per generation): %e",get_sampling(lt_rate));
            printf("\n\t\t-Gene conversion rate prior (per generation): %e",get_sampling(lgc_rate));
            printf("\n\t\t-Minimum number of leaves: %u\n\t\t-Minimum number of leaves from different species: %u\n",min_lleaves,min_lsleaves);
        }
    }
    printf("\n\t-Gene trees: %d multilocus coalescent simulations\n",ng_trees);
    printf("\nParameters:\n\t-Haploid efective population size: %e",get_sampling(Ne));
    printf("\n\t-Global substitution rate: %e",get_sampling(mu));
    sprintf(buffer, "%e",get_sampling(alpha_s));
    printf("\n\t-Substitution rate heterogeneities\n\t\t-Lineage (species) specific rate heterogeneity gamma shape: %s",get_sampling(alpha_s)==0?"No heterogeneity":buffer);
    sprintf(buffer, "%e",get_sampling(alpha_l));
    printf("\n\t\t-Gene family (locus tree) specific rate heterogeneity gamma shape: %s",get_sampling(alpha_l)==0?"No heterogeneity":buffer);
    sprintf(buffer, "%e",get_sampling(salpha_g));
    printf("\n\t\t-Gene tree branch specific rate heterogeneity gamma shape: HyperPrior %s",get_sampling(salpha_g)==0?"No heterogeneity":buffer);
    printf("\n\t-Individuals per species: %e \n\nSimulation:\n-----------\n",get_sampling(ind_per_sp));
    fflush(stdout);
    
    free(buffer);
    
}


/**
 * Writes the config of the program in stdout for each sample (locus tree).
 *
 * \param curr_ltree
 *  Locus tree id.
 * \param b_rate
 *  Birth rate of the locus tree simulation.
 * \param d_rate
 *  Death rate of the locus tree simulation.
 * \param t_rate
 *  Transfer rate of the locus tree simulation.
 * \param gc_rate
 *  Gene conversion rate of the locus tree simulation.
 * \param gamma_l
 *  Locus tree rate heterogeneity multiplier. Gene family specific rate heterogeneity.
 * \param lalpha_g
 *  Hyperparemeter (locus-tree related) for the Alpha parameter of the gamma distribution to sample rate heterogeneity mulipliers. Gene tree branch specific rate heterogeneity.
 *******************************************************************************/
void PrintSettingsLloop(int curr_ltree,sampling_unit lb_rate,sampling_unit ld_rate,sampling_unit lt_rate,sampling_unit lgc_rate,double gamma_l,sampling_unit lalpha_g)
{
    if (get_sampling(lb_rate)!=0 || get_sampling(ld_rate)!=0 || get_sampling(lt_rate)!=0 || get_sampling(lgc_rate)!=0)
        printf("\n\tLocus tree %d parameters:\n\t-----------------------\n\t-Birth rate (per generation): %e\n\t-Death rate (per generation): %e\n\t-Transference rate (per generation): %e\n\t-Gene conversion rate (per generation): %e\n\t-Gene family rate multiplier %e\n\t-Gene tree branch specific heterogeneity hyperparameter %e\n\n",curr_ltree,get_sampling(lb_rate),get_sampling(ld_rate),get_sampling(lt_rate),get_sampling(lgc_rate),gamma_l,get_sampling(lalpha_g));
    else
        printf("\n\tLocus tree %d parameters:\n\t-----------------------\n\t-Fixed\n\t-Gene family rate multiplier %e\n\t-Gene tree branch specific heterogeneity hyperparameter %e\n\n",curr_ltree,gamma_l,get_sampling(lalpha_g));
    
}

/**
 * Writes the config of the program in stdout for each sample (gene tree).
 *
 * \param curr_gtree
 *  Gene tree id.
 * \param alpha_g
 *  Paremeter (gene-tree related) for the Alpha parameter of the gamma distribution to sample rate heterogeneity mulipliers. Gene tree branch specific rate heterogeneity.
 *******************************************************************************/
void PrintSettingsGloop(int curr_gtree,sampling_unit alpha_g)
{
    printf("\n\t\tGene tree %d parameter:\n\t\t--------------------\n\t\t-Gene tree branch specific heterogeneity parameter %e\n",curr_gtree,get_sampling(alpha_g));
}

long int CheckSampledSettingsSloop(sampling_unit bds_leaves, sampling_unit bds_length, sampling_unit sb_rate, sampling_unit sd_rate, sampling_unit outgroup, sampling_unit ind_per_sp, sampling_unit nl_trees, sampling_unit b_rate, sampling_unit d_rate, sampling_unit t_rate, sampling_unit gc_rate,sampling_unit Ne, sampling_unit alpha_s, sampling_unit alpha_l, sampling_unit salpha_g, sampling_unit mu, sampling_unit gen_time, int min_lleaves)
{
    int is_error=0;
    
    if (is_sampling_set(bds_leaves)&& (get_sampling(bds_leaves)<=0))
    {
        fprintf(stderr,"\n\tImproper value sampling the parameter -Sl, Number of species tree leaves. Please, check your sampling settings and try again\n");
        is_error=1;
    }
    if (is_sampling_set(bds_length)&& (get_sampling(bds_length)<=1))
    {
        fprintf(stderr,"\n\tImproper value sampling the parameter -St, Species tree length. Please, check your sampling settings and try again\n");
        is_error=1;
    }
    if (is_sampling_set(sb_rate)&& (get_sampling(sb_rate)<=0))
    {
        fprintf(stderr,"\n\tImproper value sampling the parameter -Sb, Species tree birth rate. Please, check your sampling settings and try again\n");
        is_error=1;
    }
    if (is_sampling_set(sd_rate)&&(get_sampling(sd_rate)>get_sampling(sb_rate)))
    {
        fprintf(stderr,"\n\tImproper value sampling the parameter -Sd, Species tree death rate. Please, check your sampling settings and try again\n");
        is_error=1;
    }
    if (is_sampling_set(outgroup)&&(get_sampling(outgroup)<=0))
    {
        fprintf(stderr,"\n\tImproper value sampling the parameter -So, Internal branch length deviation from the half of the height of the ingroup. Please, check your sampling settings and try again\n");
        is_error=1;
    }
    if (get_sampling(ind_per_sp)<1)
    {
        fprintf(stderr,"\n\tImproper value sampling the parameter -I, Number of individuals per species. Please, check your sampling settings and try again\n");
        is_error=1;
    }
    if (get_sampling(nl_trees)<1)
    {
        fprintf(stderr,"\n\tImproper value sampling the parameter -Rl, Number of locus trees per species tree. Please, check your sampling settings and try again\n");
        is_error=1;
    }
    if (is_sampling_set(b_rate)&&(get_sampling(b_rate)<0))
    {
        fprintf(stderr,"\n\tImproper value sampling the parameter -Pb, Locus tree birth rate prior. Please, check your sampling settings and try again\n");
        is_error=1;
    }
    if (is_sampling_set(d_rate)&&(get_sampling(d_rate)>get_sampling(b_rate)))
    {
        fprintf(stderr,"\n\tImproper value sampling the parameter -Pd, Locus tree death rate prior. Please, check your sampling settings and try again\n");
        is_error=1;
    }
    if (is_sampling_set(t_rate)&&(get_sampling(t_rate)<0))
    {
        fprintf(stderr,"\n\tImproper value sampling the parameter -Pt, Locus tree transfer rate prior. Please, check your sampling settings and try again\n");
        is_error=1;
    }
    if (is_sampling_set(gc_rate)&&(get_sampling(gc_rate)<0))
    {
        fprintf(stderr,"\n\tImproper value sampling the parameter -Pg, Locus tree gene conversion rate prior. Please, check your sampling settings and try again\n");
        is_error=1;
    }
    if (get_sampling(Ne)<2)
    {
        fprintf(stderr,"\n\tImproper value sampling the parameter -P, Haploid population size. Please, check your sampling settings and try again\n");
        is_error=1;
    }
    if (is_sampling_set(alpha_s)&&(get_sampling(alpha_s)<0))
    {
        fprintf(stderr,"\n\tImproper value sampling the parameter -Hs, Alpha parameter for the gamma species-specific heterogeneity. Please, check your sampling settings and try again\n");
        is_error=1;
    }
    if (is_sampling_set(alpha_l)&&(get_sampling(alpha_l)<0))
    {
        fprintf(stderr,"\n\tImproper value sampling the parameter -Hl, Alpha parameter for the gamma locus-specific heterogeneity. Please, check your sampling settings and try again\n");
        is_error=1;
    }
    if (is_sampling_set(salpha_g)&&(get_sampling(salpha_g)<0))
    {
        fprintf(stderr,"\n\tImproper value sampling the parameter -PP, Alpha hyperprior for the gamma gene tree branch-specific heterogeneity. Please, check your sampling settings and try again\n");
        is_error=1;
    }
    if (is_sampling_set(mu)&&(get_sampling(mu)<=0))
    {
        fprintf(stderr,"\n\tImproper value sampling the parameter -U, Substitution rate. Please, check your sampling settings and try again\n");
        is_error=1;
    }
    if (is_sampling_set(gen_time)&&(get_sampling(gen_time)<=0))
    {
        fprintf(stderr,"\n\tImproper value sampling the parameter -G, Generation time. Please, check your sampling settings and try again\n");
        is_error=1;
    }
    
    if((get_sampling(sb_rate)<get_sampling(sd_rate))&&(get_sampling(bds_leaves)>0))
    {
        fprintf(stderr,"\n\tERROR!!! The BDSA algorithm of the conditioned birth death simulation process of the species tree does not allow birth rates less than the death rate\n");
        
        is_error=1;
    }
    else if((get_sampling(bds_leaves)>0) && (get_sampling(bds_leaves)<min_lleaves))
    {
        fprintf(stderr,"\n\tERROR:The minimum number of locus tree nodes is bigger than the number of desired nodes of the species trees. Please, check the -Ll and -Sl parameters.\n");
        fflush(stderr);
        is_error=1;
    }
    
    if (is_error)
    {
        fflush(stderr);
        return SETTINGS_ERROR;
    }
    else
        return NO_ERROR;
}

long int CheckSampledSettingsLloop(sampling_unit lb_rate, sampling_unit ld_rate, sampling_unit lt_rate, sampling_unit lgc_rate, sampling_unit lalpha_g)
{
    int is_error=0;
    
    if (is_sampling_set(lb_rate)&&(get_sampling(lb_rate)<0))
    {
        fprintf(stderr,"\n\tImproper value sampling the parameter -Lb, Locus tree birth rate. Please, check your sampling settings and try again\n");
        is_error=1;
    }
    if (is_sampling_set(ld_rate)&&(get_sampling(ld_rate)>get_sampling(lb_rate)))
    {
        fprintf(stderr,"\n\tImproper value sampling the parameter -Ld, Locus tree death rate. Please, check your sampling settings and try again\n");
        is_error=1;
    }
    if (is_sampling_set(lt_rate)&&(get_sampling(lt_rate)<0))
    {
        fprintf(stderr,"\n\tImproper value sampling the parameter -Lt, Locus tree transfer rate. Please, check your sampling settings and try again\n");
        is_error=1;
    }
    if (is_sampling_set(lgc_rate)&&(get_sampling(lgc_rate)<0))
    {
        fprintf(stderr,"\n\tImproper value sampling the parameter -Lg, Locus tree gene conversion rate. Please, check your sampling settings and try again\n");
        is_error=1;
    }
    if (is_sampling_set(lalpha_g)&&(get_sampling(lalpha_g)<0))
    {
        fprintf(stderr,"\n\tImproper value sampling the parameter -Ph, Alpha hyperparameter for the gamma gene tree branch-specific heterogeneity. Please, check your sampling settings and try again\n");
        is_error=1;
    }
    
    if (is_error)
    {
        fflush(stderr);
        return SETTINGS_ERROR;
    }
    else
        return NO_ERROR;
}

long int CheckSampledSettingsGloop(sampling_unit alpha_g)
{
    int is_error=0;
    
    if (is_sampling_set(alpha_g)&&(get_sampling(alpha_g)<0))
    {
        fprintf(stderr,"\n\tImproper value sampling the parameter -HG, Alpha parameter for the gamma gene tree branch-specific heterogeneity. Please, check your sampling settings and try again\n");
        is_error=1;
    }
    if (is_error)
    {
        fflush(stderr);
        return SETTINGS_ERROR;
    }
    else
        return NO_ERROR;
}

///@}

void PrintXStringError(char **string, int x, char *errormsg)
{
    int i;
    printf("\n\nSETTINGS ERROR\n");
    for (i=1; i<x; ++i)
    {
        printf("%s ", *(string+i));
    }
    printf(" %s \n", errormsg);
    fflush(stdout);
    
}