#include "genpari.h"

#define big 10000000
#define giganto 9223372036854775807LL
#define S4    (A1*A1*A1*A1-4*A1*A1*A2 + 2* A2*A2 +4*A1*A3 - 4 * A4)
#define S3    (-A1*A1*A1+3*A1*A2-3*A3)
#define S2    (A1*A1-2*A2)
#define S5    (-A1*A1*A1*A1*A1 + 5*A1*A1*A1*A2 - 5*A1*A2*A2 - \
 5*A1*A1*A3 + 5*A2*A3 + 5*A1*A4 - 5*A5)
#define S6    (A1*A1*A1*A1*A1*A1 - 6*A1*A1*A1*A1*A2 + 9*A1*A1*A2*A2 -\
 2*A2*A2*A2 + 6*A1*A1*A1*A3 - 12*A1*A2*A3 + 3*A3*A3 - 6*A1*A1*A4 + \
 6*A2*A4 + 6*A1*A5 - 6*A6)
#define thef(m,i,x,K)  ((i)*fpow(K,(1.0*(m))/(i))*fpow(x,(i-DEG)*(1.0*(m))/(i))+(DEG-i)*fpow(x,1.0*(m)))
#define mylabs(a)   (a<0 ? -1*a : a)
/* #define t2tester      1.5*    */
#define t2tester    /* 1.1*   */

#define Dmaxroot1  (0.000001 + Dmaxroot + A1*A1/(1.0*DEG))

/* Convert bounds on A[i] to bounds on counter */
#define A2c_upper(B,iii,mmm)  floor(((B)-shift[iii])/mmm)
#define A2c_lower(B,iii,mmm)  ceil(((B)-shift[iii])/mmm)

/* This is really for DEC alphas, but it shouldn't hurt anyone.  */
#define fpow(a,b) pow((double) a,(double) b)

     /* These make the convention that centers are divided by -i, so that
        S_i/(-i) = A_i - centers[i].  B is a bound on -S_i.  */
     /* For reference, s[1]/-1 = a[1],
        s[2]/-2 = a[2] - a[1]^2/2,
        s[3]/-3 = a[3] - (a[1] a[2] - a[1]^3/3),
        s[4]/-4 = a[4] - (a[1]^4/4 - a[1]^2 a[2] + a[2]^2/2 + a[1] a[3])
        s[5]/-5 = a[5] - (a[1] a[4] + a[2] a[3] +a[1]^3 a[2] - a[1] a[2]^2
                   - a[1]^2 a[3] - a[1]^5/5)
        s[6]/-6 = a[6] - (a[1]^6/6 - a[1]^4 a[2] + (3/2) a[1]^2 a[2]^2
                   - a[2]^3/3 + a[1]^3 a[3] - 2 a[1] a[2] a[3]
                   + a[3]^2/2 - a[1]^2 a[4] + a[2] a[4] + a[1] a[5])
     */

/* Convert bounds on S[i] to bounds on counters. */

#define S2c_upper(B,iii)  A2c_upper((B)/(1.0*iii)+centers[iii],iii,arc[iii])
#define S2c_lower(B,iii)  A2c_lower((B)/(1.0*iii)+centers[iii],iii,arc[iii])

#if DEG==3
#define An    A3
#else
#if DEG==4
#define An    A4
#else
#if DEG==5
#define An    A5
#else
#if DEG==6
#define An    A6
#else
#if DEG==7
#define An    A7
#else
#if DEG==8
#define An    A8
#endif
#endif
#endif
#endif
#endif
#endif



     /*
			           ----- GLOBAL VARIABLES -----
		 */


		 /* New data structures first */

#define MAX_NUM_INDIV_CASES 1000000
#define MAX_NUM_CASES 50
#define COMMENT_CHAR '%'
     /* At most 7 primes, 9 coeffs.
        Entry #0 is c s.t. these coefs are modulo p^c.
        Entry #10 is for the input-file-code.
        The last entry (number 11), is the expo for the discrim. */

typedef struct {
  long discexp;
  long start;
  long end;
  long (*next_func)();
  int (*init_func)();
} case_data;

extern case_data cases[15][MAX_NUM_CASES];

extern long  caselist[MAX_NUM_INDIV_CASES][15];
extern long  casepntr;
/* extern long  casepointer[7]; */

/*
typedef struct 
{
  long (*next_func)();
  int (*init_func)();
  char name[10];
} registry_data;

extern registry_data registry[15];
extern int reg_cnt;
*/

extern long pdata[10][10][10];
extern long pflags[][8]; /* first by prime, then by i of a_i
                            0 = call init and next
                            -1 = call init only
                            1 = don't call anyone */

#define CALL_BOTH 0
#define CALL_INIT -1
#define CALL_NONE 1

extern long (*(f_next_list[]))();
extern int (*(f_init_list[]))();
extern long (*(pf_next_list[]))();
extern int (*(pf_init_list[]))();



extern double Dmaxroot;
extern GEN  coefs[10], td, tdtmp, gprimes[15], t2, c, a, sqrtdust;
extern GEN  gbig, gbigbig;  /* For conversions with long longs */
extern FILE  *logfile;
extern long  vc, maxexpo[10], ell;
extern long  polcount, postdisc, discfcount, t2count; /* Counters */
extern long  p_inps[15][25];  /* Which cases to search */
extern long  plist[15], maxp[15], totprimes, totprimesplus;
extern int   DO_ALL_CASES, ALTERNATING;
extern long  verbosity;

/*
  plist contains the primes
  maxp is a list of the number cases for each prime
  totprimes is the number of primes for screening
  totprimesplus also contains primes added because of the degree
*/

extern char outfilename[250];
extern long long A1, A2, A3, A4, A5, A6, A7, A8, shift[10], arc[10];
extern long slop;
extern int  do_slop;
extern GEN  slop_primes[40];
extern double centers[10];

/* A1..A6 are the coeffs.  shift and arc hold the data for this sub-case.
   centers are used to shift S_i bounds to a_i bounds */
extern double bounds[10], pbounds[10], mygamma[10];
extern int trueprimes[40];
/*
  Just a list the small primes
*/

/*
									 ----- FUNCTIONS -----
*/

void myinitprimes(void);
int checkprime(long);
long myreadint(void);
void getprimelist(int, char * []);
void getsearchcases(void);
void initmygens(void);
void spitsplittingmenu(void);
int iswild(long, long);
void initcase(long *);
void plint(long long int );
void plintf(FILE *, long long int);
long long intpow(long long, long);
long long truemod(long long, long long);
long long intpowmod(long long, long, long long);
long long gcdtrick(long long, long long);
long myval(long long, long);
double newt(double, long, double);
void myaffsi(long long int, GEN);
void showorefactors(void);
void welcome(void);
long pohst_thing(void);
void makeoutfilename(char *, long *);
void getlogfilename(char *);
void go_fetch_input(long, long, char *);
void splittingmenustring(char *, long);
long fmyreadint(FILE *);
long long minv(long long, long long);

/*  Functions for generating congruences */
#if DEG==7
int t7init(long, long, long);    long t7next(long, long);
int t61init(long, long, long);    long t61next(long, long);
int t52init(long, long, long);   long t52next(long, long);
int t511init(long, long, long);   long t511next(long, long);
int t43init(long, long, long);   long t43next(long, long);
int t421init(long, long, long);   long t421next(long, long);
int t4111init(long, long, long);  long t4111next(long, long);
int t331init(long, long, long);  long t331next(long, long);
int t322init(long, long, long);  long t322next(long, long);
int t3211init(long, long, long);  long t3211next(long, long);
int t31111init(long, long, long);  long t31111next(long, long);
int t2221init(long, long, long);  long t2221next(long, long);
int t22111init(long, long, long);  long t22111next(long, long);
int t211111init(long, long, long);  long t211111next(long, long);
#else
#if DEG==6
int t6init(long, long, long);    long t6next(long, long);
int t51init(long, long, long);   long t51next(long, long);
int t42init(long, long, long);   long t42next(long, long);
int t411init(long, long, long);  long t411next(long, long);
int t33init(long, long, long);  long t33next(long, long);
int t321init(long, long, long);  long t321next(long, long);
int t3111init(long, long, long);  long t3111next(long, long);
int t222init(long, long, long);  long t222next(long, long);
int t2211init(long, long, long);  long t2211next(long, long);
int t21111init(long, long, long);  long t21111next(long, long);
#else
#if DEG==5
int t5init(long, long, long);    long t5next(long, long);
int t41init(long, long, long);   long t41next(long, long);
int t32init(long, long, long);   long t32next(long, long);
int t311init(long, long, long);  long t311next(long, long);
int t221init(long, long, long);  long t221next(long, long);
int t2111init(long, long, long);  long t2111next(long, long);
#else
#if DEG==4
int t4init(long, long, long);    long t4next(long, long);
int t31init(long, long, long);    long t31next(long, long);
int t211init(long, long, long);   long t211next(long, long);
int t22init(long, long, long);   long t22next(long, long);
#else
#if DEG==3
int t3init(long, long, long);    long t3next(long, long);
int t21init(long, long, long);   long t21next(long, long);
#endif
#endif
#endif
#endif
#endif
int w_file_init(long, long, long);  long w_file_next(long, long);

long next_all_primes(long , long);
long init_all_primes(long , long);

int is_in_list(long, long []);

/*
int lookup_registry(char *);
void register_known_funcs();
*/
