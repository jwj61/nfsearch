#include <stdio.h>
#include <time.h>
#include "c-target.h"

static long loadcasepntr=0;
extern GEN bdisc, slopgen;

GEN bdisc, slopgen;


/*
		     	----- GLOBAL VARIABLES -----
*/

long pdata[10][10][10];
long pflags[10][8];
long (*(pf_next_list[10]))();
int (*(pf_init_list[10]))();
#if DEG==6
long (*(f_next_list[]))() = {t6next, t51next, t42next, t411next, t33next,
                             t321next, t3111next, t222next, t2211next,
                             t21111next,w_file_next};
int (*(f_init_list[]))() = {t6init, t51init, t42init, t411init, t33init,
                            t321init, t3111init, t222init, t2211init,
                            t21111init, w_file_init};
#else
#if DEG==8
long (*(f_next_list[]))() = {w_file_next};
int (*(f_init_list[]))() = {w_file_init};
#else
#if DEG==7
long (*(f_next_list[]))() = {t7next, t61next, t52next, t511next, t43next,
                             t421next, t4111next, t331next,t322next,
                             t3211next, t31111next, t2221next, t22111next,
                             t211111next,w_file_next};
int (*(f_init_list[]))() = {t7init, t61init, t52init, t511init, t43init,
                            t421init, t4111init, t331init,t322init,
                            t3211init, t31111init, t2221init, t22111init,
                            t211111init, w_file_init};
#else
#if DEG==5
long (*(f_next_list[]))() = {t5next, t41next, t32next, t311next, 
                             t221next, t2111next,w_file_next};
int (*(f_init_list[]))() = {t5init, t41init, t32init, t311init,
                            t221init, t2111init, w_file_init};
#else
#if DEG==4
long (*(f_next_list[]))() = {t4next, t31next, t22next, t211next, 
                             w_file_next};
int (*(f_init_list[]))() =  {t4init, t31init, t22init, t211init,
                             w_file_init};
#else
#if DEG==3
long (*(f_next_list[]))() = {t3next, t21next, w_file_next};
int (*(f_init_list[]))() = {t3init, t21init,  w_file_init};
#endif
#endif
#endif
#endif
#endif
#endif

double Dmaxroot;
GEN  coefs[10], td, tdtmp, gprimes[15], t2, c, a, sqrtdust;
GEN  gbig, gbigbig;  /* For conversions with long longs */
FILE  *logfile;
long  vc, maxexpo[10], ell;
long  polcount, postdisc, discfcount, t2count; /* Counters */
long  p_inps[15][25];  /* Which cases to search, just used for output */
long  plist[15], maxp[15], totprimes;
long  ore[10][10];
case_data cases[15][MAX_NUM_CASES];
long  caselist[MAX_NUM_INDIV_CASES][15];
long  casepntr;
/* long  casepointer[7]; */

/* How much screen output should we have.
   0 = none
   1 = Just end of sub-case summary
   2 = start and end of cases
   3 = loop lines
*/
long  verbosity=2;

/*
  plist contains the primes
  maxp is a list of the number cases for each prime
  totprimes is the number of primes for screening
*/

char outfilename[250];
int  DO_ALL_CASES, ALTERNATING=0;
long long A1, A2, A3, A4, A5, A6, A7, A8, shift[10], arc[10];
/* Extra factor for discrim. */
long   slop = 1;
int    do_slop = 1;
GEN    slop_primes[40];

double centers[10];
double bounds[10], pbounds[10], mygamma[10];
int    trueprimes[40];
/*
  Just a list the small primes
*/

void myinitprimes()
{
  trueprimes[0]=2; trueprimes[1]=3; trueprimes[2]=5;
  trueprimes[3]=7; trueprimes[4]=11; trueprimes[5]=13;
  trueprimes[6]=17; trueprimes[7]=19; trueprimes[8]=23;
  trueprimes[9]=29; trueprimes[10]=31; trueprimes[11]=37;
  trueprimes[12]=41; trueprimes[13]=43; trueprimes[14]=47;
  trueprimes[15]=53; trueprimes[16]=59; trueprimes[17]=61;
  trueprimes[18]=67; trueprimes[19]=71; trueprimes[20]=73;
  trueprimes[21]=79; trueprimes[22]=83; trueprimes[23]=89;
  trueprimes[24]=97; trueprimes[25]=101; trueprimes[26]=103;
  trueprimes[27]=107; trueprimes[28]=109; trueprimes[29]=113;
  trueprimes[30]=127; trueprimes[31]=131; trueprimes[32]=137;
  trueprimes[33]=139; trueprimes[34]=149; trueprimes[35]=151;
  trueprimes[36]=157; trueprimes[37]=163; trueprimes[38]=167;
  trueprimes[39]=173;
  
  /* We cheat and init the mygammas too */
  mygamma[1] = 1.0; mygamma[2] = 4.0/3.0; mygamma[3] = 2.0;
  mygamma[4] = 4.0; mygamma[5] = 8.0; mygamma[6] = 64.0/3.0;
  mygamma[7] = 64.0;
}

int checkprime(long inp)
{
  int jj, found=0;

  if(inp==0) return(0);
  if(inp==1) inp = 10;
  for(jj=0; jj<40 && !found; jj++)
    if((inp % trueprimes[jj]) == 0) {
      found=1;
      inp = inp/trueprimes[jj];
    }
  if(found)
    return(inp==1);
  else return(1);
}

long myreadint()
{
  long num=0, tmpc;

  while(!isdigit(tmpc=getchar())) ; /* Find a digit */
  num = tmpc - '0';
  while(isdigit(tmpc=getchar()))
    num = 10*num + tmpc - '0';
  return(num);
}

long fmyreadint(FILE *inf)
{
  long num=0, tmpc;

  while(!isdigit(tmpc=getc(inf))) ; /* Find a digit */
  num = tmpc - '0';
  while(isdigit(tmpc=getc(inf)))
    num = 10*num + tmpc - '0';
  return(num);
}

long fmyreadintnice(FILE *inf)
{
  long num=0, tmpc;

  while(!isdigit(tmpc=getc(inf))) ; /* Find a digit */
  num = tmpc - '0';
  while(isdigit(tmpc=getc(inf)))
    num = 10*num + tmpc - '0';
  ungetc(tmpc, inf);
  return(num);
}

long long
intpow(long long base, long exp)
{
  long long result;
  long cntr;

  result = 1;
  for(cntr=0; cntr< exp; cntr++) result *= base;
  return(result);
}

/* Loads plist with a list of desired primes, either from the command
	 line, or by asking the user. */

void getprimelist(int argc, char * argv[] )
{
  int jj, ii;

  DO_ALL_CASES = (argv[0][strlen(argv[0])-1] == 'a');
  ALTERNATING = (argv[0][0] == 'a');
  if(argc > 1) {
    jj=1;
    for(ii=1; ii<argc; ii++, jj++) {
      plist[jj-1] = atoi(argv[ii]);
      if(!checkprime(plist[jj-1])) {
        if(plist[jj-1]==0) {
          jj--;
          do_slop = 1;
        } else {
          printf("I do not think that %1ld is prime!\n\n", plist[jj-1]);
          exit(1);
        }
      }
    }
    totprimes = jj-1;
  } else {
    printf("\nHow many primes are we talking about? ");
    fflush(stdout);
    totprimes = myreadint();
    for(jj=0; jj<totprimes; jj++) {
      printf("What is prime number %1d? ", jj+1);
      plist[jj] = myreadint();
      if(!checkprime(plist[jj])) {
        printf("I do not think that %1ld is prime!\n\n", plist[jj]);
        exit(1);
      }
    }
  }
  if(do_slop) {
    printf("\nWhat is the slop factor? ");
    fflush(stdout);
    slop = (long) myreadint();
	/* New meaning of slop -- it is artificial upper bound for T2 bound */
	/*
    if(slop > trueprimes[39]) {
      printf("\nSlop is too big.\n\n");
      exit(1);
    }
	*/
  }
}

void getsearchcases()
{
  long j;
  /*  long linecount; */
  int ii;
  char t[50], outs[50];

  getlogfilename(t);
  if((logfile = fopen(t, "a")) == NULL) {
    printf("Cannot open %s\n",t);
    exit(1);
  }
  outfilename[0] = '\0';
  casepntr=0;
  
  if(DO_ALL_CASES) { /* We are not ready for this anyway! */
    for(ii=0; ii<totprimes; ii++) {
      /* if(plist[ii]>DEG)
         sprintf(s, "%1d-facpat.tame",DEG);
         else */
      /*				sprintf(s, "input-%1d_%1ld", DEG, plist[ii]);
                                        if((fp = fopen(s, "r")) == NULL) {
                                        printf("Cannot open %s\n", s);
                                        exit(1);
                                        }
                                        linecount=0;
                                        while(getc(fp) != '-') {
                                        while((j=getc(fp)) != '\n')
                                        ;
                                        linecount++;
                                        }
                                        for(j=0; j<linecount; j++)
                                        p_inps[ii][j] = j+1;
                                        maxp[ii] = linecount-1; */
    }
  } else { /* We are not doing all cases here */
    for(ii=0; ii<totprimes; ii++) {
      printf("              For p = %1ld\n", plist[ii]);
      spitsplittingmenu();
      while(!(isalpha(j=getchar()) || j=='*')) ;
      ungetc(j, stdin);
      /* maxp is the number of cases for a given prime, and p_inps tells us
         what the cases mean. */
      maxp[ii]= -1;
      if(ii>0) strcat(outfilename, ".");
      sprintf(outfilename, "%s%1ld", outfilename, plist[ii]);
      while((j=getchar()) != '\n')
        if((j >= 'a' && j <='z') || j == '*') {
          go_fetch_input(ii, j - 'a', outs);
          /* This will have to be fixed to adjust for the name thing */
          sprintf(outfilename, "%s%s", outfilename, outs);
/*           casepointer[ii+1] = casepointer[ii]+maxp[ii]+1; */
        }
    }
  }
  /* We are preparing for several types of splittings.  Now, each
     splitting is potentially subdivided, and we would have to ask
     about those too. */
  for(ii=totprimes; ii<10; ii++)
    maxp[ii] = 0;
  printf("I have to search %1ld ", maxp[0]+1);
  for(ii=1; ii<totprimes; ii++)
    printf("x %1ld ", maxp[ii]+1);
  printf("cases\n");
}

void initmygens()
{
  int i;

  for(i=0; i<totprimes; i++) gprimes[i] = stoi(plist[i]);
  i=0;
  /* Another change for slop */
  /*
  while(trueprimes[i]<=slop) {
    slop_primes[i] = stoi(trueprimes[i]);
    i++;
  }
  */
  gbig = stoi(big);
  c = gp_read_str("x");
  vc = varn(c);
  a = cgetg(DEG+3, t_POL); /*a[1]=evalsigne(1)+evallgef(DEG+3);setvarn(a,vc); */
  setsigne(a, 1);
  setvarn(a,vc);
  for(i=0; i<=DEG; i++) coefs[i] = cgeti(12);
  for(i=0; i<=DEG; i++) affsi((long) 1, coefs[i]);
  for(i=0; i<=DEG; i++) a[DEG+2-i] = (long) coefs[i];
  /* td = cgeti(32);   The poly discr */
  /* tdtmp = cgeti(32);  The poly discr */
  sqrtdust = cgeti(32);  /* The non-six part of the poly discr */
}


/* This should do the CRT thing, and load the result into arc[i] and shift[i]
	 For now, all arc[i] will be equal. */
void initcase(long *casevec)
{
  long long maxdisc=1;
  int i;

  for(i=1; i<=DEG; i++) centers[i] = 0;

  for(i=0; i<totprimes; i++) {
    pf_next_list[i] = cases[i][casevec[i]].next_func;
    pf_init_list[i] = cases[i][casevec[i]].init_func;
    pf_init_list[i](i, 0, casevec[i]);
  }
  /*  init_all_primes(0,0);
   replaced with direct calls above*/
  
  /* This does shift and arc */
  init_all_primes(1,0);

  for(i=0; i<totprimes; i++) {
    maxdisc *= intpow((long long) plist[i],
                      (long long) cases[i][casevec[i]].discexp);
  }
  
  /* Dmaxroot = 1.0 *maxdisc/(1.0* DEG ); */
  bdisc = stoi(maxdisc);
  slopgen = stoi(slop);
  bdisc = gmul(bdisc,slopgen);
  Dmaxroot = 1.0 *slop* maxdisc/(1.0* DEG );  /* Old slop usage */
  Dmaxroot *= mygamma[DEG-1];
  /*  End of stuff I probably need */
  /* This needs mygamma values put in as well */
  Dmaxroot = t2tester fpow(Dmaxroot, 1.0/(DEG-1));
  if(verbosity>1) {
    printf("\nMaximum T_2 is %2.3f\n\n", Dmaxroot);
    fflush(stdout);
  }
  
}

int is_in_list(long num, long list[]) {
  int j=0;

  if(list[0] == 0)
    return 1;
  
  while((list[j] != 0) && (list[j] != num))
    j++;
  return (list[j] == num);
}


void
go_fetch_input(long p_index, long wch, char * infilestr) {
  char buf[100], sms[10];
  FILE *infile;
  int nc, c, k, discexp, old_c, old_d, j, wildfile;
  long disclist[30];
  
  wildfile = 0;
  if((wch+'a') != '*') {
    splittingmenustring(sms, wch);
    sprintf(buf, "%s/input%1d_%s.%1ld", INFILES, DEG, sms, plist[p_index]);
    infilestr[0] = 'a'+wch;
    infilestr[1] = '\0';
    /* Look for discrim list next */
    while((c=getchar()) == ' ')
      ; /* Eat spaces, c now contains first non-space */
    j=0;
    ungetc(c, stdin); /* Undo extra reading */
    while(isdigit(c)){ /* Hey, it was a digit */
      disclist[j++] = fmyreadintnice(stdin);
      while(((c=getchar())==' ') || c == ',')
        ; /* ate more spaces and commas*/
      ungetc(c, stdin); /* Undo extra reading */
    }
    disclist[j] = 0;
    k=0;
    if(plist[p_index]<=DEG) {
      while(k<j) { /* Now adjust infilestr if there are special discs */
        sprintf(infilestr, "%s_%ld", infilestr, disclist[k]);
        k++;          
      }
    } else { /*Must be tame, so forget about disclist */
      j=0;
      disclist[j]=0;
    }
    
  } else {
    wildfile=1;
    printf("\nName your input file: ");
    while(isspace(c=getchar())) ;
    k=0;
    infilestr[k++] = c;
    while(!isspace(c=getchar())) {
      infilestr[k++] = c;
    }
    ungetc(c, stdin);
    infilestr[k] = '\0';
    /* infilestr gets appended to output file name (it is the letters part) */
    /*
      This might be a good place to look for numbers.
       Ignore them if p>=deg,
      add them to infilestr
      Look for comma separated list of numbers
      We will have to hold them in an array, (terminate with a 0)
      Need a function to see if a number is in the array
    */
    sprintf(buf, "input%1d_%s.%1ld", DEG, infilestr, plist[p_index]);
  }

  if((!iswild(wch, plist[p_index])) && wildfile==0)  { /* Tame, use function */
    (maxp[p_index])++;
    printf("Tame function, choosing function %ld\n", wch);
    cases[p_index][maxp[p_index]].next_func = f_next_list[wch];
    cases[p_index][maxp[p_index]].init_func = f_init_list[wch];
    cases[p_index][maxp[p_index]].init_func(p_index,0,0);
    cases[p_index][maxp[p_index]].discexp = pdata[p_index][0][0];
  } else { /* Use file */
    if((infile = fopen(buf, "r")) == NULL) {
      printf("\nCannot open file \"%s\".  I quit.\n\n", buf);
      exit(1);
    }
    (maxp[p_index])++;
    cases[p_index][maxp[p_index]].next_func = w_file_next;
    cases[p_index][maxp[p_index]].init_func = w_file_init;
    cases[p_index][maxp[p_index]].start = loadcasepntr;
    
    c = 1;       old_c = 1;
    discexp = 0; old_d = 0;
    nc=getc(infile);
    while(nc != EOF) {
      switch (nc) {
      case COMMENT_CHAR: /* Skip comment lines. */
        while(((nc=getc(infile)) != EOF) && (nc != '\n')) ;
        break;
      case EOF: break;
      case 'c': old_c = c;
        c=fmyreadint(infile);
        if(c!= old_c) { /* Start a new case? */
          if(loadcasepntr != cases[p_index][maxp[p_index]].start) {
            /* Yes, start new case */
            cases[p_index][maxp[p_index]].end = loadcasepntr;
            printf("Built a wild case with %ld values, d=%d and c=%d\n",
                   cases[p_index][maxp[p_index]].end-
                   cases[p_index][maxp[p_index]].start,discexp, old_c);
            (maxp[p_index])++;
            cases[p_index][maxp[p_index]].next_func = w_file_next;
            cases[p_index][maxp[p_index]].init_func = w_file_init;
            cases[p_index][maxp[p_index]].start = loadcasepntr;
            cases[p_index][maxp[p_index]].discexp = discexp;
          } 
          /* continue in this case */
        }
        break;
      case 'd': old_d = discexp;
        discexp=fmyreadint(infile);
        if(discexp!= old_d) { /* Start a new case? */
          if(loadcasepntr != cases[p_index][maxp[p_index]].start) {
            /* Yes, start new case */
            cases[p_index][maxp[p_index]].end = loadcasepntr;
            printf("Built a wild case with %ld values, d=%d and c=%d\n",
                   cases[p_index][maxp[p_index]].end-
                   cases[p_index][maxp[p_index]].start, old_d, c);
            (maxp[p_index])++;
            cases[p_index][maxp[p_index]].next_func = w_file_next;
            cases[p_index][maxp[p_index]].init_func = w_file_init;
            cases[p_index][maxp[p_index]].start = loadcasepntr;
          }
          cases[p_index][maxp[p_index]].discexp = discexp;
          /* continue in this case */
        }
        /* Skip over wrong discs */
        if((disclist[0]>0) && (plist[p_index]<=DEG) &&
           !is_in_list(discexp, disclist)) {
          while(((nc=getc(infile)) != EOF) &&  (nc != 'd')) {
            if(nc == COMMENT_CHAR) {
              while(((nc=getc(infile)) != EOF) && (nc != '\n'))
                ;
              ungetc(nc, infile);
            }
            if(nc == 'c') {
              c=fmyreadint(infile);
            }
          }
          ungetc(nc,infile);
        }
      case '\n': break;
      default:
        break;
      }
      if(isdigit(nc)) { /* line starting with a digit, load it */
        if(!ALTERNATING || ((discexp % 2) ==0 )) {
          /* loadcasepntr++; */
          ungetc(nc, infile);
          for(k=1; k<= DEG; k++)
            caselist[loadcasepntr][k]
                    = fmyreadint(infile);
          
          caselist[loadcasepntr][0] = c;
          caselist[loadcasepntr][11] = discexp;
          caselist[loadcasepntr][10] = wch;
          if(loadcasepntr> MAX_NUM_INDIV_CASES) {
            printf("\nToo many cases! %ld\n\n",loadcasepntr);
            exit(1);
          }
          loadcasepntr++; 
        }
      }
      if(nc != EOF) nc = getc(infile);
    }
    cases[p_index][maxp[p_index]].end = loadcasepntr;
    if(loadcasepntr != cases[p_index][maxp[p_index]].start) {
      printf("Built a wild case with %ld values, d=%d and c=%d\n",
             cases[p_index][maxp[p_index]].end-
             cases[p_index][maxp[p_index]].start, discexp, c);
    } else { /* Last case had no entries */
      maxp[p_index]--;
    }
  } 
}


/* For printing out long longs */
void plint(long long nn)
{

  long long int  bigbig;
  long long int  n1, n2, n3;

  /* The compiler won't allow a constant this big!  In fact, since
     big is a constant, it has type int and we have to cast its type */
  bigbig = ((long long int) big) * big;

  /* Take care of negative numbers */
  if(nn<0) {printf("-"); nn = -nn;}
  if(nn<0) {
    printf("Cannot correct nn!\n");
  }
  

  /* Chop it into long size pieces */
  n1 = nn/bigbig;
  n2 = (nn % bigbig)/big;
  n3 = nn % big;

  if(n1 == 0) { if(n2 == 0) printf("%ld", (long) n3); else
    printf("%ld%.7ld", (long) n2, (long) n3);} else
      printf("%ld%.7ld%.7ld", (long) n1,(long) n2,(long) n3);
}

void plintf(FILE *outfile, long long int nn)
{

#define big 10000000

  long long int  bigbig;
  long long int  n1, n2, n3;

  /* The compiler won't allow a constant this big!  In fact, since
     big is a constant, it has type int and we have to cast its type */
  bigbig = ((long long int) big) * big;

  /* Take care of negative numbers */
  if(nn<0) {fprintf(outfile, "-"); nn = -nn;}

  /* Chop it into long size pieces */
  n1 = nn/bigbig;
  n2 = (nn % bigbig)/big;
  n3 = nn % big;

  if(n1 == 0) { if(n2 == 0) fprintf(outfile, "%ld", (long) n3); else
    fprintf(outfile, "%ld%.7ld", (long) n2, (long) n3);} else
      fprintf(outfile, "%ld%.7ld%.7ld", (long) n1,(long) n2,(long) n3);
}

/* I don't trust the mod operator in C to return non-negative
   results  */

long long truemod(long long a, long long m)
{
  long long ans = a % m;
  if(ans<0)
    ans += m;
  return ans;
}

/* Multiplicative inverse of a modulo b.
   Assumes that 0<a<b and that the inverse exists
*/
long long minv(long long a, long long b) {
  long long c[2], d[2], r, q, m = b;
  int j=0;

  a = a % b;
  c[0] = 1;
  d[0] = 0;
  c[1] = a/b;
  d[1] = 1;
  r = a % b;
  while(r > 0) {
    q = a/b;
    c[j] = q*c[1-j]+c[j];
    d[j] = q*d[1-j]+d[j];
    a = b;
    b = r;
    r = a % b;
    j = 1-j;
  }
  
  return (j ? c[1-j] : m-c[1-j]);
}

long long
intpowmod(long long bse, long expo, long long modulus)
{
  long j;
  long long ans=1;

  for(j=1; j<=expo; j++) ans = (bse * ans) % modulus;
  return(ans);
}

/* Assumes aa and bb>0, and relatively prime.  Returns the
   coefficient of aa */
long long
gcdtrick(long long aa, long long bb)
{
  long long r;

  r = aa % bb;
  if(r == 0) return(0);
  else return((1-gcdtrick(bb, r)*bb)/r);
}

/* Computes the valuation.  No big deal. */
long myval(long long n, long p)
{
  if((n % p) == 0) return(1+myval(n/p, p)); else return(0);
}

/* The parameters should be set for this application of
   Newton's method. */

double newt(double inval, long tt, double Kpow)
{
  double ypow,opow,oopow;

  opow = fpow(inval, ((2.0*DEG)/tt) - 2.0);
  oopow = opow*inval;
  ypow = oopow*inval;
  /* printf("    In newt: top = %f, bottom = %f\n",Kpow*tt + (6.0-tt)*ypow - 
     Dmaxroot*opow, (12/tt)*(6.0-tt)*opow*inval- Dmaxroot*((12/tt)-2)*opow/inval); */

  return( (Kpow*tt + (DEG-tt)*ypow - Dmaxroot1*opow)/
          ( ((2.0*DEG)/tt)*(DEG-tt)*oopow -
            Dmaxroot1*(((2.0*DEG)/tt)-2.0)*opow/inval));
}


void myaffsi(long long int innum, GEN whereto)
{
  long long int temp;
  long sign, ltop;
	
  ltop = avma;
  if(innum<0) {sign = -1; innum = -1*innum;} else sign = 1;
	
  temp=(innum/1073741824);
	
  affsi((long) temp/1073741824, whereto);
  gshiftz(whereto, 30, whereto);
  temp = temp % 1073741824;
  gaddz(stoi((long) temp), whereto, whereto);
  gshiftz(whereto, 30, whereto);
  gaddz(stoi((long) (innum % 1073741824)), whereto, whereto);
  if(sign == -1) gmulz(stoi((long) -1), whereto, whereto);
  avma = ltop;
}

void showorefactors()
{
  int j;

  if(verbosity>1) {
    printf("Local conditions are:\n");
    for(j=1; j<=DEG; j++) {
      printf("a_%1d == %1ld (mod %1ld)\n", j, (long) shift[j], (long) arc[j]);
    }
  }
}
