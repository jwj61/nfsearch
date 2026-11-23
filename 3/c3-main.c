#include <stdio.h>
#ifdef ALPHA
#include <sys/time.h>
#else
/* #include <time.h> */
#include <sys/time.h>
#include <unistd.h>
#endif
#include "c-target.h"
#include "c-doout.h"

/* aim for at least this many groups of A2 congruence classes (for each A1)
   to share out among parallel workers */
#define PAR_GROUPS3 1024


int main(int argc, char * argv [])
{
  struct        timeval time1, time2;
  struct        timezone tzn;
  FILE          *afile;
  unsigned long timetot;
  long          sgn_n, j_p[10];
  register long i, count[10];


  char      s[512];
  long long high[9], low[9];
  long      ltop, vp, hgrp;
  double    boundn;

  /* End of variable declarations */

  char *pargs[20];
  int pall, pjobs, pnargs;

  pari_init(2000000,500000);  /* take some memory for the stack */

  welcome();
  myinitprimes();
  pnargs = par_options(argc, argv, pargs, &pall, &pjobs);   /* -a, -jN */
  getprimelist(pnargs, pargs);
  initmygens();
  if(pall) par_all_types();
  getsearchcases();

  makeoutfilename(s, j_p);

  /* Make a little file which says what case we are currently doing */
  afile = fopen("Current", "w");
  fprintf(afile, "Working on %s.\n", outfilename);
  fclose(afile);
  
  afile = fopen(s, "a");

  polcount = 0;
  postdisc = 0;
  t2count = 0;
  discfcount = 0;
  gettimeofday(&time1, &tzn);
  /* Loop through cases */

  if(par_start(pjobs, afile)) {   /* serial run or parallel worker */
  for(j_p[0]=0; j_p[0]<=maxp[0]; ++j_p[0]) {
    /*    printf("\n%1ld", j_p[0]+1); fflush(stdout); */
    for(j_p[1]=0; j_p[1]<=maxp[1]; ++j_p[1]) {
      /* printf("."); fflush(stdout); */
      for(j_p[2]=0; j_p[2]<=maxp[2]; ++j_p[2]) {
        for(j_p[3]=0; j_p[3]<=maxp[3]; ++j_p[3]) {
          for(j_p[4]=0; j_p[4]<=maxp[4]; ++j_p[4]) {
            for(j_p[5]=0; j_p[5]<=maxp[5]; ++j_p[5]) {
              for(j_p[6]=0; j_p[6]<=maxp[6]; ++j_p[6]) {
                for(j_p[7]=0; j_p[7]<=maxp[7]; ++j_p[7]) {
                  for(j_p[8]=0; j_p[8]<=maxp[8]; ++j_p[8]) {

                  if(verbosity>1) {
                    printf("\nThis is case number %ld", j_p[0]+1);
                    for(i=1; i<totprimes; i++)
                      printf(" x %ld", j_p[i]+1);
                    printf(" of %1ld", maxp[0]+1);
                    for(i=1; i<totprimes; i++)
                      printf(" x %ld", maxp[i]+1);
                    printf("\n");
                  }
                  
                  
                  initcase(j_p); /* This will do the init 0 things */
                  qf_setup(j_p);  /* modular discriminant filter */
                  /*                  showorefactors(); */  /* Does this still make sense? */

                  /*    makeoutfilename(s, j_p);
                  afile = fopen(s, "a"); */


                  /*
                    gettimeofday(&time1, &tzn);
                  */
                  ltop = avma;

                  /* ------------ Counter 1 ---------- */
                  do {
                  for(count[1] = 0;
                      2*count[1]* arc[1] <= (DEG-2*shift[1]); count[1]++) {
                    A1 = arc[1]*count[1] + shift[1];
                    myaffsi(A1, coefs[1]);

                    /* This is |S_i| */
                    for(i=3;i<=DEG; i++) {
                      bounds[i] =  fpow(Dmaxroot + A1*A1/(1.0*DEG), i/2.0);
                      /* printf("\nBounds[%1d] = %6.2f ", i, bounds[i]); */
                    }

                    boundn = fpow(Dmaxroot1/DEG, DEG/2.0);
                    centers[2] = A1*A1*(DEG-1)/(2.0*DEG);

                    if(init_all_primes(2,A1)) {

                    /* Parallel runs share out the congruence classes for A2.
                       In degree 3 with tame primes there are about prod(p) of
                       them for each A1, and most contain no A2 in range, so
                       nearly all of the time goes into stepping through them
                       and combining the congruences by CRT.  The classes are
                       taken in groups (the classes at the primes j >= hgrp
                       fixed, see next_all_primes_range); a worker claims whole
                       groups, steps through and CRTs only the classes of its
                       own groups, and jumps over the others in one step.  The
                       order of the classes is the same as with
                       next_all_primes. */
                    hgrp = group_prime(2, PAR_GROUPS3);
                    do {                  /* groups */
                    if(par_claim()) {
                    do {                  /* classes in the group */
                    docrt(2);     /* arc[2], shift[2] for this class */

                    high[2] = (long long) S2c_upper(Dmaxroot, 2);
                    low[2] = (long long) S2c_lower(- Dmaxroot, 2);

                    /*  ------------  Counter 2 ---------- */

                    for(count[2] = low[2]; count[2] <= high[2]; count[2]++) {
                      A2 = arc[2]*count[2] + shift[2];
                      myaffsi( A2, coefs[2]);

                      /* This is generic */
                      if(init_all_primes(DEG,A2)) {
                      do {
                      
                      high[DEG] = (long long) A2c_upper(boundn, DEG, arc[DEG]);
                      /* New adjustment for deg odd! */
                      if(A1==0)
                        low[DEG] = 0;
                      else
                        low[DEG] = (long long)
                          A2c_lower(-boundn, DEG, arc[DEG]);
                      

                      centers[3] = A1*A2 - A1*A1*A1/3.0;

                      if(verbosity>2) {
                        printf("\nLoop at (%1ld,%1ld) ---> (%1d,",
                               count[1], count[2],
                               (int) ((DEG-2*shift[1])/(2*arc[1])));
                        plint(high[2]);
                        printf("), %1d range: ", DEG);
                        plint(low[DEG]);
                        printf(" to ");
                        plint(high[DEG]);
                        printf("\n");
                        fflush(stdout);
                      }
                      
		
                      /*  ------------  Counter 3 ---------- */

                      qf_newloop();  /* new innermost loop for the filter */
                      for(count[DEG] =(long) low[DEG];
                          count[DEG]<= high[DEG]; count[DEG]++) {
                        An = arc[DEG]*count[DEG] + shift[DEG];

                        /*		printf("\nA5 = "); plint(An); */
                        if(( An != 0) && (mylabs(An) <= boundn)) {
                          
                          /* It mostly does side effects, and captures
                             the sign of An (whoppee!!). */
                          
                                  
                                  /*  ------------  The actual testing   ---------- */
                                  
                                  
                                  polcount++;
                                  
                                  
                                  /* Cheap necessary test before the exact discriminant */
                                  
                                  
                                  if(!qf_test()) { avma = ltop; continue; }
                                  myaffsi(An, coefs[DEG]);  /* only for polynomials that pass the filter */
                                  /* take absolute value on the way */
                                  td=absi(discsr(a));
                                  
                                  if(!gcmp0(td)) {
                                    tdtmp = gcopy(td);

                                    if(!ALTERNATING || Z_issquare(td)) {
                                    
                                    for(i=0; i<totprimes; i++) {
                                      vp = Z_pvalrem(td, gprimes[i], &tdtmp);
                                      gaffect(tdtmp, td);
                                    }
                                    if(gcmp1(td)) {
                                      if(polisirreducible(a)) {
                                        pari_fprintf(afile, "%Ps\n", polredabs0(a,0));
                                        fflush(afile);
                                      }
                                    } else {
                                      if(Z_issquare(td)) {
                                        postdisc++;
#ifdef TTWO
                                        if(gcmpsg((long)
                                                  ceil(Dmaxroot),
                                                  gnorml2(roots(a,3))) > 0) {
                                          t2count++;
#endif
                                          if(polisirreducible(a)) {
                                            td=absi( discf(a));
                                            discfcount++;
                                            for(i=0; i<totprimes; i++) {
                                              vp = Z_pvalrem(td, gprimes[i], &tdtmp);
                                              gaffect(tdtmp, td);
                                            }
                                            if(gcmp1(td)) {
                                              pari_fprintf(afile, "%Ps\n", polredabs0(a,0));
                                              fflush(afile);
                                            } /* Discf check */
                                          } /* Irreducible */
#ifdef TTWO
                                        } /* T2 check */
#endif
                                      } /* square times ramifying primes */
                                    }
                                    } /* Alternating group check */
                                    
                                  }
                                  
                                  avma = ltop;
                        } /* Is An != 0 */
                      } /* Counter 3 */
                      } while(next_all_primes(3,A2));
                      }
                    } /* Counter 2 */
                    } while(next_all_primes_range(2, A1, 0, hgrp));          /* next class in the group */
                    } /* claimed group */
                    } while(next_all_primes_range(2, A1, hgrp, totprimes));  /* next group */
                    }
                  } /* Counter 1 */
                  }  while (next_all_primes(1,0));
                  if(verbosity>2) {
                    if(polcount!= 0) {
                      printf("\nSearched ");
                      plint(polcount);
                      printf(" polynomials for %1ld^%1ld", plist[0],
                             cases[0][j_p[0]].discexp);
                      for(i=1; i<totprimes; i++)
                        printf(" %1ld^%1ld", plist[i],
                               cases[i][j_p[i]].discexp);
                      printf(".  ");
                      printf("Poly disc check %1ld times, ", postdisc);
                      printf("Field disc check %1ld times.\n",
                             discfcount);
                    } else printf("0");
                    
                  }
                  /*                    
                  gettimeofday(&time2, &tzn);
                  timetot = time2.tv_sec - time1.tv_sec;
                  fprintf(logfile,
                          "(%1ld:%d",plist[0], (int) j_p[0]-1+'a');
                  for(i=1; i<totprimes; i++)
                    fprintf(logfile, ",%1ld:%d", plist[i], (int) j_p[i]-1+'a');
                  fprintf(logfile,"), %1ld polynomials, %1ld seconds, \
%3.2f T2, %1ld past disc, %1ld past T2, %1ld discf\n",
                          polcount, timetot, Dmaxroot, postdisc,
                          t2count, discfcount);
                  fflush(logfile);
                  */
                }
                }
              }
            }
          }
        }
      }
    }
  }
  par_finish(afile);               /* a worker exits here */
  }
  gettimeofday(&time2, &tzn);
  timetot = time2.tv_sec - time1.tv_sec;
  fprintf(logfile, "(%s",outfilename);
  fprintf(logfile,"), %1ld polynomials, %1ld seconds, \
%3.2f T2, %1ld past disc, %1ld past T2, %1ld discf\n",
          polcount, timetot, Dmaxroot, postdisc, t2count, discfcount);
  fflush(logfile);
  fclose(afile);
  printf("\nModular filter: %ld polynomials tested in %ld innermost loops, %ld rejected.", qf_tested, qf_loops, qf_rejected);
  printf("\nDone.");
  printf("\n");
  return 0;
}

