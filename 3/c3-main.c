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
  long      ltop, vp;
  double    boundn;

  /* End of variable declarations */

  pari_init(2000000,500000);  /* take some memory for the stack */

  welcome();
  myinitprimes();
  getprimelist(argc, argv);
  initmygens();
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
                      bounds[i] =  fpow(Dmaxroot + A1*A1/DEG*1.0, i/2.0);
                      /* printf("\nBounds[%1d] = %6.2f ", i, bounds[i]); */
                    }

                    boundn = fpow(Dmaxroot1/DEG, DEG/2.0);
                    centers[2] = A1*A1*(DEG-1)/(2.0*DEG);

                    if(init_all_primes(2,A1)) {
                    
                    do {
                    /* This is generic */

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

                      for(count[DEG] =(long) low[DEG];
                          count[DEG]<= high[DEG]; count[DEG]++) {
                        An = arc[DEG]*count[DEG] + shift[DEG];
                        myaffsi(An, coefs[DEG]);

                        /*		printf("\nA5 = "); plint(An); */
                        if(( An != 0) && (mylabs(An) <= boundn)) {
                          
                          /* It mostly does side effects, and captures
                             the sign of An (whoppee!!). */
                          
                                  
                                  /*  ------------  The actual testing   ---------- */
                                  
                                  
                                  polcount++;
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
                                      if(isirreducible(a)) {
                                        DUMPIT
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
                                          if(isirreducible(a)) {
                                            td=absi( discf(a));
                                            discfcount++;
                                            for(i=0; i<totprimes; i++) {
                                              vp = Z_pvalrem(td, gprimes[i], &tdtmp);
                                              gaffect(tdtmp, td);
                                            }
                                            if(gcmp1(td)) {
                                              DUMPIT
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
                    } while(next_all_primes(2,A1));
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
  gettimeofday(&time2, &tzn);
  timetot = time2.tv_sec - time1.tv_sec;
  fprintf(logfile, "(%s",outfilename);
  fprintf(logfile,"), %1ld polynomials, %1ld seconds, \
%3.2f T2, %1ld past disc, %1ld past T2, %1ld discf\n",
          polcount, timetot, Dmaxroot, postdisc, t2count, discfcount);
  fflush(logfile);
  fclose(afile);
  printf("\nDone.");
  printf("\n");
  return 0;
}

