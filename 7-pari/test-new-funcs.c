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
  long          min_disc;
  register long i, count[10], ii=1;


  char      s[512];
  long long high[9], low[9], llval;
  long      ltop, vp;
  double    btmp, boundn;
  double    newval, lval;

  long long bigmodulus;


  /* End of variable declarations */

  init(2000000,500000);  /* take some memory for the stack */

  verbosity = 0;
  
  /*  welcome(); */
  myinitprimes();
  getprimelist(argc, argv);
  initmygens();
  getsearchcases();

  makeoutfilename(s, j_p);

  /* Make a little file which says what case we are currently doing */
  
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

                  if(verbosity>10) {
                    printf("\nThis is case number %ld", j_p[0]+1);
                    for(i=1; i<totprimes; i++)
                      printf(" x %ld", j_p[i]+1);
                    printf(" of %1ld", maxp[0]+1);
                    for(i=1; i<totprimes; i++)
                      printf(" x %ld", maxp[i]+1);
                    printf("\n");
                    fflush(stdout);
                  }
                  
                  
                  initcase(j_p); /* This will do the init 0 things */
                  /*                  showorefactors(); */  /* Does this still make sense? */

                  /*    makeoutfilename(s, j_p);
                        afile = fopen(s, "a"); */


                  /*
                    gettimeofday(&time1, &tzn);
                  */
                  bigmodulus=1;
                  for(i=0; i<totprimes; i++) {
                    bigmodulus *= intpow(plist[i],
                                         /*    caselist[cases[i][j_p[i]].start][0] */
                                         1
                                         );
                  }
                  printf("\nBig modulus is %lld\n", bigmodulus); 
                  
                  ltop = avma;

                  /* Do new initialization stuff */

                  /* ------------ Counter 1 ---------- */
                  do {
                    /*                    printf("shift[1] = %lld;arc[1]=%lld", shift[1], arc[1]); */
                    for(count[1] = 0;
                      2*count[1]* arc[1] <= (DEG-2*shift[1]); count[1]++) {
                    A1 = arc[1]*count[1] + shift[1];
                    myaffsi(A1, coefs[1]);
                    
                    if(init_all_primes(2,A1)) {

                    /* This is |S_i| */
                    for(i=3;i<=DEG; i++) {
                      bounds[i] =  fpow(Dmaxroot + A1*A1/DEG*1.0, i/2.0);
                      /* printf("\nBounds[%1d] = %6.2f ", i, bounds[i]); */
                    }

                    boundn = fpow(Dmaxroot1/DEG, DEG/2.0);
                    centers[2] = A1*A1* (DEG-1)/(2.0*DEG);

                    do {
                    /* This is generic */

                      
                    /*  ------------  Counter 2 ---------- */

                    for(count[2] = 0; count[2] < bigmodulus/arc[2]; count[2]++) {
                      A2 = arc[2]*count[2] + shift[2];
                      myaffsi( A2, coefs[2]);

                      /* This is generic */

                      if(init_all_primes(DEG,A2)) {
                        
                      do {
                      
                      
		
                      /*  ------------  Counter 7 ---------- */

                      for(count[DEG] =0;
                          count[DEG]<  bigmodulus/arc[DEG]; count[DEG]++) {
                        An = arc[DEG]*count[DEG] + shift[DEG];
                        myaffsi(An, coefs[DEG]);

                        if(init_all_primes(3,An)) {
                        do {

                        /*		printf("\nA6 = "); plint(An); */
                        if(1) {
                          
                          /* It mostly does side effects, and captures
                             the sign of An (whoppee!!). */
                          
                          /*  ------------  Counter 3 ---------- */


                          for(count[3] =0 ; count[3]< bigmodulus/arc[3];
                              count[3]++) {
                            A3 = arc[3] * count[3] + shift[3];
                            myaffsi(A3, coefs[3]);

                            if(init_all_primes(4,A3)) {
                            do {
                            


                            /*  ------------  Counter 4 ---------- */

                            for(count[4] = 0; count[4]<bigmodulus/arc[4];
                                count[4]++) {
                              A4 = arc[4] * count[4] + shift[4];
                              myaffsi(A4, coefs[4]);

                              if(init_all_primes(5,A4)) {
                              do {

                                /*  ------------  Counter 5 ---------- */
                                
                                for(count[5] =0; count[5]<bigmodulus/arc[5];
                                    count[5]++) {
                                  A5 = arc[5] * count[5] + shift[5];
                                  myaffsi(A5, coefs[5]);

                                  if(init_all_primes(6,A5)) {
                              do {


                                    for(count[6] =0; count[6]<bigmodulus/arc[6];
                                        count[6]++) {
                                      A6 = arc[6] * count[6] + shift[6];
                                      myaffsi(A6, coefs[6]);
                                      
                                    
                                      /*  ------  The actual testing  --- */
                                      if(ii) {
                                        ii=0;
                                        printf("a={");
                                      } else {
                                        printf(",");
                                      }
                                      
                                      
                                      printf("{%lld, %lld, %lld, %lld, %lld, %lld,%lld}",
                                             A1, A2, A3, A4, A5, A6,A7);

                                      
                                      avma = ltop;
                                    } /* Counter 6 */
                              } while(next_all_primes(6,A5));
                                  } /* init*/
                                } /* Counter 5 */
                                /* }   Safety check */
                                } while(next_all_primes(5,A4));
                                  } /* init*/
                            } /* Counter 4 */
                            } while(next_all_primes(4,A3));
                                  } /* init*/
                          } /* Counter 3 */
                          /* } A6 satisfying a basic congruence */
                        } /* Is A6 != 0 */
                        } while(next_all_primes(3,An));
                                  } /* init*/
                      } /* Counter 6 */
                      } while(next_all_primes(DEG,A2));
                                  } /* init*/
                    } /* Counter 2 */
                    } while(next_all_primes(2,A1));
                    } /* init*/
                    } /* Counter 1 */
                  }   while (next_all_primes(1,0));
                  if(verbosity>2) {
                    if(polcount!= 0) {printf("\nSearched ");
                    plint(polcount);
                    printf(" polynomials for %1ld^%1ld", plist[0],
                           cases[0][j_p[0]].discexp);
                    for(i=1; i<totprimes; i++)
                      printf(" %1ld^%1ld", plist[i],
                             cases[i][j_p[i]].discexp);
                    printf(".");
                    printf("Poly disc check %1ld times,", postdisc);
                    printf("Field disc check %1ld times.\n", discfcount);

                    } else printf("0");
                    
                  }
                  /*gettimeofday(&time2, &tzn);
                    timetot = time2.tv_sec - time1.tv_sec;
                    fprintf(logfile, "(%1ld:%d",plist[0], (int) j_p[0]-1+'a');
                    for(i=1; i<totprimes; i++)
                    fprintf(logfile, ",%1ld:%d", plist[i], (int) j_p[i]-1+'a');
                    fprintf(logfile,"), %1ld polynomials, %1ld seconds, \
                    %3.2f T2, %1ld past disc, %1ld past T2, %1ld discf\n",
                    polcount, timetot, Dmaxroot, postdisc, t2count, discfcount);
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
  printf("};\n");
  
  return 1;
}
