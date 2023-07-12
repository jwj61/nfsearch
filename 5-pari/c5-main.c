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
  long long high[9], low[9], llval, newbnd;
  long      ltop, vp;
  double    boundn, btmp, savelval, lval, newval;

  /* End of variable declarations */

  pari_init(2000000,500000);  /* take some memory for the stack */

  welcome();
  myinitprimes();
  getprimelist(argc, argv);
  initmygens();
  getsearchcases();

  verbosity = 3;
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
                      
		
                      /*  ------------  Counter 5 ---------- */

                      for(count[DEG] =(long) low[DEG];
                          count[DEG]<= high[DEG]; count[DEG]++) {
                        An = arc[DEG]*count[DEG] + shift[DEG];
                        myaffsi(An, coefs[DEG]);

/*
if(A1 == 0 && A2 == 70 && An == 1364) {
printf("testing ");plint(A1);printf(",");plint(A2);printf(",");plint(An);
printf("\n");
}
*/
                        if(init_all_primes(3,An)) {

                        do {
                        /*		printf("\nA5 = "); plint(An); */
                        if(( An != 0) && (mylabs(An) <= boundn)) {
                          
                          /* It mostly does side effects, and captures
                             the sign of An (whoppee!!). */
                          
                          sgn_n = pohst_thing();
                          
                          low[3] =  (long long) S2c_lower(-pbounds[3],3);
                          high[3] = (long long) S2c_upper(pbounds[3],3);
                          

                          /* Inherited from An bound */
                          savelval =  -An*pbounds[2]/2.0;
                          if(An>0) {
                            newbnd = (long long) ceil((savelval - shift[3])
                                                      /(1.0*arc[3]));
                            if(low[3]<newbnd) low[3] = newbnd;
                          } else {
                            newbnd = (long long) floor((savelval - shift[3])
                                                       /(1.0*arc[3]));
                            if(high[3]>newbnd) high[3] = newbnd;
                          }
                          
                          
                          /*  ------------  Counter 3 ---------- */

                          for(count[3] = low[3]; count[3]<= high[3];
                              count[3]++) {
                            A3 = arc[3]*count[3] + shift[3];
                            myaffsi(A3, coefs[3]);

                            if(init_all_primes(4,A3)) {
                            do {

                            centers[4] = A1*A1*A1*A1/4.0 -A1*A1*A2 +
                              A2*A2/2.0 +A1*A3;

                            if(((2.0*A3)/An+pbounds[2]) < 0) {
                              printf("Funny square root: %3.2f\n",
                                     ((2.0*A3)/An+pbounds[2]));
                              printf("%3.2f is the bound, An = ", savelval);
                              plint(An);
                              printf(" and A3 = ");
                              plint(A3);
                              printf("\n");
                              printf("Is %3.2f plus %3.2f positive?\n",
                                       (2.0 * A3)*(An),An*An*pbounds[2]);
                              } else {
                                btmp = (Dmaxroot1- S2)*(Dmaxroot1-S2) -
                                  2.0*S3*S3/(S2+Dmaxroot1);

                                if(btmp > pbounds[4])  btmp = pbounds[4];
                                low[4] =  (long long) S2c_lower(-pbounds[4],4);
                                high[4] = (long long) S2c_upper(btmp,4);

                                /* |A4| <= lval */
                                lval = sqrt(((2.0*A3)/An+pbounds[2]))*An*sgn_n;
                                /* |A5| <= newval */
                                newval = pbounds[1]*mylabs(An);
                                if(lval > newval) lval = newval;

                                llval = (long long) A2c_upper(lval,4,arc[4]);
                                if(llval< high[4]) high[4] = llval;
			
                                llval = (long long) A2c_lower(-lval,4,arc[4]);
                                if(llval> low[4]) low[4] = llval;



                                /*  ------------  Counter 4 ---------- */

                                for(count[4] = low[4]; count[4]<= high[4];
                                    count[4]++) {
                                  A4 = arc[4] * count[4] + shift[4];
                                  myaffsi(A4, coefs[4]);
                                  
                                  /*  ------------  The actual testing   ---------- */
                                  
                                  
                                  polcount++;
                                  /* take absolute value on the way */
                                  td=absi(discsr(a));
                                  
                                  if(!gcmp0(td)) {
                                    tdtmp = gcopy(td);

/* Note, sqrtdust is not used, and current function sets it improperly */
                                    if(!ALTERNATING || Z_issquare(td)) {
                                    
                                    for(i=0; i<totprimes; i++) {
                                      vp = Z_pvalrem(td, gprimes[i], &tdtmp);
                                      gaffect(tdtmp, td);
                                    }
                                    if(gcmp1(td)) {
                                      if(isirreducible(a)) {
									    pari_fprintf(afile, "%Ps\n", polredabs0(a,0));
                                        fflush(afile);
									  }
                                    } else {
                                      if(Z_issquare(td)) {
                                        postdisc++;
/* outbrute(td);printf(",\n"); 
printf("%l\n", gissquarerem(td, &sqrtdust)); */
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
                                } /* Counter 4 */
                              } /* square root check */
                            } while(next_all_primes(4,A3));
                            }
                          } /* Counter 3 */
                        } /* Is A5 != 0 */
                        } while(next_all_primes(3,An));
                        }
                      } /* Counter 5 */
                      } while(next_all_primes(5,A2));
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

