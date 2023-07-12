/*
  To Do list:
    o  Check on the myaffsi's to see if things are really long or long long.
    o  Add "functions"
        - currently, we initialize by passing j_p, which gives pointers to the congruences
        - instead, we would have a short global list of the current vectors
        - for each j_p, we first initialize the current function
            then we loop until we have no entries left
            each time through the loop, we have to pass info to the function so it can recreate
            the previous state, compute the next vector, save the new state, load the info into
            some vector, and pass back a key.  The key either signals no more entries, or it will
            be used to pass back to the function so it can find its state.
        - functions need to keep a static array with states from primes which have called it.
            Maybe this should be a list over all possible p.
        - the standard function is just a front end to the list of entries from the input file.  When
           a prime calls, we set a pointer to the beginning of the list.  The static state info is where
           we currently are in the list.  We return done when this list has been exhausted.
*/


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
  register long i, count[10], ii;


  char      s[512];
  long long high[9], low[9], llval, newbnd;
  long      ltop, vp;
  double    btmp, boundn;
  double    newval, lval, savelval;

  long long bigmodulus;
  

  /* End of variable declarations */

  init(2000000,500000);  /* take some memory for the stack */

  welcome();
  myinitprimes();
  getprimelist(argc, argv);
  initmygens();
  /*  register_known_funcs(); */

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

                  if(verbosity>1) {
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

                  bigmodulus=1;
                  for(i=0; i<totprimes; i++) {
                    bigmodulus *= intpow(plist[i],
                                         caselist[cases[i][j_p[i]].start][0]);
                  }
                  printf("\nBig modulus is %lld\n", bigmodulus);
                  

                  /*
                    gettimeofday(&time1, &tzn);
                  */
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


                    do {
                    /* This is generic */


                    /*
                      if(count[1]==3) {
                      printf(" dmax = %7.2f low|high = %lld|%lld; arc=%lld; shift=%lld,centers=%7.2f\n",
                             Dmaxroot,low[2],high[2],arc[2],shift[2], centers[2]);
                    }
                    */
                    
                    /*  ------------  Counter 2 ---------- */

                    for(count[2] = 0; count[2] < bigmodulus/arc[2]; count[2]++) {
                      A2 = arc[2]*count[2] + shift[2];
                      myaffsi( A2, coefs[2]);

                      /* This is generic */

                      if(init_all_primes(DEG,A2)) {
                        

                      do {
                      
                      /*  ------------  Counter 6 ---------- */

                      for(count[DEG] =(long) 0;
                          count[DEG]< bigmodulus/arc[DEG]; count[DEG]++) {
                        An = arc[DEG]*count[DEG] + shift[DEG];

                        myaffsi(An, coefs[DEG]);
                        
                        if(init_all_primes(3,An)) {
                        do {

                        /*		printf("\nA6 = "); plint(An); */
                        if(1) {
                          
                          /* It mostly does side effects, and captures
                             the sign of An (whoppee!!). */
                          
			
                          /*  ------------  Counter 3 ---------- */


                          for(count[3] = 0; count[3]< bigmodulus/arc[3];
                              count[3]++) {
                            A3 = arc[3] * count[3] + shift[3];
                            myaffsi(A3, coefs[3]);
                            
                            if(init_all_primes(4,A3)) {
                            do {

                            centers[4] = A1*A1*A1*A1/4.0 -A1*A1*A2 +
                              A2*A2/2.0 +A1*A3;



                            /* Old inequality for s3  -
                               now lower bound for s4 */
                            btmp = (Dmaxroot1- S2)*(Dmaxroot1-S2) -
                              2.0*S3*S3/(S2+Dmaxroot1);

                            if(btmp > pbounds[4]) btmp = pbounds[4];


                            low[4] =  (long long) S2c_lower(-pbounds[4],4);
                            high[4] = (long long) S2c_upper(btmp,4);

                            /* Inherited from An bound */
                            savelval =  -An*pbounds[2]/2.0;
                            if(An>0) {
                              newbnd = (long long) ceil((savelval - shift[4])
                                                        /(1.0*arc[4]));
                              if(low[4]<newbnd) low[4] = newbnd;
                            } else {
                              newbnd = (long long) floor((savelval - shift[4])
                                                         /(1.0*arc[4]));
                              if(high[4]>newbnd) high[4] = newbnd;
                            }


                            /*  ------------  Counter 4 ---------- */
                            for(count[4] = 0; count[4]< bigmodulus/arc[4];
                                count[4]++) {
                              A4 = arc[4] * count[4] + shift[4];
                              myaffsi(A4, coefs[4]);
                              if(init_all_primes(5,A4)) {
                              do {

                              centers[5] = A1*A4 + A2*A3 - A1*A1*A3 - A1*A2*A2
                                + A1*A1*A1*A2 - A1*A1*A1*A1*A1/5.0;
                              
                              if(0) {
                                printf("Funny square root: %3.2f\n",
                                       ((2.0*A4)/An+pbounds[2]));
                                printf("%3.2f is the bound, An = ", savelval);
                                plint(An);
                                printf(" and A4 = ");
                                plint(A4);
                                printf("\n");
                                printf("Is %3.2f plus %3.2f positive?\n",
                                       (2.0 * A4)*(An),An*An*pbounds[2]);
                              } else {
                                
                                low[5] =  (long long) S2c_lower(-pbounds[5],5);
                                high[5] = (long long) S2c_upper(pbounds[5],5);

                                /* |A5| <= lval */
                                lval = sqrt(((2.0*A4)/An+pbounds[2]))*An*sgn_n;
                                /* |A5| <= newval */
                                newval = pbounds[1]*mylabs(An);
                                if(lval > newval) lval = newval;

                                llval = (long long) A2c_upper(lval,5,arc[5]);
                                if(llval< high[5]) high[5] = llval;
			
                                llval = (long long) A2c_lower(-lval,5,arc[5]);
                                if(llval> low[5]) low[5] = llval;

                                /*
                                  if((2.0*A4)/(1.0*A6) > pbounds[2]) {
                                  Comment				 OK to hit these values, but not in between 
                                  bot5 = (long long) floor((((-sqrt((2.0*A4*A6)-pbounds[2]*A6*A6)-
                                  shift[5])/arc[5])-invertme5)/divmod5);
                                  top5 = (long long) ceil((((sqrt((2.0*A4*A6)-pbounds[2]*A6*A6)-shift[5])
                                  /arc[5])-invertme5)/divmod5);
                                  if((bot5 < low[5]) && (top5 > high[5])) printf("Savings of %1ld\n",
                                  (long)		(high[5]-low[5]+1));
                                  else if((top5 <= low[5]) || (bot5 >= high[5])) ;
                                  else if((bot5>=low[5]) && (top5<=high[5]))
                                  printf("Savings of %1ld\n", (long) (top5-bot5-1));
                                  else if(top5<high[5]) printf("Savings of %1ld\n", (long) (top5-low[5]));
                                  else if(bot5>low[5]) printf("Savings of %1ld\n", (long) (high[5]-bot5));
                                  else printf("Confused\n");
                                  }
                                */

                                /*  ------------  Counter 5 ---------- */
                                for(count[5] =0; count[5]< bigmodulus/arc[5];
                                    count[5]++) {
                                  A5 = arc[5] * count[5] + shift[5];
                                  myaffsi(A5, coefs[5]);

                                  /*  ------  The actual testing  --- */

                                  printf("{%lld, %lld, %lld, %lld, %lld, %lld},", A1, A2, A3, A4, A5, A6);
                                  
                                  

                                  avma = ltop;
                                } /* Counter 5 */
                              }  /* Safety check */
                              } while(next_all_primes(5,A4));
                              } /* init */
                            } /* Counter 4 */
                            } while(next_all_primes(4,A3));
                            } /* init */
                          } /* Counter 3 */
                          /* } A6 satisfying a basic congruence */
                        } /* Is A6 != 0 */
                        } while(next_all_primes(3,A6));
                        } /* init 6 */
                      } /* Counter 6 */
                      } while(next_all_primes(6,A2));
                      } /* init 2 */
                    } /* Counter 2 */
                    } while(next_all_primes(2,A1));
                    } /* if init 1 */
                    } /* Counter 1 */
                  }  while (next_all_primes(1,0));
                  
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
  gettimeofday(&time2, &tzn);
  timetot = time2.tv_sec - time1.tv_sec;

  
  /*  system("cat .message | mail -s \"Sextic searches\" jj@asu.edu"); */
  
  printf("\nDone.");
  printf("\n");
  return 0;
}

