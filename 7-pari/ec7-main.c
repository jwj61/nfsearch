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
  register long i, count[10], ii;


  char      s[512];
  long long high[9], low[9], llval;
  long      ltop, vp;
  double    btmp, boundn;
  double    newval, lval;

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

                    /* This is |S_i| */
                    for(i=3;i<=DEG; i++) {
                      bounds[i] =  fpow(Dmaxroot + A1*A1/DEG*1.0, i/2.0);
                      /* printf("\nBounds[%1d] = %6.2f ", i, bounds[i]); */
                    }

                    boundn = fpow(Dmaxroot1/DEG, DEG/2.0);
                    centers[2] = A1*A1* (DEG-1)/(2.0*DEG);

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
                      centers[3] = A1*A2 - A1*A1*A1/3.0;

                      if(init_all_primes(DEG,A2)) {
                      do {
                        
                      high[DEG] = (long long) A2c_upper(boundn, DEG, arc[DEG]);
                      low[DEG] = (long long) A2c_lower(-boundn, DEG, arc[DEG]);

                      /* Since DEG is odd, we normalize here. */
                      if(A1==0) low[DEG] = 0;
                     verbosity=5; 
                      if(verbosity>2) {
                        printf("\nLoop at (%1ld,%1ld) ---> (%1d,%lld), %1d range: %lld to %lld mults [%lld,%lld,%lld]\n",
                               count[1], count[2],
                               (int) ((DEG-2*shift[1])/(2*arc[1])),
                          high[2], DEG,low[DEG],high[DEG],arc[1],arc[2],arc[DEG]);
                        fflush(stdout);
                      }
                      
		
                      /*  ------------  Counter 7 ---------- */

                      for(count[DEG] =(long) low[DEG];
                          count[DEG]<= high[DEG]; count[DEG]++) {
                        An = arc[DEG]*count[DEG] + shift[DEG];
                        myaffsi(An, coefs[DEG]);

                        if(init_all_primes(3,An)) {
                        do {

                        /*		printf("\nA6 = "); plint(An); */
                        if(( An != 0) && (mylabs(An) <= boundn)) {
                          
                          /* It mostly does side effects, and captures
                             the sign of An (whoppee!!). */
                          sgn_n = pohst_thing();
			
                          
                          low[3] =  (long long) S2c_lower(-pbounds[3],3);
                          high[3] = (long long) S2c_upper(pbounds[3],3);
                          
                          /*  ------------  Counter 3 ---------- */


/* printf("[%lld|%lld]", low[3],high[3]);fflush(stdout); */
                          for(count[3] = low[3]; count[3]<= high[3];
                              count[3]++) {
                            A3 = arc[3] * count[3] + shift[3];
                            myaffsi(A3, coefs[3]);

                            if(init_all_primes(4,A3)) {
                            do {
                            
/* printf(".");fflush(stdout); */
                            centers[4] = A1*A1*A1*A1/4.0 -A1*A1*A2 +
                              A2*A2/2.0 +A1*A3;

                            /* Old inequality for s3  -
                               now lower bound for s4 */
                            btmp = (Dmaxroot1- S2)*(Dmaxroot1-S2) -
                              2.0*S3*S3/(S2+Dmaxroot1);

                            if(btmp > pbounds[4]) btmp = pbounds[4];


                            low[4] =  (long long) S2c_lower(-pbounds[4],4);
                            high[4] = (long long) S2c_upper(btmp,4);


                            /*  ------------  Counter 4 ---------- */

/* printf("/%lld|%lld/", low[4],high[4]);fflush(stdout); */
                            for(count[4] = low[4]; count[4]<= high[4];
                                count[4]++) {
                              A4 = arc[4] * count[4] + shift[4];
                              myaffsi(A4, coefs[4]);

                              if(init_all_primes(5,A4)) {
                              do {

/* printf("+");fflush(stdout); */
                              centers[5] = A1*A4 + A2*A3 - A1*A1*A3 - A1*A2*A2
                                + A1*A1*A1*A2 - A1*A1*A1*A1*A1/5.0;
                              
                              
                                
                              low[5] =  (long long) S2c_lower(-pbounds[5],5);
                              high[5] = (long long) S2c_upper(pbounds[5],5);

/* Inherited from An bound */
                              /*
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
                              */


                                
                                /*  ------------  Counter 5 ---------- */
                                
                                for(count[5] =low[5]; count[5]<= high[5];
                                    count[5]++) {
                                  A5 = arc[5] * count[5] + shift[5];
                                  myaffsi(A5, coefs[5]);

                                  if(init_all_primes(6,A5)) {
                              do {

                                  centers[6] = A1*A5 - A1*A1*A1*A1*A2 +
                                    1.5*A1*A1*A2*A2 - A2*A2*A2/3.0 +
                                    A1*A1*A1*A3 - 2*A1*A2*A3
                                    + A3*A3/2.0 - A1*A1*A4 + A2*A4
                                    + A1*A1*A1*A1*A1*A1/6.0;

                                  low[6] =  (long long)
                                    S2c_lower(-pbounds[6],6);
                                  high[6] = (long long)
                                    S2c_upper(pbounds[6],6);

                                  /* Sanity check before taking a square root */
                                  if((lval = (2.0*A5)/An+pbounds[2]) >= 0) {
                                    lval = sqrt(lval)*An*sgn_n;
                                  
                                    /* |A6| <= newval = | S_{-1} An| */
                                    newval = pbounds[1]*mylabs(An);
                                    if(lval > newval) lval = newval;

                                    llval = (long long) A2c_upper(lval,6,arc[6]);
                                    if(llval< high[6]) high[6] = llval;
                                    
                                    llval = (long long) A2c_lower(-lval,6,arc[6]);
                                    if(llval> low[6]) low[6] = llval;
                                    

                                    for(count[6] =low[6]; count[6]<= high[6];
                                        count[6]++) {
                                      A6 = arc[6] * count[6] + shift[6];
                                      myaffsi(A6, coefs[6]);
                                      
                                    
                                      /*  ------  The actual testing  --- */
                                      

                                      polcount++;
                                      /* take absolute value on the way */
                                      td = absi(discsr(a));
                                      if(!gcmp0(td) && Z_issquare(td)) {

                                        tdtmp = gcopy(td);
                                        for(i=0; i<totprimes; i++) {
                                          vp = pvaluation(td, gprimes[i], &tdtmp);
                                          gaffect(tdtmp, td);
                                        }
										/* Slop change */
										/*                                        if(gcmpgs(td, slop) <= 0) { */
                                        if(gcmpgs(td, 1) <= 0) {
                                          if(isirreducible(a)) {
											pari_fprintf(afile, "%Ps\n", polredabs0(a,0));
                                              }
                                        } else {
                                          min_disc = 1;
                                          i=0;
                                          while(min_disc <= /* slop */1 && trueprimes[i] <= /* slop */ 1) {
                                            vp = pvaluation(td, slop_primes[i], &tdtmp);
                                            gaffect(tdtmp, td);
                                            if((vp % 2) != 0)
                                              min_disc *= trueprimes[i];
                                            i++;
                                          }
                                          if(min_disc <= slop) {
                                            
                                            if(Z_issquarerem(td, &tdtmp)) {
                                              postdisc++;
#ifdef TTWO
                                              if(gcmpsg((long) ceil(Dmaxroot),
                                                        gnorml2(roots(a,3))) > 0) {
                                                t2count++;
#endif
                                                if(isirreducible(a)) {
                                                  td= absi(discf(a));
                                                  discfcount++;
                                                  for(i=0; i<totprimes; i++) {
                                                    vp = pvaluation(td, gprimes[i], &tdtmp);
                                                    gaffect(tdtmp, td);
                                                  }
                                              if(gcmpgs(td, slop) <= 0) {
												/* DUMPIT */
												pari_fprintf(afile, "%Ps\n", polredabs0(a,0));
                                                  } /* Discf check */
                                                } /* Irreducible */
#ifdef TTWO
                                              } /* T2 check */
#endif
                                            } /* Is there a square left over? */
                                          } /* can we get below slop with no small primes? */
                                        }
                                        
                                      }
                                      
                                      avma = ltop;
                                    } /* Counter 6 */
                                  } /* Sanity check */
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
  gettimeofday(&time2, &tzn);
  timetot = time2.tv_sec - time1.tv_sec;
  fprintf(logfile, "(%s",outfilename);
  fprintf(logfile,"), %1ld polynomials, %1ld seconds, \
%3.2f T2, %1ld past disc, %1ld past T2, %1ld discf\n",
          polcount, timetot, Dmaxroot, postdisc, t2count, discfcount);
  fflush(logfile);
  fclose(afile);

  afile = fopen("Current", "a");
  fprintf(afile, "Done.\n");
  fclose(afile);

  
  afile = fopen(".message", "w");
  fprintf(afile, "\nHi John,\n\nI have finished septic search %s.\n",
          outfilename);
  fclose(afile);
  /*  system("cat .message | mail -s \"Septic searches\" jj@asu.edu");*/
  
  printf("\nDone.");
  printf("\n");
  return 1;
}
