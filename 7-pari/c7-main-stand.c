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
  long          sgn_n;
  register long i, count[10];


  char      s[512];
  long long high[9], low[9], llval, bigB, bigB2;
  long      ltop;
  double    btmp, boundn;
  double    newval, lval;
  GEN		bigbgen;
  char outfilename[250];
  FILE      *logfile;
  logfile = fopen("eSeptic.log", "a");

  /* End of variable declarations */

  pari_init(2000000,500000);  /* take some memory for the stack */

  welcome();
  myinitprimes();
  totprimes = 0;
  initmygens();

  bigB = atoi(argv[1]);
  bigbgen = stoi((long) bigB);
  if(argc>2) {
    bigB2 = atoi(argv[2]);
	printf("Cutoff is %Ld\n", bigB2);
    bigbgen = stoi((long) bigB2);
  }
  verbosity = 3;

  sprintf(s, "StdSearch_%Ld", bigB);

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

  Dmaxroot = 1.0 * bigB/(1.0* DEG );
  Dmaxroot *= mygamma[DEG-1];
  Dmaxroot = t2tester fpow(Dmaxroot, 1.0/(DEG-1));
  Dmaxroot += 9./7.;
  printf("Max T2 for %Ld is %5.3f\n", bigB, Dmaxroot);

  for(i=0; i<8; i++) {
    arc[i] = 1;
    shift[i] = 0;
  }

                  ltop = avma;

                  /* Do new initialization stuff */

                  /* ------------ Counter 1 ---------- */
					for(count[1] = 0;
                      2*count[1] <= DEG; count[1]++) {
                    A1 = count[1];
                    myaffsi(A1, coefs[1]);

                    /* This is |S_i| */
                    for(i=3;i<=DEG; i++) {
                      bounds[i] =  fpow(Dmaxroot + A1*A1/DEG*1.0, i/2.0);
                      /* printf("\nBounds[%1d] = %6.2f ", i, bounds[i]); */
                    }

                    boundn = fpow(Dmaxroot1/DEG, DEG/2.0);
                    centers[2] = A1*A1* (DEG-1)/(2.0*DEG);

                    /* This is generic */

                    high[2] = (long long) S2c_upper(Dmaxroot, 2);
                    low[2] = (long long) S2c_lower(- Dmaxroot, 2);

                    /*  ------------  Counter 2 ---------- */

                    for(count[2] = low[2]; count[2] <= high[2]; count[2]++) {
                      A2 = count[2];
                      myaffsi( A2, coefs[2]);

                      /* This is generic */
                      centers[3] = A1*A2 - A1*A1*A1/3.0;

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
                        An = count[DEG];
                        myaffsi(An, coefs[DEG]);

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
                            A3 = count[3];
                            myaffsi(A3, coefs[3]);

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
                              A4 = count[4];
                              myaffsi(A4, coefs[4]);

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
                                  A5 = count[5];
                                  myaffsi(A5, coefs[5]);

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
                                      A6 = count[6];
                                      myaffsi(A6, coefs[6]);
                                      
                                    
                                      /*  ------  The actual testing  --- */
                                      

                                      polcount++;
									  td= discsr(a);
									  if(!gcmp0(td) /*  &&  Z_issquare(td) */ ) {
									    postdisc++;
									    if(isirreducible(a)) {
									      td=absi( discf(a)); 
									      if( gcmp(td, bigbgen)<=0 ) {
										  
									     /* if( 1 ) { */
									        discfcount++;
									        /* DUMPIT */
										    pari_fprintf(afile, "%Ps\n", polredabs0(a,0));
										    fflush(afile);
										  } /* small disc*/
									    } /* Irreducible */
									  } /* even */

                                      avma = ltop;
                                    } /* Counter 6 */
                                  } /* Sanity check */
                                } /* Counter 5 */
                            } /* Counter 4 */
                          } /* Counter 3 */
                        } /* Is A6 != 0 */
                      } /* Counter 6 */
                    } /* Counter 2 */
                    } /* Counter 1 */
                  if(verbosity>2) {
                    if(polcount!= 0) {printf("\nSearched ");
                    plint(polcount);
                    printf(" polynomials \n");
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
  gettimeofday(&time2, &tzn);
  timetot = time2.tv_sec - time1.tv_sec;
  fprintf(logfile, "(%Ld",bigB);
  fprintf(logfile,"), %1ld polynomials, %1ld seconds, \
%3.2f T2, %1ld past disc, %1ld past T2, %1ld discf\n",
          polcount, timetot, Dmaxroot, postdisc, t2count, discfcount);
  fflush(logfile);
  fclose(afile);

  afile = fopen("Current", "a");
  fprintf(afile, "Done.\n");
  fclose(afile);
  
  printf("\nDone.");
  printf("\n");
  return 1;
}
