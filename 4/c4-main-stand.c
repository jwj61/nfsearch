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
  long long high[9], low[9], llval, newbnd, bigB;
  long      ltop, vp;
  double    btmp, boundn;
  double    newval, lval, rdin;
  GEN		bigbgen;

  /* End of variable declarations */

  pari_init(2000000,500000);  /* take some memory for the stack */

  welcome();
  myinitprimes();
  initmygens();

  bigB = atoi(argv[1]);
  bigbgen = stoi((long) bigB);
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
                    for(count[1] = 0; 2*count[1] <= DEG; count[1]++) {
                    A1 = count[1];
                    myaffsi(A1, coefs[1]);

                    /* This is |S_i| */
                    for(i=3;i<=DEG; i++) {
                      bounds[i] =  fpow(Dmaxroot + A1*A1/DEG*1.0, i/2.0);
                      /* printf("\nBounds[%1d] = %6.2f ", i, bounds[i]); */
                    }

                    boundn = fpow(Dmaxroot1/DEG, DEG/2.0);
                    centers[2] = A1*A1*(DEG-1)/(2.0*DEG);

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

                      if(verbosity>2) {
                        printf("\nLoop at (%1ld,%1ld) ---> (%1d,%lld), %1d range: %lld to %lld\n",
                               count[1], count[2],
                               (int) ((DEG-2*shift[1])/(2*arc[1])),
                          high[2], DEG,low[DEG],high[DEG]);
                        fflush(stdout);
                      }
                      
                      /*  ------------  Counter 4 ---------- */

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

                          if(A1==0) low[3] = 0;
                          
                          /*  ------------  Counter 3 ---------- */


                          for(count[3] = low[3]; count[3]<= high[3];
                              count[3]++) {
                            A3 = arc[3] * count[3] + shift[3];
                            myaffsi(A3, coefs[3]);

                            
                                  /*  ------  The actual testing  --- */


                                  polcount++;
                                  
			if(polisirreducible(a)) {
				td=absi( discf(a));
				if(gcmp(td, bigbgen)<=0) {
					  discfcount++;
					  pari_fprintf(afile, "%Ps\n", polredabs0(a,0));
				} /* small disc*/

								  } /* Irreducible */
                                  avma = ltop;
                                  
                          } /* Counter 3 */
                        } /* An != 0 */
                        
                          
                        } /* Counter 4 */
                    } /* Counter 2 */
                    } /* Counter 1 */
                  if(verbosity>2) {
                    if(polcount!= 0) {printf("\nSearched ");
                    plint(polcount);
                    printf(" polynomials.  ");
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
  /*fprintf(logfile, "(%s",outfilename);
  fprintf(logfile,"), %1ld polynomials, %1ld seconds, \
%3.2f T2, %1ld past disc, %1ld past T2, %1ld discf\n",
          polcount, timetot, Dmaxroot, postdisc, t2count, discfcount);
  fflush(logfile); */
  fclose(afile);

  afile = fopen("Current", "a");
  fprintf(afile, "Done.\n");
  fclose(afile);

  
  afile = fopen(".message", "w");
  fprintf(afile, "\nHi John,\n\nI have finished sextic search %s on hobbes.\n",
          outfilename);
  fclose(afile);
  /*  system("cat .message | mail -s \"Sextic searches\" jj@asu.edu"); */
  
  printf("\nDone.");
  printf("\n");
  return 0;
}

