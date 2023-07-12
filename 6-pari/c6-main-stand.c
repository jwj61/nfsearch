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
  long          sgn_n;
  register long i, count[10];


  char      s[512];
  long long high[9], low[9], llval, newbnd, bigB;
  long      ltop;
  long		a1set=0,a2set=0,a1isset=0,a2isset=0,a3set=0,a3isset=0,justcount=0;
  double    btmp, boundn;
  double    newval, lval, savelval;
  GEN       bigbgen;

  /* End of variable declarations */

  pari_init(2000000,500000);  /* take some memory for the stack */

  /* welcome(); */
  myinitprimes();
  totprimes = 0;
  initmygens();

  bigB = atoi(argv[1]);
  bigbgen = stoi((long) bigB);
  sprintf(s, "StdSearch_%Ld", bigB);
  verbosity = 3;
  if(argc>2) {
	a1isset=1;
	a1set=atoi(argv[2]);
    sprintf(s, "StdSearch_%Ld.%ld", bigB, a1set);
	if(a1set<0) {
	  justcount=1;
	  verbosity=0;
	}
  }
  if(argc>3 && justcount==0) {
	a2isset=1;
	a2set=atoi(argv[3]);
    sprintf(s, "StdSearch_%Ld.%ld.%ld", bigB, a1set, a2set);
  }
  if(argc>4 && justcount==0) {
	a3isset=1;
	a3set=atoi(argv[4]);
    sprintf(s, "StdSearch_%Ld.%ld.%ld.%ld", bigB, a1set, a2set,a3set);
  }
  if(justcount) {
    a2isset=0; a1isset=0;
  }


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
  Dmaxroot += 1.5; /* a_1^2/n */
  if(justcount==0) printf("Max T2 for %Ld is %5.3f\n", bigB, Dmaxroot);

  for(i=0; i<8; i++) {
    arc[i] = 1;
    shift[i] = 0;
  }


  ltop = avma;
  
  /* Do new initialization stuff */
  if(a1isset) {
	low[1]=a1set;
	high[1]=a1set;
  } else {
    low[1]=0;
    high[1]= (int) (DEG/2);
  }
  
  /* ------------ Counter 1 ---------- */
  /* for(count[1] = 0; 2*count[1]* arc[1] <= DEG; count[1]++) { */
  for(count[1] = low[1]; count[1] <= high[1]; count[1]++) {
	A1 = count[1];
	myaffsi(A1, coefs[1]);

	for(i=3;i<=DEG; i++) {
	  bounds[i] =  fpow(Dmaxroot + A1*A1/DEG*1.0, i/2.0);
	  /* printf("\nBounds[%1d] = %6.2f ", i, bounds[i]); */
	}

	boundn = fpow(Dmaxroot1/DEG, DEG/2.0);
	centers[2] = A1*A1* (DEG-1)/(2.0*DEG);

	/* This is generic */
                      
  if(a2isset) {
	low[2]=a2set;
	high[2]=a2set;
  } else {
	high[2] = (long long) S2c_upper(Dmaxroot, 2);
	low[2] = (long long) S2c_lower(- Dmaxroot, 2);
  }


	/*
	  if(count[1]==3) {
	  printf(" dmax = %7.2f low|high = %lld|%lld; arc=%lld; shift=%lld,centers=%7.2f\n",
	  Dmaxroot,low[2],high[2],arc[2],shift[2], centers[2]);
	  }
	*/
    
	/*  ------------  Counter 2 ---------- */
	
	for(count[2] = low[2]; count[2] <= high[2]; count[2]++) {
	  A2 = count[2];
	  myaffsi( A2, coefs[2]);
      
	  /* This is generic */
	  
	  centers[3] = A1*A2 - A1*A1*A1/3.0;
	  
                        
	  
	  /*    if(A1==1 && A2 == -3 && shift[DEG] == 232)
			printf("[%5.3f, %lld, %5.3f: %5.2f, %lld, %lld] ",
			boundn, arc[DEG],
			((-boundn)-shift[DEG])/arc[DEG],
			A2c_lower(-boundn, DEG, arc[DEG]),
			(long long) A2c_lower(-boundn, DEG, arc[DEG]),
			low[DEG]);
	  */
                             
  if(a3isset) {
	low[DEG]=a3set;
	high[DEG]=a3set;
  } else {
	  high[DEG] = (long long) A2c_upper(boundn, DEG, arc[DEG]);
	  low[DEG] = (long long) A2c_lower(-boundn, DEG, arc[DEG]);
  }
		  if(justcount) {
			printf("[%ld, %ld, %Ld, %Ld]\n", count[1], count[2], low[DEG], high[DEG]);
		    low[DEG]=10; high[DEG]=-10;
		  }
                      
	  if(verbosity>2) {
		printf("\nLoop at (%1ld,%1ld) ---> (%1d,%lld), %1d range: %lld to %lld\n",
			   count[1], count[2],
			   (int) ((DEG-2*shift[1])/(2*arc[1])),
			   high[2], DEG,low[DEG],high[DEG]);
		fflush(stdout);
	  }
                      
	  /*  ------------  Counter 6 ---------- */

	  for(count[DEG] =(long) low[DEG];
		  count[DEG]<= high[DEG]; count[DEG]++) {
		An = arc[DEG]*count[DEG] + shift[DEG];
		myaffsi(An, coefs[DEG]);
		
	  if(verbosity>2) {
		printf(".");fflush(stdout);
	  }
		/*		printf("\nA6 = "); plint(An); */
		if(( An != 0) && (mylabs(An) <= boundn)) {
		  
		  /* It mostly does side effects, and captures
			 the sign of An (whoppee!!). */
		  
		  sgn_n = pohst_thing();
		  
          
		  low[3] =  (long long) S2c_lower(-pbounds[3],3);
		  high[3] = (long long) S2c_upper(pbounds[3],3);
		  
		  if(A1==0) low[3] = 0;
          
		  /*  ------------  Counter 3 ---------- */
		  
		  for(count[3] = low[3]; count[3]<= high[3]; count[3]++) {
			A3 = count[3];
			myaffsi(A3, coefs[3]);
			
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
			for(count[4] = low[4]; count[4]<= high[4]; count[4]++) {
			  A4 = count[4];
			  myaffsi(A4, coefs[4]);

			  centers[5] = A1*A4 + A2*A3 - A1*A1*A3 - A1*A2*A2
				+ A1*A1*A1*A2 - A1*A1*A1*A1*A1/5.0;
			  
			  if(((2.0*A4)/An+pbounds[2]) < 0) {
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
				for(count[5] =low[5]; count[5]<= high[5]; count[5]++) {
				  A5 = count[5];
				  myaffsi(A5, coefs[5]);
				  
				  /*  ------  The actual testing  --- */
				  
				  
				  polcount++;
                  
				  /* take absolute value on the way */
				  if(polisirreducible(a)) {
					td=absi( discf(a));
					if(1 || gcmp(td, bigbgen)<=0) {
					  discfcount++;
					  pari_fprintf(afile, "%Ps\n", polredabs0(a,0));
					  fflush(afile);
						} /* small disc*/
				  } /* Irreducible */
                 
				  avma = ltop;
				} /* Counter 5 */
			  }  /* Safety check */
			} /* Counter 4 */
		  } /* Counter 3 */
		  /* } A6 satisfying a basic congruence */
		} /* Is A6 != 0 */
	  } /* Counter 6 */
	} /* Counter 2 */
  } /* Counter 1 */
  if(verbosity>2) {
	if(polcount!= 0) {printf("\nSearched ");
	  plint(polcount);
	  printf(" polynomials ");
	  printf(".");
	  printf("Poly disc check %1ld times,", postdisc);
	  printf("Field disc check %1ld times.\n", discfcount);
	  
	} else printf("0");
    
  }

  gettimeofday(&time2, &tzn);
  timetot = time2.tv_sec - time1.tv_sec;
  timetot = timetot+1-1;
  
  char t[50];
  getlogfilename(t);
  if((logfile = fopen(t, "a")) == NULL) {
    printf("Cannot open %s\n",t);
    exit(1);
  }
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
  fprintf(afile, "\nHi John,\n\nI have finished sextic search %s on hobbes.\n",
          outfilename);
  fclose(afile);
  /*  system("cat .message | mail -s \"Sextic searches\" jj@asu.edu"); */
  
  if(verbosity>2) {
    printf("\nDone.");
    printf("\n");
  }
  return 0;
}

